"""Verified union of many solids, one at a time.

A Boolean union over many solids often fails on the whole selection because
of one pair with flush faces, and when it does return a solid, "it returned
a solid" is not "it worked": the call can hand back one input with the
other silently dropped. The manual fix is to union one piece at a time
until the culprit shows up, nudge it a hair so the faces stop being
coincident, and union again. union_all() runs that routine and checks every
step:

1. Every input must be a valid closed solid with a volume. Otherwise
   nothing is united and the offenders are named.
2. The largest piece (by volume) seeds the result, and the others join it
   largest first, each one unioned into the result it touches.
3. A piece that refuses is set aside and retried plainly once everything
   else is in (it often unions fine once its neighbours have merged).
4. Only then is it nudged a small distance along +Z, -Z, +X, -X, +Y, -Y in
   turn, against the results it failed on. Still failing: it is left out
   and named.
5. Results that grew into each other without noticing (each only ever met
   the piece being placed) are merged pairwise at the end.

A union is accepted only when one valid closed solid comes back that still
covers both inputs' bounding boxes and whose volume adds up: at least as big
as either input and no bigger than both together. When one input adds
almost nothing (enclosed by the other, or dropped), the volume of their
common part decides which, or, when there is no common part to compute (a
Boolean intersection of a piece fully inside another can return nothing),
whether the piece lies inside the other. Pieces that touch nothing stay separate solids;
whether two pieces touch is decided by their intersection curves.

The geometry calls are passed in as an ops object (see Ops), so this module
needs no Rhino; meshwright.rhino.union provides the Brep version.
"""

from __future__ import absolute_import, division

from .report import Report

# Nudge directions in the order they are tried. Z first: flush faces are
# mostly parts sitting on top of each other.
NUDGES = (
    ("Z+", (0.0, 0.0, 1.0)),
    ("Z-", (0.0, 0.0, -1.0)),
    ("X+", (1.0, 0.0, 0.0)),
    ("X-", (-1.0, 0.0, 0.0)),
    ("Y+", (0.0, 1.0, 0.0)),
    ("Y-", (0.0, -1.0, 0.0)),
)

NUDGE_PER_TOLERANCE = 10.0  # default nudge: 0.01 at the usual tolerance of 0.001
VOLUME_SLACK = 0.01  # volumes must add up within 1% of the smaller input

MERGED, DISJOINT, FAILED = "merged", "disjoint", "failed"


class Ops(object):
    """The geometry calls union_all() needs. Subclass it for your own solids.

    problem(shape)          None for a valid closed solid, else a short reason ("open")
    volume(shape)           its volume, 0 when it cannot be computed
    box(shape)              (min corner, max corner) of its bounding box
    union(a, b)             list of solids the union gave back, or None when it failed
    intersection_volume(a, b)  volume of the common part, None when that failed
    touch(a, b)             True when they share surface or volume (intersection curves)
    encloses(outer, inner)  True when inner lies entirely inside outer
    moved(shape, vector)    a moved copy, or None when it could not be moved
    """

    def problem(self, shape):
        raise NotImplementedError

    def volume(self, shape):
        raise NotImplementedError

    def box(self, shape):
        raise NotImplementedError

    def union(self, a, b):
        raise NotImplementedError

    def intersection_volume(self, a, b):
        raise NotImplementedError

    def touch(self, a, b):
        raise NotImplementedError

    def encloses(self, outer, inner):
        raise NotImplementedError

    def moved(self, shape, vector):
        raise NotImplementedError


class Island(object):
    """A growing union result: starts as one input and absorbs the ones that unite with it."""

    def __init__(self, index, shape, volume):
        self.shape = shape
        self.volume = volume
        self.members = [index]  # input indices, in the order they joined


class UnionResult(object):
    """The outcome of union_all().

    islands  list of Island, one per output solid (seed order)
    failed   input indices that could not be united even after nudging
    nudged   {input index: direction label} for pieces that needed a nudge
    report   the Report
    """

    def __init__(self, report):
        self.report = report
        self.islands = []
        self.failed = []
        self.nudged = {}

    @property
    def status(self):
        return self.report.status

    @property
    def solids(self):
        return [island.shape for island in self.islands]

    def untouched(self, island):
        """True for an island that is one input exactly as given: it met nothing."""
        return len(island.members) == 1 and island.members[0] not in self.nudged


def union_all(shapes, ops, tolerance=0.001, nudge=None):
    """Unite shapes (a list of solids) through ops. Returns a UnionResult.

    nudge is the distance a refusing piece is moved (default 10 x tolerance).
    """
    report = Report("verified union")
    res = UnionResult(report)
    if not tolerance > 0.0:
        report.fail("the tolerance must be more than zero")
        return res
    if nudge is None:
        nudge = NUDGE_PER_TOLERANCE * tolerance
    if not shapes:
        report.fail("no solids to unite")
        return res

    volumes = []
    problems = []
    for i, shape in enumerate(shapes):
        if shape is None:
            problems.append("input %d is empty" % i)
            volumes.append(0.0)
            continue
        why = ops.problem(shape)
        if why:
            problems.append("input %d is %s" % (i, why))
        v = ops.volume(shape) if not why else 0.0
        volumes.append(v)
        if not why and not v > 0.0:
            problems.append("input %d has no volume that can be computed" % i)
    if problems:
        report.fail("nothing united: " + "; ".join(problems))
        return res
    if len(shapes) == 1:
        res.islands.append(Island(0, shapes[0], volumes[0]))
        report.note("one solid, nothing to unite")
        return res

    plan = _Plan(shapes, volumes, ops, tolerance, nudge, res)
    # Largest first: the main body seeds the result and the small parts join it.
    # Real volume, not bounding-box volume: long thin arms have the biggest
    # boxes, but the part they all meet in is the one each unites with cleanly.
    order = sorted(range(len(shapes)), key=lambda k: -volumes[k])
    res.islands.append(Island(order[0], shapes[order[0]], volumes[order[0]]))
    deferred = [i for i in order[1:] if not plan.place(i, False)]
    late = 0
    for i in deferred:
        if plan.place(i, True):
            late += 1
        else:
            res.failed.append(i)
    plan.merge_islands()
    _summarize(res, len(shapes), late, nudge)
    return res


class _Plan(object):
    def __init__(self, shapes, volumes, ops, tolerance, nudge, res):
        self.shapes = shapes
        self.volumes = volumes
        self.ops = ops
        self.tol = tolerance
        self.nudge = nudge
        self.res = res

    def place(self, i, last_chance):
        """Unite input i with an island it meets. True when it ended up in an
        island (merged, or touching nothing: its own island). False when it met
        an island and every attempt failed; nudges are tried on the last chance only."""
        shape, volume = self.shapes[i], self.volumes[i]
        reach = _inflate(self.ops.box(shape), 2.0 * self.tol)
        failed_on = []
        for island in [isl for isl in self.res.islands if _overlap(self.ops.box(isl.shape), reach)]:
            attempt = self.try_union(island, shape, volume)
            if attempt == MERGED:
                island.members.append(i)
                return True
            if attempt == FAILED:
                failed_on.append(island)
        if not failed_on:
            # Bounding boxes overlapped (or not) but the solids never met.
            self.res.islands.append(Island(i, shape, volume))
            return True
        if not last_chance:
            return False
        # Shift it a hair so flush faces stop being coincident, and try the islands it failed on.
        for label, d in NUDGES:
            moved = self.ops.moved(shape, (d[0] * self.nudge, d[1] * self.nudge, d[2] * self.nudge))
            if moved is None:
                continue
            for island in failed_on:
                if self.try_union(island, moved, volume) == MERGED:
                    island.members.append(i)
                    self.res.nudged[i] = label
                    return True
        return False

    def try_union(self, island, shape, volume):
        """One checked union attempt; on MERGED the island holds the result."""
        ops = self.ops
        out = ops.union(island.shape, shape)
        if (out and len(out) == 1 and out[0] is not None and not ops.problem(out[0])
                and self.covers(out[0], island.shape) and self.covers(out[0], shape)):
            u = ops.volume(out[0])
            if not self.adds_up(island.shape, island.volume, shape, volume, u):
                return FAILED
            island.shape = out[0]
            island.volume = u
            return MERGED
        if out and len(out) == 2 and not ops.touch(island.shape, shape):
            return DISJOINT
        return FAILED

    def covers(self, result, given):
        """The result's box (with some slack) still holds the input's box. The
        cheap tell of an input silently dropped."""
        lo, hi = _inflate(self.ops.box(result), 10.0 * self.tol)
        glo, ghi = self.ops.box(given)
        return all(lo[k] <= glo[k] and ghi[k] <= hi[k] for k in range(3))

    def adds_up(self, a, va, b, vb, union):
        """The union is at least as big as either input, no bigger than both,
        and each input adds some volume. An input that adds (almost) nothing is
        either enclosed by the other (fine) or was dropped; the common part
        tells which (union = a + b - common). A dropped piece can sit inside
        the other's bounding box, so covers() alone misses it. A piece fully
        inside the other can leave no common part to compute; then it must
        lie inside the other."""
        if not (va > 0.0 and vb > 0.0 and union is not None and union > 0.0):
            return False
        slack = VOLUME_SLACK * min(va, vb)
        if union < max(va, vb) - slack or union > va + vb + slack:
            return False
        if union - va >= slack and union - vb >= slack:
            return True
        common = self.ops.intersection_volume(a, b)
        if common is not None:
            return abs(va + vb - common - union) <= slack
        if union - va < slack:
            return self.ops.encloses(a, b)
        return self.ops.encloses(b, a)

    def merge_islands(self):
        """Merge islands that overlap, pairwise through the same checked union.
        One union call over many islands can drop a piece silently. A pair that
        will not merge stays as two sound solids. Passes repeat until none
        merges: an island that just grew may now reach one examined before."""
        islands = self.res.islands
        merged = True
        while merged and len(islands) > 1:
            merged = False
            for i in range(1, len(islands)):
                small = islands[i]
                reach = _inflate(self.ops.box(small.shape), 2.0 * self.tol)
                for big in islands[:i]:
                    if not _overlap(self.ops.box(big.shape), reach):
                        continue
                    if self.try_union(big, small.shape, small.volume) != MERGED:
                        continue
                    big.members.extend(small.members)
                    del islands[i]
                    merged = True
                    break
                if merged:
                    break


def _summarize(res, total, late, nudge):
    report = res.report
    lone = len([isl for isl in res.islands if res.untouched(isl)])
    united = total - len(res.failed) - lone
    solids = len(res.islands) - lone
    report.note("%d of %d inputs united into %d solid%s" % (
        united, total, solids, "" if solids == 1 else "s"))
    for k, isl in enumerate(res.islands):
        if not res.untouched(isl):
            report.note("solid %d: inputs %s" % (k, _list(sorted(isl.members))))
    if late:
        report.note("%d input%s united on the second pass" % (late, "" if late == 1 else "s"))
    if lone:
        report.warn("%d input%s nothing, left as they are: %s" % (
            lone, " touches" if lone == 1 else "s touch",
            _list(sorted(isl.members[0] for isl in res.islands if res.untouched(isl)))))
    if res.nudged:
        report.warn("moved %g to unite: %s" % (
            nudge, ", ".join("input %d (%s)" % (i, res.nudged[i]) for i in sorted(res.nudged))))
    if res.failed:
        report.messy("%d input%s could not be united even after nudging, left out: %s" % (
            len(res.failed), "" if len(res.failed) == 1 else "s", _list(sorted(res.failed))))


def _list(indices):
    return ", ".join(str(i) for i in indices)


def _inflate(box, d):
    lo, hi = box
    return (lo[0] - d, lo[1] - d, lo[2] - d), (hi[0] + d, hi[1] + d, hi[2] + d)


def _overlap(a, b):
    return all(a[0][k] <= b[1][k] and b[0][k] <= a[1][k] for k in range(3))
