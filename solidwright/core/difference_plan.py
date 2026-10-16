"""Verified Boolean difference: one solid cut by many cutters, every result checked.

A Boolean difference with many cutters has the same two ways to go wrong as
a union. The loud one: one awkward cutter is enough for the call over all
cutters to fail (Rhino returns nothing at all for a cutter wholly inside the
solid). The quiet one: the call returns something that looks fine and is
not: an open Brep (cutters touching each other along a line), or a closed
solid with a cutter skipped, so a hole is missing. cut() runs the routine
that finds and works around the bad cutters, and checks every result before
keeping it:

1. The solid must be a valid closed solid with a volume, or nothing is cut.
   A cutter that is not is left out and named (unless it misses the solid
   anyway).
2. Cutters whose bounding box misses the solid are skipped, and so are
   cutters that reach its box but never touch it.
3. All cutters are tried in one call, at the tolerance, then twice it, then
   a tenth of it.
4. When that fails, the cutters are split in halves, and halves again,
   until the ones that break the call are isolated; the others still go
   through in bulk.
5. A cutter that fails on its own goes through a rescue ladder: a sealed
   hollow when it lies wholly inside the solid, twice and a tenth of the
   tolerance, the job scaled up 10 times and back, the cutter grown by a
   few tolerances (over-cuts a hair, clears faces and lines where things
   only touch), and last a loose tolerance. Still failing: it is left out
   and named.

A result is accepted only when every piece is a valid closed solid, the
volume went down by no more than the cutters hold, the pieces stay inside
the solid's bounding box, and two sampled tests pass: points just inside
each cutter's surface that were inside the solid are no longer inside the
result (a skipped cutter fails this), and points just inside the result's
surface are not inside a cutter.

The geometry calls are passed in as an ops object (see Ops), so this module
needs no Rhino; solidwright.rhino.difference provides the Brep version.
"""

from __future__ import absolute_import, division

from . import vec
from .report import Report

VOLUME_SLACK = 0.001  # volumes are compared within 0.1% of the solid's volume
DEPTH = 10.0  # sample points are moved this many tolerances off a surface
BOX_PAD = 10.0  # tolerances of slack on bounding box tests
PROBES = 8  # probe points kept per cutter
SAMPLES = 60  # sample points of a solid for the "swallowed" test
GROW = 5.0  # tolerances a cutter is grown by in the over-cut rescue
SCALE_UP = 10.0  # factor of the scaled-up rescue
LOOSE = 10.0  # tolerance factor of the last rescue


class Ops(object):
    """The geometry calls cut() needs. Subclass it for your own solids.

    problem(shape)          None for a valid closed solid, else a short reason ("open")
    volume(shape)           its volume, 0 when it cannot be computed
    box(shape)              (min corner, max corner) of its bounding box
    difference(solids, cutters, tolerance)
                            list of solids left when every cutter is taken
                            away from every solid, or None when it failed
    touch(a, b)             True when their surfaces meet (intersection curves)
    inside(shape, point)    True when point lies inside the solid, off its surface
    samples(shape)          list of (point, outward unit normal) spread over its surface
    scaled(shape, centre, factors)
                            a copy scaled by factors (x, y, z) about centre, or None
    hollowed(shape, cutter) a copy of shape with cutter, which lies wholly
                            inside it, as a sealed hollow; or None
    """

    def problem(self, shape):
        raise NotImplementedError

    def volume(self, shape):
        raise NotImplementedError

    def box(self, shape):
        raise NotImplementedError

    def difference(self, solids, cutters, tolerance):
        raise NotImplementedError

    def touch(self, a, b):
        raise NotImplementedError

    def inside(self, shape, point):
        raise NotImplementedError

    def samples(self, shape):
        raise NotImplementedError

    def scaled(self, shape, centre, factors):
        raise NotImplementedError

    def hollowed(self, shape, cutter):
        raise NotImplementedError


class DifferenceResult(object):
    """The outcome of cut().

    solids    list of what is left of the solid (several when it was cut apart)
    applied   cutter indices that are in the result
    failed    {cutter index: [reasons]} for cutters left out
    missed    cutter indices that do not reach the solid (nothing to cut)
    rescued   {cutter index: rescue name} for cutters that needed a rescue
    report    the Report
    """

    def __init__(self, report):
        self.report = report
        self.solids = []
        self.applied = []
        self.failed = {}
        self.missed = []
        self.rescued = {}

    @property
    def status(self):
        return self.report.status


class _Cutter(object):
    def __init__(self, index, shape, volume, box):
        self.i = index
        self.shape = shape
        self.volume = volume
        self.box = box
        self.probes = []

    def variant(self, shape, volume, box):
        """The same cutter with other geometry (grown): keeps its index and probes."""
        other = _Cutter(self.i, shape, volume, box)
        other.probes = self.probes
        return other


def cut(solid, cutters, ops, tolerance=0.001):
    """Take cutters (a list of solids) away from solid through ops. Returns a DifferenceResult."""
    report = Report("verified difference")
    res = DifferenceResult(report)
    if not tolerance > 0.0:
        report.fail("the tolerance must be more than zero")
        return res
    if solid is None:
        report.fail("no solid to cut")
        return res
    why = ops.problem(solid)
    volume = ops.volume(solid) if not why else 0.0
    if not why and not volume > 0.0:
        why = "without a volume that can be computed"
    if why:
        report.fail("nothing cut: the solid is %s" % why)
        return res

    plan = _Plan(ops, tolerance, res)
    pad = BOX_PAD * tolerance
    solid_box = ops.box(solid)
    active, empty, unusable = [], [], []
    for i, shape in enumerate(cutters or []):
        if shape is None:
            empty.append(i)
            continue
        why = ops.problem(shape)
        v = ops.volume(shape) if not why else 0.0
        if not why and not v > 0.0:
            why = "without a volume that can be computed"
        box = ops.box(shape)
        if not _overlap(box, _inflate(solid_box, pad)):
            res.missed.append(i)
            continue
        if why:
            res.failed[i] = ["the cutter is %s" % why]
            unusable.append(i)
            continue
        active.append(_Cutter(i, shape, v, box))

    pieces = [solid]
    swallowed = None
    for c in active:
        c.probes = plan.probes(c, solid)
    reaching = []
    for c in active:
        if c.probes or ops.touch(solid, c.shape):
            reaching.append(c)
        elif plan.encloses(c, solid):
            swallowed = c.i
        else:
            res.missed.append(c.i)
    if swallowed is not None:
        pieces = []
        res.applied = [c.i for c in reaching] + [swallowed]
    elif reaching:
        pieces = plan.run(pieces, reaching)
        if plan.steps > 1:
            # Built in several steps: no later step may have brought back an earlier cut.
            for c in plan.applied:
                if plan.left_inside(c, pieces):
                    res.failed[c.i] = ["its cut was undone by a later step"]
        res.applied = [c.i for c in plan.applied if c.i not in res.failed]
    res.solids = pieces
    _summarize(res, len(cutters or []), empty, unusable, swallowed, volume, plan)
    return res


class _Plan(object):
    def __init__(self, ops, tolerance, res):
        self.ops = ops
        self.tol = tolerance
        self.depth = DEPTH * tolerance
        self.pad = BOX_PAD * tolerance
        self.res = res
        self.applied = []  # _Cutter objects in the result
        self.steps = 0  # results accepted
        self.calls = 0  # difference calls made
        self.notes = []  # (severity, text) from rescues

    # ---- sample points

    def probes(self, c, solid):
        """Points just inside the cutter's surface that are inside the solid.
        After the cut none of them may be inside the result."""
        found = []
        sbox = self.ops.box(solid)
        for p, n in self.ops.samples(c.shape):
            q = vec.sub(p, vec.scale(n, self.depth))
            if _holds(sbox, q) and self.ops.inside(c.shape, q) and self.ops.inside(solid, q):
                found.append(q)
        if len(found) > PROBES:
            step = len(found) / float(PROBES)
            found = [found[int(k * step)] for k in range(PROBES)]
        return found

    def inner_points(self, shape, limit=None):
        """Points just inside the surface of shape."""
        pts = [vec.sub(p, vec.scale(n, self.depth)) for p, n in self.ops.samples(shape)]
        if limit is not None and len(pts) > limit:
            step = len(pts) / float(limit)
            pts = [pts[int(k * step)] for k in range(limit)]
        return pts

    def encloses(self, c, solid):
        """True when the whole solid lies inside cutter c."""
        pts = self.inner_points(solid, SAMPLES)
        return bool(pts) and all(_holds(c.box, p) and self.ops.inside(c.shape, p) for p in pts)

    def left_inside(self, c, pieces):
        """True when a probe of cutter c is still inside one of the pieces."""
        for q in c.probes:
            for p in pieces:
                if _holds(self.ops.box(p), q) and self.ops.inside(p, q):
                    return True
        return False

    # ---- the checks

    def accept(self, out, before, cutters):
        """None when out is a sound result of cutting before with cutters, else the reason."""
        ops = self.ops
        if out is None:
            return "the difference failed"
        if any(p is None for p in out):
            return "the difference returned an empty piece"
        for p in out:
            why = ops.problem(p)
            if why:
                return "the difference returned a piece that is %s" % why
        if not out:
            if all(any(_holds(c.box, q) and ops.inside(c.shape, q) for c in cutters)
                   for q in self.inner_points_of(before)):
                return None  # the cutters swallow these pieces whole
            return "the difference returned nothing"
        v0 = sum(ops.volume(b) for b in before)
        v1 = sum(ops.volume(p) for p in out)
        vc = sum(c.volume for c in cutters)
        slack = VOLUME_SLACK * v0
        if v1 > v0 + slack:
            return "the volume grew (%.6g to %.6g)" % (v0, v1)
        if v1 < v0 - vc - slack:
            return "more volume went than the cutters hold (%.6g to %.6g)" % (v0, v1)
        if not _inside_box(_box_union([ops.box(p) for p in out]),
                           _inflate(_box_union([ops.box(b) for b in before]), self.pad)):
            return "the result reaches outside the solid's bounding box"
        for c in cutters:
            if self.left_inside(c, out):
                return "cutter %d was skipped: its material is still inside the result" % c.i
        reach = _inflate(_box_union([c.box for c in cutters]), self.pad)
        for p in out:
            for q in self.inner_points(p):
                if not _holds(reach, q):
                    continue
                for c in cutters:
                    if _holds(c.box, q) and ops.inside(c.shape, q):
                        return "part of the result lies inside cutter %d" % c.i
        return None

    def inner_points_of(self, shapes):
        pts = []
        for s in shapes:
            pts.extend(self.inner_points(s, SAMPLES))
        return pts

    # ---- attempts: each returns (new pieces, None) or (None, reason)

    def targets(self, pieces, cutters):
        """The pieces the cutters can reach, and the ones they cannot."""
        reach = _inflate(_box_union([c.box for c in cutters]), self.pad)
        hit = [p for p in pieces if _overlap(self.ops.box(p), reach)]
        idle = [p for p in pieces if not _overlap(self.ops.box(p), reach)]
        return hit, idle

    def plain(self, pieces, cutters, tolerance):
        hit, idle = self.targets(pieces, cutters)
        if not hit:
            return pieces, None
        self.calls += 1
        out = self.ops.difference(hit, [c.shape for c in cutters], tolerance)
        out = None if out is None else list(out)
        return self.finish(out, hit, idle, cutters)

    def finish(self, out, hit, idle, cutters):
        why = self.accept(out, hit, cutters)
        if why:
            return None, why
        self.steps += 1
        return out + idle, None

    def scaled_up(self, pieces, c):
        hit, idle = self.targets(pieces, [c])
        if not hit:
            return pieces, None
        centre = _centre(_box_union([self.ops.box(p) for p in hit]))
        up = (SCALE_UP, SCALE_UP, SCALE_UP)
        down = (1.0 / SCALE_UP, 1.0 / SCALE_UP, 1.0 / SCALE_UP)
        big = [self.ops.scaled(p, centre, up) for p in hit]
        big_cutter = self.ops.scaled(c.shape, centre, up)
        if big_cutter is None or any(b is None for b in big):
            return None, "could not scale the geometry"
        self.calls += 1
        out = self.ops.difference(big, [big_cutter], self.tol)
        if out is not None:
            out = [self.ops.scaled(p, centre, down) for p in out]
        return self.finish(out, hit, idle, [c])

    def hollow(self, pieces, c):
        """A cutter wholly inside a piece, its surface touching nothing, leaves a sealed hollow."""
        hit, idle = self.targets(pieces, [c])
        if any(self.ops.touch(p, c.shape) for p in hit):
            return None, "the cutter meets the surface"
        rim = [p for p, _n in self.ops.samples(c.shape)]
        holders = [p for p in hit if rim and all(self.ops.inside(p, q) for q in rim)]
        if not holders:
            return None, "the cutter is not wholly inside"
        out = self.ops.hollowed(holders[0], c.shape)
        rest = [p for p in hit if p is not holders[0]]
        return self.finish(None if out is None else [out] + rest, hit, idle, [c])

    def grown(self, pieces, c, uniform):
        """The cutter grown about its box centre: per axis, every side moves
        GROW tolerances; uniform, the thinnest direction moves that much."""
        lo, hi = c.box
        half = [(hi[k] - lo[k]) / 2.0 for k in range(3)]
        if min(half) <= 0.0:
            return None, "the cutter is flat"
        d = GROW * self.tol
        if uniform:
            s = 1.0 + d / min(half)
            factors = (s, s, s)
        else:
            factors = tuple(1.0 + d / h for h in half)
        shape = self.ops.scaled(c.shape, _centre(c.box), factors)
        if shape is None:
            return None, "could not grow the cutter"
        bigger = c.variant(shape, self.ops.volume(shape), self.ops.box(shape))
        res, why = self.plain(pieces, [bigger], self.tol)
        if res is None:
            return None, why
        # How far the farthest corner moved: the most the cut can overshoot.
        return res, vec.length([(factors[k] - 1.0) * half[k] for k in range(3)])

    # ---- control

    def run(self, pieces, cutters):
        if len(cutters) > 1:
            first = None
            for factor in (1.0, 2.0, 0.1):
                res, why = self.plain(pieces, cutters, factor * self.tol)
                if res is not None:
                    self.applied.extend(cutters)
                    return res
                if first is None:
                    first = why
            return self.bisect(pieces, cutters, first)
        return self.bisect(pieces, cutters)

    def bisect(self, pieces, cutters, known_failure=None):
        why = known_failure
        if why is None:
            res, why = self.plain(pieces, cutters, self.tol)
            if res is not None:
                self.applied.extend(cutters)
                return res
        if len(cutters) == 1:
            return self.rescue(pieces, cutters[0], why)
        mid = len(cutters) // 2
        pieces = self.bisect(pieces, cutters[:mid])
        return self.bisect(pieces, cutters[mid:])

    def rescue(self, pieces, c, first_why):
        tol = self.tol
        reasons = ["at the tolerance: " + first_why]
        ladder = (
            ("as a sealed hollow", lambda: self.hollow(pieces, c)),
            ("twice the tolerance", lambda: self.plain(pieces, [c], 2.0 * tol)),
            ("a tenth of the tolerance", lambda: self.plain(pieces, [c], tol / 10.0)),
            ("scaled up %g times" % SCALE_UP, lambda: self.scaled_up(pieces, c)),
            ("grown per axis", lambda: self.grown(pieces, c, False)),
            ("grown uniformly", lambda: self.grown(pieces, c, True)),
            ("%g times the tolerance" % LOOSE, lambda: self.plain(pieces, [c], LOOSE * tol)),
        )
        for name, attempt in ladder:
            res, extra = attempt()
            if res is not None:
                self.applied.append(c)
                self.res.rescued[c.i] = name
                if name == "as a sealed hollow":
                    self.notes.append((1, "cutter %d lies wholly inside the solid: it leaves a sealed hollow" % c.i))
                elif name.startswith("grown"):
                    self.notes.append((1, "cutter %d cut only when %s, over-cuts up to %.3g"
                                       % (c.i, name, extra)))
                elif name.startswith("%g" % LOOSE):
                    self.notes.append((1, "cutter %d cut only at %s (%g)" % (c.i, name, LOOSE * tol)))
                else:
                    self.notes.append((0, "cutter %d cut at %s" % (c.i, name)))
                return res
            reasons.append("%s: %s" % (name, extra))
        self.res.failed[c.i] = reasons
        return pieces


def _summarize(res, total, empty, unusable, swallowed, volume, plan):
    report = res.report
    ops = plan.ops
    report.note("%d of %d cutters cut the solid" % (len(res.applied), total))
    if swallowed is not None:
        report.warn("cutter %d swallows the whole solid: nothing is left" % swallowed)
    elif res.solids:
        left = sum(ops.volume(p) for p in res.solids)
        report.note("left: %d solid%s, volume %.6g of %.6g" % (
            len(res.solids), "" if len(res.solids) == 1 else "s", left, volume))
    else:
        report.warn("the cutters swallow the whole solid: nothing is left")
    if plan.calls:
        report.note("%d Boolean difference call%s" % (plan.calls, "" if plan.calls == 1 else "s"))
    if res.missed:
        report.note("do not reach the solid: cutter%s %s" % (
            "" if len(res.missed) == 1 else "s", _list(sorted(res.missed))))
    if empty:
        report.warn("empty input%s skipped: %s" % ("" if len(empty) == 1 else "s", _list(empty)))
    for severity, text in plan.notes:
        report.add(severity, text)
    if res.failed:
        report.messy("%d cutter%s could not be applied, the result is missing %s: %s" % (
            len(res.failed), "" if len(res.failed) == 1 else "s",
            "its cut" if len(res.failed) == 1 else "their cuts", _list(sorted(res.failed))))
        for i in sorted(res.failed):
            for why in res.failed[i]:
                report.note("cutter %d, %s" % (i, why))


def _list(indices):
    return ", ".join(str(i) for i in indices)


def _inflate(box, d):
    lo, hi = box
    return (lo[0] - d, lo[1] - d, lo[2] - d), (hi[0] + d, hi[1] + d, hi[2] + d)


def _overlap(a, b):
    return all(a[0][k] <= b[1][k] and b[0][k] <= a[1][k] for k in range(3))


def _holds(box, p):
    return all(box[0][k] <= p[k] <= box[1][k] for k in range(3))


def _inside_box(inner, outer):
    return all(outer[0][k] <= inner[0][k] and inner[1][k] <= outer[1][k] for k in range(3))


def _box_union(boxes):
    return (tuple(min(b[0][k] for b in boxes) for k in range(3)),
            tuple(max(b[1][k] for b in boxes) for k in range(3)))


def _centre(box):
    return tuple((box[0][k] + box[1][k]) / 2.0 for k in range(3))
