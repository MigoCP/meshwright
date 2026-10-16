"""Grouping parts by size.

Each part comes in as one measure (a volume, an area or a length) and is
compared by its linear size, measure ** power (1/3 for volumes, 1/2 for
areas), so a part 5% longer counts as 5% bigger whatever was measured.
Three ways to group, in this order of priority:

1. references: one measure per size; each part goes to the nearest one on a
   log scale (so 10% too big and 10% too small are equally far).
2. count: cut the sorted sizes at the count - 1 largest relative gaps.
3. neither: every relative gap bigger than the tolerance starts a new group.

Groups are numbered smallest first, and with references every reference
gets a group even when no part lands in it, so group k always means the
k-th smallest size.
"""

from __future__ import absolute_import, division

import math

from .report import Report

DEFAULT_TOLERANCE = 0.05


class SizeGroups(object):
    """index: group per input (-1 when skipped); sizes: linear size per group
    (from the mean measure of its members, None when empty); count: number of groups."""

    def __init__(self, report):
        self.report = report
        self.index = []
        self.sizes = []
        self.count = 0

    @property
    def status(self):
        return self.report.status

    def members(self):
        """Input indices per group."""
        out = [[] for _ in range(self.count)]
        for i, g in enumerate(self.index):
            if g >= 0:
                out[g].append(i)
        return out


def group_by_size(values, tolerance=DEFAULT_TOLERANCE, count=None, references=None, power=1.0):
    """SizeGroups for a list of measures (None or <= 0 are skipped and named)."""
    report = Report("size groups")
    res = SizeGroups(report)
    res.index = [-1] * len(values)
    if not tolerance > 0.0:
        report.fail("the tolerance must be more than zero")
        return res
    if count is not None and count < 1:
        count = None
    good = [i for i, v in enumerate(values) if _usable(v)]
    kept = set(good)
    bad = [i for i in range(len(values)) if i not in kept]
    if bad:
        report.warn("skipped %d input%s with no size: %s" % (len(bad), _s(bad), _list(bad)))
    if not good:
        report.fail("no inputs with a size")
        return res
    vals = [float(values[i]) for i in good]

    if references:
        refs = [float(r) for r in references if _usable(r)]
        dropped = len(references) - len(refs)
        if dropped:
            report.warn("skipped %d reference%s with no size" % (dropped, "" if dropped == 1 else "s"))
        if not refs:
            report.fail("no references with a size")
            return res
        if count is not None and count != len(refs):
            report.warn("the count is ignored when references are given (%d references)" % len(refs))
        groups, n = _by_references(vals, refs)
        how = "nearest of %d references" % len(refs)
    else:
        groups, n = _by_gaps(vals, power, tolerance, count, report)
        how = "%d asked for" % count if count else "gaps over %g%%" % (100.0 * tolerance)

    for k, g in zip(good, groups):
        res.index[k] = g
    res.count = n
    for g in range(n):
        members = [v for v, gg in zip(vals, groups) if gg == g]
        res.sizes.append((sum(members) / len(members)) ** power if members else None)
    empty = [g for g in range(n) if res.sizes[g] is None]
    report.note("%d group%s (%s): %s" % (n, "" if n == 1 else "s", how, ", ".join(
        "-" if s is None else "%.4g" % s for s in res.sizes)))
    if empty:
        report.warn("empty group%s %s: no input is nearest to %s reference" % (
            _s(empty), _list(empty), "its" if len(empty) == 1 else "their"))
    return res


def _by_references(vals, refs):
    """Nearest reference on a log scale; references ranked by size."""
    rank = sorted(range(len(refs)), key=lambda j: refs[j])
    rank_of = dict((j, r) for r, j in enumerate(rank))
    logs = [math.log(r) for r in refs]
    out = []
    for v in vals:
        lv = math.log(v)
        j = min(range(len(refs)), key=lambda j: abs(lv - logs[j]))
        out.append(rank_of[j])
    return out, len(refs)


def _by_gaps(vals, power, tolerance, count, report):
    """Cut the sorted values at the largest relative gaps in linear size."""
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    gaps = []
    for k in range(1, len(order)):
        a, b = vals[order[k - 1]], vals[order[k]]
        gaps.append(((b / a) ** power - 1.0, k))
    real = [g for g in gaps if g[0] > tolerance]
    if count:
        if count - 1 > len(gaps):
            report.warn("asked for %d groups but there are only %d inputs" % (count, len(vals)))
        elif count - 1 > len(real):
            report.warn("only %d sizes differ by more than the tolerance; some of the %d groups split "
                        "near-identical inputs" % (len(real) + 1, count))
        elif count - 1 < len(real):
            report.warn("%d distinct sizes found but %d asked for; some sizes were merged" % (
                len(real) + 1, count))
        chosen = sorted(gaps, reverse=True)[:count - 1]
    else:
        chosen = real
    cuts = set(k for _, k in chosen)
    out = [0] * len(vals)
    g = 0
    for r, i in enumerate(order):
        if r in cuts:
            g += 1
        out[i] = g
    return out, len(cuts) + 1


def _usable(v):
    try:
        return v is not None and float(v) > 0.0 and not math.isinf(float(v))
    except (TypeError, ValueError):
        return False


def _s(items):
    return "" if len(items) == 1 else "s"


def _list(indices):
    return ", ".join(str(i) for i in indices)
