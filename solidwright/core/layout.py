"""Stacking parts along an axis with a fixed gap.

Works on extents only: each part is its (low, high) interval along the
axis, so the geometry side just measures and moves.
"""

from __future__ import absolute_import, division

import math

from .report import Report


class Stack(object):
    """shifts: distance to move each input along the axis (None when skipped);
    order: input indices in stacking order, first one staying put."""

    def __init__(self, report):
        self.report = report
        self.shifts = []
        self.order = []

    @property
    def status(self):
        return self.report.status


def stack(extents, gap, keep_order=False):
    """Stack intervals end to end with gap between neighbours.

    extents is a list of (low, high) along the axis, or None for a missing
    part. The parts go lowest first (by their low end) unless keep_order,
    then in input order. The first part stays where it is.
    """
    report = Report("stack")
    res = Stack(report)
    res.shifts = [None] * len(extents)
    try:
        gap = float(gap)
    except (TypeError, ValueError):
        gap = float("nan")
    if math.isnan(gap) or math.isinf(gap):
        report.fail("the gap must be a number")
        return res
    present = [i for i, e in enumerate(extents) if e is not None]
    missing = [i for i, e in enumerate(extents) if e is None]
    if missing:
        report.warn("skipped %d empty input%s: %s" % (
            len(missing), "" if len(missing) == 1 else "s", ", ".join(str(i) for i in missing)))
    if not present:
        report.fail("nothing to stack")
        return res
    if gap < 0.0:
        report.warn("the gap is negative: neighbours overlap by %g" % -gap)
    order = present if keep_order else sorted(present, key=lambda i: extents[i][0])
    top = None
    for i in order:
        low, high = extents[i]
        if top is None:
            res.shifts[i] = 0.0
            top = high
        else:
            res.shifts[i] = top + gap - low
            top = top + gap + (high - low)
    res.order = order
    report.note("%d part%s stacked, %g long" % (
        len(order), "" if len(order) == 1 else "s", top - extents[order[0]][0]))
    return res
