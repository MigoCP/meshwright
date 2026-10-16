"""SizeGroups: sort parts into size groups, smallest first.

Each part is measured by its volume (any orientation works), or by its area
when some parts are open, and compared by linear size, so 5% means 5%
longer. Grouping, in this order of priority:

1. Refs connected: each part goes to the nearest reference part's size.
   Every reference gets a branch, even an empty one.
2. N given: the sizes are cut at the N - 1 largest gaps.
3. Neither: every gap bigger than Tol starts a new group.

Inputs:
    G     list of Mesh or Brep, list access. The parts, in any order.
    Refs  list of Mesh or Brep, list access, optional. One part per size.
    N     int, item access, optional. Number of groups when no Refs are given.
    Tol   float, item access, default 0.05. Relative size gap that starts a new group (5%).

Outputs:
    Groups  tree. Branch {k} holds the parts of group k, smallest group first.
    Index   list of int. Group of each input (-1: skipped), in input order.
    Sizes   list of float. Linear size of each group (cube root of the mean
            volume, or square root of the mean area).
    Status  int. 0 good, 1 warning, 2 messy, 3 failed.
    OK      bool. True for status 0 or 1.
    Info    str. Every reason, one per line.
"""

from __future__ import absolute_import, division

from meshwright.core import cluster1d
from meshwright.core.report import Report
from meshwright.gh.io import read_input, show_report, tree
from meshwright.rhino import measure

_ns = globals()
_parts = list(read_input(_ns, "G", []))
_refs = [r for r in read_input(_ns, "Refs", []) if r is not None]
_count = read_input(_ns, "N")
_tol = float(read_input(_ns, "Tol", cluster1d.DEFAULT_TOLERANCE))

Groups, Index, Sizes = tree([]), [], []
_report = Report("size groups")
if not _parts:
    _report.fail("no parts connected")
else:
    _measures, _power, _measured = measure.sizes(_parts + _refs)
    _report.merge(_measured)
    _result = cluster1d.group_by_size(_measures[:len(_parts)], _tol, None if _count is None else int(_count),
                                      _measures[len(_parts):] or None, _power)
    _report.merge(_result.report)
    Groups = tree([[_parts[i] for i in m] for m in _result.members()])
    Index = _result.index
    Sizes = _result.sizes

Status = _report.status
OK = _report.ok
Info = str(_report)
show_report(_ns.get("ghenv"), _report)
