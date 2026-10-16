"""StackAlongAxis: line parts up end to end along an axis with a fixed gap.

The lowest part (along the axis) stays put and the others follow it, lowest
first, each moved along the axis only, so it sits Gap past the previous
one. Extents are measured in the axis direction, so any axis works, not
only X, Y and Z.

Inputs:
    G     list of Geometry, list access. Breps, meshes, curves...
    Gap   float, item access, default 0. Space between neighbours (negative: overlap).
    A     Vector3d or Plane, item access, optional. The axis (a plane means its
          Z axis). Default: world Z.
    Keep  bool, item access, default False. Stack in input order instead of lowest first.

Outputs:
    Moved   list of Geometry. The moved copies, in input order.
    X       list of Transform. The move for each input, to apply to related geometry.
    Order   list of int. Input indices in stacking order.
    Status  int. 0 good, 1 warning, 2 messy, 3 failed.
    OK      bool. True for status 0 or 1.
    Info    str. Every reason, one per line.
"""

from __future__ import absolute_import, division

import Rhino.Geometry as rg

from solidwright.core import layout, vec
from solidwright.core.report import Report
from solidwright.gh.io import read_input, show_report
from solidwright.rhino import convert, measure

_ns = globals()
_items = list(read_input(_ns, "G", []))
_gap = read_input(_ns, "Gap", 0.0)
_axis = convert.from_vector(read_input(_ns, "A")) or (0.0, 0.0, 1.0)
_keep = bool(read_input(_ns, "Keep", False))

Moved, X, Order = [], [], []
_report = Report("stack along axis")
_unit = vec.unit(_axis)
if _unit is None:
    _report.fail("the axis is a zero vector")
else:
    _result = layout.stack([measure.extent(g, _unit) for g in _items], _gap, _keep)
    _report.merge(_result.report)
    for _g, _d in zip(_items, _result.shifts):
        if _d is None:
            Moved.append(None)
            X.append(None)
            continue
        _x = rg.Transform.Translation(rg.Vector3d(_unit[0] * _d, _unit[1] * _d, _unit[2] * _d))
        _copy = _g.Duplicate()
        _copy.Transform(_x)
        Moved.append(_copy)
        X.append(_x)
    Order = _result.order

Status = _report.status
OK = _report.ok
Info = str(_report)
show_report(_ns.get("ghenv"), _report)
