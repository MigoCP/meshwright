"""Deform: compress, curl or radially fit parts above a start height along an axis.

The axis is the Z axis of P; the start height is a fraction S of the parts'
combined height along it (0 at the lowest point, 1 at the top). Below the
start nothing moves.

Modes:
    0 compress  lengths above the start times K (K < 1 squeezes, K > 1 stretches)
    1 curl      everything above the start rolls Ang degrees in towards the axis
                by the top, about a hinge circle of radius R around the axis
                (default: the mean distance from the axis above the start).
                The turn eases in, so there is no crease at the hinge.
    2 fit       radial scaling about the axis so the From outline lands on the
                To outline: fully at and below the start, fading out over Fade above it.

Inputs:
    G     list of Brep or Mesh, list access. The parts.
    Mode  int, item access, default 0. 0 compress, 1 curl, 2 fit.
    P     Plane, item access, optional. Axis frame: its Z axis is the axis and
          its origin lies on it. Default: world XY.
    S     float, item access, default 0. Start height as a fraction of the height.
    K     float, item access, default 1. Compress factor.
    Ang   float, item access, default 0. Curl angle in degrees.
    R     float, item access, optional. Curl hinge radius.
    From  Curve, item access. Fit: the outline as it is (closed, around the axis).
    To    Curve, item access. Fit: the outline it should become.
    Fade  float, item access, optional. Fit: length above the start over which
          the fit fades out (default: half the height above the start).
    Tol   float, item access, optional. Morph tolerance (default: the document's).

Outputs:
    Out     list of Brep or Mesh. The deformed copies, in input order.
    Hinge   Circle. Curl: the hinge circle.
    Status  int. 0 good, 1 warning, 2 messy, 3 failed.
    OK      bool. True for status 0 or 1.
    Info    str. Every reason, one per line.
"""

from __future__ import absolute_import, division

import math

import Rhino.Geometry as rg

from meshwright.core import deform as core
from meshwright.core.report import Report
from meshwright.gh.io import read_input, show_report
from meshwright.rhino import convert, deform

_ns = globals()
_items = [g for g in read_input(_ns, "G", []) if g is not None]
_mode = int(read_input(_ns, "Mode", 0))
_plane = read_input(_ns, "P", rg.Plane.WorldXY)
_s = float(read_input(_ns, "S", 0.0))
_tol = float(read_input(_ns, "Tol", convert.document_tolerance()))

Out, Hinge = [], None
_report = Report("deform")
_span = deform.extent(_items, _plane) if _items else None
if not _items:
    _report.fail("no parts connected")
elif _mode not in (0, 1, 2):
    _report.fail("Mode must be 0 (compress), 1 (curl) or 2 (fit), got %d" % _mode)
elif not 0.0 <= _s <= 1.0:
    _report.fail("S must be between 0 and 1, got %g" % _s)
else:
    _low, _high = _span
    _start = _low + _s * (_high - _low)
    _map, _affine = None, None
    if _mode == 0:
        _k = float(read_input(_ns, "K", 1.0))
        _map, _made = core.compress(_start, _k)
        if _map is not None and _s == 0.0:
            # Starting at the very bottom it is an exact scale along the axis.
            _base = rg.Plane(_plane.PointAt(0.0, 0.0, _start), _plane.XAxis, _plane.YAxis)
            _affine = rg.Transform.Scale(_base, 1.0, 1.0, _k)
    elif _mode == 1:
        _r = read_input(_ns, "R")
        _map, _made = core.curl(_start, _high, math.radians(float(read_input(_ns, "Ang", 0.0))),
                                None if _r is None else float(_r),
                                deform.local_points(_items, _plane) if _r is None else None)
        if _map is not None:
            Hinge = rg.Circle(rg.Plane(_plane.PointAt(0.0, 0.0, _start), _plane.XAxis, _plane.YAxis), _map.radius)
    else:
        _from, _to = read_input(_ns, "From"), read_input(_ns, "To")
        _fade = read_input(_ns, "Fade")
        if _from is None or _to is None:
            _map, _made = None, Report("radial fit").fail("connect both From and To outlines")
        else:
            _map, _made = core.radial_fit(deform.outline(_from, _plane), deform.outline(_to, _plane), _start,
                                          0.5 * (_high - _start) if _fade is None else float(_fade))
    _report.merge(_made)
    if _map is not None:
        Out, _done = deform.morph(_items, _map, _plane, _tol, _affine)
        _report.merge(_done)

Status = _report.status
OK = _report.ok
Info = str(_report)
show_report(_ns.get("ghenv"), _report)
