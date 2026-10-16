"""RobustSweep: a one-rail sweep that tries hard and says how it went.

Cleans the curves, finds where each section sits on the rail, sorts and
aligns the sections (direction and seam, compared in rotation-minimizing
frames), then tries solvers in order and keeps the first clean result:
Sweep1 (and the other miter types when the rail has kinks), Sweep1 with
sections refit, Sweep1 on a rebuilt rail, its own frames + loft (kinks
mitered, closed-rail twist spread evenly), Sweep1 at 10 x tolerance,
CreateFromSweepSegmented, CreateFromSweep. Info lists every attempt.

Several rails at once: graft Rail and Sec so each branch holds one rail and
its own sections.

Inputs:
    Rail   Curve, item access. The rail.
    Sec    list of Curve, list access. One or more sections, in any order.
    Miter  int, item access, default 1. 0 none, 1 trimmed, 2 untrimmed (Grasshopper's "rotated").
    Cap    bool, item access, default False. Cap planar ends when the rail is open.
    Tol    float, item access, optional. Sweep tolerance (default: the document's).

Outputs:
    Brep    list of Brep. The sweep (one Brep when it went well).
    Frames  list of Plane. Rail frames where each section attaches.
    Status  int. 0 good, 1 warning (a fallback method was needed), 2 messy, 3 failed.
    OK      bool. True for status 0 or 1.
    Info    str. Every reason, one per line.
"""

from __future__ import absolute_import, division

from solidwright.gh.io import read_input, show_report
from solidwright.rhino import convert, sweep

_ns = globals()
_rail = read_input(_ns, "Rail")
_secs = read_input(_ns, "Sec", [])
if not isinstance(_secs, (list, tuple)):
    _secs = [_secs]  # Sec left on item access
_secs = [s for s in _secs if s is not None]
_miter = int(read_input(_ns, "Miter", 1))
_cap = bool(read_input(_ns, "Cap", False))
_tol = float(read_input(_ns, "Tol", convert.document_tolerance()))

_result = sweep.robust_sweep(_rail, _secs, _miter, _cap, _tol, convert.document_angle_tolerance())
_report = _result.report

Brep = _result.breps
Frames = _result.frames
Status = _report.status
OK = _report.ok
Info = str(_report)
show_report(_ns.get("ghenv"), _report)
