"""DomeCap: rounded caps on flat outlines, or a part cut off and capped.

A dome rises from a closed planar outline to an apex above its centroid
through shrinking copies of the outline (a quarter ellipse in profile), so
it meets the outline's plane square on. Its height is H times the radius of
a circle with the outline's area: H = 1 on a circle gives a half ball.

Two ways to use it:
    C connected: a dome on each outline (open surfaces, their edge is the outline).
    B and P connected: each part is cut by P, the side above P (along its Z)
        is dropped, every cut section gets a dome, and the lot is joined into
        a closed Brep. A hollow cut (a section with a hole) is refused.

Inputs:
    C     list of Curve, list access, optional. Closed planar outlines.
    D     Vector3d, item access, optional. Side the domes on C rise to
          (default: towards world +Z).
    B     list of Brep, list access, optional. Parts to cut and cap.
    P     Plane, item access. The cut (needed with B); domes rise along its Z.
    H     float, item access, default 1. Dome height ratio (more than 0).
    N     int, item access, default 6. Sections between the outline and the apex.
    Tol   float, item access, optional. Tolerance (default: the document's).

Outputs:
    Brep      list of Brep. The capped parts (with B), or the domes (with C only).
    Domes     list of Brep. The domes alone.
    Sections  list of Curve. The outlines that were capped.
    Status    int. 0 good, 1 warning, 2 messy, 3 failed (worst of all inputs).
    OK        bool. True for status 0 or 1.
    Info      str. Every reason, one per line.
"""

from __future__ import absolute_import, division

from meshwright.core.report import Report
from meshwright.gh.io import read_input, show_report
from meshwright.rhino import convert, dome

_ns = globals()
_curves = [c for c in read_input(_ns, "C", []) if c is not None]
_up = read_input(_ns, "D")
_parts = [b for b in read_input(_ns, "B", []) if b is not None]
_plane = read_input(_ns, "P")
_ratio = float(read_input(_ns, "H", 1.0))
_steps = int(read_input(_ns, "N", 6))
_tol = float(read_input(_ns, "Tol", convert.document_tolerance()))

Brep, Domes, Sections = [], [], []
_report = Report("dome cap")
if not _ratio > 0.0:
    _report.fail("H must be more than zero, got %g" % _ratio)
elif _steps < 1:
    _report.fail("N must be at least 1, got %d" % _steps)
elif not _curves and not _parts:
    _report.fail("connect outlines (C), or parts (B) and a cutting plane (P)")
elif _parts and _plane is None:
    _report.fail("connect a cutting plane (P) for the parts")
else:
    for _k, _b in enumerate(_parts):
        _breps, _domes, _outlines, _done = dome.cap_part(_b, _plane, _ratio, _steps, _tol)
        _report.merge(_done, prefix="part %d: " % _k if len(_parts) > 1 else "")
        Brep.extend(_breps)
        Domes.extend(_domes)
        Sections.extend(_outlines)
    for _k, _c in enumerate(_curves):
        _made, _height, _why = dome.dome(_c, _ratio, _steps, _tol, _up)
        if _made is None:
            _report.messy("outline %d: %s" % (_k, _why))
            continue
        _report.note("outline %d: dome %.4g high" % (_k, _height))
        Domes.extend(_made)
        Sections.append(_c)
        if not _parts:
            Brep.extend(_made)

Status = _report.status
OK = _report.ok
Info = str(_report)
show_report(_ns.get("ghenv"), _report)
