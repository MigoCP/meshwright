"""VerifiedUnion: Boolean union of many solids, one at a time, every step checked.

For the many-solids case where one Boolean union over everything fails (one
pair with flush faces is enough). Starts from the largest solid and unites
the others into it one by one, largest first. A solid that refuses is
retried once everything else is in, and only then nudged a hair (Z+, Z-,
X+, X-, Y+, Y-) so flush faces stop being coincident. A union is kept only
when it gives one valid closed solid that covers both inputs and whose
volume adds up, so a piece the Boolean silently dropped is caught.

Solids that touch nothing come out as they went in. Open or invalid inputs
stop everything before any union, and are named in Info.

Inputs:
    B      list of Brep, list access. Closed solids.
    Tol    float, item access, optional. Boolean tolerance (default: the
           document's absolute tolerance).
    Nudge  float, item access, optional. How far a refusing solid is moved
           (default 10 x Tol).

Outputs:
    Solids  list of Brep. The united solids, plus the ones that touch nothing.
    Failed  list of int. Indices of inputs that could not be united, left out of Solids.
    Status  int. 0 good, 1 warning, 2 messy, 3 failed.
    OK      bool. True for status 0 or 1.
    Info    str. Every reason, one per line.
"""

from __future__ import absolute_import, division

from solidwright.core import union_plan
from solidwright.gh.io import read_input, show_report
from solidwright.rhino import convert
from solidwright.rhino.union import BrepOps

_ns = globals()
_breps = list(read_input(_ns, "B", []))
_tol = float(read_input(_ns, "Tol", convert.document_tolerance()))
_nudge = read_input(_ns, "Nudge")

_result = union_plan.union_all(_breps, BrepOps(_tol), _tol, None if _nudge is None else float(_nudge))
_report = _result.report

Solids = _result.solids
Failed = sorted(_result.failed)
Status = _report.status
OK = _report.ok
Info = str(_report)
show_report(_ns.get("ghenv"), _report)
