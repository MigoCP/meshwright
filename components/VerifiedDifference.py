"""VerifiedDifference: Boolean difference with many cutters, every result checked.

Cuts solid A with every cutter in B. All cutters go in one Boolean call
first; when that fails, the cutters are split in halves until the ones that
break it are found, and each of those goes through a rescue ladder (other
tolerances, the job scaled up 10 times, the cutter grown a hair, a loose
tolerance). A result is kept only when every piece is a valid closed solid,
the volume adds up, and sampled points show that no cutter was skipped and
no part of the result lies inside a cutter.

A cutter that still fails is left out and named in Failed and Info, so a
missing cut never goes unnoticed. Cutters that do not reach A are skipped.

Several solids: A is item access, so Grasshopper cuts each A with the whole
B list of its branch; put each solid's cutters in a matching branch.

Inputs:
    A      Brep, item access. The closed solid to cut.
    B      list of Brep, list access. Closed solids to cut away.
    Tol    float, item access, optional. Boolean tolerance (default: the
           document's absolute tolerance).

Outputs:
    Solids  list of Brep. What is left of A (several when the cuts split it).
    Failed  list of int. Indices of cutters that could not be applied.
    Status  int. 0 good, 1 warning, 2 messy, 3 failed.
    OK      bool. True for status 0 or 1.
    Info    str. Every reason, one per line.
"""

from __future__ import absolute_import, division

from solidwright.core import difference_plan
from solidwright.gh.io import read_input, show_report
from solidwright.rhino import convert
from solidwright.rhino.difference import BrepOps

_ns = globals()
_solid = read_input(_ns, "A")
_cutters = list(read_input(_ns, "B", []))
_tol = float(read_input(_ns, "Tol", convert.document_tolerance()))

_result = difference_plan.cut(_solid, _cutters, BrepOps(_tol), _tol)
_report = _result.report

Solids = _result.solids
Failed = sorted(_result.failed)
Status = _report.status
OK = _report.ok
Info = str(_report)
show_report(_ns.get("ghenv"), _report)
