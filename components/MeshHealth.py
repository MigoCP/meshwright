"""MeshHealth: check a mesh and clean what can be cleaned safely.

Welds coincident vertices, removes zero-area and duplicate faces, drops loose
junk, makes face directions consistent and points closed pieces outward. Then
grades what is left: holes, extra pieces, non-manifold edges.

Inputs:
    M      Mesh, item access. The mesh to check.
    Tol    float, item access, default 0. Weld vertices closer than this
           (0 welds only exact duplicates).
    Junk   float, item access, default 0.01. Loose pieces smaller than this
           share of the largest piece's area are dropped.

Outputs:
    Clean   Mesh. The cleaned mesh (empty when the check failed).
    Status  int. 0 good, 1 warning, 2 messy, 3 failed.
    OK      bool. True when the mesh can be used as is (status 0 or 1).
    Info    str. Every reason, one per line.
"""

from __future__ import absolute_import, division

from meshwright.core import topology
from meshwright.gh.io import read_input, show_report
from meshwright.rhino import convert

_ns = globals()
_mesh = read_input(_ns, "M")
_tol = float(read_input(_ns, "Tol", 0.0))
_junk = float(read_input(_ns, "Junk", topology.JUNK_RATIO))

if _mesh is None:
    Clean, Status, OK, Info = None, 3, False, "no mesh connected"
else:
    _vertices, _faces = convert.mesh_to_lists(_mesh)
    _health = topology.health(_vertices, _faces, tolerance=_tol, junk_ratio=_junk)
    Clean = convert.lists_to_mesh(_health.vertices, _health.faces) if _health.faces else None
    Status = _health.status
    OK = _health.report.ok
    Info = str(_health.report)
    show_report(_ns.get("ghenv"), _health.report)
