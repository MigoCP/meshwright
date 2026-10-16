"""OffsetSolid: a closed Brep around a faceted convex part, every facet pushed out by a clearance.

Checks and cleans the mesh first (as MeshHealth does), groups its triangles
into flat facets, moves each facet plane out until it touches the part and
then out by Clr, and intersects all of them. Corners that the offset splits
into tiny clusters are merged only when the merged solid is still closed,
flat and clear of the part.

With an axis in A, the clearance is applied square to it only: faces square
to the axis stay where they are, faces along it move the full Clr. Feed it
from NaturalAxis to follow the part's own axis.

The solid is convex. A part that is not convex gets a convex solid that
bridges its hollows, and Status says so.

Inputs:
    M     Mesh, item access. The part.
    Clr   float, item access, default 0. Clearance (0: the solid touches the part).
    A     Vector3d or Plane, item access, optional. Apply the clearance square
          to this axis only (a plane means its Z axis). Leave empty for the
          same clearance in every direction.
    Tol   float, item access, optional. Modelling tolerance (default: the
          document's absolute tolerance).

Outputs:
    Solid   Brep. The closed solid (empty when it failed).
    Status  int. 0 good, 1 warning, 2 messy, 3 failed (worst of mesh, solid and Brep).
    OK      bool. True for status 0 or 1.
    Info    str. Every reason, one per line.
"""

from __future__ import absolute_import, division

from solidwright.core import polyhedron, topology
from solidwright.core.report import Report
from solidwright.gh.io import read_input, show_report
from solidwright.rhino import build, convert

_ns = globals()
_mesh = read_input(_ns, "M")
_clr = float(read_input(_ns, "Clr", 0.0))
_axis = convert.from_vector(read_input(_ns, "A"))
_tol = float(read_input(_ns, "Tol", convert.document_tolerance()))

Solid = None
_report = Report("offset solid")
if _mesh is None:
    _report.fail("no mesh connected")
else:
    _vertices, _faces = convert.mesh_to_lists(_mesh)
    _health = topology.health(_vertices, _faces)
    _report.merge(_health.report, prefix="mesh: ")
    if not _health.report.failed:
        _result = polyhedron.offset_solid(_health.vertices, _health.faces, clearance=_clr,
                                          axis=_axis, tolerance=_tol)
        _report.merge(_result.report, prefix="solid: ")
        if _result.solid is not None:
            Solid, _built = build.polyhedron_to_brep(_result.solid, _tol)
            _report.merge(_built, prefix="brep: ")

Status = _report.status
OK = _report.ok
Info = str(_report)
show_report(_ns.get("ghenv"), _report)
