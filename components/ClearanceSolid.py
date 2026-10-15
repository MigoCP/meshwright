"""ClearanceSolid: a closed Brep around any part, its sections kept a clearance away sideways.

Checks and cleans the mesh first (as MeshHealth does), then describes the part
by cross sections square to its axis: the widest band, a section wherever
the profile bends, both ends (flat, or tapering to a point or a ridge), and
extra sections wherever a vertex came out closer than the clearance. Each
section is the part's convex outline at that height, pushed out by Clr.
The solid joins the sections with flat faces, so it is exactly what was
checked.

Works for curved and scanned parts too (OffsetSolid needs flat facets). Each
section is convex: hollows within a section are bridged, but the profile may
go in and out along the axis.

Inputs:
    M     Mesh, item access. The part (closed).
    Clr   float, item access, default 0. Sideways clearance.
    A     Vector3d or Plane, item access, optional. The axis (a plane means
          its Z axis). Leave empty to detect the part's natural axis.
    Tol   float, item access, optional. Modelling tolerance (default: the
          document's absolute tolerance).
    Dirs  int, item access, default 96. Support directions for sections that
          are not polygons (more: closer to round, more faces).

Outputs:
    Solid     Brep. The closed solid (empty when it failed).
    Sections  list of Polyline. The section at every level, closed.
    Status    int. 0 good, 1 warning, 2 messy, 3 failed (worst of every step).
    OK        bool. True for status 0 or 1.
    Info      str. Every reason, one per line.
"""

from __future__ import absolute_import, division

from meshwright.core import orient, sections, topology
from meshwright.core.report import Report
from meshwright.gh.io import read_input, show_report
from meshwright.rhino import build, convert

_ns = globals()
_mesh = read_input(_ns, "M")
_clr = float(read_input(_ns, "Clr", 0.0))
_axis = convert.from_vector(read_input(_ns, "A"))
_tol = float(read_input(_ns, "Tol", convert.document_tolerance()))
_dirs = int(read_input(_ns, "Dirs", sections.DIRECTIONS))

Solid, Sections = None, []
_report = Report("clearance solid")
if _mesh is None:
    _report.fail("no mesh connected")
else:
    _vertices, _faces = convert.mesh_to_lists(_mesh)
    _health = topology.health(_vertices, _faces)
    _report.merge(_health.report, prefix="mesh: ")
    if not _health.report.failed:
        if _axis is None:
            _found = orient.natural_axis(_health.vertices, _health.faces)
            _report.merge(_found.report, prefix="axis: ")
            _axis = _found.axis
        if _axis is not None:
            _result = sections.clearance_solid(_health.vertices, _health.faces, _axis, clearance=_clr,
                                               tolerance=_tol, directions=_dirs)
            _report.merge(_result.report, prefix="solid: ")
            if _result.model is not None and _result.model.levels:
                for _z, _hs in _result.model.levels:
                    _pts = [_result.to_world((x, y, _z)) for x, y in _result.model.corners(_hs)]
                    Sections.append(convert.to_polyline(_pts, closed=True))
            if _result.solid is not None:
                Solid, _built = build.polyhedron_to_brep(_result.solid, _tol)
                _report.merge(_built, prefix="brep: ")

Status = _report.status
OK = _report.ok
Info = str(_report)
show_report(_ns.get("ghenv"), _report)
