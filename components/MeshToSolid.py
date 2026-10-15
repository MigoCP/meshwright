"""MeshToSolid: a closed Brep with one exact flat face per planar region of a mesh.

Checks and cleans the mesh first (as MeshHealth does), groups its triangles
into flat regions, solves every corner from the planes that meet there and
builds the solid. Works for parts that are not convex and for faces with
holes. Corners closer than Tol are merged when the solid stays closed and
flat.

Meant for meshes of parts made of flat faces (STL exports, scans of machined
parts). A curved mesh still gives a valid solid, with one face per flat
patch, and Status says it looks curved.

Inputs:
    M     Mesh, item access. The part (closed).
    Tol   float, item access, optional. Modelling tolerance (default: the
          document's absolute tolerance). Faces are flat within half of it;
          a mesh noisier than that keeps faces split, and Status says so.

Outputs:
    Solid   Brep. The closed solid (empty when it failed).
    Status  int. 0 good, 1 warning, 2 messy, 3 failed (worst of mesh, rebuild and Brep).
    OK      bool. True for status 0 or 1.
    Info    str. Every reason, one per line.
"""

from __future__ import absolute_import, division

from meshwright.core import rebuild, topology
from meshwright.core.report import Report
from meshwright.gh.io import read_input, show_report
from meshwright.rhino import build, convert

_ns = globals()
_mesh = read_input(_ns, "M")
_tol = float(read_input(_ns, "Tol", convert.document_tolerance()))

Solid = None
_report = Report("mesh to solid")
if _mesh is None:
    _report.fail("no mesh connected")
else:
    _vertices, _faces = convert.mesh_to_lists(_mesh)
    _health = topology.health(_vertices, _faces)
    _report.merge(_health.report, prefix="mesh: ")
    if not _health.report.failed:
        _result = rebuild.rebuild(_health.vertices, _health.faces, tolerance=_tol)
        _report.merge(_result.report, prefix="rebuild: ")
        if _result.solid is not None:
            Solid, _built = build.polyhedron_to_brep(_result.solid, _tol)
            _report.merge(_built, prefix="brep: ")

Status = _report.status
OK = _report.ok
Info = str(_report)
show_report(_ns.get("ghenv"), _report)
