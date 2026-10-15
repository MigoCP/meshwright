"""NaturalAxis: the direction a part obviously sits along, found from its mesh.

Checks and cleans the mesh first (as MeshHealth does), then finds the axis
from the part's flat facets and outline. The axis points to the part's top:
toward its stronger flat end, or away from its tapering end. Give an axis in
A to have it checked against the mesh instead of detected.

Inputs:
    M   Mesh, item access. The part.
    A   Vector3d or Plane, item access, optional. An axis to check (a plane
        means its Z axis). Leave empty to detect the axis.

Outputs:
    Axis    Vector3d. Unit axis pointing to the top.
    Frame   Plane. Origin at the part's centre, Z along Axis, X as close to
            world X as possible.
    Status  int. 0 good, 1 warning, 2 messy, 3 failed (worst of mesh and axis).
    OK      bool. True for status 0 or 1.
    Info    str. Every reason, one per line.
"""

from __future__ import absolute_import, division

from meshwright.core import orient, topology, vec
from meshwright.core.report import Report
from meshwright.gh.io import read_input, show_report
from meshwright.rhino import convert

_ns = globals()
_mesh = read_input(_ns, "M")
_given = convert.from_vector(read_input(_ns, "A"))

Axis, Frame = None, None
_report = Report("natural axis")
if _mesh is None:
    _report.fail("no mesh connected")
else:
    _vertices, _faces = convert.mesh_to_lists(_mesh)
    _health = topology.health(_vertices, _faces)
    _report.merge(_health.report, prefix="mesh: ")
    if not _health.report.failed:
        if _given is None:
            _result = orient.natural_axis(_health.vertices, _health.faces)
        else:
            _result = orient.check_axis(_health.vertices, _health.faces, _given)
        _report.merge(_result.report, prefix="axis: ")
        if _result.axis is not None:
            Axis = convert.to_vector(_result.axis)
            Frame = convert.to_plane(vec.mean(_health.vertices), _result.frame)

Status = _report.status
OK = _report.ok
Info = str(_report)
show_report(_ns.get("ghenv"), _report)
