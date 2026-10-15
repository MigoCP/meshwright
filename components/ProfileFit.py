"""ProfileFit: a part's side outline and where it bends, as a three-piece broken line.

Cuts the part with a few planes that contain the axis, keeps the outermost
point at every height (so noisy or coarse meshes still give one clean
outline), and fits three straight pieces to it by least squares, trying
every pair of break points. The two knees, where neighbouring pieces meet,
are where the outline bends: for a part that tapers in, runs straight and
tapers out again, the ends of its straight band.

Inputs:
    M   Mesh, item access. The part.
    A   Vector3d or Plane, item access, optional. The axis the outline runs
        along (a plane means its Z axis). Leave empty to detect it.
    D   Vector3d, item access, optional. The side to look at (default: the
        X axis of the part's frame).

Outputs:
    Profile  Polyline. The side outline.
    Lines    list of Line. The three fitted pieces.
    Knees    list of Point3d. The two bend points.
    Status   int. 0 good, 1 warning, 2 messy, 3 failed.
    OK       bool. True for status 0 or 1.
    Info     str. Every reason, one per line.
"""

from __future__ import absolute_import, division

import Rhino.Geometry as rg

from meshwright.core import orient, sections, topology, vec
from meshwright.core.report import Report
from meshwright.gh.io import read_input, show_report
from meshwright.rhino import convert

_ns = globals()
_mesh = read_input(_ns, "M")
_axis = convert.from_vector(read_input(_ns, "A"))
_side = convert.from_vector(read_input(_ns, "D"))

Profile, Lines, Knees = None, [], []
_report = Report("profile fit")
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
            if _side is None:
                _side = vec.frame_from_axis(_axis)[0]
            _result = sections.side_profile(_health.vertices, _health.faces, _axis, _side)
            _report.merge(_result.report, prefix="profile: ")
            if _result.points:
                Profile = convert.to_polyline([_result.to_world(x, z) for z, x in _result.points])
            if _result.fit is not None:
                _zs = [p[0] for p in _result.points]
                _spans = ((_zs[0], _result.fit.knees[0][0]), (_result.fit.knees[0][0], _result.fit.knees[1][0]),
                          (_result.fit.knees[1][0], _zs[-1]))
                for _line, (_z0, _z1) in zip(_result.fit.lines, _spans):
                    Lines.append(rg.Line(convert.to_point(_result.to_world(_line.at(_z0), _z0)),
                                         convert.to_point(_result.to_world(_line.at(_z1), _z1))))
                Knees = [convert.to_point(_result.to_world(x, z)) for z, x in _result.fit.knees]

Status = _report.status
OK = _report.ok
Info = str(_report)
show_report(_ns.get("ghenv"), _report)
