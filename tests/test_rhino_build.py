import Rhino.Geometry as rg

from solidwright.core import polyhedron as ph
from solidwright.core import rebuild
from solidwright.rhino import build

import shapes


def offset(mesh, clearance=0.05):
    v, f = mesh.as_lists()
    return ph.offset_solid(v, f, clearance=clearance, axis=(0.0, 0.0, 1.0)).solid


def test_planar_faces_join_into_one_solid():
    poly = offset(shapes.faceted_part())
    brep, report = build.polyhedron_to_brep(poly, 0.001)
    assert report.status == 0
    assert brep.IsSolid and brep.SolidOrientation == rg.BrepSolidOrientation.Outward
    assert brep.Faces.Count == len(poly.faces)


def test_falls_back_to_a_welded_mesh(monkeypatch):
    monkeypatch.setattr(rg.Brep, "fail_planar", True)
    poly = offset(shapes.faceted_part())
    brep, report = build.polyhedron_to_brep(poly, 0.001)
    assert report.status == 1 and "through a welded mesh" in str(report)
    assert brep.IsSolid
    assert brep.Faces.Count == len(poly.faces)  # coplanar triangles merged back


def test_inside_out_faces_are_flipped():
    poly = offset(shapes.box())
    inside_out = ph.Polyhedron(poly.vertices, [list(reversed(f)) for f in poly.faces])
    brep, report = build.polyhedron_to_brep(inside_out, 0.001)
    assert brep.SolidOrientation == rg.BrepSolidOrientation.Outward


def test_open_faces_fail(monkeypatch):
    monkeypatch.setattr(rg.Brep, "fail_planar", True)
    poly = offset(shapes.box())
    opened = ph.Polyhedron(poly.vertices, poly.faces[:-1])
    brep, report = build.polyhedron_to_brep(opened, 0.001)
    assert brep is None and report.failed
    assert "could not make a closed solid from the 5 faces" in str(report)


def rebuilt(mesh):
    v, f = mesh.as_lists()
    return rebuild.rebuild(v, f, 0.001).solid


def test_faces_with_holes_and_non_convex_faces(monkeypatch):
    for mesh in (shapes.frame(), shapes.plate_with_round_hole(), shapes.u_block()):
        poly = rebuilt(mesh)
        brep, report = build.polyhedron_to_brep(poly, 0.001)
        assert report.status == 0 and brep.IsSolid
        assert brep.Faces.Count == len(poly.faces)
    # Through the welded mesh, the holed and L-shaped faces are triangulated and merged back.
    monkeypatch.setattr(rg.Brep, "fail_planar", True)
    for mesh in (shapes.frame(), shapes.u_block()):
        poly = rebuilt(mesh)
        brep, report = build.polyhedron_to_brep(poly, 0.001)
        assert report.status == 1 and brep.IsSolid
        assert brep.Faces.Count == len(poly.faces)

