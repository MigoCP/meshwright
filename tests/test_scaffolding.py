"""The test tools themselves: procedural shapes and the Rhino mock.

If these are wrong, every other test is meaningless, so they get checked too.
"""

import math

import numpy as np
import pytest

import shapes

import Rhino.Geometry as rg
from Rhino.Geometry.Intersect import Intersection


# ---------------------------------------------------------------- shapes

@pytest.mark.parametrize("make", [
    lambda: shapes.faceted_part(),
    lambda: shapes.faceted_part(sides=7, width=3.0),
    lambda: shapes.box((1.0, 2.0, 3.0)),
    lambda: shapes.l_block(),
    lambda: shapes.stacked_loops([(0, 8, 1.0, 0.5), (1, 8, 1.0, 0.5)]),
])
def test_clean_shapes_are_closed_and_outward(make):
    m = make()
    assert m.euler_characteristic() == 2
    # Closed and consistently oriented: every directed edge appears once, reversed once.
    directed = set()
    for a, b, c in m.faces:
        for e in ((a, b), (b, c), (c, a)):
            assert e not in directed
            directed.add(e)
    assert all((j, i) in directed for i, j in directed)
    assert m.volume() > 0


def test_known_volumes():
    assert abs(shapes.box((1.0, 2.0, 3.0)).volume() - 6.0) < 1e-12
    assert abs(shapes.l_block(size=2.0, notch=1.0, height=1.0).volume() - 3.0) < 1e-12


def test_placement_keeps_volume_and_shape():
    g = shapes.rng(0)
    m = shapes.faceted_part()
    placed, rot = shapes.randomly_placed(m, g)
    assert abs(placed.volume() - m.volume()) < 1e-9
    assert np.allclose(rot.dot(rot.T), np.eye(3))


def test_damage_helpers():
    g = shapes.rng(1)
    m = shapes.faceted_part()
    assert len(shapes.with_holes(m, g, 2).faces) == len(m.faces) - 2
    assert len(shapes.with_duplicate_faces(m, g, 3).faces) == len(m.faces) + 3
    flipped = shapes.with_flipped_faces(m, g, 1)
    assert abs(flipped.volume() - m.volume()) > 1e-6
    junk = shapes.with_junk(m, g, count=2)
    assert len(junk.faces) == len(m.faces) + 8
    f32 = shapes.as_float32(shapes.moved(m, translation=(1000.1, 0, 0)))
    assert np.abs(f32.vertices - shapes.moved(m, translation=(1000.1, 0, 0)).vertices).max() > 0
    u = shapes.unwelded(m)
    assert len(u.vertices) == 3 * len(m.faces)
    assert abs(u.volume() - m.volume()) < 1e-12


# ---------------------------------------------------------------- Rhino mock

def test_mock_is_installed():
    import Rhino
    assert getattr(Rhino, "IS_MOCK", False)


def test_topology_vertices_weld_unwelded_mesh():
    m = rg.Mesh.from_trimesh(shapes.unwelded(shapes.box()))
    assert m.Vertices.Count == 36
    top = m.TopologyVertices
    assert top.Count == 8
    for t in range(top.Count):
        members = top.MeshVertexIndices(t)
        assert all(top.TopologyVertexIndex(i) == t for i in members)


def test_mesh_plane_section_of_box_is_one_closed_square():
    m = rg.Mesh.from_trimesh(shapes.box((2.0, 2.0, 2.0)))
    plane = rg.Plane(rg.Point3d(0, 0, 0.7), rg.Vector3d(0, 0, 1))
    lines = Intersection.MeshPlane(m, plane)
    assert len(lines) == 1
    line = lines[0]
    assert line.IsClosed
    pts = np.array([[p.X, p.Y, p.Z] for p in line[:-1]])
    assert np.allclose(pts[:, 2], 0.7)
    x, y = pts[:, 0], pts[:, 1]
    area = 0.5 * abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))
    assert abs(area - 4.0) < 1e-12


def test_mesh_plane_two_pieces_and_a_miss():
    two = shapes.combined(shapes.box(), shapes.moved(shapes.box(), translation=(5, 0, 0)))
    m = rg.Mesh.from_trimesh(two)
    lines = Intersection.MeshPlane(m, rg.Plane(rg.Point3d(0, 0, 0.5), rg.Vector3d(0, 0, 1)))
    assert len(lines) == 2 and all(l.IsClosed for l in lines)
    assert Intersection.MeshPlane(m, rg.Plane(rg.Point3d(0, 0, 9), rg.Vector3d(0, 0, 1))) is None


def test_mesh_plane_open_mesh_gives_open_chain():
    box = shapes.box((1.0, 1.0, 1.0))
    # Remove the two side triangles of the x = 0 wall: the section becomes an open U.
    side = [i for i, f in enumerate(box.faces) if np.allclose(box.vertices[f][:, 0], 0.0)]
    open_box = shapes.TriMesh(box.vertices, np.delete(box.faces, side, axis=0))
    lines = Intersection.MeshPlane(rg.Mesh.from_trimesh(open_box),
                                   rg.Plane(rg.Point3d(0, 0, 0.5), rg.Vector3d(0, 0, 1)))
    assert len(lines) == 1 and not lines[0].IsClosed


def test_plane_to_plane_transform():
    a = rg.Plane.WorldXY
    b = rg.Plane(rg.Point3d(1, 2, 3), rg.Vector3d(1, 0, 0), rg.Vector3d(0, 0, 1))
    x = rg.Transform.PlaneToPlane(a, b)
    p = rg.Point3d(1, 0, 0)
    p.Transform(x)
    assert p.DistanceTo(rg.Point3d(2, 2, 3)) < 1e-12
    ok, inv = x.TryGetInverse()
    assert ok
    p.Transform(inv)
    assert p.DistanceTo(rg.Point3d(1, 0, 0)) < 1e-12


def test_point3f_rounds_to_single_precision():
    p = rg.Point3f(0.1, 1000.1, 0.0)
    assert p.X == float(np.float32(0.1))
    assert abs(p.Y - 1000.1) > 1e-6
    assert math.isfinite(p.Z)
