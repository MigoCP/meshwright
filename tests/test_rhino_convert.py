import Rhino.Geometry as rg

from meshwright.core import topology
from meshwright.rhino import convert

import shapes


def test_round_trip_keeps_geometry():
    m = shapes.faceted_part()
    v, f = convert.mesh_to_lists(convert.lists_to_mesh(*m.as_lists()))
    assert len(v) == len(m.vertices) and len(f) == len(m.faces)
    assert abs(topology.signed_volume(v, f) - m.volume()) < 1e-14  # double precision kept


def test_single_precision_meshes_still_read():
    mesh = rg.Mesh()
    mesh.Vertices.UseDoublePrecisionVertices = False
    for p in [(1.0 / 3.0, 0, 0), (1, 0, 0), (0, 1, 0)]:
        mesh.Vertices.Add(*p)
    mesh.Faces.AddFace(0, 1, 2)
    v, _ = convert.mesh_to_lists(mesh)
    assert v[0][0] != 1.0 / 3.0 and abs(v[0][0] - 1.0 / 3.0) < 1e-7
    exact = convert.lists_to_mesh([(1.0 / 3.0, 0, 0), (1, 0, 0), (0, 1, 0)], [(0, 1, 2)])
    assert convert.mesh_to_lists(exact)[0][0][0] == 1.0 / 3.0


def test_split_vertices_arrive_welded():
    rhino_mesh = rg.Mesh.from_trimesh(shapes.unwelded(shapes.box()))
    v, f = convert.mesh_to_lists(rhino_mesh)
    assert len(v) == 8
    assert topology.EdgeCensus(f).closed


def test_quads_split_on_the_shorter_diagonal():
    m = rg.Mesh()
    for p in [(0, 0, 0), (4, 0, 0), (4, 1, 0), (0, 1, 0)]:
        m.Vertices.Add(*p)
    m.Faces.AddFace(0, 1, 2, 3)                 # diagonals 0-2 and 1-3 have equal length
    v, f = convert.mesh_to_lists(m)
    assert len(f) == 2
    m2 = rg.Mesh()
    for p in [(0, 0, 0), (4, 0, 0), (5, 1, 0), (0, 1, 0)]:
        m2.Vertices.Add(*p)
    m2.Faces.AddFace(0, 1, 2, 3)                # 1-3 (length ~4.1) is shorter than 0-2 (~5.1)
    v, f = convert.mesh_to_lists(m2)
    assert sorted(map(sorted, f)) == [[0, 1, 3], [1, 2, 3]]
    assert abs(topology.surface_area(v, f) - 4.5) < 1e-9


def test_document_tolerance_outside_rhino():
    assert convert.document_tolerance(0.01) == 0.01


def test_from_vector_reads_vectors_planes_and_sequences():
    # A Vector3d also has (static) XAxis/YAxis/ZAxis under pythonnet: it must not be read as a plane.
    assert convert.from_vector(rg.Vector3d(1.0, 2.0, 3.0)) == (1.0, 2.0, 3.0)
    plane = rg.Plane(rg.Point3d(5, 5, 5), rg.Vector3d(1.0, 0.0, 0.0))
    assert convert.from_vector(plane) == (1.0, 0.0, 0.0)
    assert convert.from_vector([0, 0, 2]) == (0.0, 0.0, 2.0)
    assert convert.from_vector(None) is None and convert.from_vector("x") is None


def test_polylines_and_points():
    line = convert.to_polyline([(0, 0, 0), (1, 0, 0), (1, 1, 0)], closed=True)
    assert line.Count == 4 and line.IsClosed
    assert convert.to_polyline([(0, 0, 0), (1, 0, 0)]).Count == 2
    p = convert.to_point((1.0, 2.0, 3.0))
    assert (p.X, p.Y, p.Z) == (1.0, 2.0, 3.0)

