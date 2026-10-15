import Rhino.Geometry as rg

from meshwright.core import topology
from meshwright.rhino import convert

import shapes


def test_round_trip_keeps_geometry():
    m = shapes.faceted_part()
    v, f = convert.mesh_to_lists(convert.lists_to_mesh(*m.as_lists()))
    assert len(v) == len(m.vertices) and len(f) == len(m.faces)
    assert abs(topology.signed_volume(v, f) - m.volume()) < 1e-5  # Point3f storage


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
