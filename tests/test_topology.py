import math

import numpy as np
import pytest

from solidwright.core import topology as topo
from solidwright.core.report import FAILED, GOOD, MESSY, WARNING

import shapes


def lists(mesh):
    return mesh.as_lists()


def check(mesh, **kw):
    v, f = lists(mesh)
    return topo.health(v, f, **kw)


def is_closed_and_outward(h):
    census = topo.EdgeCensus(h.faces)
    return census.closed and not census.flipped and topo.signed_volume(h.vertices, h.faces) > 0


def euler(h):
    census = topo.EdgeCensus(h.faces)
    return len(h.vertices) - census.edge_count + len(h.faces)


# ---------------------------------------------------------------- building blocks

def test_union_find_groups():
    uf = topo.UnionFind(6)
    uf.union(4, 1)
    uf.union(1, 3)
    uf.union(5, 2)
    assert uf.groups() == [[0], [1, 3, 4], [2, 5]]
    assert not uf.union(3, 4)


def test_measures_of_a_box():
    v, f = lists(shapes.box((1.0, 2.0, 3.0)))
    assert abs(topo.signed_volume(v, f) - 6.0) < 1e-12
    assert abs(topo.surface_area(v, f) - 22.0) < 1e-12
    assert abs(topo.bounding_size(v) - math.sqrt(14.0)) < 1e-12


def test_signed_volume_far_from_origin():
    # Summing tetrahedra against the mesh centre keeps precision at large offsets.
    m = shapes.moved(shapes.box((1.0, 1.0, 1.0)), translation=(1e6, -2e6, 3e6))
    v, f = lists(m)
    assert abs(topo.signed_volume(v, f) - 1.0) < 1e-6


def test_weld_exact_and_with_tolerance():
    box = shapes.box()
    v, f = lists(shapes.unwelded(box))
    wv, wf, remap = topo.weld(v, f)
    assert len(wv) == 8
    assert topo.EdgeCensus(wf).closed
    g = shapes.rng(0)
    noisy = shapes.noisy(shapes.unwelded(box), g, 1e-6)
    v, f = lists(noisy)
    assert len(topo.weld(v, f)[0]) == len(v)          # nothing equal exactly
    wv, wf, _ = topo.weld(v, f, tolerance=1e-4)
    assert len(wv) == 8 and topo.EdgeCensus(wf).closed


def test_weld_is_transitive_and_drops_unused():
    v = [(0, 0, 0), (0.6, 0, 0), (1.2, 0, 0), (9, 9, 9), (5, 5, 5)]
    f = [(0, 1, 2), (2, 3, 0)]
    wv, wf, remap = topo.weld(v, f, tolerance=0.7)
    assert remap[0] == remap[1] == remap[2]
    assert len(wv) == 2                                   # vertex 4 is unused
    assert remap[4] == -1


def test_census_counts():
    v, f = lists(shapes.box())
    c = topo.EdgeCensus(f)
    assert c.edge_count == 18 and c.closed and not c.flipped
    f2 = list(f)
    f2[0] = (f2[0][0], f2[0][2], f2[0][1])
    assert len(topo.EdgeCensus(f2).flipped) == 3


def test_pieces_touching_at_one_vertex_are_two_pieces():
    a = shapes.box()
    b = shapes.moved(shapes.box(), translation=(1.0, 1.0, 1.0))  # shares only corner (1,1,1)
    v, f = lists(shapes.combined(a, b))
    wv, wf, _ = topo.weld(v, f)
    assert len(wv) == 15                                  # the shared corner welded
    assert len(topo.pieces(wf)) == 2


def test_boundary_loops_of_one_removed_face():
    v, f = lists(shapes.box())
    loops = topo.boundary_loops(f[1:])
    assert len(loops) == 1 and len(loops[0]) == 3


# ---------------------------------------------------------------- health: clean input

@pytest.mark.parametrize("make", [
    lambda: shapes.faceted_part(),
    lambda: shapes.box((3.0, 1.0, 0.2)),
    lambda: shapes.l_block(),
])
def test_clean_shapes_are_good(make):
    m = make()
    h = check(m)
    assert h.status == GOOD, str(h.report)
    assert h.closed and is_closed_and_outward(h)
    assert euler(h) == 2
    assert abs(h.stats["volume"] - m.volume()) < 1e-9


def test_status_does_not_depend_on_placement():
    g = shapes.rng(1)
    base = shapes.faceted_part()
    for _ in range(10):
        placed, _ = shapes.randomly_placed(base, g, spread=1000.0)
        h = check(shapes.as_float32(placed))
        assert h.status == GOOD, str(h.report)
        assert abs(h.stats["volume"] - base.volume()) < 1e-4


def test_unwelded_input_is_welded():
    h = check(shapes.unwelded(shapes.faceted_part()))
    assert h.status == GOOD and h.closed
    assert h.stats["welded_vertices"] > 0


# ---------------------------------------------------------------- health: damage

def test_small_hole_is_a_warning():
    g = shapes.rng(2)
    h = check(shapes.with_holes(shapes.faceted_part(sides=32), g, 1))
    assert h.status == WARNING, str(h.report)
    assert h.stats["holes"] == 1 and not h.closed
    assert h.stats["volume"] is None


def test_big_hole_is_messy():
    m = shapes.faceted_part(sides=16)
    top = m.vertices[:, 2].max()
    keep = [f for f in m.faces if not np.allclose(m.vertices[f][:, 2], top)]  # remove the flat end
    h = check(shapes.TriMesh(m.vertices, keep))
    assert h.status == MESSY, str(h.report)
    assert h.stats["holes"] == 1


def test_flipped_faces_are_fixed():
    g = shapes.rng(3)
    base = shapes.faceted_part()
    h = check(shapes.with_flipped_faces(base, g, 3))
    assert h.status == WARNING, str(h.report)
    assert h.stats["faces_reoriented"] == 3
    assert is_closed_and_outward(h)
    assert abs(h.stats["volume"] - base.volume()) < 1e-9


def test_inside_out_mesh_is_reversed():
    base = shapes.faceted_part()
    inside_out = shapes.TriMesh(base.vertices, base.faces[:, ::-1])
    h = check(inside_out)
    assert h.status == WARNING
    assert "inward" in h.report.summary()
    assert is_closed_and_outward(h)


def test_duplicate_faces_are_removed():
    g = shapes.rng(4)
    h = check(shapes.with_duplicate_faces(shapes.box(), g, 2))
    assert h.status == WARNING
    assert h.stats["duplicate_faces"] == 2
    assert h.closed and len(h.faces) == 12


def test_junk_is_dropped():
    g = shapes.rng(5)
    base = shapes.faceted_part()
    h = check(shapes.with_junk(base, g, count=3))
    assert h.status == WARNING, str(h.report)
    assert h.stats["junk_dropped"] == 3 and h.stats["pieces"] == 1
    assert abs(h.stats["volume"] - base.volume()) < 1e-9


def test_junk_is_judged_piece_by_piece():
    # One sizeable extra piece must not stop the crumbs from being dropped.
    g = shapes.rng(6)
    base = shapes.faceted_part()
    second = shapes.moved(shapes.box((0.3, 0.3, 0.3)), translation=(5, 0, 0))
    h = check(shapes.with_junk(shapes.combined(base, second), g, count=2))
    assert h.stats["junk_dropped"] == 2
    assert h.stats["pieces"] == 2
    assert h.status == MESSY


def test_two_parts_are_messy():
    two = shapes.combined(shapes.box(), shapes.moved(shapes.box(), translation=(3, 0, 0)))
    h = check(two)
    assert h.status == MESSY and h.stats["pieces"] == 2


def test_zero_area_faces_are_removed():
    v, f = lists(shapes.box())
    v = v + [(0.5, 0.0, 0.0)]                      # on the edge from (0,0,0) to (1,0,0)
    f = f + [(0, 0, 1), (0, 8, 4)]                 # repeated index, and a collinear sliver
    h = topo.health(v, f)
    assert h.stats["degenerate_faces"] == 2
    assert h.closed and len(h.faces) == 12


def test_non_manifold_fin():
    v, f = lists(shapes.box())
    a, b = f[0][0], f[0][1]
    v = v + [(5.0, 5.0, 5.0)]
    f = f + [(a, b, len(v) - 1)]                   # a fin hanging off an existing edge
    h = topo.health(v, f)
    assert h.stats["non_manifold_edges"] == 1
    assert h.status >= WARNING
    assert "non-manifold" in str(h.report)


def test_orientation_conflict_is_messy():
    # A strip whose ends are joined with a half twist cannot be oriented.
    n = 12
    v, f = [], []
    for k in range(n):
        t = 2 * math.pi * k / n
        w = 0.3 * math.cos(t / 2)
        h_ = 0.3 * math.sin(t / 2)
        r = 2.0
        v.append(((r + w) * math.cos(t), (r + w) * math.sin(t), h_))
        v.append(((r - w) * math.cos(t), (r - w) * math.sin(t), -h_))
    for k in range(n):
        a, b = 2 * k, 2 * k + 1
        if k < n - 1:
            c, d = 2 * k + 2, 2 * k + 3
        else:
            c, d = 1, 0                              # the half twist
        f += [(a, b, d), (a, d, c)]
    hh = topo.health(v, f)
    assert hh.stats["orientation_conflicts"] > 0
    assert hh.status == MESSY


# ---------------------------------------------------------------- health: unusable input

def test_empty_mesh_fails():
    assert topo.health([], []).status == FAILED


def test_bad_index_fails():
    h = topo.health([(0, 0, 0), (1, 0, 0), (0, 1, 0)], [(0, 1, 5)])
    assert h.status == FAILED and "do not exist" in h.report.summary()


def test_non_finite_fails():
    h = topo.health([(0, 0, 0), (1, 0, 0), (0, float("nan"), 0)], [(0, 1, 2)])
    assert h.status == FAILED


def test_all_degenerate_fails():
    h = topo.health([(0, 0, 0), (1, 0, 0), (2, 0, 0)], [(0, 1, 2)])
    assert h.status == FAILED
