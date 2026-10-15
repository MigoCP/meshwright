import math

import pytest

from meshwright.core import polyhedron as ph
from meshwright.core import rebuild as rb
from meshwright.core import vec
from meshwright.core.topology import EdgeCensus, signed_volume

import shapes


def run(mesh, tolerance=0.001):
    v, f = mesh.as_lists()
    return rb.rebuild(v, f, tolerance)


def assert_solid(poly):
    assert ph.closed_problems(poly) == []
    v, tris = poly.triangles()
    assert EdgeCensus(tris).closed
    assert abs(signed_volume(v, tris) - poly.volume()) < 1e-9 * abs(poly.volume())


# ---------------------------------------------------------------- corner solve

def test_three_planes_meet_exactly():
    planes = [((1.0, 0.0, 0.0), 2.0), ((0.0, 1.0, 0.0), 3.0), (vec.unit((0.0, 1.0, 1.0)), 5.0)]
    p = rb.solve_corner(planes, (0.0, 0.0, 0.0))
    for n, d in planes:
        assert abs(vec.dot(n, p) - d) < 1e-12


def test_four_planes_least_squares():
    # x = 1 and x = -1 cannot both hold: the least squares point sits between them.
    planes = [((1.0, 0.0, 0.0), 1.0), ((-1.0, 0.0, 0.0), 1.0),
              ((0.0, 1.0, 0.0), 0.0), ((0.0, 0.0, 1.0), 0.0)]
    p = rb.solve_corner(planes, (5.0, 5.0, 5.0))
    assert max(abs(c) for c in p) < 1e-12


def test_directions_the_planes_do_not_fix_keep_the_start():
    # Two planes meet in a line along Z: the corner keeps its own z.
    planes = [((1.0, 0.0, 0.0), 1.0), ((0.0, 1.0, 0.0), 2.0)]
    p = rb.solve_corner(planes, (0.3, 0.4, 7.0))
    assert max(abs(a - b) for a, b in zip(p, (1.0, 2.0, 7.0))) < 1e-12
    # Two planes 0.5 degrees apart barely fix the direction across them: their exact
    # crease is at x = 0.0044, but the corner stays near its own x = 0.2.
    tilt = vec.unit((math.sin(math.radians(0.5)), 0.0, math.cos(math.radians(0.5))))
    p = rb.solve_corner([((0.0, 0.0, 1.0), 1.0), (tilt, 1.0)], (0.2, 0.0, 1.0))
    assert abs(p[0] - 0.2) < 1e-4 and abs(p[2] - 1.0) < 1e-3


# ---------------------------------------------------------------- clean meshes

@pytest.mark.parametrize("mesh, faces, corners, holed", [
    (shapes.box((2.0, 3.0, 1.0)), 6, 8, 0),
    (shapes.l_block(), 8, 12, 0),
    (shapes.u_block(), 10, 16, 0),
    (shapes.faceted_part(), 65, 49, 0),
    (shapes.pyramid(), 5, 5, 0),
    (shapes.frame(), 10, 16, 2),
    (shapes.plate_with_round_hole(), 30, 56, 2),
], ids=["box", "L", "U", "faceted part", "pyramid", "frame", "round hole"])
def test_clean_meshes_rebuild_exactly(mesh, faces, corners, holed):
    for m in (mesh, shapes.subdivided(mesh, 2)):  # interior and crease vertices are dropped
        res = run(m)
        assert res.status == 0, str(res.report)
        s = res.solid
        assert len(s.faces) == faces and len(s.vertices) == corners
        assert sum(1 for h in s.holes if h) == holed
        assert abs(s.volume() - mesh.volume()) < 1e-9
        assert_solid(s)


def test_coplanar_faces_that_do_not_touch_stay_separate():
    # The ends of the U's two arms lie on one plane but are two faces.
    res = run(shapes.u_block())
    facing = [k for k in range(len(res.solid.faces)) if res.solid.plane(k)[0][1] > 0.999]
    offsets = sorted(round(res.solid.plane(k)[1], 9) for k in facing)
    assert offsets == [1.0, 2.0, 2.0]  # the bottom of the slot, then both arm ends


def test_a_sliver_below_the_tolerance_joins_its_neighbour():
    # A box corner cut off by 1e-4: the sliver lies within the tolerance of a side.
    e = 1e-4
    pts = [(x, y, z) for x in (0, 1) for y in (0, 1) for z in (0, 1) if (x, y, z) != (1, 1, 1)]
    pts += [(1 - e, 1, 1), (1, 1 - e, 1), (1, 1, 1 - e)]
    res = run(shapes.convex_hull(pts))
    assert res.status == 0
    assert len(res.solid.faces) == 6 and len(res.solid.vertices) == 8
    assert abs(res.solid.volume() - 1.0) < 1e-12  # the joined sliver does not tilt the face
    # With a tolerance below the cut, the sliver is a face of its own.
    res = run(shapes.convex_hull(pts), tolerance=1e-5)
    assert res.status == 0 and len(res.solid.faces) == 7


def test_a_step_below_the_tolerance_is_flattened():
    # One side of this block steps out by 2e-4: within half the tolerance, so one face.
    step = shapes.prism([(0, 0), (2, 0), (2, 1), (1, 1), (1, 1.0002), (0, 1.0002)], 1.0)
    res = run(step)
    assert res.status == 0 and len(res.solid.faces) == 6
    assert abs(res.solid.volume() - step.volume()) < 1e-6
    # With a tolerance below the step, the step is kept.
    assert len(run(step, tolerance=1e-4).solid.faces) == 8


def test_corners_closer_than_the_tolerance_are_merged():
    # A coarse tolerance (0.2 on a part 1 across) joins faces into few big ones,
    # and two of the corners left end up closer than the tolerance.
    res = run(shapes.faceted_part(), tolerance=0.2)
    assert res.status == 0 and res.merged
    assert "corners closer than the tolerance merged" in str(res.report)
    assert min(res.solid.edge_lengths()) >= 0.2
    # Faces are only flat within 0.1 here, so check closure, not volume agreement.
    assert ph.closed_problems(res.solid) == []
    assert EdgeCensus(res.solid.triangles()[1]).closed


# ---------------------------------------------------------------- placed, single precision, noisy

@pytest.mark.parametrize("sigma", [0.0, 1e-5, 1e-4])
def test_placed_noisy_float32_meshes(sigma):
    g = shapes.rng(41)
    for mesh, faces in ((shapes.box((2.0, 3.0, 1.0)), 6), (shapes.l_block(), 8), (shapes.faceted_part(), 65),
                        (shapes.frame(), 10), (shapes.plate_with_round_hole(), 30)):
        for _ in range(3):
            placed, _ = shapes.randomly_placed(shapes.subdivided(mesh, 2), g, spread=10.0)
            if sigma:
                placed = shapes.noisy(placed, g, sigma)
            res = run(shapes.as_float32(placed))
            assert res.status == 0, str(res.report)
            assert len(res.solid.faces) == faces
            assert abs(res.solid.volume() - mesh.volume()) < 2e-3 * mesh.volume()
            assert_solid(res.solid)


def test_noise_beyond_the_tolerance_is_never_called_good():
    # Noise of 0.3 x the tolerance: faces can no longer always be joined. Whatever
    # comes out, a result with the wrong faces must not be graded good.
    g = shapes.rng(42)
    for mesh, faces in ((shapes.faceted_part(), 65), (shapes.plate_with_round_hole(), 30), (shapes.frame(), 10)):
        for _ in range(4):
            placed, _ = shapes.randomly_placed(shapes.subdivided(mesh, 2), g, spread=10.0)
            res = run(shapes.as_float32(shapes.noisy(placed, g, 3e-4)))
            if res.solid is None or len(res.solid.faces) != faces:
                assert res.status >= 1
                assert "noisier than the tolerance" in str(res.report)


# ---------------------------------------------------------------- refusals

def test_refuses_meshes_it_cannot_trust():
    v, f = shapes.box().as_lists()
    assert "not closed" in str(rb.rebuild(v, f[:-1]).report)
    flipped = f[:-1] + [tuple(reversed(f[-1]))]
    assert "disagree" in str(rb.rebuild(v, flipped).report)
    inside_out = [tuple(reversed(t)) for t in f]
    assert "inside out" in str(rb.rebuild(v, inside_out).report)
    assert "no faces" in str(rb.rebuild([], []).report)
    assert "more than zero" in str(rb.rebuild(v, f, tolerance=0.0).report)
