import math

import numpy as np

from meshwright.core import facets as fac

import shapes


def angle(a, b):
    return math.degrees(math.acos(min(1.0, abs(float(np.dot(a, b))))))


def test_box_has_six_facets_with_exact_planes():
    v, f = shapes.box((1.0, 2.0, 3.0)).as_lists()
    result = fac.cluster_planes(v, f)
    assert len(result) == 6
    assert [round(fc.area, 12) for fc in result] == [6.0, 6.0, 3.0, 3.0, 2.0, 2.0]
    for fc in result:
        assert len(fc.faces) == 2
        assert fc.rms < 1e-12
        # Normals point outward: the box centre is behind every facet.
        assert fc.distance((0.5, 1.0, 1.5)) < 0


def test_subdivided_flat_regions_stay_one_facet_each():
    m = shapes.subdivided(shapes.box(), 3)
    v, f = m.as_lists()
    result = fac.cluster_planes(v, f)
    assert len(result) == 6
    assert sum(len(fc.faces) for fc in result) == len(f)


def test_reversed_triangles_join_their_facet():
    g = shapes.rng(0)
    m = shapes.with_flipped_faces(shapes.subdivided(shapes.box(), 2), g, 10)
    v, f = m.as_lists()
    result = fac.cluster_planes(v, f)
    assert len(result) == 6
    for fc in result:
        assert fc.distance((0.5, 0.5, 0.5)) < 0   # direction follows the majority


def test_planes_are_refitted_and_rotation_invariant():
    g = shapes.rng(1)
    base = shapes.subdivided(shapes.faceted_part(sides=12), 2)
    top_normal = (0.0, 0.0, 1.0)
    for _ in range(5):
        placed, rot = shapes.randomly_placed(base, g)
        v, f = shapes.as_float32(placed).as_lists()
        result = fac.cluster_planes(v, f)
        tops = [fc for fc in result if angle(fc.normal, rot.dot(top_normal)) < 0.01]
        assert len(tops) == 1
        assert np.dot(tops[0].normal, rot.dot(top_normal)) > 0


def test_curved_surface_does_not_creep_into_one_facet():
    # A fine 64-sided band: neighbours differ by 5.6 degrees, so no merging.
    v, f = shapes.stacked_loops([(0.0, 64, 1.0), (0.5, 64, 1.0)]).as_lists()
    result = fac.cluster_planes(v, f, angle_deg=1.0)
    sides = [fc for fc in result if abs(fc.normal[2]) < 0.5]
    assert len(sides) == 64


def test_gentle_curve_is_not_chained_together():
    # Normals turn 0.5 degree per step: each step is within tolerance of its
    # neighbour, but membership is tested against the seed, so a facet spans
    # only a few steps instead of the whole arc.
    n = 60
    v, f = [], []
    for k in range(n + 1):
        a = math.radians(0.5 * k)
        x, z = math.sin(a) * 50.0, (1.0 - math.cos(a)) * 50.0
        v += [(x, 0.0, z), (x, 1.0, z)]
    for k in range(n):
        a, b, c, d = 2 * k, 2 * k + 1, 2 * k + 2, 2 * k + 3
        f += [(a, c, b), (b, c, d)]
    result = fac.cluster_planes(v, f, angle_deg=1.0)
    assert len(result) >= 15
    assert max(len(fc.faces) for fc in result) <= 10


def test_offset_tolerance_separates_parallel_steps():
    # Two parallel flat steps 0.01 apart: one facet each with a tight tolerance.
    v = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
         (2, 0, 0.01), (3, 0, 0.01), (3, 1, 0.01), (2, 1, 0.01)]
    f = [(0, 1, 2), (0, 2, 3), (4, 5, 6), (4, 6, 7)]
    assert len(fac.cluster_planes(v, f, offset_tol=0.001)) == 2
    assert len(fac.cluster_planes(v, f, offset_tol=0.05)) == 1


def test_zero_area_faces_are_skipped():
    v = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (2, 0, 0)]
    f = [(0, 1, 2), (0, 1, 3)]
    result = fac.cluster_planes(v, f)
    assert len(result) == 1 and result[0].faces == [0]
