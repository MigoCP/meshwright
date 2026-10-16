import math

import numpy as np
import pytest

from solidwright.core import orient
from solidwright.core.facets import cluster_planes
from solidwright.core.report import FAILED, GOOD, MESSY, WARNING

import shapes

Z = (0.0, 0.0, 1.0)


def angle(a, b):
    """Angle between two lines (sign ignored), degrees."""
    c = abs(float(np.dot(a, b)) / (np.linalg.norm(a) * np.linalg.norm(b)))
    return math.degrees(math.acos(min(1.0, c)))


def points_up(axis, truth):
    return float(np.dot(axis, truth)) > 0.0


def tall_tip_part():
    """Flat top, widest band near it, then a long point: its side view is
    larger than its top view, so outline area alone would pick a side view."""
    return shapes.stacked_loops([
        (0.16, 16, 0.275, 0.275, math.pi / 16),
        (0.02, 16, 0.5),
        (0.0, 16, 0.5),
        (-2.0, 1, 0.0),
    ])


def bicone():
    """No flat end at all: a short point up, a long point down."""
    return shapes.stacked_loops([(0.3, 1, 0.0), (0.0, 32, 0.5), (-0.9, 1, 0.0)])


# ---------------------------------------------------------------- accuracy

@pytest.mark.parametrize("make", [shapes.faceted_part, tall_tip_part, bicone])
def test_finds_axis_and_top_at_any_orientation(make):
    g = shapes.rng(10)
    base = make()
    for _ in range(8):
        placed, rot = shapes.randomly_placed(base, g, spread=100.0)
        v, f = shapes.as_float32(placed).as_lists()
        res = orient.natural_axis(v, f)
        truth = rot.dot(Z)
        assert res.status == GOOD, res.report.summary()
        assert angle(res.axis, truth) < 0.01
        assert points_up(res.axis, truth)


def test_frame_is_right_handed_around_the_axis():
    v, f = shapes.faceted_part().as_lists()
    res = orient.natural_axis(v, f)
    x, y, z = res.frame
    assert np.allclose(np.cross(x, y), z) and np.allclose(z, res.axis)


@pytest.mark.parametrize("sigma,limit", [(0.001, 0.3), (0.004, 0.5), (0.01, 1.0)])
def test_dense_noisy_mesh(sigma, limit):
    g = shapes.rng(11)
    base = shapes.subdivided(shapes.faceted_part(sides=32), 2)
    for _ in range(4):
        placed, rot = shapes.randomly_placed(base, g)
        v, f = shapes.noisy(placed, g, sigma).as_lists()
        res = orient.natural_axis(v, f)
        truth = rot.dot(Z)
        assert res.report.ok, res.report.summary()
        assert angle(res.axis, truth) < limit
        assert points_up(res.axis, truth)


def test_noisy_flat_end_falls_back_to_the_outline_and_says_so():
    g = shapes.rng(13)
    base = shapes.subdivided(shapes.faceted_part(sides=32), 3)
    top = base.vertices[:, 2].max()
    on_top = np.where(np.isclose(base.vertices[:, 2], top))[0]
    rough = base.copy()
    pick = g.choice(on_top, len(on_top) // 5, replace=False)
    rough.vertices[pick, 2] += g.normal(0.0, 0.003, len(pick))
    placed, rot = shapes.randomly_placed(rough, g)
    v, f = placed.as_lists()
    res = orient.natural_axis(v, f)
    assert res.status == WARNING and "noisy" in res.report.summary()
    assert angle(res.axis, rot.dot(Z)) < 0.2 and points_up(res.axis, rot.dot(Z))


# ---------------------------------------------------------------- the mirror-plane trap

def test_mirror_plane_view_looks_level_but_is_not_chosen():
    v, f = tall_tip_part().as_lists()
    facets = cluster_planes(v, f)
    size = orient.bounding_size(v)
    points = orient._sample(v, f)
    # The view square to a mirror plane (here the XZ plane, seen along Y) IS level...
    side = orient.evaluate(points, (0.0, 1.0, 0.0), facets, size)
    top = orient.evaluate(points, Z, facets, size)
    assert side.level and top.level
    # ...and shows MORE outline than the true axis, so area alone would be fooled.
    assert side.outline_area > top.outline_area
    # The flat end is what decides it.
    assert top.strong_end is not None and side.strong_end is None
    res = orient.natural_axis(v, f)
    assert angle(res.axis, Z) < 1e-6 and res.status == GOOD


def test_true_tie_is_reported_not_guessed():
    # A tall square bar: both side views are equally good answers.
    v, f = shapes.box((1.0, 1.0, 3.0)).as_lists()
    res = orient.natural_axis(v, f)
    assert res.status == MESSY
    assert "equally well" in res.report.summary()


# ---------------------------------------------------------------- other shapes

def test_plate_axis_is_its_normal():
    v, f = shapes.box((4.0, 3.0, 0.2)).as_lists()
    res = orient.natural_axis(v, f)
    assert angle(res.axis, Z) < 1e-9 and res.status == GOOD
    assert "arbitrary" in str(res.report)          # both faces alike


def test_prism_axis_runs_along_its_walls():
    v, f = shapes.l_block(size=2.0, notch=1.0, height=1.0).as_lists()
    res = orient.natural_axis(v, f)
    assert angle(res.axis, Z) < 1e-9 and res.status == GOOD


def test_no_flat_end_points_away_from_the_taper():
    v, f = bicone().as_lists()
    res = orient.natural_axis(v, f)
    assert angle(res.axis, Z) < 0.01 and res.axis[2] > 0
    assert "tapering" in str(res.report)


def test_no_level_direction_is_messy():
    g = shapes.rng(12)
    blob = shapes.convex_hull(g.normal(size=(60, 3)) * (1.0, 1.6, 0.7))
    v, f = blob.as_lists()
    res = orient.natural_axis(v, f)
    assert res.axis is not None
    assert res.status == MESSY


def test_empty_mesh_fails():
    assert orient.natural_axis([], []).status == FAILED


# ---------------------------------------------------------------- checking a given axis

def test_check_axis_grades_the_difference():
    v, f = shapes.faceted_part().as_lists()
    tilted = lambda deg: (math.sin(math.radians(deg)), 0.0, math.cos(math.radians(deg)))
    assert orient.check_axis(v, f, tilted(0.2)).status == GOOD
    assert orient.check_axis(v, f, tilted(3.0)).status == WARNING
    assert orient.check_axis(v, f, tilted(20.0)).status == MESSY


def test_check_axis_keeps_the_given_axis_signed_to_the_top():
    v, f = shapes.faceted_part().as_lists()
    res = orient.check_axis(v, f, (0.0, 0.0, -5.0))
    assert np.allclose(res.axis, Z)
    assert orient.check_axis(v, f, (0.0, 0.0, 0.0)).status == FAILED


def test_clean_dense_flat_end_next_to_a_slope_is_not_called_noisy():
    # Thin sloped triangles just below a flat end sit within the outline band;
    # they must not count against the flat end (a false "noisy" warning).
    v, f = shapes.subdivided(shapes.faceted_part(sides=48), 3).as_lists()
    res = orient.natural_axis(v, f)
    assert res.status == GOOD, res.report.summary()
    assert angle(res.axis, Z) < 1e-6
