import math

import pytest

from solidwright.core import poly2d
from solidwright.core import polyhedron as ph
from solidwright.core import sections as sc
from solidwright.core import vec

import shapes

UP = (0.0, 0.0, 1.0)


def solid(mesh, clearance=0.0, axis=UP, tolerance=0.001, **kw):
    v, f = mesh.as_lists()
    return sc.clearance_solid(v, f, axis, clearance, tolerance, **kw)


def local(res, p):
    d = vec.sub(p, res.origin)
    return tuple(vec.dot(d, a) for a in res.frame)


def worst_slack(res, mesh):
    v, f = mesh.as_lists()
    return min(res.model.slack(local(res, p)) for p in v if res.model.slack(local(res, p)) is not None)


def notes(res):
    return " | ".join(t for _, t in res.report.reasons)


# ---------------------------------------------------------------- slicing and the model

def test_slices_of_a_box():
    v, f = shapes.box((2.0, 3.0, 1.0)).as_lists()
    hulls = sc.section_hulls(v, f, [0.5, 0.25, 2.0])
    assert abs(poly2d.area(hulls[0]) - 6.0) < 1e-12 and abs(poly2d.area(hulls[1]) - 6.0) < 1e-12
    assert hulls[2] is None
    # A plane through the top face: the face itself is skipped, its sides still give the rim.
    top = sc.section_hulls(v, f, [1.0])[0]
    assert abs(poly2d.area(top) - 6.0) < 1e-12


def test_model_corners_offsets_and_slack():
    m = sc.SectionModel(sc.fan(4), (0.0, 0.0), 0.0)
    square = m.corners([1.0, 1.0, 1.0, 1.0])
    assert abs(poly2d.area(square) - 4.0) < 1e-12
    m.levels = [(0.0, [1.0] * 4), (2.0, [2.0] * 4)]
    assert m.offsets_at(1.0) == [1.5] * 4 and m.offsets_at(3.0) is None
    assert abs(m.slack((0.5, 0.0, 1.0)) - 1.0) < 1e-12
    assert m.slack((2.0, 0.0, 1.0)) < 0.0
    assert m.clean([1.0] * 4, 0.0) and not m.clean([1.0, 1.0, -1.5, 1.0], 0.0)


def test_pinch_and_safe_blend():
    m = sc.SectionModel(sc.fan(4), (0.0, 0.0), 0.0)
    # Widths 4 at t = 0 and 2 at t = 1 reach 0.5 at t = 1.75.
    t = sc.pinch(m, [2.0] * 4, [1.0] * 4, 0.5)
    assert abs(t - 1.75) < 1e-12
    assert sc.pinch(m, [1.0] * 4, [2.0] * 4, 0.5) is None  # widening: never pinches
    hs, used = sc.safe_blend(m, [2.0] * 4, [1.0] * 4, 3.0, 0.1, 0.0)
    assert used < 3.0 and min(hs) * 2 >= 0.1  # backed off before turning inside out


# ---------------------------------------------------------------- clearance solids

@pytest.mark.parametrize("clearance", [0.0, 0.05])
def test_box_is_exact(clearance):
    res = solid(shapes.box((2.0, 3.0, 1.0)), clearance)
    assert res.status == 0
    assert "4-sided polygon" in notes(res)  # sections follow the box's own edges: sharp corners
    c = clearance
    assert abs(res.solid.volume() - (2 + 2 * c) * (3 + 2 * c) * 1.0) < 1e-9
    assert len(res.solid.faces) == 6 and ph.closed_problems(res.solid) == []


@pytest.mark.parametrize("name, mesh, end", [
    ("faceted part", shapes.faceted_part(), "the bottom tapers to a point"),
    ("pyramid", shapes.pyramid(), "the top tapers to a point"),
    ("wedge", shapes.wedge(), "the bottom tapers to a ridge"),
    ("oval part", shapes.oval_part(), "the bottom tapers to a point"),
    ("cylinder", shapes.cylinder(), None),
    ("ball", shapes.ball(), None),
    ("hourglass", shapes.hourglass(), None),
])
@pytest.mark.parametrize("clearance", [0.0, 0.03])
def test_parts_get_a_closed_solid_that_keeps_the_clearance(name, mesh, end, clearance):
    res = solid(mesh, clearance)
    assert res.status == 0, str(res.report)
    assert ph.closed_problems(res.solid) == []
    if end:
        assert end in notes(res)
    # Every vertex, each at its own height: at least 90% of the clearance (the tolerance when 0).
    assert worst_slack(res, mesh) >= sc.KEEP_CLEARANCE * clearance - 0.001
    assert res.solid.volume() >= mesh.volume() - 1e-3


def test_round_parts_use_a_fan_and_stay_close():
    res = solid(shapes.cylinder(sides=96), 0.05)
    assert "96 support directions" in notes(res)
    # A 96-direction fan around a 96-gon: within 0.1% of the true offset cylinder.
    want = math.pi * 1.05 ** 2 * 2.0
    assert abs(res.solid.volume() - want) < 1e-3 * want


def test_refinement_follows_a_curved_profile():
    res = solid(shapes.ball(), 0.0)
    assert res.added > 0
    assert worst_slack(res, shapes.ball()) >= -0.001


def test_the_profile_may_go_in_and_out():
    # An hourglass is narrow in the middle: its sections follow the waist.
    res = solid(shapes.hourglass(), 0.02)
    waist = res.model.offsets_at(0.0)  # local z = 0 is the middle height
    assert max(waist) < 0.45
    assert res.solid.volume() < 1.1 * shapes.hourglass().volume()


@pytest.mark.parametrize("seed", range(3))
def test_placed_noisy_single_precision(seed):
    g = shapes.rng(seed)
    for mesh in (shapes.faceted_part(), shapes.wedge(), shapes.ball(), shapes.oval_part()):
        placed, rot = shapes.randomly_placed(shapes.subdivided(mesh, 2), g, spread=10.0)
        placed = shapes.as_float32(shapes.noisy(placed, g, 1e-4))
        axis = tuple(rot.dot(UP))
        res = solid(placed, 0.03, axis=axis)
        assert res.report.ok, str(res.report)
        assert ph.closed_problems(res.solid) == []
        assert res.min_clearance >= 0.9 * 0.03 - 0.001


def test_tilted_axis_gives_the_same_solid():
    base = solid(shapes.faceted_part(), 0.03).solid.volume()
    g = shapes.rng(4)
    placed, rot = shapes.randomly_placed(shapes.faceted_part(), g)
    res = solid(placed, 0.03, axis=tuple(rot.dot(UP)))
    assert abs(res.solid.volume() - base) < 1e-6 * base


@pytest.mark.parametrize("kwargs, text", [
    ({"clearance": -1.0}, "zero or more"),
    ({"tolerance": 0.0}, "more than zero"),
    ({"axis": (0.0, 0.0, 0.0)}, "zero vector"),
])
def test_bad_inputs(kwargs, text):
    args = dict(clearance=0.0, axis=UP, tolerance=0.001)
    args.update(kwargs)
    v, f = shapes.box().as_lists()
    res = sc.clearance_solid(v, f, args["axis"], args["clearance"], args["tolerance"])
    assert res.status == 3 and text in str(res.report)


def test_no_faces_and_flat_parts():
    assert "no faces" in str(sc.clearance_solid([], [], UP).report)
    flat = [(0, 0, 0), (1, 0, 0), (0, 1, 0)]
    assert "flat along the axis" in str(sc.clearance_solid(flat, [(0, 1, 2)], UP).report)


# ---------------------------------------------------------------- side profile

def test_profile_knees_land_on_the_band_edges():
    m = shapes.faceted_part()
    v, f = m.as_lists()
    res = sc.side_profile(v, f, UP, (1.0, 0.0, 0.0))
    assert res.status == 0
    (z0, x0), (z1, x1) = res.fit.knees
    top = res.to_world(x1, z1)
    bottom = res.to_world(x0, z0)
    assert abs(bottom[2] - 0.0) < 1e-6 and abs(top[2] - 0.02) < 1e-6  # the widest band, z 0 to 0.02
    assert abs(bottom[0] - 0.5) < 1e-6


def test_profile_of_a_noisy_dense_mesh():
    g = shapes.rng(6)
    placed, rot = shapes.randomly_placed(shapes.subdivided(shapes.faceted_part(), 2), g)
    v, f = shapes.noisy(placed, g, 1e-4).as_lists()
    res = sc.side_profile(v, f, tuple(rot.dot(UP)), tuple(rot.dot((1.0, 0.0, 0.0))))
    assert res.status == 0
    (z0, _), (z1, _) = res.fit.knees
    assert abs((z1 - z0) - 0.02) < 2e-3


def test_profile_reports_bad_inputs_and_slope_limits():
    v, f = shapes.faceted_part().as_lists()
    assert "to a side" in str(sc.side_profile(v, f, UP, UP).report)
    res = sc.side_profile(v, f, UP, (1.0, 0.0, 0.0), slope_limits=[(0.0, 0.1), None, None])
    assert "hit its slope limit" in str(res.report)
