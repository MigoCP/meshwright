import math

import numpy as np
import pytest

from meshwright.core import vec

import shapes


def close(a, b, tol=1e-9):
    return all(abs(x - y) <= tol for x, y in zip(a, b))


def degrees_between_axes(a, b):
    """Angle between two lines (sign of the direction ignored), in degrees."""
    c = abs(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))
    return math.degrees(math.acos(min(1.0, c)))


# ---------------------------------------------------------------- basics

def test_basic_arithmetic():
    a, b = (1.0, 2.0, 3.0), (4.0, -5.0, 6.0)
    assert vec.add(a, b) == (5.0, -3.0, 9.0)
    assert vec.sub(a, b) == (-3.0, 7.0, -3.0)
    assert vec.scale(a, 2.0) == (2.0, 4.0, 6.0)
    assert vec.dot(a, b) == 12.0
    assert close(vec.cross(a, b), np.cross(a, b))
    assert vec.length((3.0, 4.0, 0.0)) == 5.0
    assert vec.distance(a, a) == 0.0
    assert vec.lerp(a, b, 0.5) == (2.5, -1.5, 4.5)


def test_integer_input_divides_as_float():
    # Python 2 would floor-divide ints without `from __future__ import division`.
    assert vec.mean([(0, 0, 0), (1, 1, 1)]) == (0.5, 0.5, 0.5)
    assert close(vec.unit((0, 0, 3)), (0.0, 0.0, 1.0))


def test_unit_of_zero_is_none():
    assert vec.unit((0.0, 0.0, 0.0)) is None
    assert vec.unit((1e-15, 0.0, 0.0)) is None


def test_mean_of_nothing_raises():
    with pytest.raises(ValueError):
        vec.mean([])


def test_angle_between():
    assert abs(vec.angle_between((1, 0, 0), (0, 1, 0)) - math.pi / 2) < 1e-12
    assert vec.angle_between((1, 0, 0), (2, 0, 0)) == 0.0
    assert abs(vec.angle_between((1, 0, 0), (-1, 1e-9, 0)) - math.pi) < 1e-6
    with pytest.raises(ValueError):
        vec.angle_between((0, 0, 0), (1, 0, 0))


def test_rotate_matches_scipy():
    g = shapes.rng(1)
    for _ in range(20):
        v = tuple(g.normal(size=3))
        axis = tuple(g.normal(size=3))
        angle = float(g.uniform(-math.pi, math.pi))
        from scipy.spatial.transform import Rotation
        expected = Rotation.from_rotvec(np.array(axis) / np.linalg.norm(axis) * angle).apply(v)
        assert close(vec.rotate(v, axis, angle), expected, 1e-9)


def test_area_vector_of_non_convex_polygon():
    # L-shape of area 3 in the XY plane, counter-clockwise: normal +Z.
    poly = [(0, 0, 0), (2, 0, 0), (2, 1, 0), (1, 1, 0), (1, 2, 0), (0, 2, 0)]
    assert close(vec.area_vector(poly), (0.0, 0.0, 3.0))
    assert close(vec.area_vector(poly[::-1]), (0.0, 0.0, -3.0))


def test_area_vector_tilted_square():
    g = shapes.rng(2)
    rot = shapes.random_rotation(g)
    square = [(0, 0, 0), (2, 0, 0), (2, 2, 0), (0, 2, 0)]
    moved = [tuple(rot.dot(p)) for p in square]
    av = vec.area_vector(moved)
    assert abs(vec.length(av) - 4.0) < 1e-9
    assert close(av, rot.dot((0, 0, 4.0)), 1e-9)


# ---------------------------------------------------------------- frames

@pytest.mark.parametrize("n", [(0, 0, 1), (0, 0, -1), (1, 0, 0), (1, 1, 1), (1e-3, 0, 1)])
def test_perpendicular_pair_is_right_handed(n):
    u, v = vec.perpendicular_pair(n)
    z = vec.unit(n)
    assert abs(vec.length(u) - 1) < 1e-12 and abs(vec.length(v) - 1) < 1e-12
    assert abs(vec.dot(u, v)) < 1e-12 and abs(vec.dot(u, z)) < 1e-12
    assert close(vec.cross(u, v), z, 1e-12)


def test_frame_from_axis_follows_hint():
    x, y, z = vec.frame_from_axis((0, 0, 5))
    assert close(x, (1, 0, 0)) and close(y, (0, 1, 0)) and close(z, (0, 0, 1))
    # Tilted axis: x stays as close to world X as possible.
    x, y, z = vec.frame_from_axis((1, 0, 1))
    assert close(z, (math.sqrt(0.5), 0, math.sqrt(0.5)))
    assert close(x, (math.sqrt(0.5), 0, -math.sqrt(0.5)))
    assert close(vec.cross(x, y), z)


def test_frame_from_axis_hint_parallel_to_axis():
    x, y, z = vec.frame_from_axis((1, 0, 0), x_hint=(3, 0, 0))
    assert abs(vec.dot(x, z)) < 1e-12
    assert close(vec.cross(x, y), z)


def test_angular_order_counter_clockwise_about_normal():
    g = shapes.rng(3)
    pts = shapes.loop(0.0, 12, 2.0, 1.0)  # ellipse, already counter-clockwise about +Z
    order = list(g.permutation(len(pts)))
    shuffled = [pts[i] for i in order]
    result = vec.angular_order(shuffled, normal=(0, 0, 1))
    restored = [order[i] for i in result]
    k = restored.index(0)
    assert restored[k:] + restored[:k] == list(range(12))
    # Seen from below, the order reverses.
    down = [order[i] for i in vec.angular_order(shuffled, normal=(0, 0, -1))]
    k = down.index(0)
    assert down[k:] + down[:k] == [0] + list(range(11, 0, -1))


# ---------------------------------------------------------------- eigen

def test_eig_sym3_matches_numpy_on_random_matrices():
    g = shapes.rng(4)
    for _ in range(200):
        a = g.normal(size=(3, 3)) * 10 ** g.uniform(-6, 6)
        m = (a + a.T) / 2
        values, vectors = vec.eig_sym3(m.tolist())
        expected = np.linalg.eigvalsh(m)
        scale = max(1.0, np.abs(expected).max())
        assert np.allclose(values, expected, atol=1e-10 * scale)
        assert values == sorted(values)
        for lam, v in zip(values, vectors):
            assert np.allclose(m.dot(v), lam * np.array(v), atol=1e-8 * scale)


def test_eig_sym3_repeated_and_diagonal():
    values, vectors = vec.eig_sym3([[2, 0, 0], [0, 2, 0], [0, 0, 5]])
    assert close(values, (2, 2, 5))
    basis = np.array(vectors)
    assert np.allclose(basis.dot(basis.T), np.eye(3), atol=1e-12)
    values, _ = vec.eig_sym3([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
    assert values == [0.0, 0.0, 0.0]


# ---------------------------------------------------------------- fits

def test_plane_fit_recovers_tilted_noisy_plane():
    g = shapes.rng(5)
    for _ in range(25):
        normal = g.normal(size=3)
        normal /= np.linalg.norm(normal)
        u, v = vec.perpendicular_pair(tuple(normal))
        origin = g.uniform(-100, 100, 3)
        pts = []
        for _ in range(200):
            a, b = g.uniform(-5, 5, 2)
            p = origin + a * np.array(u) + b * np.array(v) + g.normal(0, 0.001) * normal
            pts.append(tuple(p))
        o, n, rms, flatness = vec.plane_fit(pts)
        assert degrees_between_axes(n, normal) < 0.05
        assert abs(rms - 0.001) < 0.0003
        assert flatness < 1e-5
        assert abs(np.dot(np.array(o) - origin, normal)) < 0.001


def test_plane_fit_flags_non_planar_points():
    g = shapes.rng(6)
    blob = [tuple(p) for p in g.normal(size=(500, 3))]
    flatness = vec.plane_fit(blob)[3]
    assert flatness > 0.5
    with pytest.raises(ValueError):
        vec.plane_fit([(0, 0, 0), (1, 0, 0)])


def test_covariance_weights():
    pts = [(0, 0, 0), (2, 0, 0)]
    c, _ = vec.covariance(pts, weights=[3.0, 1.0])
    assert close(c, (0.5, 0, 0))
    with pytest.raises(ValueError):
        vec.covariance(pts, weights=[0.0, 0.0])


def test_least_variance_direction_finds_rod_axis():
    # Side normals of a tilted 24-sided rod are all square to its axis.
    g = shapes.rng(7)
    rot = shapes.random_rotation(g)
    axis = rot.dot((0, 0, 1))
    normals = [tuple(rot.dot((math.cos(a), math.sin(a), 0.0)))
               for a in np.linspace(0, 2 * math.pi, 24, endpoint=False)]
    d, ratio = vec.least_variance_direction(normals)
    assert degrees_between_axes(d, axis) < 1e-6
    assert ratio < 1e-9


def test_least_variance_direction_reports_ball_like_input():
    g = shapes.rng(8)
    normals = [tuple(n / np.linalg.norm(n)) for n in g.normal(size=(2000, 3))]
    _, ratio = vec.least_variance_direction(normals)
    assert ratio > 0.8
