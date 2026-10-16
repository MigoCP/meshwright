import math

import numpy as np

from solidwright.core import fit

import shapes


def test_solve_matches_numpy_and_flags_singular():
    g = shapes.rng(0)
    for n in (1, 2, 3, 5):
        a = g.normal(size=(n, n)) + n * np.eye(n)
        b = g.normal(size=n)
        assert np.allclose(fit.solve(a.tolist(), b.tolist()), np.linalg.solve(a, b))
    assert fit.solve([[1, 2], [2, 4]], [1, 2]) is None
    assert fit.solve([[0, 1], [1, 0]], [2, 3]) == [3.0, 2.0]   # needs pivoting


def test_line_fit_exact_noisy_and_degenerate():
    xs = [0.0, 1.0, 2.0, 3.0]
    line = fit.line_fit(xs, [1.0 + 2.0 * x for x in xs])
    assert abs(line.slope - 2) < 1e-12 and abs(line.intercept - 1) < 1e-12 and line.rms < 1e-12
    assert abs(line.root() + 0.5) < 1e-12
    g = shapes.rng(1)
    xs = list(np.linspace(1000.0, 1010.0, 200))
    ys = [3.0 - 0.5 * x + g.normal(0, 0.01) for x in xs]
    line = fit.line_fit(xs, ys)
    assert abs(line.slope + 0.5) < 0.002 and abs(line.rms - 0.01) < 0.002
    assert fit.line_fit([1.0], [2.0]) is None
    assert fit.line_fit([2.0, 2.0, 2.0], [1.0, 2.0, 3.0]) is None


def test_weighted_line_fit_ignores_zero_weight_outlier():
    xs = [0.0, 1.0, 2.0, 3.0]
    ys = [0.0, 1.0, 2.0, 100.0]
    line = fit.line_fit(xs, ys, weights=[1, 1, 1, 0])
    assert abs(line.slope - 1.0) < 1e-12


def test_parabola_fit_far_from_origin():
    xs = list(np.linspace(5000.0, 5002.0, 50))
    ys = [0.3 * x * x - 2999.0 * x + 7.0 for x in xs]
    a, b, c, rms = fit.parabola_fit(xs, ys)
    assert abs(a - 0.3) < 1e-6
    assert rms < 1e-6 * max(abs(y) for y in ys)
    assert fit.parabola_fit([1.0, 1.0, 2.0], [1.0, 2.0, 3.0]) is None


def profile(knots, n=120, noise=0.0, seed=0):
    """Samples of a piecewise-linear profile through knots [(x, y), ...]."""
    g = shapes.rng(seed)
    xs = list(np.linspace(knots[0][0], knots[-1][0], n))
    ys = [float(np.interp(x, [k[0] for k in knots], [k[1] for k in knots])) + g.normal(0, noise)
          for x in xs]
    return xs, ys


def test_breakpoints_find_the_bends():
    knots = [(0.0, 0.0), (2.0, 1.0), (3.0, 1.0), (5.0, -2.0)]
    xs, ys = profile(knots, noise=0.002)
    cuts = fit.breakpoints(xs, ys, tolerance=0.01)
    assert len(cuts) == 4
    for cut, (kx, _) in zip(cuts[1:-1], knots[1:-1]):
        assert abs(xs[cut] - kx) < 0.1


def test_breakpoints_stop_at_max_segments_and_on_straight_data():
    xs, ys = profile([(0.0, 0.0), (1.0, 1.0)])
    assert fit.breakpoints(xs, ys, 1e-6) == [0, len(xs) - 1]
    xs, ys = profile([(0, 0), (1, 1), (2, 0), (3, 1), (4, 0)])
    assert len(fit.breakpoints(xs, ys, 1e-6, max_segments=2)) == 3


def test_three_segment_fit_recovers_knees():
    # Rising, nearly flat, falling: the classic waist band of a side profile.
    knots = [(-0.4, 0.0), (0.0, 0.5), (0.03, 0.5), (0.2, 0.3)]
    xs, ys = profile(knots, n=150, noise=0.001, seed=2)
    result = fit.three_segment_fit(xs, ys)
    (x1, y1), (x2, y2) = result.knees
    assert abs(x1 - 0.0) < 0.01 and abs(x2 - 0.03) < 0.01
    assert abs(y1 - 0.5) < 0.005 and abs(y2 - 0.5) < 0.005
    assert result.rms < 0.002


def test_three_segment_fit_respects_slope_limits():
    knots = [(0.0, 0.0), (1.0, 1.0), (2.0, 1.2), (3.0, 0.0)]
    xs, ys = profile(knots, n=90)
    free = fit.three_segment_fit(xs, ys)
    assert abs(free.lines[1].slope - 0.2) < 0.01
    flat = fit.three_segment_fit(xs, ys, slope_limits=[None, (-0.05, 0.05), None])
    assert abs(flat.lines[1].slope) <= 0.05 + 1e-12
    assert flat.rms > free.rms
    assert fit.three_segment_fit([0.0, 1.0, 2.0], [0.0, 1.0, 2.0]) is None


def test_taper_apex_point_end():
    # A cone: radius shrinks linearly to zero at height -2.
    hs = list(np.linspace(-1.5, -0.5, 12))
    areas = [math.pi * (0.4 * (h + 2.0)) ** 2 for h in hs]
    centres = [(1.0 + 0.1 * h, -2.0) for h in hs]
    apex = fit.taper_apex(hs, areas, centres)
    assert apex.ok and apex.kind == "point"
    assert abs(apex.height + 2.0) < 1e-9
    assert np.allclose(apex.position, (0.8, -2.0))


def test_taper_apex_ridge_end():
    # A wedge: constant length 3, width shrinking linearly to zero at height 1.
    hs = list(np.linspace(0.2, 0.8, 10))
    areas = [3.0 * 0.5 * (1.0 - h) for h in hs]
    apex = fit.taper_apex(hs, areas)
    assert apex.ok and apex.kind == "ridge" and abs(apex.height - 1.0) < 1e-9


def test_taper_apex_with_noise_and_failures():
    g = shapes.rng(3)
    hs = list(np.linspace(-1.5, -0.5, 20))
    areas = [math.pi * (0.4 * (h + 2.0)) ** 2 * (1 + g.normal(0, 0.01)) for h in hs]
    apex = fit.taper_apex(hs, areas)
    assert apex.kind == "point" and abs(apex.height + 2.0) < 0.05
    assert not fit.taper_apex([0.0, 1.0], [1.0, 2.0]).ok
    flat = fit.taper_apex([0.0, 1.0, 2.0], [1.0, 1.0, 1.0])
    assert not flat.ok and "taper" in flat.reason
