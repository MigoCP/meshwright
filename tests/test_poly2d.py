import math

import pytest

import shapes

from meshwright.core import poly2d


def test_hull_drops_inner_and_collinear_points():
    pts = [(0, 0), (2, 0), (1, 0), (2, 2), (0, 2), (1, 1), (1, 2)]
    hull = poly2d.convex_hull(pts)
    assert sorted(hull) == [(0.0, 0.0), (0.0, 2.0), (2.0, 0.0), (2.0, 2.0)]
    assert poly2d.signed_area(hull) > 0                      # counter-clockwise


def test_hull_of_degenerate_input():
    assert poly2d.convex_hull([(1, 1), (1, 1)]) == [(1.0, 1.0)]
    assert len(poly2d.convex_hull([(0, 0), (1, 1), (2, 2)])) < 3


def test_areas():
    square = [(0, 0), (2, 0), (2, 2), (0, 2)]
    assert poly2d.signed_area(square) == 4.0
    assert poly2d.signed_area(square[::-1]) == -4.0
    assert poly2d.area(square[::-1]) == 4.0
    circle = [(math.cos(a), math.sin(a)) for a in [2 * math.pi * k / 400 for k in range(400)]]
    assert abs(poly2d.area(circle) - math.pi) < 1e-3


def test_contains_convex():
    square = [(0, 0), (2, 0), (2, 2), (0, 2)]
    assert poly2d.contains_convex(square, (1, 1))
    assert poly2d.contains_convex(square, (2, 1))           # on the edge
    assert not poly2d.contains_convex(square, (2.01, 1))
    assert poly2d.contains_convex(square, (2.01, 1), eps=0.02)
    assert not poly2d.contains_convex([(0, 0), (1, 1)], (0.5, 0.5))


def ngon(n, rx=1.0, ry=None, centre=(0.0, 0.0), rotation=0.0, phase=0.0):
    ry = rx if ry is None else ry
    out = []
    for k in range(n):
        a = phase + 2 * math.pi * k / n
        x, y = rx * math.cos(a), ry * math.sin(a)
        c, s = math.cos(rotation), math.sin(rotation)
        out.append((centre[0] + x * c - y * s, centre[1] + x * s + y * c))
    return out


def close(a, b, tol=1e-9):
    return all(abs(x - y) <= tol for x, y in zip(a, b))


def test_centroid_and_perimeter():
    rect = [(1, 1), (5, 1), (5, 3), (1, 3)]
    assert close(poly2d.centroid(rect), (3.0, 2.0))
    assert poly2d.perimeter(rect) == 12.0
    l_shape = [(0, 0), (2, 0), (2, 1), (1, 1), (1, 2), (0, 2)]
    assert close(poly2d.centroid(l_shape), (5.0 / 6.0, 5.0 / 6.0))
    assert close(poly2d.centroid([(0, 0), (2, 0), (4, 0)]), (2.0, 0.0))  # no area: vertex average


def test_resample_is_even_and_starts_at_the_first_point():
    rect = [(0, 0), (4, 0), (4, 2), (0, 2)]
    pts = poly2d.resample(rect, 12)
    assert len(pts) == 12 and pts[0] == (0.0, 0.0)
    gaps = [math.hypot(pts[(i + 1) % 12][0] - pts[i][0], pts[(i + 1) % 12][1] - pts[i][1])
            for i in range(12)]
    assert max(gaps) - min(gaps) < 1e-9 or all(g <= 1.0 + 1e-9 for g in gaps)
    assert (4.0, 1.0) in [tuple(round(c, 9) for c in p) for p in pts]


def test_corners_of_polygons_are_exact_and_round_shapes_have_none():
    tri = ngon(3, centre=(3.0, -1.0), rotation=0.4)
    found = poly2d.corners(tri)
    assert len(found) == 3
    assert sorted(p for p, _ in found) == sorted(tri)
    assert all(abs(ang - 60.0) < 5.0 for _, ang in found)
    assert poly2d.corners(ngon(200)) == []


def test_corners_ignore_small_serrations():
    # A square whose edges carry tiny steps: still four corners.
    pts = []
    for (x0, y0), (x1, y1) in zip([(0, 0), (4, 0), (4, 4), (0, 4)], [(4, 0), (4, 4), (0, 4), (0, 0)]):
        for k in range(40):
            t = k / 40.0
            bump = 0.01 if k % 2 else 0.0
            nx, ny = (y1 - y0) / 4.0, -(x1 - x0) / 4.0
            pts.append((x0 + (x1 - x0) * t + nx * bump, y0 + (y1 - y0) * t + ny * bump))
    assert len(poly2d.corners(pts)) == 4


def test_min_area_rect_of_a_rotated_offset_rectangle():
    rect = [(0, 0), (3, 0), (3, 1), (0, 1)]
    c, s = math.cos(0.6), math.sin(0.6)
    moved = [(10 + x * c - y * s, -5 + x * s + y * c) for x, y in rect]
    centre, angle, half_long, half_short = poly2d.min_area_rect(moved + [(10.5, -4.0)])
    assert close(centre, (10 + 1.5 * c - 0.5 * s, -5 + 1.5 * s + 0.5 * c))
    assert abs(angle - 0.6) < 1e-9
    assert abs(half_long - 1.5) < 1e-9 and abs(half_short - 0.5) < 1e-9


def test_inertia_ellipse_is_exact_for_ellipses_anywhere():
    poly = ngon(720, rx=2.0, ry=0.7, centre=(10.0, 4.0), rotation=0.7)
    centre, angle, a, b = poly2d.inertia_ellipse(poly)
    assert close(centre, (10.0, 4.0), 1e-9)
    assert abs(angle - 0.7) < 1e-6 and abs(a - 2.0) < 1e-3 and abs(b - 0.7) < 1e-3
    assert close(poly2d.inertia_ellipse(poly[::-1])[0], centre)        # clockwise input too


def test_classify():
    cases = [
        (ngon(200), "circle"),
        (ngon(200, centre=(5, -3)), "circle"),
        (ngon(16), "circle"),                                  # within 2%
        (ngon(8), "polygon"),
        (ngon(200, rx=2, ry=1, centre=(10, 4), rotation=0.7), "ellipse"),
        (ngon(4, phase=math.pi / 4), "polygon"),
        ([(0, 0), (3, 0), (3, 1), (0, 1)], "polygon"),
        (ngon(3), "polygon"),
        ([(math.cos(a) * (1 + 0.4 * math.sin(a)), 0.8 * math.sin(a))
          for a in [2 * math.pi * k / 200 for k in range(200)]], "free"),
    ]
    for poly, kind in cases:
        assert poly2d.classify(poly).kind == kind, (kind, poly2d.classify(poly))
    off = poly2d.classify(ngon(200, rx=2, ry=1, centre=(10, 4), rotation=0.7))
    assert close(off.centre, (10.0, 4.0), 1e-6) and abs(off.angle - 0.7) < 1e-3
    square = poly2d.classify(ngon(4, centre=(2, 2), phase=math.pi / 4))
    assert len(square.corners) == 4 and close(square.centre, (2.0, 2.0))


def test_line_angle_folds_directions():
    assert abs(poly2d.line_angle(math.pi) - 0.0) < 1e-12
    assert abs(poly2d.line_angle(-math.pi / 2) - math.pi / 2) < 1e-12
    assert abs(poly2d.line_angle(2.0) - (2.0 - math.pi)) < 1e-12


@pytest.mark.parametrize("base,kind", [
    (ngon(300), "circle"),
    (ngon(300, rx=2.0, ry=1.0, rotation=1.0), "ellipse"),
    (ngon(3), "polygon"),
    (ngon(6), "polygon"),
    (ngon(4, phase=math.pi / 4), "polygon"),
    ([(0, 0), (3, 0), (3, 1), (0, 1)], "polygon"),
])
def test_classify_survives_scan_noise(base, kind):
    # Noise of 1% of the size on every point: one stray point must not decide.
    g = shapes.rng(21)
    for _ in range(10):
        pts = [(x + g.normal(0, 0.01), y + g.normal(0, 0.01)) for x, y in poly2d.resample(base, 300)]
        assert poly2d.classify(pts).kind == kind


def test_noisy_polygon_corners_land_on_the_true_corners():
    g = shapes.rng(22)
    tri = ngon(3, centre=(4.0, 2.0), rotation=0.3)
    pts = [(x + g.normal(0, 0.01), y + g.normal(0, 0.01)) for x, y in poly2d.resample(tri, 300)]
    found = poly2d.classify(pts).corners
    assert len(found) == 3
    for c in found:
        assert min(math.hypot(c[0] - t[0], c[1] - t[1]) for t in tri) < 0.02


def test_prune_corners_drops_points_on_straight_edges():
    square = [(0, 0), (1, 0.001), (2, 0), (2, 2), (0, 2)]
    assert poly2d.prune_corners(square, 0.01) == [(0, 0), (2, 0), (2, 2), (0, 2)]
    assert poly2d.prune_corners(square, 0.0001) == square
