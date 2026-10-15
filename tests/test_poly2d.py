import math

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
