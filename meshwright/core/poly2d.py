"""2D outline tools on plain (x, y) tuples.

Polygons are lists of points, counter-clockwise unless stated otherwise, with
no repeated closing point.
"""

from __future__ import absolute_import, division


def cross(o, a, b):
    """z of (a - o) x (b - o): positive when o, a, b turn counter-clockwise."""
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def convex_hull(points):
    """Convex hull, counter-clockwise, collinear points dropped (monotone chain).

    Returns fewer than 3 points when the input is degenerate (all equal or on
    one line).
    """
    pts = sorted(set((float(p[0]), float(p[1])) for p in points))
    if len(pts) < 3:
        return pts
    lower = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0.0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0.0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


def signed_area(polygon):
    """Shoelace area: positive for counter-clockwise polygons."""
    total = 0.0
    n = len(polygon)
    for i in range(n):
        x0, y0 = polygon[i]
        x1, y1 = polygon[(i + 1) % n]
        total += x0 * y1 - x1 * y0
    return 0.5 * total


def area(polygon):
    return abs(signed_area(polygon))


def contains_convex(hull, point, eps=0.0):
    """True when point lies inside (or within eps of) a counter-clockwise convex polygon."""
    n = len(hull)
    if n < 3:
        return False
    for i in range(n):
        a, b = hull[i], hull[(i + 1) % n]
        ex, ey = b[0] - a[0], b[1] - a[1]
        length = (ex * ex + ey * ey) ** 0.5
        if length == 0.0:
            continue
        if cross(a, b, point) / length < -eps:
            return False
    return True
