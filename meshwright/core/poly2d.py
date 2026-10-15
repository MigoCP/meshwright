"""2D outline tools on plain (x, y) tuples.

Polygons are lists of points, counter-clockwise unless stated otherwise, with
no repeated closing point.
"""

from __future__ import absolute_import, division

import math


def cross(o, a, b):
    """z of (a - o) x (b - o): positive when o, a, b turn counter-clockwise."""
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def _dist(a, b):
    return ((b[0] - a[0]) ** 2 + (b[1] - a[1]) ** 2) ** 0.5


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


def perimeter(polygon):
    n = len(polygon)
    return sum(_dist(polygon[i], polygon[(i + 1) % n]) for i in range(n))


def centroid(polygon):
    """Area centroid; the plain vertex average when the polygon has no area."""
    n = len(polygon)
    if n == 0:
        raise ValueError("centroid of an empty polygon")
    a = signed_area(polygon)
    if abs(a) <= 1e-15 * max(1.0, perimeter(polygon) ** 2):
        return (sum(p[0] for p in polygon) / n, sum(p[1] for p in polygon) / n)
    cx = cy = 0.0
    for i in range(n):
        x0, y0 = polygon[i]
        x1, y1 = polygon[(i + 1) % n]
        w = x0 * y1 - x1 * y0
        cx += (x0 + x1) * w
        cy += (y0 + y1) * w
    return (cx / (6.0 * a), cy / (6.0 * a))


def contains_convex(hull, point, eps=0.0):
    """True when point lies inside (or within eps of) a counter-clockwise convex polygon."""
    n = len(hull)
    if n < 3:
        return False
    for i in range(n):
        a, b = hull[i], hull[(i + 1) % n]
        length = _dist(a, b)
        if length == 0.0:
            continue
        if cross(a, b, point) / length < -eps:
            return False
    return True


def distance_to_polygon(point, polygon):
    """Distance from a point to the nearest edge of a closed polygon."""
    best = None
    n = len(polygon)
    for i in range(n):
        a, b = polygon[i], polygon[(i + 1) % n]
        ex, ey = b[0] - a[0], b[1] - a[1]
        l2 = ex * ex + ey * ey
        t = 0.0
        if l2 > 0.0:
            t = max(0.0, min(1.0, ((point[0] - a[0]) * ex + (point[1] - a[1]) * ey) / l2))
        d = _dist(point, (a[0] + ex * t, a[1] + ey * t))
        if best is None or d < best:
            best = d
    return best


def resample(polygon, count):
    """count points evenly spaced along a closed polygon, starting at its first point."""
    n = len(polygon)
    if n == 0 or count <= 0:
        return []
    total = perimeter(polygon)
    if total <= 0.0:
        return [tuple(polygon[0])] * count
    step = total / count
    out = []
    i, walked = 0, 0.0
    seg = _dist(polygon[0], polygon[1 % n])
    for k in range(count):
        target = k * step
        while walked + seg < target and i < n - 1:
            walked += seg
            i += 1
            seg = _dist(polygon[i], polygon[(i + 1) % n])
        a, b = polygon[i], polygon[(i + 1) % n]
        t = min(1.0, max(0.0, (target - walked) / seg)) if seg > 0.0 else 0.0
        out.append((a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t))
    return out


def corners(polygon, max_angle_deg=150.0, span=0.04, samples=256):
    """Sharp corners of a closed outline, ignoring small serrations.

    The outline is resampled evenly, then each point's angle is measured
    between neighbours `span` of the perimeter away on either side: a wide
    stencil, so small steps along an edge do not count. Points sharper than
    max_angle_deg are candidates; the sharpest point of each run of
    candidates is the corner, snapped to the outline's own vertex when one is
    within a sample spacing (so a true polygon's corners come back exactly).
    Returns [(point, angle_deg)] in outline order.
    """
    pts = resample(polygon, samples)
    spacing = perimeter(polygon) / samples
    k = max(1, int(round(span * samples)))
    angles = []
    for i in range(samples):
        p, a, b = pts[i], pts[(i - k) % samples], pts[(i + k) % samples]
        u = (a[0] - p[0], a[1] - p[1])
        v = (b[0] - p[0], b[1] - p[1])
        nu, nv = _dist((0.0, 0.0), u), _dist((0.0, 0.0), v)
        if nu == 0.0 or nv == 0.0:
            angles.append(180.0)
            continue
        c = max(-1.0, min(1.0, (u[0] * v[0] + u[1] * v[1]) / (nu * nv)))
        angles.append(math.degrees(math.acos(c)))
    sharp = [ang <= max_angle_deg for ang in angles]
    if all(sharp) or not any(sharp):
        return []
    start = sharp.index(False)  # walk runs of sharp points from a calm one
    found, run = [], []
    for step in range(1, samples + 1):
        i = (start + step) % samples
        if sharp[i]:
            run.append(i)
        elif run:
            best = min(run, key=lambda j: angles[j])
            p = pts[best]
            nearest = min(polygon, key=lambda q: _dist(q, p))
            if _dist(nearest, p) <= spacing:
                p = (float(nearest[0]), float(nearest[1]))
            found.append((best, p, angles[best]))
            run = []
    found.sort()
    return [(p, ang) for _, p, ang in found]


def min_area_rect(points):
    """Smallest-area rectangle around points: (centre, angle, half_long, half_short).

    angle (radians, in (-pi/2, pi/2]) is the direction of the long side. Every hull
    edge direction is tried, which always includes the optimum.
    """
    hull = convex_hull(points)
    if len(hull) < 3:
        raise ValueError("points span no area")
    best = None
    n = len(hull)
    for i in range(n):
        a, b = hull[i], hull[(i + 1) % n]
        length = _dist(a, b)
        if length == 0.0:
            continue
        ux, uy = (b[0] - a[0]) / length, (b[1] - a[1]) / length
        us = [p[0] * ux + p[1] * uy for p in hull]
        vs = [-p[0] * uy + p[1] * ux for p in hull]
        rect_area = (max(us) - min(us)) * (max(vs) - min(vs))
        if best is None or rect_area < best[0]:
            best = (rect_area, ux, uy, min(us), max(us), min(vs), max(vs))
    _, ux, uy, u0, u1, v0, v1 = best
    cu, cv = (u0 + u1) / 2.0, (v0 + v1) / 2.0
    centre = (cu * ux - cv * uy, cu * uy + cv * ux)
    hu, hv = (u1 - u0) / 2.0, (v1 - v0) / 2.0
    angle = math.atan2(uy, ux)
    if hv > hu:
        hu, hv = hv, hu
        angle += math.pi / 2.0
    return centre, line_angle(angle), hu, hv


def line_angle(angle):
    """A direction angle folded to (-pi/2, pi/2]: lines have no sign."""
    while angle > math.pi / 2.0:
        angle -= math.pi
    while angle <= -math.pi / 2.0:
        angle += math.pi
    return angle


def inertia_ellipse(polygon):
    """The ellipse with the polygon's area centroid and second moments.

    Returns (centre, angle, a, b): a >= b are the semi-axes and angle
    (radians) is the direction of a. For an elliptic outline this is that
    ellipse exactly, wherever it sits.
    """
    if signed_area(polygon) < 0.0:
        polygon = polygon[::-1]
    poly_area = signed_area(polygon)
    if poly_area <= 0.0:
        raise ValueError("polygon has no area")
    cx, cy = centroid(polygon)
    sxx = syy = sxy = 0.0
    n = len(polygon)
    for i in range(n):
        x0, y0 = polygon[i][0] - cx, polygon[i][1] - cy
        x1, y1 = polygon[(i + 1) % n][0] - cx, polygon[(i + 1) % n][1] - cy
        w = x0 * y1 - x1 * y0
        sxx += (x0 * x0 + x0 * x1 + x1 * x1) * w
        syy += (y0 * y0 + y0 * y1 + y1 * y1) * w
        sxy += (x0 * y1 + 2.0 * x0 * y0 + 2.0 * x1 * y1 + x1 * y0) * w
    sxx /= 12.0 * poly_area
    syy /= 12.0 * poly_area
    sxy /= 24.0 * poly_area
    mean = (sxx + syy) / 2.0
    root = (((sxx - syy) / 2.0) ** 2 + sxy * sxy) ** 0.5
    big, small = mean + root, max(mean - root, 0.0)
    angle = 0.5 * math.atan2(2.0 * sxy, sxx - syy)
    # A filled ellipse with semi-axis a has variance a^2 / 4 along it.
    return (cx, cy), angle, 2.0 * big ** 0.5, 2.0 * small ** 0.5


def prune_corners(points, tolerance):
    """Drop corners that barely bend the polygon.

    Repeatedly removes the corner closest to the straight line between its
    two neighbours, while that distance is within tolerance. Noise can fake
    a corner along a straight edge; such a corner sits almost on that line.
    """
    pts = list(points)
    while len(pts) > 3:
        best, at = None, None
        n = len(pts)
        for i in range(n):
            d = distance_to_polygon(pts[i], [pts[i - 1], pts[(i + 1) % n]])
            if best is None or d < best:
                best, at = d, i
        if best > tolerance:
            break
        pts.pop(at)
    return pts


def refine_corners(corners_, outline, margin=0.15):
    """Move each corner to where the lines fitted to its two edges meet.

    outline is the closed outline, evenly sampled, in the same order as the
    corners. Each edge line is a total least-squares fit through the outline
    points between two corners, leaving out the stretch within margin of
    either corner (taken by position along the outline, so a sharp corner
    never pulls in points of the neighbouring edge). On a noisy outline this
    places corners far better than any single outline point can. A corner
    whose edges are near parallel keeps its position.
    """
    n = len(corners_)
    m = len(outline)
    if n < 3 or m < 3 * n:
        return list(corners_)
    at = [min(range(m), key=lambda k: _dist(outline[k], c)) for c in corners_]
    lines = []
    for i in range(n):
        start, stop = at[i], at[(i + 1) % n]
        span = (stop - start) % m
        lo, hi = int(start + margin * span), int(start + (1.0 - margin) * span)
        along = [outline[k % m] for k in range(lo, hi + 1)]
        lines.append(_line_through(along) if len(along) >= 2 else None)
    out = []
    for i in range(n):
        before, after = lines[i - 1], lines[i]
        hit = _meet(before, after) if before and after else None
        out.append(hit if hit is not None else corners_[i])
    return out


def _line_through(points):
    """(centre, unit direction) of the total least-squares line through points."""
    m = len(points)
    cx = sum(p[0] for p in points) / m
    cy = sum(p[1] for p in points) / m
    sxx = sum((p[0] - cx) ** 2 for p in points)
    syy = sum((p[1] - cy) ** 2 for p in points)
    sxy = sum((p[0] - cx) * (p[1] - cy) for p in points)
    angle = 0.5 * math.atan2(2.0 * sxy, sxx - syy)
    return (cx, cy), (math.cos(angle), math.sin(angle))


def _meet(l1, l2):
    (p, d), (q, e) = l1, l2
    det = d[0] * e[1] - d[1] * e[0]
    if abs(det) < 1e-6:
        return None
    t = ((q[0] - p[0]) * e[1] - (q[1] - p[1]) * e[0]) / det
    return (p[0] + d[0] * t, p[1] + d[1] * t)


class Outline(object):
    """Result of classify().

    kind     "circle", "ellipse", "polygon" or "free"
    centre   centre of the fitted shape
    angle    direction of the long axis (radians), 0 for a circle
    a, b     semi-axes of the outline's inertia ellipse (equal for a circle)
    corners  corner points for the polygon kind, else []
    error    how far the outline strays from the fitted shape (90th
             percentile over evenly spaced outline points, so one noisy spike
             does not decide it): relative radius error for curves, distance
             to the corner polygon over the outline's reach for polygons
    """

    def __init__(self, kind, centre, angle, a, b, corners, error):
        self.kind = kind
        self.centre = centre
        self.angle = angle
        self.a = a
        self.b = b
        self.corners = corners
        self.error = error

    def __repr__(self):
        return "Outline(%s, a=%.4g, b=%.4g, corners=%d, error=%.3g)" % (
            self.kind, self.a, self.b, len(self.corners), self.error)


CIRCLE_TOL = 0.02     # circle: ellipse fit good and axes within 2% of each other
ELLIPSE_TOL = 0.03    # ellipse: outline within 3% of the fitted ellipse (90th percentile)
POLYGON_TOL = 0.025   # polygon: outline within 2.5% of the polygon through its corners (90th percentile)


def classify(polygon, circle_tol=CIRCLE_TOL, ellipse_tol=ELLIPSE_TOL, polygon_tol=POLYGON_TOL,
             max_angle_deg=150.0, span=0.04):
    """Name the shape of a closed outline: circle, ellipse, polygon or free.

    The ellipse tried is the outline's own inertia ellipse (same area centroid
    and second moments): no iterative fitting, and it stays centred on an
    off-centre outline. A polygon is accepted when the outline follows the
    straight edges between its detected corners closely.
    """
    centre, angle, a, b = inertia_ellipse(polygon)
    pts = resample(polygon, 256)
    ca, sa = math.cos(angle), math.sin(angle)
    errors = []
    for x, y in pts:
        dx, dy = x - centre[0], y - centre[1]
        u, v = dx * ca + dy * sa, -dx * sa + dy * ca
        r = ((u / a) ** 2 + (v / b) ** 2) ** 0.5 if b > 0.0 else 1e300
        errors.append(abs(r - 1.0))
    worst = _percentile(errors, 0.9)
    if worst <= ellipse_tol:
        if (a - b) / a <= circle_tol:
            radius = (a + b) / 2.0
            return Outline("circle", centre, 0.0, radius, radius, [], worst)
        return Outline("ellipse", centre, angle, a, b, [], worst)
    found = corners(polygon, max_angle_deg, span)
    reach = max(_dist(p, centre) for p in pts)
    poly = prune_corners([p for p, _ in found], 0.5 * polygon_tol * reach)
    if len(poly) >= 3:
        poly = refine_corners(poly, pts)
        err = _percentile([distance_to_polygon(p, poly) for p in pts], 0.9) / reach
        if err <= polygon_tol:
            return Outline("polygon", centroid(poly), angle, a, b, poly, err)
    return Outline("free", centre, angle, a, b, [], worst)


def _percentile(values, q):
    ordered = sorted(values)
    return ordered[min(len(ordered) - 1, int(q * len(ordered)))]
