"""Least-squares fits for measured profiles.

Small, dependency-free fitting tools: a line, a parabola, the break points of
a piecewise-linear profile, an optimal three-segment broken line with slope
limits, and a taper model that extrapolates cross-section areas to the apex
of a tapering end (a point or a ridge).

Inputs are plain sequences of floats. Fits return None (or an object whose
`ok` is False) instead of raising when the data cannot support them.
"""

from __future__ import absolute_import, division

import math


def solve(matrix, rhs):
    """Solve a small dense linear system (Gaussian elimination, partial pivoting).

    Returns the solution list, or None when the system is singular.
    """
    n = len(rhs)
    a = [list(map(float, row)) + [float(rhs[i])] for i, row in enumerate(matrix)]
    scale = max(abs(v) for row in a for v in row[:n]) or 1.0
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(a[r][col]))
        if abs(a[pivot][col]) <= 1e-12 * scale:
            return None
        a[col], a[pivot] = a[pivot], a[col]
        for r in range(col + 1, n):
            f = a[r][col] / a[col][col]
            if f != 0.0:
                for c in range(col, n + 1):
                    a[r][c] -= f * a[col][c]
    x = [0.0] * n
    for r in range(n - 1, -1, -1):
        x[r] = (a[r][n] - sum(a[r][c] * x[c] for c in range(r + 1, n))) / a[r][r]
    return x


class Line(object):
    """y = slope * x + intercept, with the rms of the fit's residuals."""

    def __init__(self, slope, intercept, rms=0.0):
        self.slope = slope
        self.intercept = intercept
        self.rms = rms

    def at(self, x):
        return self.slope * x + self.intercept

    def root(self):
        """x where the line crosses zero, or None for a flat line."""
        if self.slope == 0.0:
            return None
        return -self.intercept / self.slope

    def __repr__(self):
        return "Line(slope=%.6g, intercept=%.6g, rms=%.3g)" % (self.slope, self.intercept, self.rms)


def line_fit(xs, ys, weights=None):
    """Least-squares line y(x). None with fewer than 2 distinct x values."""
    n = len(xs)
    if n < 2 or len(ys) != n:
        return None
    if weights is None:
        weights = [1.0] * n
    sw = float(sum(weights))
    if sw <= 0.0:
        return None
    mx = sum(w * x for w, x in zip(weights, xs)) / sw
    my = sum(w * y for w, y in zip(weights, ys)) / sw
    sxx = sum(w * (x - mx) ** 2 for w, x in zip(weights, xs))
    if sxx <= 1e-24 * max(1.0, mx * mx) * sw:
        return None
    sxy = sum(w * (x - mx) * (y - my) for w, x, y in zip(weights, xs, ys))
    slope = sxy / sxx
    intercept = my - slope * mx
    sse = sum(w * (y - slope * x - intercept) ** 2 for w, x, y in zip(weights, xs, ys))
    return Line(slope, intercept, math.sqrt(sse / sw))


def parabola_fit(xs, ys):
    """Least-squares y = a x^2 + b x + c. Returns (a, b, c, rms) or None.

    x is centred and scaled before solving, which keeps the normal equations
    well conditioned for data far from the origin.
    """
    n = len(xs)
    if n < 3 or len(set(xs)) < 3:
        return None
    mx = sum(xs) / float(n)
    sx = max(abs(x - mx) for x in xs)
    t = [(x - mx) / sx for x in xs]
    s = [float(n), sum(t), sum(v * v for v in t), sum(v ** 3 for v in t), sum(v ** 4 for v in t)]
    rhs = [sum(y * v * v for y, v in zip(ys, t)), sum(y * v for y, v in zip(ys, t)), float(sum(ys))]
    sol = solve([[s[4], s[3], s[2]], [s[3], s[2], s[1]], [s[2], s[1], s[0]]], rhs)
    if sol is None:
        return None
    p, q, r = sol
    # Back to x: y = p ((x - mx)/sx)^2 + q (x - mx)/sx + r
    a = p / (sx * sx)
    b = q / sx - 2.0 * p * mx / (sx * sx)
    c = p * mx * mx / (sx * sx) - q * mx / sx + r
    rms = math.sqrt(sum((y - (a * x * x + b * x + c)) ** 2 for x, y in zip(xs, ys)) / n)
    return a, b, c, rms


# ---------------------------------------------------------------- piecewise-linear profiles

def breakpoints(xs, ys, tolerance, max_segments=8):
    """Where a profile bends: indices splitting it into near-straight pieces.

    Top-down: a piece fits when its least-squares line passes within
    tolerance of every inner point. A piece that does not is split where it strays
    furthest from the chord joining its ends (that is where it bends; the
    least-squares line strays most near the ends instead). Worst piece
    first, until every piece fits or max_segments is reached. Returns sorted
    indices including the first and the last point.
    """
    n = len(xs)
    if n < 3:
        return list(range(n))
    cuts = [0, n - 1]

    def worst(i, j):
        if j - i < 2:
            return (0.0, None)
        line = line_fit(xs[i:j + 1], ys[i:j + 1])
        # Judged on inner points: a piece's ends are shared with its neighbours,
        # and a bend between two samples leaves either end slightly off one line.
        misfit = max(abs(ys[k] - line.at(xs[k])) for k in range(i + 1, j)) if line else 0.0
        if misfit <= tolerance:
            return (misfit, None)
        x0, y0, x1, y1 = xs[i], ys[i], xs[j], ys[j]
        slope = (y1 - y0) / (x1 - x0) if x1 != x0 else 0.0
        at = max(range(i + 1, j), key=lambda k: abs(ys[k] - (y0 + slope * (xs[k] - x0))))
        return (misfit, at)

    while len(cuts) - 1 < max_segments:
        pieces = [(worst(cuts[m], cuts[m + 1]), m) for m in range(len(cuts) - 1)]
        (dev, at), _ = max(pieces, key=lambda item: item[0][0])
        if at is None or dev <= tolerance:
            break
        cuts.append(at)
        cuts.sort()
    return cuts


class BrokenLine(object):
    """Result of three_segment_fit()."""

    def __init__(self, breaks, lines, rms, knees):
        self.breaks = breaks  # (i, j): segments are [0, i), [i, j), [j, n)
        self.lines = lines    # three Line objects
        self.rms = rms
        self.knees = knees    # two (x, y) points where neighbouring lines meet

    def __repr__(self):
        return "BrokenLine(breaks=%s, knees=%s, rms=%.3g)" % (self.breaks, self.knees, self.rms)


def three_segment_fit(xs, ys, slope_limits=None, min_points=2):
    """Best three-piece broken line through a profile, by exhaustive search.

    xs must be increasing. Each piece gets its own least-squares line; with
    slope_limits = [(lo, hi), (lo, hi), (lo, hi)] (None entries for no limit)
    a piece's slope is clamped into its range and its intercept refitted.
    Every split (i, j) is tried; prefix sums make each try O(1), so the whole
    search is O(n^2). The knees are where neighbouring lines meet, kept within
    the span of the two pieces they join. None when there are too few points.
    """
    n = len(xs)
    if n < 3 * min_points:
        return None
    if slope_limits is None:
        slope_limits = [None, None, None]
    # Prefix sums of 1, x, y, xx, xy, yy (centred for precision).
    mx, my = sum(xs) / float(n), sum(ys) / float(n)
    cx = [x - mx for x in xs]
    cy = [y - my for y in ys]
    sums = [[0.0] * 6]
    for x, y in zip(cx, cy):
        last = sums[-1]
        sums.append([last[0] + 1.0, last[1] + x, last[2] + y,
                     last[3] + x * x, last[4] + x * y, last[5] + y * y])

    def piece(i, j, limit):
        a, b = sums[j], sums[i]
        k, sx, sy, sxx, sxy, syy = [a[t] - b[t] for t in range(6)]
        var = sxx - sx * sx / k
        slope = (sxy - sx * sy / k) / var if var > 0.0 else 0.0
        if limit is not None:
            lo, hi = limit
            slope = max(lo, min(hi, slope))
        icpt = (sy - slope * sx) / k
        sse = (syy - 2.0 * slope * sxy - 2.0 * icpt * sy + slope * slope * sxx
               + 2.0 * slope * icpt * sx + k * icpt * icpt)
        return max(sse, 0.0), slope, icpt

    best = None
    for i in range(min_points, n - 2 * min_points + 1):
        e0 = piece(0, i, slope_limits[0])
        for j in range(i + min_points, n - min_points + 1):
            e1 = piece(i, j, slope_limits[1])
            e2 = piece(j, n, slope_limits[2])
            total = e0[0] + e1[0] + e2[0]
            if best is None or total < best[0]:
                best = (total, i, j, e0, e1, e2)
    total, i, j, e0, e1, e2 = best
    # Undo the centring: y - my = m (x - mx) + c  ->  y = m x + (c + my - m mx).
    lines = [Line(e[1], e[2] + my - e[1] * mx) for e in (e0, e1, e2)]
    knees = [_knee(lines[0], lines[1], xs[0], xs[j - 1]), _knee(lines[1], lines[2], xs[i], xs[n - 1])]
    return BrokenLine((i, j), lines, math.sqrt(total / n), knees)


def _knee(a, b, lo, hi):
    """Where two lines meet, kept within [lo, hi]."""
    if abs(a.slope - b.slope) <= 1e-12:
        x = (lo + hi) / 2.0
    else:
        x = max(lo, min(hi, (b.intercept - a.intercept) / (a.slope - b.slope)))
    return (x, (a.at(x) + b.at(x)) / 2.0)


# ---------------------------------------------------------------- tapering ends

class Apex(object):
    """Result of taper_apex().

    ok        False when the sections could not support a fit (see reason)
    height    where the cross-section area reaches zero
    kind      "point" (a cone-like end: sqrt(area) falls linearly) or
              "ridge" (a wedge-like end: area falls linearly)
    misfit    relative misfit of the chosen model and of the other one
    position  (x, y) of the apex when section centres were given, else None
    """

    def __init__(self):
        self.ok = False
        self.reason = ""
        self.height = None
        self.kind = None
        self.misfit = None
        self.other_misfit = None
        self.position = None

    def __repr__(self):
        if not self.ok:
            return "Apex(failed: %s)" % self.reason
        return "Apex(%s at %.6g, misfit=%.3g vs %.3g)" % (
            self.kind, self.height, self.misfit, self.other_misfit)


def taper_apex(heights, areas, centres=None):
    """Extrapolate cross sections of a tapering end to its apex.

    A cone-like end has sections similar about a point, so sqrt(area) falls
    linearly with height; a wedge-like end keeps one dimension while the other
    shrinks, so the area itself falls linearly. Both models are fitted and the
    one with the smaller relative misfit wins. With section centres, a line
    through them gives the apex position at the apex height.
    """
    res = Apex()
    pts = [(h, a) for h, a in zip(heights, areas) if a > 0.0]
    if len(pts) < 3:
        res.reason = "fewer than 3 sections with area"
        return res
    hs = [h for h, _ in pts]
    point_fit = line_fit(hs, [math.sqrt(a) for _, a in pts])
    ridge_fit = line_fit(hs, [a for _, a in pts])
    if point_fit is None or ridge_fit is None or point_fit.slope == 0.0 or ridge_fit.slope == 0.0:
        res.reason = "the sections do not taper"
        return res
    mean_root = sum(math.sqrt(a) for _, a in pts) / len(pts)
    mean_area = sum(a for _, a in pts) / len(pts)
    point_misfit = point_fit.rms / mean_root
    ridge_misfit = ridge_fit.rms / mean_area
    if point_misfit <= ridge_misfit:
        res.kind, chosen, res.misfit, res.other_misfit = "point", point_fit, point_misfit, ridge_misfit
    else:
        res.kind, chosen, res.misfit, res.other_misfit = "ridge", ridge_fit, ridge_misfit, point_misfit
    res.height = chosen.root()
    res.ok = True
    if centres is not None:
        cs = [c for c, a in zip(centres, areas) if a > 0.0]
        fx = line_fit(hs, [c[0] for c in cs])
        fy = line_fit(hs, [c[1] for c in cs])
        if fx is not None and fy is not None:
            res.position = (fx.at(res.height), fy.at(res.height))
    return res
