"""Point maps for deforming parts along an axis.

Every map works in a local frame: z runs along the axis, x and y across it,
and the axis passes through x = y = 0. Each map takes a point (x, y, z) and
returns the moved point; meshwright.rhino.deform applies one to Breps and
meshes as a Rhino SpaceMorph.

    Compress  squeeze (or stretch) everything above a start height along the axis
    Curl      roll everything above a start height in towards the axis (or out)
    RadialFit scale radially so one outline lands on another, fading out above a height

dome_height() and dome_levels() give the height of a rounded cap on a flat
outline and the scale and height of its sections.
"""

from __future__ import absolute_import, division

import math

from .report import Report

FIT_BINS = 360  # angular resolution of a radial fit (fine enough for sharp corners to track)


class Compress(object):
    """z above start becomes start + (z - start) * factor; below start nothing moves."""

    def __init__(self, start, factor):
        self.start = float(start)
        self.factor = float(factor)

    def __call__(self, p):
        x, y, z = p
        if z <= self.start:
            return (x, y, z)
        return (x, y, self.start + (z - self.start) * self.factor)


class Curl(object):
    """Roll the zone above start (length long) about a hinge circle of the
    given radius around the axis, turning up to angle radians at the top.

    The turn grows with the square of the height above the hinge (an ease-in:
    no crease at the hinge). A positive angle rolls the top in towards the
    axis. Points on the axis do not move.
    """

    def __init__(self, start, length, radius, angle):
        self.start = float(start)
        self.length = float(length)
        self.radius = float(radius)
        self.angle = float(angle)

    def __call__(self, p):
        x, y, z = p
        s = z - self.start
        if s <= 0.0 or self.length <= 0.0:
            return (x, y, z)
        r = math.sqrt(x * x + y * y)
        if r < 1e-12:
            return (x, y, z)
        t = min(s / self.length, 1.0)
        phi = self.angle * t * t
        a = r - self.radius
        c, si = math.cos(phi), math.sin(phi)
        rho = self.radius + a * c - s * si
        return (x / r * rho, y / r * rho, self.start + a * si + s * c)


class RadialFit(object):
    """Scale each point radially by factors[angle bin] (interpolated between
    bins), fully at and below start, fading smoothly to no change at start + fade."""

    def __init__(self, factors, start, fade):
        self.factors = list(factors)
        self.start = float(start)
        self.fade = float(fade)

    def factor_at(self, angle):
        n = len(self.factors)
        fb = (angle + math.pi) / (2.0 * math.pi) * n - 0.5
        i0 = int(math.floor(fb)) % n
        u = fb - math.floor(fb)
        return self.factors[i0] * (1.0 - u) + self.factors[(i0 + 1) % n] * u

    def __call__(self, p):
        x, y, z = p
        r = math.sqrt(x * x + y * y)
        if r < 1e-12:
            return (x, y, z)
        f = self.factor_at(math.atan2(y, x))
        if z > self.start:
            h = min((z - self.start) / self.fade, 1.0) if self.fade > 0.0 else 1.0
            f = 1.0 + (f - 1.0) * (1.0 - h * h * (3.0 - 2.0 * h))
        return (x * f, y * f, z)


def compress(start, factor):
    """(Compress map or None, report). factor must be more than zero (1: no change)."""
    report = Report("compress")
    if not factor > 0.0:
        report.fail("the factor must be more than zero, got %r" % (factor,))
        return None, report
    report.note("above %g along the axis, lengths times %g" % (start, factor))
    return Compress(start, factor), report


def curl(start, top, angle, radius=None, points=None):
    """(Curl map or None, report) for the zone from start up to top.

    radius is the hinge circle's radius; without one it is the mean distance
    from the axis of the points above start (a list of local points).
    """
    report = Report("curl")
    length = top - start
    if not length > 0.0:
        report.fail("nothing above the start height to curl")
        return None, report
    if radius is None:
        radius = mean_radius(points or [], start)
        if radius is None:
            report.fail("no points above the start height to find the hinge radius from")
            return None, report
        report.note("hinge radius %.4g (mean distance from the axis above the start)" % radius)
    if not radius > 0.0:
        report.fail("the hinge radius must be more than zero")
        return None, report
    report.note("curled %.4g degrees over %.4g" % (math.degrees(angle), length))
    return Curl(start, length, radius, angle), report


def mean_radius(points, start):
    """Mean distance from the axis of the points above start, or None."""
    rs = [math.sqrt(x * x + y * y) for x, y, z in points if z > start]
    return sum(rs) / len(rs) if rs else None


def radial_fit(source, target, start, fade, bins=FIT_BINS):
    """(RadialFit map or None, report) that moves the 2D outline source onto
    target (both closed polygons around the axis, as (x, y) lists).

    Each bin's factor is target radius / source radius along the bin's middle
    ray (the nearest crossing). An outline the axis is outside of has rays
    that miss it; those bins keep factor 1 and are counted.
    """
    report = Report("radial fit")
    if len(source) < 3 or len(target) < 3:
        report.fail("both outlines need at least 3 points")
        return None, report
    if not fade >= 0.0:
        report.fail("the fade length must be zero or more")
        return None, report
    factors, missed = [], 0
    for i in range(bins):
        th = -math.pi + 2.0 * math.pi * (i + 0.5) / bins
        d = (math.cos(th), math.sin(th))
        rs, rt = ray_hit(source, d), ray_hit(target, d)
        if rs is None or rt is None or rs <= 1e-12:
            factors.append(1.0)
            missed += 1
        else:
            factors.append(rt / rs)
    if missed == bins:
        report.fail("the axis is outside both outlines: no ray from it meets them")
        return None, report
    if missed:
        report.messy("%d of %d directions from the axis miss an outline (is the axis inside both?): "
                     "left unscaled there" % (missed, bins))
    report.note("radial factors %.4g to %.4g" % (min(factors), max(factors)))
    return RadialFit(factors, start, fade), report


def ray_hit(polygon, direction):
    """Distance from the origin along direction to the nearest crossing of a closed polygon, or None."""
    dx, dy = direction
    best = None
    n = len(polygon)
    for k in range(n):
        ax, ay = polygon[k]
        bx, by = polygon[(k + 1) % n]
        ex, ey = bx - ax, by - ay
        den = dx * ey - dy * ex
        if abs(den) < 1e-15:
            continue
        t = (ax * ey - ay * ex) / den  # along the ray
        u = (ax * dy - ay * dx) / den  # along the edge
        if t > 0.0 and -1e-12 <= u <= 1.0 + 1e-12 and (best is None or t < best):
            best = t
    return best


def dome_height(area, ratio):
    """Height of a cap over an outline of this area: ratio times the radius of a
    circle with the same area (ratio 1 over a circle gives a half ball)."""
    return ratio * math.sqrt(area / math.pi)


def dome_levels(height, steps=6):
    """[(scale, z)] for the sections of a quarter-ellipse cap over a flat
    outline: the outline scaled about its centre and lifted, from just above
    the base to just below the apex (which is at z = height)."""
    out = []
    for k in range(1, steps + 1):
        t = 0.5 * math.pi * k / (steps + 1.0)
        out.append((math.cos(t), height * math.sin(t)))
    return out
