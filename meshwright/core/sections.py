"""A solid described by cross sections, with sideways clearance.

Works in a local frame whose z runs along the part's axis. Every section of
the solid is a convex polygon given by its support offsets along a fixed set
of 2D directions n_i: the region {p : n_i . (p - c) <= h_i}. For a measured
section, h_i is how far the part's section reaches along n_i, plus the
clearance; adding the clearance to every offset moves every edge out by
exactly that much (a radius would under-offset the flanks of a long part).

Between two levels each offset is interpolated linearly. That is exact across
any flat side facet: in cross section a facet is a line whose offset slides
linearly with height. So the model is exact once there is a level at every
height where some side of the part bends, and the solid made from it (loops
of section corners joined by flat faces) is exactly the model.

Levels
    1. Widest band: a scan of section areas; the run of heights where the
       area stays within BAND_RATIO of its largest value.
    2. Bends: where sqrt(area) stops being straight, above and below the band
       (fit.breakpoints). Straight runs of sqrt(area) mean similar sections.
    3. Ends: each end is measured END_MARGIN inside and then extrapolated. A
       thin end tapers to an apex, a point or a ridge (pinch()); a thick end
       is flat and capped at the part's extent. Extrapolation backs off
       whenever a section would turn inside out: stopping early only extends
       the measured section straight, which makes the solid bigger, never
       smaller.
    4. Refinement: any part vertex that comes out closer than KEEP_CLEARANCE
       of the clearance asks for a measured level at its own height (the
       worst one per span), until none are left. This is what catches sides
       that bend at heights the overall area averages away.

The clearance is applied square to the axis only: the flat ends of the
solid sit at the part's own extent.
"""

from __future__ import absolute_import, division

import bisect
import math

from . import fit, poly2d, vec
from .polyhedron import FLAT_SHARE, Polyhedron, merge_coplanar, merge_corners, validate
from .report import Report

DIRECTIONS = 96        # support directions for outlines that are not polygons
SLICES = 200           # sections in the area scan
BAND_RATIO = 0.995     # widest band: area within this share of the largest
MIN_BAND = 0.006       # ...and at least this share of the height thick
BAND_STEPS = 9         # sections across the band, merged into its outline
END_MARGIN = 0.02      # ends are measured this share of the height inside
BEND_TOL = 0.004       # a bend: sqrt(area) strays this x the reach from a straight run
MAX_BENDS = 4          # bends looked for on each side of the band
TAPER_RATIO = 0.15     # an end thinner than this x the reach tapers to an apex
POINT_RATIO = 0.05     # an apex section smaller than this x the reach is a point
MAX_OVERSHOOT = 0.35   # an apex may lie at most this x the height past the end
STUB = 0.01            # a ridge apex is carried this x the height past the end
SAFE_T = 3.0           # extrapolate at most this many times the baseline
KEEP_CLEARANCE = 0.9   # refinement target: this share of the clearance
REFINE_PASSES = 8
MAX_CHECK = 5000       # vertices used for refinement and checks (evenly thinned)


# ---------------------------------------------------------------- slicing

def slice_segments(points, faces, heights, axis=2):
    """Where a triangle mesh crosses the planes coordinate[axis] = h, per height.

    heights must be sorted ascending. Returns one list per height of
    segments (pairs of points), one segment for each triangle crossing that
    plane. Each triangle is only tested against the heights within its own
    span, so many sections of a large mesh stay cheap.
    """
    out = [[] for _ in heights]
    for f in faces:
        a, b, c = points[f[0]], points[f[1]], points[f[2]]
        za, zb, zc = a[axis], b[axis], c[axis]
        lo = min(za, zb, zc)
        hi = max(za, zb, zc)
        if hi <= lo:
            continue  # lies in a plane square to the axis: its rim is found by its neighbours
        for k in range(bisect.bisect_left(heights, lo), bisect.bisect_right(heights, hi)):
            z = heights[k]
            ends = []
            for p, q, zp, zq in ((a, b, za, zb), (b, c, zb, zc), (c, a, zc, za)):
                if zp == z:
                    ends.append(p)
                elif (zp - z) * (zq - z) < 0.0:
                    ends.append(vec.lerp(p, q, (z - zp) / (zq - zp)))
            if len(ends) >= 2:
                out[k].append((ends[0], ends[1]))
    return out


def section_hulls(points, faces, heights):
    """Convex hull (2D, x and y) of the mesh section at each height."""
    order = sorted(range(len(heights)), key=lambda k: heights[k])
    sorted_h = [heights[k] for k in order]
    segs = slice_segments(points, faces, sorted_h)
    hulls = [None] * len(heights)
    for k, s in zip(order, segs):
        pts = [(p[0], p[1]) for seg in s for p in seg]
        hull = poly2d.convex_hull(pts) if len(pts) >= 3 else []
        hulls[k] = hull if len(hull) >= 3 else None
    return hulls


# ---------------------------------------------------------------- the model

class SectionModel(object):
    """Levels of support offsets; see the module docs.

    directions  unit 2D directions, counter-clockwise
    centre      the 2D point offsets are measured from
    levels      [(z, offsets)], bottom to top
    tips        {-1: bottom point or None, +1: top point or None}, local (x, y, z)
    """

    def __init__(self, directions, centre, clearance):
        self.directions = directions
        self.centre = centre
        self.clearance = clearance
        self.levels = []
        self.tips = {-1: None, 1: None}
        n = len(directions)
        self.opposite = []
        for i in range(n):
            ax, ay = directions[i]
            self.opposite.append(min(range(n), key=lambda j: ax * directions[j][0] + ay * directions[j][1]))

    def offsets(self, polygon):
        """Support offsets of a 2D point set, clearance included."""
        cx, cy = self.centre
        return [max((x - cx) * nx + (y - cy) * ny for x, y in polygon) + self.clearance
                for nx, ny in self.directions]

    def corners(self, offsets):
        """Corners of the section {p : n_i . (p - c) <= h_i}: where each support
        line meets the next. None when two neighbouring directions are parallel."""
        cx, cy = self.centre
        n = len(self.directions)
        out = []
        for i in range(n):
            ax, ay = self.directions[i]
            bx, by = self.directions[(i + 1) % n]
            det = ax * by - ay * bx
            if abs(det) < 1e-12:
                return None
            h1, h2 = offsets[i], offsets[(i + 1) % n]
            out.append((cx + (h1 * by - h2 * ay) / det, cy + (h2 * ax - h1 * bx) / det))
        return out

    def clean(self, offsets, slack):
        """True when every support line still bounds the section in order (no
        edge running backwards by more than slack)."""
        pts = self.corners(offsets)
        if pts is None or poly2d.signed_area(pts) <= 0.0:
            return False
        for i in range(len(pts)):
            nx, ny = self.directions[i]
            ax, ay = pts[i - 1]
            bx, by = pts[i]
            if (bx - ax) * -ny + (by - ay) * nx < -slack:
                return False
        return True

    def bottom(self):
        tip = self.tips[-1]
        return tip[2] if tip is not None else self.levels[0][0]

    def top(self):
        tip = self.tips[1]
        return tip[2] if tip is not None else self.levels[-1][0]

    def offsets_at(self, z):
        """Offsets of the solid's section at height z, or None outside it."""
        lv = self.levels
        if z < lv[0][0]:
            return self._to_tip(-1, z)
        if z > lv[-1][0]:
            return self._to_tip(1, z)
        k = bisect.bisect_right([l[0] for l in lv], z) - 1
        k = max(0, min(k, len(lv) - 2))
        z0, h0 = lv[k]
        z1, h1 = lv[k + 1]
        t = (z - z0) / (z1 - z0) if z1 > z0 else 0.0
        return [a + (b - a) * t for a, b in zip(h0, h1)]

    def _to_tip(self, side, z):
        tip = self.tips[side]
        if tip is None:
            return None
        end_z, end_h = self.levels[0] if side < 0 else self.levels[-1]
        if (z - tip[2]) * (end_z - tip[2]) <= 0.0:
            return None
        cx, cy = self.centre
        at_tip = [(tip[0] - cx) * nx + (tip[1] - cy) * ny for nx, ny in self.directions]
        t = (z - tip[2]) / (end_z - tip[2])
        return [a + (b - a) * t for a, b in zip(at_tip, end_h)]

    def slack(self, p):
        """Sideways clearance of a local point (x, y, z) inside the solid:
        negative when it lies outside. None when it is beyond an end."""
        hs = self.offsets_at(p[2])
        if hs is None:
            return None
        dx, dy = p[0] - self.centre[0], p[1] - self.centre[1]
        return min(h - (dx * nx + dy * ny) for h, (nx, ny) in zip(hs, self.directions))

    def polyhedron(self, to_world):
        """The solid: loops of section corners joined by flat faces, in world
        coordinates (to_world maps a local point)."""
        vertices, faces = [], []
        loops = []
        for z, hs in self.levels:
            loop = []
            for x, y in self.corners(hs):
                loop.append(len(vertices))
                vertices.append(to_world((x, y, z)))
            loops.append(loop)
        n = len(self.directions)
        for k in range(len(loops) - 1):
            a, b = loops[k], loops[k + 1]
            for j in range(n):
                j2 = (j + 1) % n
                faces.append([a[j], a[j2], b[j2], b[j]])
        for side, loop in ((-1, loops[0]), (1, loops[-1])):
            tip = self.tips[side]
            if tip is None:
                faces.append(list(reversed(loop)) if side < 0 else list(loop))
                continue
            t = len(vertices)
            vertices.append(to_world(tip))
            for j in range(n):
                j2 = (j + 1) % n
                faces.append([t, loop[j2], loop[j]] if side < 0 else [loop[j], loop[j2], t])
        return Polyhedron(vertices, faces)


def fan(count, start=0.0):
    """count unit directions evenly around the circle, counter-clockwise from
    the angle start (radians)."""
    return [(math.cos(start + 2.0 * math.pi * i / count), math.sin(start + 2.0 * math.pi * i / count))
            for i in range(count)]


def edge_normals(polygon):
    """Outward unit normals of a counter-clockwise polygon's edges."""
    if poly2d.signed_area(polygon) < 0.0:
        polygon = list(reversed(polygon))
    out = []
    n = len(polygon)
    for i in range(n):
        x1, y1 = polygon[i]
        x2, y2 = polygon[(i + 1) % n]
        length = math.hypot(x2 - x1, y2 - y1)
        if length > 1e-12:
            out.append(((y2 - y1) / length, -(x2 - x1) / length))
    return out


def pinch(model, h0, h1, floor):
    """The smallest t > 1 at which the section blended from h0 (t = 0) to h1
    (t = 1) thins to floor in some direction, or None if it never does
    within 8 baselines.

    A section's width along a direction is the sum of two opposite offsets,
    so it is linear in t. When no width reaches the floor there is no honest
    answer, and blending on to some cap would only turn the section inside out.
    """
    best = None
    for i, j in enumerate(model.opposite):
        w0 = h0[i] + h0[j]
        w1 = h1[i] + h1[j]
        slope = w1 - w0
        if slope < -1e-12:
            t = (floor - w0) / slope
            if 1.0 < t <= 8.0 and (best is None or t < best):
                best = t
    return best


def safe_blend(model, h0, h1, t, floor, slack):
    """Extrapolate to t, backing off towards t = 1 while the section is too thin
    somewhere or turns inside out. Returns (offsets, t used)."""
    while t > 1.0 + 1e-3:
        hs = [a + (b - a) * t for a, b in zip(h0, h1)]
        thin = any(hs[i] + hs[j] < floor for i, j in enumerate(model.opposite))
        if not thin and model.clean(hs, slack):
            return hs, t
        t = 1.0 + (t - 1.0) * 0.5
    return list(h1), 1.0


# ---------------------------------------------------------------- building the model

class ClearanceResult(object):
    """solid: Polyhedron (world coordinates) or None. model: the SectionModel
    (local frame). frame: (x, y, z) unit axes, origin: the local frame's origin.
    band: (bottom, top) heights of the widest band, local. min_clearance:
    smallest sideways clearance of the part's vertices."""

    def __init__(self):
        self.solid = None
        self.model = None
        self.frame = None
        self.origin = None
        self.band = None
        self.min_clearance = None
        self.added = 0
        self.report = Report("clearance solid")

    @property
    def status(self):
        return self.report.status

    def to_world(self, p):
        x, y, z = self.frame
        o = self.origin
        return (o[0] + x[0] * p[0] + y[0] * p[1] + z[0] * p[2],
                o[1] + x[1] * p[0] + y[1] * p[1] + z[1] * p[2],
                o[2] + x[2] * p[0] + y[2] * p[1] + z[2] * p[2])


def clearance_solid(vertices, faces, axis, clearance=0.0, tolerance=0.001, directions=DIRECTIONS,
                    slices=SLICES, max_bends=MAX_BENDS, x_hint=(1.0, 0.0, 0.0)):
    """A closed solid around a part, its sections kept clearance away sideways.

    vertices, faces  the part (a cleaned, closed triangle mesh)
    axis             the part's axis: sections are square to it
    clearance        sideways clearance (0: the solid touches the part)
    tolerance        modelling tolerance: corners closer than this are merged
    directions       support directions for outlines that are not polygons
                     (a polygonal widest band uses its own edge directions)
    Returns a ClearanceResult. See the module docs for how it is built.
    """
    res = ClearanceResult()
    r = res.report
    if clearance < 0.0:
        r.fail("the clearance must be zero or more, got %.4g" % clearance)
        return res
    if tolerance <= 0.0:
        r.fail("the tolerance must be more than zero, got %.4g" % tolerance)
        return res
    if not faces:
        r.fail("the mesh has no faces")
        return res
    z_axis = vec.unit(axis)
    if z_axis is None:
        r.fail("the axis is a zero vector")
        return res
    used = sorted(set(i for f in faces for i in f))
    res.frame = vec.frame_from_axis(z_axis, x_hint)
    res.origin = vec.mean([vertices[i] for i in used])
    fx, fy, fz = res.frame
    local = {}
    for i in used:
        d = vec.sub(vertices[i], res.origin)
        local[i] = (vec.dot(d, fx), vec.dot(d, fy), vec.dot(d, fz))
    points = [local.get(i, (0.0, 0.0, 0.0)) for i in range(len(vertices))]
    zs = [local[i][2] for i in used]
    zmin, zmax = min(zs), max(zs)
    height = zmax - zmin
    if height <= 0.0:
        r.fail("the part is flat along the axis")
        return res

    # 1. Widest band from an area scan.
    scan_z = [zmin + (k + 0.5) / slices * height for k in range(slices)]
    hulls = section_hulls(points, faces, scan_z)
    areas = [poly2d.area(h) if h else 0.0 for h in hulls]
    top_k = max(range(slices), key=lambda k: areas[k])
    if areas[top_k] <= 0.0:
        r.fail("no cross section of the part has any area")
        return res
    lo = hi = top_k
    while lo > 0 and areas[lo - 1] >= BAND_RATIO * areas[top_k]:
        lo -= 1
    while hi < slices - 1 and areas[hi + 1] >= BAND_RATIO * areas[top_k]:
        hi += 1
    z_lo, z_hi = scan_z[lo], scan_z[hi]
    if z_hi - z_lo < MIN_BAND * height:
        mid = 0.5 * (z_lo + z_hi)
        z_lo, z_hi = mid - 0.5 * MIN_BAND * height, mid + 0.5 * MIN_BAND * height
    z_lo, z_hi = max(z_lo, zmin), min(z_hi, zmax)
    band_z = [z_lo + (z_hi - z_lo) * k / (BAND_STEPS - 1.0) for k in range(BAND_STEPS)]
    band_pts = [p for h in section_hulls(points, faces, band_z) if h for p in h]
    band = poly2d.convex_hull(band_pts)
    res.band = (z_lo, z_hi)
    centre = poly2d.centroid(band)
    reach = max(math.hypot(x - centre[0], y - centre[1]) for x, y in band)

    # Directions: a polygonal band keeps its own edges, so corners stay sharp.
    shape = poly2d.classify(band)
    dirs = None
    if shape.kind == "polygon" and 3 <= len(shape.corners) <= 16:
        dirs = edge_normals(shape.corners)
        r.note("the widest band is a %d-sided polygon; sections follow its edges" % len(dirs))
    if not dirs or len(dirs) < 3:
        # Start the fan at the band's furthest point, so the solid does not
        # depend on how the part happens to sit in the frame.
        far = max(band, key=lambda p: (p[0] - centre[0]) ** 2 + (p[1] - centre[1]) ** 2)
        dirs = fan(max(8, int(directions)), math.atan2(far[1] - centre[1], far[0] - centre[0]))
    model = SectionModel(dirs, centre, clearance)
    res.model = model
    slack_tol = 2e-4 * reach

    def measure(z):
        e = 1e-5 * height
        pts = [p for h in section_hulls(points, faces, [z - e, z + e]) if h for p in h]
        return model.offsets(pts) if len(pts) >= 3 else None

    band_h = model.offsets(band)
    root = [math.sqrt(a) for a in areas]

    def bends(z_from, z_to):
        """Heights strictly between z_from and z_to where sqrt(area) bends."""
        a_, b_ = min(z_from, z_to), max(z_from, z_to)
        idx = [k for k in range(slices) if a_ <= scan_z[k] <= b_ and areas[k] > 0.0]
        if len(idx) < 5 or max_bends < 1:
            return []
        xs = [scan_z[k] for k in idx]
        ys = [root[k] for k in idx]
        cuts = fit.breakpoints(xs, ys, BEND_TOL * reach, max_bends + 1)
        return [xs[c] for c in cuts[1:-1]]

    # 2 and 3. Each side of the band: bends, then the end.
    levels = [(z_lo, band_h), (z_hi, band_h)]
    for side, band_edge, end in ((1, z_hi, zmax), (-1, z_lo, zmin)):
        z_end_in = end - side * END_MARGIN * height
        if (z_end_in - band_edge) * side <= 0.0:
            levels.append((end, list(band_h)))  # the widest band reaches this end: a flat cap
            continue
        heights = sorted(bends(band_edge, z_end_in), key=lambda z: side * z) + [z_end_in]
        measured = []
        for z in heights:
            hs = measure(z)
            measured.append((z, hs if hs is not None else list(band_h)))
        levels.extend(measured)
        _close_end(model, levels, measured, band_edge, band_h, end, side, height, reach, slack_tol, r)

    levels.sort(key=lambda lv: lv[0])
    model.levels = _dedup(levels, max(1e-4 * height, 2.0 * tolerance))

    # 4. Refinement where vertices come out short.
    check = [points[i] for i in used[::max(1, int(math.ceil(len(used) / float(MAX_CHECK))))]]
    target = KEEP_CLEARANCE * clearance if clearance > 0.0 else -tolerance
    res.added = _refine(model, check, target, measure, res.band, height, tolerance)

    model.levels = _prune(model.levels, 1e-9 * max(height, reach))

    # The solid, merged within the tolerance where that keeps it closed and flat.
    # Corners that coincide to rounding (a support line touching the section at
    # one point on two levels) are always joined: that only removes faces with no area.
    exact = merge_corners(model.polyhedron(res.to_world), 1e-9 * max(height, reach))
    merged = merge_corners(exact, tolerance)
    problems, _ = validate(merged, FLAT_SHARE * tolerance)
    if not problems:
        res.solid = merged
    else:
        problems, _ = validate(exact, FLAT_SHARE * tolerance)
        if problems:
            r.fail("the solid fails its own checks: " + "; ".join(problems))
            return res
        res.solid = exact
        r.warn("short edges kept: merging corners within %.3g would break the solid" % tolerance)
    # A straight side spanning several sections is one face.
    res.solid = merge_coplanar(res.solid, 1e-12)
    _grade(res, model, [points[i] for i in used], clearance, tolerance, reach)
    r.note("%d sections (%d added to follow the profile), %d support directions" % (
        len(model.levels), res.added, len(dirs)))
    return res


def _close_end(model, levels, measured, band_edge, band_h, end, side, height, reach, slack, report):
    """Finish one end: a flat cap at the part's extent, or an apex (point or ridge)."""
    z1, h1 = measured[-1]
    # Extrapolation needs a baseline long enough to mean something.
    base = max(END_MARGIN * height, abs(z1 - band_edge) * 0.2)
    z0, h0 = band_edge, band_h
    for z, hs in measured[:-1]:
        if abs(z1 - z) >= base:
            z0, h0 = z, hs
    span = z1 - z0
    if abs(span) <= 1e-12:
        levels.append((end, list(h1)))
        return
    width = min(h1[i] + h1[j] for i, j in enumerate(model.opposite)) - 2.0 * model.clearance
    word = "top" if side > 0 else "bottom"
    if width < TAPER_RATIO * reach:
        got = _apex(model, h0, h1, reach, slack)
        if got is not None:
            t_end, hs = got
            z_apex = z0 + span * t_end
            beyond = (z_apex - end) * side
            # A hair inside the end (no clearance) still counts: the apex is then
            # placed at the end itself, which only widens the taper above it.
            if beyond >= -STUB * height and beyond <= MAX_OVERSHOOT * height:
                beyond = max(beyond, 0.0)
                z_apex = end + side * beyond
                corners = model.corners(hs)
                # The part's own extent there: the section is the clearance wider all round.
                own = _span(corners) - 2.0 * model.clearance
                if own < POINT_RATIO * reach:
                    t0 = pinch(model, h0, h1, 0.0)
                    z_tip = z0 + span * (t0 if t0 is not None else t_end)
                    # Further out only widens the taper above it: the safe side.
                    z_tip = end + side * min(max((z_tip - end) * side, 0.002 * height), MAX_OVERSHOOT * height)
                    cx, cy = poly2d.centroid(corners)
                    model.tips[side] = (cx, cy, z_tip)
                    report.note("the %s tapers to a point" % word)
                    return
                levels.append((z_apex, hs))
                levels.append((z_apex + side * STUB * height, list(hs)))
                report.note("the %s tapers to a ridge %.4g long" % (word, own))
                return
        report.note("the %s tapers, but no apex could be extrapolated safely; it is "
                    "carried on straight" % word)
        levels.append((end + side * STUB * height, list(h1)))
        return
    t_end = (end - z0) / span
    hs, t = safe_blend(model, h0, h1, max(1.0, min(t_end, SAFE_T)), 0.004 * reach, slack)
    levels.append((end, hs))
    if t < t_end - 1e-6 and t_end > 1.0 + 1e-6:
        report.note("the %s end is carried on straight (extrapolating it further would "
                    "turn its section inside out)" % word)


def _apex(model, h0, h1, reach, slack):
    """The first clean section where the extrapolated taper thins to a floor,
    backing the floor off until one comes clean. (t, offsets) or None."""
    floor = 0.002 * reach
    for _ in range(14):
        t = pinch(model, h0, h1, floor)
        if t is not None:
            hs = [a + (b - a) * t for a, b in zip(h0, h1)]
            if model.clean(hs, slack):
                return t, hs
        floor *= 1.6
    return None


def _span(points):
    best = 0.0
    for i in range(len(points)):
        for j in range(i + 1, len(points)):
            best = max(best, math.hypot(points[i][0] - points[j][0], points[i][1] - points[j][1]))
    return best


def _dedup(levels, gap, keep=()):
    """Levels sorted by height with any closer than gap to the previous dropped.
    Of two levels too close together, one whose height is in keep wins."""
    out = [levels[0]]
    for z, hs in levels[1:]:
        if z - out[-1][0] > gap:
            out.append((z, hs))
        elif z in keep and out[-1][0] not in keep and len(out) > 1:
            out[-1] = (z, hs)
    return out


def _prune(levels, eps):
    """Drop levels that lie on the straight blend of their neighbours: they add
    no shape, only faces (a straight side would be split at every level)."""
    out = list(levels)
    k = 1
    while k < len(out) - 1:
        (z0, h0), (z1, h1), (z2, h2) = out[k - 1], out[k], out[k + 1]
        t = (z1 - z0) / (z2 - z0)
        if all(abs(b - (a + (c - a) * t)) <= eps for a, b, c in zip(h0, h1, h2)):
            del out[k]
        else:
            k += 1
    return out


def _refine(model, points, target, measure, band, height, tolerance):
    """Add measured levels where part vertices come out closer than target."""
    gap = max(1e-4 * height, 2.0 * tolerance)
    added = 0
    for _ in range(REFINE_PASSES):
        lv = model.levels
        zs = [l[0] for l in lv]
        worst = {}
        for p in points:
            z = p[2]
            if z <= zs[0] + gap or z >= zs[-1] - gap or band[0] <= z <= band[1]:
                continue
            s = model.slack(p)
            if s is None or s >= target:
                continue
            k = bisect.bisect_right(zs, z) - 1
            if k not in worst or s < worst[k][0]:
                worst[k] = (s, z)
        new = []
        for s, z in worst.values():
            hs = measure(z)
            if hs is not None:
                new.append((z, hs))
        if not new:
            break
        fresh = set(z for z, _ in new)
        merged = _dedup(sorted(lv + new, key=lambda l: l[0]), gap, fresh)
        if not fresh.intersection(z for z, _ in merged):
            break
        added += len(fresh.intersection(z for z, _ in merged))
        model.levels = merged
    return added


def _grade(res, model, points, clearance, tolerance, reach):
    r = res.report
    worst, outside, beyond = None, 0, 0
    for p in points:
        s = model.slack(p)
        if s is None:
            if p[2] < model.bottom() - tolerance or p[2] > model.top() + tolerance:
                beyond += 1
            continue
        if worst is None or s < worst:
            worst = s
        if s < -tolerance:
            outside += 1
    res.min_clearance = worst
    if beyond or outside:
        r.fail("the part pokes out of the solid at %d vertices" % (beyond + outside))
    elif clearance > 0.0 and worst is not None and worst < 0.5 * clearance:
        r.warn("the sideways clearance drops to %.4g (asked for %.4g): the sections miss part "
               "of the profile even after refinement" % (worst, clearance))
    if worst is not None:
        r.note("smallest sideways clearance %.4g" % worst)


# ---------------------------------------------------------------- side profile

PROFILE_SAMPLES = 400   # profile height bins across the part's width
PROFILE_PLANES = 3      # cutting planes from the centre out to the touching vertex
PROFILE_RMS = 0.02      # a broken line further off than this x the width is a poor fit


class ProfileResult(object):
    """points: the side profile as (z, x), z increasing (x: reach toward the
    direction). fit: fit.BrokenLine or None. frame: (d, e, u) unit axes and
    origin: where profile coordinates start; plane_offset: the e coordinate
    of the plane the profile is drawn on (through the touching vertex)."""

    def __init__(self):
        self.points = []
        self.fit = None
        self.frame = None
        self.origin = None
        self.plane_offset = 0.0
        self.report = Report("profile fit")

    @property
    def status(self):
        return self.report.status

    def to_world(self, x, z):
        d, e, u = self.frame
        o, y = self.origin, self.plane_offset
        return tuple(o[k] + d[k] * x + e[k] * y + u[k] * z for k in range(3))


def side_profile(vertices, faces, axis, direction, samples=PROFILE_SAMPLES, planes=PROFILE_PLANES,
                 slope_limits=None):
    """The part's outline seen from the side, and a three-piece broken line fitted to it.

    axis       the part's axis (profile heights run along it)
    direction  the side to look at: the profile is how far the part reaches
               toward it, square to the axis
    The part is cut by `planes` planes containing the axis direction and the
    side direction, spread from the centre out to the vertex that reaches
    furthest; only the outermost point per height bin is kept, so a noisy
    or coarse mesh still gives one clean outline. The broken line
    (fit.three_segment_fit) finds where the outline bends: its two knees.
    slope_limits are passed on to the fit (dx/dz per piece, None for free).
    """
    res = ProfileResult()
    r = res.report
    if not faces:
        r.fail("the mesh has no faces")
        return res
    u = vec.unit(axis)
    if u is None:
        r.fail("the axis is a zero vector")
        return res
    d = vec.unit(vec.sub(direction, vec.scale(u, vec.dot(direction, u)))) if direction is not None else None
    if d is None:
        r.fail("the direction must point to a side, not along the axis")
        return res
    e = vec.cross(u, d)
    res.frame = (d, e, u)
    used = sorted(set(i for f in faces for i in f))
    origin = vec.mean([vertices[i] for i in used])
    res.origin = origin
    local = [None] * len(vertices)
    for i in used:
        p = vec.sub(vertices[i], origin)
        local[i] = (vec.dot(p, d), vec.dot(p, e), vec.dot(p, u))
    xs = [local[i][0] for i in used]
    width = max(xs) - min(xs)
    if width <= 0.0:
        r.fail("the part has no width toward the direction")
        return res
    touch = max(used, key=lambda i: local[i][0])
    y_touch = local[touch][1]
    res.plane_offset = y_touch
    step = width / float(samples)
    keep_from = min(xs) + 0.75 * width  # only the side facing the direction
    offsets = sorted(y_touch * k / max(1.0, planes - 1.0) for k in range(max(1, planes)))
    pts = [p if p is not None else (0.0, 0.0, 0.0) for p in local]
    bins = {}
    for segs in slice_segments(pts, faces, offsets, axis=1):
        for a, b in segs:
            _bin_segment(bins, (a[0], a[2]), (b[0], b[2]), step, keep_from)
    res.points = sorted(bins.values())
    if len(res.points) < 8:
        r.fail("the side profile has only %d points" % len(res.points))
        return res
    zs = [p[0] for p in res.points]
    xs_ = [p[1] for p in res.points]
    res.fit = fit.three_segment_fit(zs, xs_, slope_limits)
    if res.fit is None:
        r.fail("too few profile points for a three-piece fit")
        return res
    if res.fit.rms > PROFILE_RMS * width:
        r.warn("the profile is %.3g off a three-piece broken line (rms); its knees are rough" % res.fit.rms)
    if slope_limits is not None:
        for k, (line, limit) in enumerate(zip(res.fit.lines, slope_limits)):
            if limit is not None and (line.slope <= limit[0] + 1e-9 or line.slope >= limit[1] - 1e-9):
                r.warn("piece %d of the broken line hit its slope limit" % (k + 1))
    (z0, x0), (z1, x1) = res.fit.knees
    r.note("knees at heights %.4g and %.4g (%.4g apart), %d profile points, rms %.3g" % (
        z0, z1, z1 - z0, len(res.points), res.fit.rms))
    return res


def _bin_segment(bins, a, b, step, keep_from):
    """Densify a 2D segment (x, z) and keep the outermost x per height bin."""
    (x0, z0), (x1, z1) = a, b
    if x0 < keep_from and x1 < keep_from:
        return
    n = max(1, int(math.ceil(math.hypot(x1 - x0, z1 - z0) / step)))
    for j in range(n + 1):
        t = j / float(n)
        x = x0 + (x1 - x0) * t
        if x < keep_from:
            continue
        z = z0 + (z1 - z0) * t
        k = int(math.floor(z / step))
        best = bins.get(k)
        if best is None or x > best[1]:
            bins[k] = (z, x)
