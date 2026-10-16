"""Natural axis of a part: the direction it "obviously" sits along.

Many parts have one: the axis of a turned or faceted body, the normal of a
plate, the length of a bar. This module finds it from the mesh alone, at any
orientation, and grades how sure it is.

Candidates
    The normals of the largest flat facets, plus the three principal axes of
    the surface (exact second moments of the triangles), each refined.

"Level"
    Look down a candidate axis and take, in each direction across it, the
    centre of whatever touches the outline (the silhouette) there. Those
    outline points form a rim: the edge of a lens, the widest band of a turned
    part, or the mid-height loop of a prism seen down its walls (a face seen
    edge-on touches the outline along a line, and its centre is what counts).
    A candidate is level when a plane fitted through its rim is square to it
    within LEVEL_DEG.

Choice
    Among level candidates, one with a strong flat end wins (a flat facet
    square to the axis, filling most of that end). Otherwise the one showing
    the largest outline wins.

The mirror-plane trap
    For a symmetric part, the view square to any mirror plane ALSO looks
    level: the silhouette is symmetric about that plane, so a plane fitted
    through it is the mirror plane itself, exactly square to the view.
    Levelness alone therefore cannot tell the true axis from a side view. The
    flat-end preference and the outline area are what separate them, and a
    remaining near tie is reported as ambiguous instead of being guessed.
"""

from __future__ import absolute_import, division

import math

from . import poly2d, vec
from .facets import cluster_planes
from .report import Report
from .topology import bounding_size, face_area

LEVEL_DEG = 1.5          # a candidate within this of square counts as level
BAND_RATIO = 0.008       # outline band: points within this x size of an extreme
DIRECTIONS = 32          # directions around the axis used to find outline points
MAX_POINTS = 2500        # points used for evaluation (evenly thinned beyond this)
MAX_FACETS = 12          # largest facets used as candidates
MIN_FACET_SHARE = 0.004  # ...when at least this share of the total area
RIM_FLATNESS = 0.05      # outline points count as a rim when this flat (plane_fit)
STRONG_COVER = 0.5       # a flat end fills at least this share of its end's outline
STRONG_CAP = 0.1         # ...and is at least this share of the whole outline
CLEAN_COVER = 0.9        # below this, the flat end is called noisy
RIVAL_DEG = 10.0         # a level candidate this far away...
RIVAL_AREA = 0.9         # ...showing this share of the chosen outline is a rival
GIVEN_WARN_DEG = 1.0     # check_axis: warn beyond this difference
GIVEN_MESSY_DEG = 5.0    # check_axis: messy beyond this difference


class Candidate(object):
    """One direction under evaluation (see evaluate())."""

    def __init__(self, axis, source):
        self.axis = axis
        self.source = source          # "facet", "surface" or "given"
        self.tilt = 90.0              # degrees from level (90 when the outline is no plane)
        self.rim_normal = None        # normal of the plane through the outline points
        self.outline_area = 0.0
        self.ends = {}                # +1 / -1 -> (facet, cover, cap share)

    @property
    def level(self):
        return self.tilt <= LEVEL_DEG

    def flat_end(self, side):
        return self.ends.get(side)

    @property
    def strong_end(self):
        """The strongest flat end as (side, facet, cover), or None."""
        best = None
        for side, (facet, cover, cap) in self.ends.items():
            if cover >= STRONG_COVER and cap >= STRONG_CAP:
                if best is None or facet.area > best[1].area:
                    best = (side, facet, cover)
        return best

    def __repr__(self):
        return "Candidate(%s, tilt=%.2f, area=%.4g, ends=%d)" % (
            self.source, self.tilt, self.outline_area, len(self.ends))


class AxisResult(object):
    """axis: unit vector, pointing to the part's top. frame: (x, y, z) around it."""

    def __init__(self):
        self.axis = None
        self.frame = None
        self.chosen = None
        self.candidates = []
        self.report = Report("natural axis")

    @property
    def status(self):
        return self.report.status


# ---------------------------------------------------------------- evaluation

def _sample(vertices, faces):
    used = sorted(set(i for f in faces for i in f))
    step = max(1, int(math.ceil(len(used) / float(MAX_POINTS))))
    return [vertices[i] for i in used[::step]]


def _line_angle(a, b):
    """Angle between two lines (direction sign ignored), in degrees."""
    c = abs(vec.dot(vec.unit(a), vec.unit(b)))
    return math.degrees(math.acos(min(1.0, c)))


def evaluate(points, axis, facets, size, source="given", mesh=None):
    """Measure one candidate axis: tilt from level, outline area, flat ends.

    mesh is (vertices, faces, face_areas, faces_of_vertex), as built by
    natural_axis(). With it, a flat end's cover is the
    share of the end's own triangles that its facet holds, exact for any end
    shape; without it, the end's outline hull stands in for those triangles.
    """
    cand = Candidate(vec.unit(axis), source)
    z = cand.axis
    u, v = vec.perpendicular_pair(z)
    local = [(vec.dot(p, u), vec.dot(p, v), vec.dot(p, z)) for p in points]
    band = BAND_RATIO * size
    cand.outline_area = poly2d.area(poly2d.convex_hull([(x, y) for x, y, _ in local]))

    hs = [h for _, _, h in local]
    h_lo, h_hi = min(hs), max(hs)

    # One outline point per direction: the centre of everything touching the
    # outline there. A face seen edge-on touches it along a whole line running
    # down the view; its centre is the point that belongs to the rim.
    rim_points = []
    for k in range(DIRECTIONS):
        a = 2.0 * math.pi * k / DIRECTIONS
        c, s = math.cos(a), math.sin(a)
        reach = [x * c + y * s for x, y, _ in local]
        m = max(reach)
        rim_points.append(vec.mean([local[i] for i, rr in enumerate(reach) if rr >= m - band]))

    _, n, _, flatness = vec.plane_fit(rim_points)
    if flatness <= RIM_FLATNESS:
        cand.tilt = _line_angle(n, (0.0, 0.0, 1.0))
        world = vec.add(vec.add(vec.scale(u, n[0]), vec.scale(v, n[1])), vec.scale(z, n[2]))
        cand.rim_normal = world if vec.dot(world, z) >= 0.0 else vec.scale(world, -1.0)

    # Flat ends: facets square to the axis, sitting at either extreme.
    cos_level = math.cos(math.radians(2.0 * LEVEL_DEG))
    heights = None
    if mesh is not None:
        heights = [vec.dot(p, z) for p in mesh[0]]
        h_lo, h_hi = min(heights), max(heights)
    for side, extreme in ((1, h_hi), (-1, h_lo)):
        best = None
        for fc in facets:
            if vec.dot(fc.normal, z) * side < cos_level:
                continue
            if abs(vec.dot(fc.centre, z) - extreme) > band:
                continue
            if best is None or fc.area > best.area:
                best = fc
        if best is None or cand.outline_area <= 0.0:
            continue
        if heights is not None:
            # The end's own triangles: at the extreme, inside the facet's outline.
            # Noise fragments of a flat end lie inside it; a sloped band around
            # the end lies outside, so it does not count against the facet.
            corners = set(i for k in best.faces for i in mesh[1][k])
            hull = poly2d.convex_hull([(vec.dot(mesh[0][i], u), vec.dot(mesh[0][i], v)) for i in corners])
            at_end = set(i for i, h in enumerate(heights) if abs(h - extreme) <= band)
            near = set(k for i in at_end for k in mesh[3][i])
            own = set(best.faces)
            end_area = best.area
            # Bounding circle of the outline: a cheap first test before the full one.
            hc = (sum(q[0] for q in hull) / len(hull), sum(q[1] for q in hull) / len(hull))
            reach_sq = max((q[0] - hc[0]) ** 2 + (q[1] - hc[1]) ** 2 for q in hull)
            for k in near:
                if k in own:
                    continue
                f = mesh[1][k]
                if f[0] in at_end and f[1] in at_end and f[2] in at_end:
                    c = vec.mean([mesh[0][i] for i in f])
                    q = (vec.dot(c, u), vec.dot(c, v))
                    if (q[0] - hc[0]) ** 2 + (q[1] - hc[1]) ** 2 > reach_sq:
                        continue
                    if poly2d.contains_convex(hull, q):
                        end_area += mesh[2][k]
        else:
            end_pts = [(x, y) for x, y, h in local if abs(h - extreme) <= band]
            end_area = poly2d.area(poly2d.convex_hull(end_pts)) if len(end_pts) >= 3 else 0.0
        if end_area > 0.0:
            cover = min(1.0, best.area / end_area)
            cand.ends[side] = (best, cover, best.area / cand.outline_area)
    return cand


def refine(points, axis, facets, size, rounds=6, mesh=None):
    """Turn a candidate toward the normal of its own rim plane until it settles."""
    cand = evaluate(points, axis, facets, size, "surface", mesh)
    for _ in range(rounds):
        if cand.rim_normal is None:
            break
        if _line_angle(cand.rim_normal, cand.axis) < 0.01:
            break
        cand = evaluate(points, cand.rim_normal, facets, size, "surface", mesh)
    return cand


def surface_axes(vertices, faces):
    """Principal axes of the surface (exact second moments of its triangles)."""
    total = 0.0
    sx = sy = sz = 0.0
    m = [[0.0] * 3 for _ in range(3)]
    for f in faces:
        a, b, c = vertices[f[0]], vertices[f[1]], vertices[f[2]]
        ar = face_area(vertices, f)
        if ar <= 0.0:
            continue
        s = vec.add(vec.add(a, b), c)
        total += ar
        sx += ar * s[0] / 3.0
        sy += ar * s[1] / 3.0
        sz += ar * s[2] / 3.0
        for p in (a, b, c, s):
            for i in range(3):
                for j in range(i, 3):
                    m[i][j] += ar * p[i] * p[j] / 12.0
    if total <= 0.0:
        return []
    centre = (sx / total, sy / total, sz / total)
    for i in range(3):
        for j in range(i, 3):
            m[i][j] = m[i][j] / total - centre[i] * centre[j]
    return vec.eig_sym3(m)[1]


# ---------------------------------------------------------------- the decision

def natural_axis(vertices, faces, facets=None):
    """Find the natural axis of a (cleaned) triangle mesh. See the module docs.

    Returns an AxisResult whose axis points to the part's top: toward its
    stronger flat end, or away from its tapering end when it has none.
    """
    res = AxisResult()
    r = res.report
    if not faces:
        r.fail("the mesh has no faces")
        return res
    size = bounding_size(vertices)
    if size <= 0.0:
        r.fail("the mesh has no extent")
        return res
    if facets is None:
        facets = cluster_planes(vertices, faces)
    points = _sample(vertices, faces)
    faces_of_vertex = [[] for _ in vertices]
    for k, f in enumerate(faces):
        for i in f:
            faces_of_vertex[i].append(k)
    mesh = (vertices, faces, [face_area(vertices, f) for f in faces], faces_of_vertex)
    total_area = sum(fc.area for fc in facets) or 1.0

    candidates = []

    def add(cand):
        # Near-duplicates: a facet candidate wins (its plane is exact), otherwise the more level one.
        for k, other in enumerate(candidates):
            if _line_angle(other.axis, cand.axis) < 0.5:
                if other.source != "facet" and (cand.source == "facet" or cand.tilt < other.tilt):
                    candidates[k] = cand
                return
        candidates.append(cand)

    for fc in facets[:MAX_FACETS]:
        if fc.area >= MIN_FACET_SHARE * total_area:
            add(evaluate(points, fc.normal, facets, size, "facet", mesh))
    for axis in surface_axes(vertices, faces):
        add(refine(points, axis, facets, size, mesh=mesh))
    res.candidates = candidates
    if not candidates:
        r.fail("no candidate direction could be measured")
        return res

    level = [c for c in candidates if c.level]
    if not level:
        best = min(candidates, key=lambda c: c.tilt)
        r.messy("no direction looks level; using the closest one (%.1f degrees off)" % best.tilt)
    else:
        flat = [c for c in level if c.strong_end is not None]
        pool = flat if flat else level
        best = max(pool, key=lambda c: c.outline_area)
        for c in level:
            if c is best or _line_angle(c.axis, best.axis) <= RIVAL_DEG:
                continue
            if (c.strong_end is not None) != (best.strong_end is not None):
                continue
            if c.outline_area >= RIVAL_AREA * best.outline_area:
                r.messy("two directions %.0f degrees apart fit about equally well" % (
                    _line_angle(c.axis, best.axis)))
                break

    axis = best.axis
    strong = best.strong_end
    # A sizeable flat end that is broken up (scan noise) cannot be trusted for
    # the exact direction: follow the outline plane instead, and say so.
    noisy = [e for e in best.ends.values() if e[2] >= STRONG_CAP and e[1] < CLEAN_COVER]
    if noisy:
        refined = refine(points, best.axis, facets, size, mesh=mesh)
        if refined.level:
            axis = refined.axis
            r.warn("the flat end is noisy; the axis follows the outline instead")
    elif best.source == "facet" and strong is None and best.rim_normal is not None:
        axis = best.rim_normal

    res.chosen = best
    res.axis = _point_to_top(points, axis, best, r)
    res.frame = vec.frame_from_axis(res.axis)
    return res


def _point_to_top(points, axis, cand, report):
    """Sign the axis: toward the stronger flat end, else away from the taper."""
    top, bottom = cand.flat_end(1), cand.flat_end(-1)
    if vec.dot(axis, cand.axis) < 0.0:  # keep "top" meaning the candidate's +1 side
        axis = vec.scale(axis, -1.0)
    a_top = top[0].area if top else 0.0
    a_bottom = bottom[0].area if bottom else 0.0
    if a_top > 0.0 or a_bottom > 0.0:
        if abs(a_top - a_bottom) <= 0.01 * max(a_top, a_bottom):
            report.note("both ends are alike; which one is the top is arbitrary")
            return _canonical(axis)
        return axis if a_top > a_bottom else vec.scale(axis, -1.0)
    # No flat end: the long thin tail of the height distribution is the taper.
    hs = [vec.dot(p, axis) for p in points]
    mean = sum(hs) / len(hs)
    third = sum((h - mean) ** 3 for h in hs) / len(hs)
    second = sum((h - mean) ** 2 for h in hs) / len(hs)
    if second <= 0.0 or abs(third) <= 1e-6 * second ** 1.5:
        report.note("no flat end and no taper; which end is the top is arbitrary")
        return _canonical(axis)
    report.note("no flat end; the axis points away from the tapering end")
    return vec.scale(axis, -1.0) if third > 0.0 else axis


def _canonical(axis):
    """Direction sign made repeatable: the largest component positive."""
    k = max(range(3), key=lambda i: abs(axis[i]))
    return axis if axis[k] >= 0.0 else vec.scale(axis, -1.0)


def check_axis(vertices, faces, axis, facets=None):
    """Grade an axis someone supplied against the detected natural axis.

    The result keeps the GIVEN axis (signed toward the detected top) and
    reports how far it is from what the mesh suggests.
    """
    detected = natural_axis(vertices, faces, facets)
    res = AxisResult()
    res.report.title = "given axis"
    given = vec.unit(axis)
    if given is None:
        res.report.fail("the given axis has zero length")
        return res
    res.candidates = detected.candidates
    res.chosen = detected.chosen
    if detected.axis is None:
        res.report.merge(detected.report, prefix="detection: ")
        res.axis = given
    else:
        off = _line_angle(given, detected.axis)
        if off > GIVEN_MESSY_DEG:
            res.report.messy("the given axis is %.1f degrees off the detected one" % off)
        elif off > GIVEN_WARN_DEG:
            res.report.warn("the given axis is %.1f degrees off the detected one" % off)
        else:
            res.report.note("the given axis agrees with the detected one (%.2f degrees)" % off)
        res.axis = given if vec.dot(given, detected.axis) >= 0.0 else vec.scale(given, -1.0)
    res.frame = vec.frame_from_axis(res.axis)
    return res
