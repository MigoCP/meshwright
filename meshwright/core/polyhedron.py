"""Convex solids from planes: half-space clipping, offsets and a validated corner merge.

A Polyhedron is a closed solid with flat polygon faces: one shared vertex
list, and faces as lists of vertex indices, counter-clockwise seen from
outside. Faces share corners by index, never by comparing coordinates, so a
cut can never leave two copies of one corner a hair apart.

Building a solid from planes
    Start from a box around the part and cut it with each plane in turn
    (clip()). Where a plane crosses an edge, the new corner is made once for
    that edge and shared by both faces along it. The hole a cut leaves is
    closed by a new face whose rim is chained from the cut faces' own edges,
    (every kept edge left without a partner), so corners on straight runs of
    the rim are kept and the rim always matches its neighbours.

Offsets
    Each facet plane is first moved out until it touches the part, so a noisy
    or merged facet can only make the solid bigger, never cut into the part.
    Then it moves out by the clearance. With an axis, the clearance is
    applied square to the axis only: a plane moves by the clearance times the
    sideways share of its normal. A flat end square to the axis stays put, a
    side parallel to the axis moves the full clearance.

Split corners and the validated merge
    Where four or more planes meet at one corner of the part, moving them by
    different amounts splits that corner into a cluster of corners very close
    together, joined by edges far shorter than any modelling tolerance.
    merge_corners() joins each cluster (union-find, so merging is transitive)
    and validate() refits every face's plane and checks that the result is
    still a closed solid with flat faces that keeps the clearance.
    offset_solid() uses the biggest merge that passes.

Only convex solids come out of this module. A part that is not convex
gets a convex solid that bridges its hollows, built from its own facet
directions only (an L-shaped block gets its bounding block), and the report
says so.
"""

from __future__ import absolute_import, division

import math

from . import poly2d, vec
from .facets import cluster_planes
from .report import Report
from .topology import UnionFind, bounding_size

MAX_PLANES = 250          # more facets than this: the part is curved, not faceted
MERGE_STEPS = (0.25, 0.1)  # corner merge distances tried, as shares of the clearance
KEEP_CLEARANCE = 0.9      # a merged solid must keep this share of the clearance
FLAT_SHARE = 0.5          # faces must stay flat within this x the tolerance
SLACK_SHARE = 0.5         # differences below this x the tolerance count as rounding
EPS_RATIO = 1e-9          # clipping: corners this x the size from a plane lie on it
ISLAND_EPS = 1000.0       # clipping: stray outside corners up to this x eps count as on it
NEAR_RATIO = 0.05         # points this x the size below a plane's support are checked
REACH_CLEARANCES = 5.0    # a corner further out than this x the clearance...
REACH_RATIO = 0.05        # ...plus this x the size points to a missing facet
CONCAVE_RMS = 4.0         # a facet crossed by the part beyond this x its rms (and the tolerance)


class Polyhedron(object):
    """A closed solid with flat faces.

    vertices  (x, y, z) tuples
    faces     lists of vertex indices, counter-clockwise seen from outside
    holes     per face, its inner loops (lists of vertex indices, clockwise
              seen from outside); empty for faces without holes
    labels    per face, what made it: the index of the cutting plane or facet,
              or None for a face left over from the starting box
    """

    def __init__(self, vertices, faces, labels=None, holes=None):
        self.vertices = list(vertices)
        self.faces = [list(f) for f in faces]
        self.labels = list(labels) if labels is not None else [None] * len(self.faces)
        if holes is None:
            self.holes = [[] for _ in self.faces]
        else:
            self.holes = [[list(h) for h in hs] for hs in holes]

    def loops(self, k):
        """Face k's outer loop followed by its holes, as vertex index lists."""
        return [self.faces[k]] + self.holes[k]

    def all_loops(self):
        for k in range(len(self.faces)):
            for loop in self.loops(k):
                yield loop

    def polygon(self, k):
        """Corner points of face k's outer loop."""
        return [self.vertices[i] for i in self.faces[k]]

    def loop_points(self, k):
        """Corner points of every loop of face k, outer loop first."""
        return [[self.vertices[i] for i in loop] for loop in self.loops(k)]

    def area_vector(self, k):
        """Face k's normal scaled by its area (holes subtracted)."""
        total = (0.0, 0.0, 0.0)
        for pts in self.loop_points(k):
            total = vec.add(total, vec.area_vector(pts))
        return total

    def plane(self, k):
        """(unit normal, offset, flatness) of face k: the Newell normal, the mean
        offset of its corners and their largest distance from that plane.
        None for a face with no area."""
        n = vec.unit(self.area_vector(k))
        if n is None:
            return None
        ds = [vec.dot(n, p) for pts in self.loop_points(k) for p in pts]
        d = sum(ds) / len(ds)
        return n, d, max(abs(x - d) for x in ds)

    def triangles(self):
        """(vertices, triangles) covering every face, holes left open.

        A convex face without holes is split into a fan; other faces are
        triangulated in their plane (poly2d.triangulate)."""
        tris = []
        for k, f in enumerate(self.faces):
            if not self.holes[k] and self._convex(k):
                for i in range(1, len(f) - 1):
                    tris.append((f[0], f[i], f[i + 1]))
                continue
            n = vec.unit(self.area_vector(k))
            if n is None:
                continue
            u, v = vec.perpendicular_pair(n)
            ids = [i for loop in self.loops(k) for i in loop]
            flat = [[(vec.dot(self.vertices[i], u), vec.dot(self.vertices[i], v)) for i in loop]
                    for loop in self.loops(k)]
            for a, b, c in poly2d.triangulate(flat[0], flat[1:]):
                tris.append((ids[a], ids[b], ids[c]))
        return list(self.vertices), tris

    def _convex(self, k):
        f = self.faces[k]
        n = self.area_vector(k)
        m = len(f)
        for i in range(m):
            a, b, c = self.vertices[f[i - 1]], self.vertices[f[i]], self.vertices[f[(i + 1) % m]]
            if vec.dot(vec.cross(vec.sub(b, a), vec.sub(c, b)), n) < 0.0:
                return False
        return True

    def volume(self):
        if not self.faces:
            return 0.0
        ref = vec.mean(self.vertices)
        total = 0.0
        for loop in self.all_loops():
            a = vec.sub(self.vertices[loop[0]], ref)
            for i in range(1, len(loop) - 1):
                b = vec.sub(self.vertices[loop[i]], ref)
                c = vec.sub(self.vertices[loop[i + 1]], ref)
                total += vec.dot(a, vec.cross(b, c))
        return total / 6.0

    def area(self):
        return sum(vec.length(self.area_vector(k)) for k in range(len(self.faces)))

    def edge_lengths(self):
        """Length of every edge, each counted once."""
        out = []
        for loop in self.all_loops():
            for k in range(len(loop)):
                a, b = loop[k], loop[(k + 1) % len(loop)]
                if a < b:
                    out.append(vec.distance(self.vertices[a], self.vertices[b]))
        return out

    def __repr__(self):
        return "Polyhedron(%d vertices, %d faces)" % (len(self.vertices), len(self.faces))


def box(lo, hi):
    """The axis-aligned box from corner lo to corner hi."""
    x0, y0, z0 = lo
    x1, y1, z1 = hi
    vertices = [(x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
                (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)]
    faces = [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]]
    return Polyhedron(vertices, faces)


# ---------------------------------------------------------------- clipping

def clip(poly, normal, offset, eps=0.0, label=None):
    """The part of a convex polyhedron (faces without holes) where dot(normal, p) <= offset.

    Corners within eps of the plane count as lying on it. Returns
    (polyhedron, problem): the polyhedron itself when nothing lies outside,
    a new one with the cut closed by a face labelled `label`, or None plus a
    readable problem when nothing would remain or the cut cannot be closed.

    Rounding near tiny edges: on an exact convex solid, the corners outside
    a plane always form one connected group. When noise splits off other
    groups (a chain of corners all within a hair of the plane), those count
    as lying on it, provided they are no further out than ISLAND_EPS x eps.
    """
    if any(poly.holes):
        raise ValueError("clip() works on convex polyhedra, which have no faces with holes")
    dist = [vec.dot(normal, p) - offset for p in poly.vertices]
    if not dist or max(dist) <= eps:
        return poly, None
    if min(dist) >= -eps:
        return None, "a plane cuts the whole solid away"
    side = [0 if abs(x) <= eps else (1 if x > 0.0 else -1) for x in dist]
    island = ISLAND_EPS * eps
    if not _one_outside_group(poly, dist, side, island):
        return None, "a plane cuts the solid in more than one place"
    vertices = list(poly.vertices)
    made = {}  # edge (i, j), i < j -> index of the corner made where the plane crosses it

    def crossing(i, j):
        key = (i, j) if i < j else (j, i)
        k = made.get(key)
        if k is None:
            a, b = key
            t = dist[a] / (dist[a] - dist[b])
            k = len(vertices)
            vertices.append(vec.lerp(poly.vertices[a], poly.vertices[b], t))
            made[key] = k
        return k

    faces, labels = [], []
    for face, lab in zip(poly.faces, poly.labels):
        kept = []
        m = len(face)
        for k in range(m):
            a, b = face[k], face[(k + 1) % m]
            if side[a] <= 0:
                kept.append(a)
            if side[a] * side[b] < 0:
                kept.append(crossing(a, b))
        if len(kept) < 3:
            continue
        if all(i >= len(side) or side[i] == 0 for i in kept):
            continue  # lies in the cutting plane: the new face replaces it
        faces.append(kept)
        labels.append(lab)

    # The new face's rim is every kept edge left without its partner running the
    # other way. Matching edges instead of testing them against the plane keeps
    # the result closed even where a cluster of corners lies within eps of it.
    edges = set()
    for f in faces:
        for k in range(len(f)):
            edges.add((f[k], f[(k + 1) % len(f)]))
    rim = {}  # corner -> next corner of the new face, counter-clockwise from outside
    for a, b in edges:
        if (b, a) not in edges:
            if b in rim:
                return None, "a cut leaves a pinched rim"
            rim[b] = a

    if len(rim) < 3:
        if max(dist) <= island:
            return poly, None  # only a corner a hair outside, with nothing to cut off
        return None, "a cut leaves no face to close it"
    start = min(rim)
    loop = [start]
    current = rim[start]
    while current != start:
        if current not in rim or len(loop) > len(rim):
            return None, "a cut leaves an open rim"
        loop.append(current)
        current = rim[current]
    if len(loop) != len(rim):
        return None, "a cut leaves more than one rim"
    faces.append(loop)
    labels.append(label)
    return _compact(vertices, faces, labels), None


def _one_outside_group(poly, dist, side, island):
    """Keep only the connected group of outside corners that holds the furthest
    one; set stray groups (none of them more than `island` out) to on-plane.
    Returns False when a stray group is further out than that."""
    outside = [i for i, s in enumerate(side) if s > 0]
    if len(outside) < 2:
        return True
    uf = UnionFind(len(side))
    for f in poly.faces:
        for k in range(len(f)):
            a, b = f[k], f[(k + 1) % len(f)]
            if side[a] > 0 and side[b] > 0:
                uf.union(a, b)
    furthest = uf.find(max(outside, key=lambda i: dist[i]))
    strays = [i for i in outside if uf.find(i) != furthest]
    if any(dist[i] > island for i in strays):
        return False
    for i in strays:
        side[i] = 0
    return True


def _compact(vertices, faces, labels, holes=None):
    """Drop unused vertices and renumber, keeping the order they are first used in."""
    new_index = {}
    out_vertices = []

    def renumber(loop):
        g = []
        for i in loop:
            k = new_index.get(i)
            if k is None:
                k = len(out_vertices)
                new_index[i] = k
                out_vertices.append(vertices[i])
            g.append(k)
        return g

    out_faces = [renumber(f) for f in faces]
    out_holes = None
    if holes is not None:
        out_holes = [[renumber(h) for h in hs] for hs in holes]
    return Polyhedron(out_vertices, out_faces, labels, out_holes)


def intersect_planes(planes, lo, hi, eps=0.0, order=None):
    """The convex solid inside every plane (normal, offset), within the box lo..hi.

    Plane k labels the face it makes. order: the plane indices in the order
    to cut with (default: as given). Returns (polyhedron, problem) like
    clip(). Faces still labelled None afterwards are box faces: the planes
    leave that direction open.
    """
    poly = box(lo, hi)
    for k in (range(len(planes)) if order is None else order):
        poly, problem = clip(poly, planes[k][0], planes[k][1], eps, label=k)
        if poly is None:
            return None, problem
    return poly, None


# ---------------------------------------------------------------- merge and checks

def merge_corners(poly, distance, place=None):
    """Join corners closer than distance (transitively).

    A merged corner goes to place(indices of the joined corners), or to
    their mean point without a place function. Faces whose outer loop
    collapses to fewer than three corners are dropped, and so are holes that
    collapse. The result is not checked here: run validate() on it.
    """
    count = len(poly.vertices)
    uf = UnionFind(count)
    if distance > 0.0:
        grid = {}
        limit = distance * distance
        for i, p in enumerate(poly.vertices):
            key = (int(math.floor(p[0] / distance)), int(math.floor(p[1] / distance)),
                   int(math.floor(p[2] / distance)))
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        for j in grid.get((key[0] + dx, key[1] + dy, key[2] + dz), ()):
                            if vec.length_sq(vec.sub(p, poly.vertices[j])) <= limit:
                                uf.union(i, j)
            grid.setdefault(key, []).append(i)
    groups = uf.groups()
    if len(groups) == count:
        return Polyhedron(poly.vertices, poly.faces, poly.labels, poly.holes)
    target = [0] * count
    vertices = []
    for g in groups:
        for i in g:
            target[i] = len(vertices)
        if len(g) == 1:
            vertices.append(poly.vertices[g[0]])
        elif place is not None:
            vertices.append(place(g))
        else:
            vertices.append(vec.mean([poly.vertices[i] for i in g]))

    def mapped(loop):
        g = []
        for i in loop:
            t = target[i]
            if not g or g[-1] != t:
                g.append(t)
        while len(g) > 1 and g[0] == g[-1]:
            g.pop()
        return g if len(g) >= 3 else None

    faces, labels, holes = [], [], []
    for k, lab in enumerate(poly.labels):
        outer = mapped(poly.faces[k])
        if outer is None:
            continue
        faces.append(outer)
        labels.append(lab)
        holes.append([h for h in (mapped(loop) for loop in poly.holes[k]) if h is not None])
    return _compact(vertices, faces, labels, holes)


def closed_problems(poly):
    """Why poly is not a closed, consistently wound solid (empty list when it is)."""
    problems = []
    if len(poly.faces) < 4:
        problems.append("fewer than 4 faces")
    edges = {}
    pinched = 0
    for k in range(len(poly.faces)):
        corners = [i for loop in poly.loops(k) for i in loop]
        if len(set(corners)) != len(corners):
            pinched += 1
        for loop in poly.loops(k):
            for j in range(len(loop)):
                e = (loop[j], loop[(j + 1) % len(loop)])
                edges[e] = edges.get(e, 0) + 1
    if pinched:
        problems.append("%d %s pinched (a corner used twice)" % (pinched, _faces_word(pinched)))
    doubled = sum(1 for c in edges.values() if c > 1)
    open_ = sum(1 for (a, b) in edges if (b, a) not in edges)
    if doubled:
        problems.append("%d edges run the same way twice" % doubled)
    if open_:
        problems.append("%d open edges" % open_)
    return problems


def side_share(normal, axis):
    """The share of a unit normal that lies square to the axis (1 without an axis)."""
    if axis is None:
        return 1.0
    c = vec.dot(normal, axis)
    return math.sqrt(max(0.0, 1.0 - c * c))


class Extent(object):
    """A part's points, set up to answer "how far out does it reach along n?".

    Built once for a set of directions (the facet normals). For each one it
    keeps the support (the largest dot(n, p)) and the points near it. A later
    query along a direction close to one of those only needs the near points:
    a point more than `band` below the support cannot become the furthest
    one while the direction turns by less than band / (2 x radius).
    """

    def __init__(self, points, normals, band=None):
        self.points = list(points)
        self.centre = vec.mean(self.points)
        self.radius = max(vec.distance(p, self.centre) for p in self.points)
        if band is None:
            band = NEAR_RATIO * bounding_size(self.points)
        self.band = band
        self.normals = list(normals)
        self.supports = []
        self.near = []
        for n in self.normals:
            dots = [vec.dot(n, p) for p in self.points]
            s = max(dots)
            self.supports.append(s)
            self.near.append([self.points[i] for i, x in enumerate(dots) if x >= s - band])

    def reach(self, normal, hint=None):
        """The largest dot(normal, p) over all points. hint: index of a direction
        from the set close to normal, to search only its near points."""
        pool = self.points
        if hint is not None:
            turn = vec.length(vec.sub(normal, self.normals[hint]))
            if 2.0 * turn * self.radius < self.band:
                pool = self.near[hint]
        return max(vec.dot(normal, p) for p in pool)

    def planes(self, clearance=0.0, axis=None):
        """One (normal, offset) per direction: touching the part, then pushed out
        by the clearance (only square to the axis when one is given)."""
        return [(n, s + clearance * side_share(n, axis)) for n, s in zip(self.normals, self.supports)]


def support_planes(points, normals, clearance=0.0, axis=None):
    """Planes with the given unit normals that touch the points from outside,
    then pushed out by the clearance. With an axis, a plane moves by the
    clearance times the sideways share of its normal. Returns (normal, offset) pairs."""
    return Extent(points, normals).planes(clearance, axis)


def validate(poly, flat_tol, extent=None, clearance=0.0, axis=None, keep=KEEP_CLEARANCE, slack=0.0):
    """Check a solid and refit its face planes. Returns (problems, planes).

    problems  readable reasons it cannot be used (empty when it passes):
              not closed, a face not flat within flat_tol, or (with an
              extent) a face closer to the part than keep x its clearance
              minus slack (with no clearance: crossing the part by more
              than slack)
    planes    per face (unit normal, offset), refitted through its corners
    """
    problems = closed_problems(poly)
    planes = []
    worst_flat = 0.0
    tight = 0
    least = None
    for k in range(len(poly.faces)):
        fitted = poly.plane(k)
        if fitted is None:
            problems.append("a face has no area")
            planes.append(None)
            continue
        n, d, flat = fitted
        planes.append((n, d))
        worst_flat = max(worst_flat, flat)
        if extent is None:
            continue
        lab = poly.labels[k]
        hint = lab if lab is not None and lab < len(extent.normals) else None
        gap = d - extent.reach(n, hint)
        need = clearance * side_share(n, axis)
        if gap < keep * need - slack:
            tight += 1
            short = gap - need
            if least is None or short < least:
                least = short
    if worst_flat > flat_tol:
        problems.append("a face is %.3g off flat (allowed %.3g)" % (worst_flat, flat_tol))
    if tight:
        problems.append("%d %s closer to the part than the clearance allows (by up to %.3g)" % (
            tight, _faces_word(tight), -least))
    return problems, planes


def _faces_word(count):
    return "face" if count == 1 else "faces"


# ---------------------------------------------------------------- the offset solid

class OffsetResult(object):
    """solid: Polyhedron or None. planes: per face (normal, offset). merged: the
    corner merge distance used (None when only corners within the tolerance were joined)."""

    def __init__(self):
        self.solid = None
        self.planes = []
        self.merged = None
        self.facets = []
        self.report = Report("offset solid")

    @property
    def status(self):
        return self.report.status


def offset_solid(vertices, faces, clearance=0.0, axis=None, tolerance=0.001, facets=None,
                 max_planes=MAX_PLANES):
    """A closed convex solid around a faceted part, every facet pushed out by the clearance.

    vertices, faces  the part (a cleaned triangle mesh)
    clearance        how far the solid stays off the part (0: it touches)
    axis             when given, the clearance is applied square to it only
                     (see the module docs); otherwise in every direction
    tolerance        modelling tolerance: no edge ends up shorter, faces are
                     flat within half of it
    facets           facets from facets.cluster_planes(), found when omitted
    Returns an OffsetResult. See the module docs for how it is built.
    """
    res = OffsetResult()
    r = res.report
    if clearance < 0.0:
        r.fail("the clearance must be zero or more, got %.4g" % clearance)
        return res
    if tolerance <= 0.0:
        r.fail("the tolerance must be more than zero, got %.4g" % tolerance)
        return res
    if axis is not None:
        axis = vec.unit(axis)
        if axis is None:
            r.fail("the axis is a zero vector")
            return res
    if not faces:
        r.fail("the mesh has no faces")
        return res
    used = sorted(set(i for f in faces for i in f))
    points = [vertices[i] for i in used]
    size = bounding_size(points)
    if size <= 0.0:
        r.fail("the mesh has no extent")
        return res
    if facets is None:
        facets = cluster_planes(vertices, faces)
    res.facets = facets
    if not facets:
        r.fail("no flat facets found")
        return res
    if len(facets) > max_planes:
        r.fail("%d facets, more than %d: the part looks curved rather than faceted" % (
            len(facets), max_planes))
        return res

    extent = Extent(points, [fc.normal for fc in facets])
    _grade_facets(facets, vertices, faces, extent, tolerance, r)

    planes = extent.planes(clearance, axis)
    # Sideways planes first, then by area: the order only affects rounding, but keeps it repeatable.
    order = sorted(range(len(planes)), key=lambda k: (-side_share(planes[k][0], axis), -facets[k].area, k))
    lo = tuple(min(p[c] for p in points) for c in range(3))
    hi = tuple(max(p[c] for p in points) for c in range(3))
    margin = size + 10.0 * clearance
    poly, problem = intersect_planes(planes, vec.sub(lo, (margin, margin, margin)),
                                     vec.add(hi, (margin, margin, margin)), EPS_RATIO * size, order)
    if poly is None:
        r.fail("the facet planes do not make a solid: " + problem)
        return res
    if any(lab is None for lab in poly.labels):
        r.fail("the facets leave the solid open in some direction (a missing facet or a hole)")
        return res

    solid, chosen, last = None, None, None
    checks = (FLAT_SHARE * tolerance, extent, clearance, axis, KEEP_CLEARANCE, SLACK_SHARE * tolerance)
    steps = [s * clearance for s in MERGE_STEPS if s * clearance > tolerance] + [tolerance]
    for distance in steps:
        merged = merge_corners(poly, distance)
        problems, fitted = validate(merged, *checks)
        if not problems:
            solid, chosen, res.planes = merged, distance, fitted
            break
        last = (distance, problems)
    if solid is None:
        problems, fitted = validate(poly, *checks)
        if problems:
            r.fail("the solid fails its own checks: " + "; ".join(problems))
            return res
        solid, res.planes = poly, fitted
        r.warn("short edges kept: merging their corners within %.3g would break the solid (%s)" % (
            last[0], last[1][0]))
    res.solid = solid
    res.merged = chosen if chosen is not None and chosen > tolerance else None

    shortest = min(solid.edge_lengths())
    if shortest < tolerance:
        r.warn("shortest edge %.3g is below the tolerance %.3g" % (shortest, tolerance))
    if res.merged is not None:
        r.note("split corners merged within %.3g" % res.merged)
    reach = _reach_beyond(solid, lo, hi)
    allowed = REACH_CLEARANCES * clearance + REACH_RATIO * size
    if reach > allowed:
        r.messy("a corner reaches %.3g beyond the part (expected up to %.3g): a facet may be missing" % (
            reach, allowed))
    unused = len(planes) - len(set(solid.labels))
    r.note("%d flat faces from %d facets%s" % (
        len(solid.faces), len(facets), ", %d facets made no face" % unused if unused else ""))
    r.note("volume %.6g" % solid.volume())
    return res


def _grade_facets(facets, vertices, faces, extent, tolerance, report):
    """Report facet planes that the part crosses (not convex) or that had to move
    out to touch it (facets not quite flat)."""
    crossed, cross_by, moved = 0, 0.0, 0.0
    for k, fc in enumerate(facets):
        own = set(i for fi in fc.faces for i in faces[fi])
        own_reach = max(vec.dot(fc.normal, vertices[i]) for i in own)
        moved = max(moved, own_reach - fc.offset)
        beyond = extent.supports[k] - own_reach
        if beyond > max(tolerance, CONCAVE_RMS * fc.rms):
            crossed += 1
            cross_by = max(cross_by, beyond)
    if crossed:
        report.messy("the part is not convex: %d facet %s %s through it (by up to %.3g), so the "
                     "solid bridges its hollows" % (
                         crossed, "plane" if crossed == 1 else "planes",
                         "passes" if crossed == 1 else "pass", cross_by))
    if moved > tolerance:
        report.warn("facet planes moved out by up to %.3g to touch the part (facets not quite flat)" % moved)


def _reach_beyond(poly, lo, hi):
    """How far the furthest corner lies outside the box lo..hi."""
    worst = 0.0
    for p in poly.vertices:
        for c in range(3):
            worst = max(worst, lo[c] - p[c], p[c] - hi[c])
    return worst
