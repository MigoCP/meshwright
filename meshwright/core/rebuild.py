"""Planar mesh to exact solid: one flat face per planar region, corners solved from planes.

A mesh of a part made of flat faces (an STL export, a scan of a machined
part) holds each face as many triangles that are never quite flat with each
other. rebuild() turns it into a Polyhedron with one exact flat face per
region. It keeps the mesh's own topology instead of intersecting
half-spaces, so it works for parts that are not convex, with holes too.

Steps
    1. Regions: triangles grouped into flat facets (facets.cluster_planes),
       each facet split into edge-connected pieces, so two coplanar faces
       that do not touch (the ends of a U's arms) stay separate faces.
       Neighbouring regions that together are flat within half the
       tolerance are joined, flattest pair first, which gathers up the
       fragments noise breaks a face into. Every region gets a plane
       refitted through its own vertices.
    2. Loops: each region's boundary, chained from the mesh edges it shares
       with other regions: one outer loop, plus one per hole.
    3. Corners: a mesh vertex where three or more regions meet is a corner.
       Its position is solved from the planes of those regions (least
       squares when more than three). A vertex where only two regions meet
       lies on their crease line between two corners, so it is dropped.
    4. Merge: corners closer than the tolerance are joined (union-find, then
       solved again from all their planes). The merge is kept only when the
       solid stays closed and every face stays flat within half the
       tolerance.

Corner solve
    The corner is the point closest to all its planes in the least squares
    sense. In a direction the planes barely pin down (two of them nearly
    parallel), it keeps the mesh vertex's own position instead of jumping
    far along a nearly parallel crease.
"""

from __future__ import absolute_import, division

import heapq
import math

from . import vec
from .facets import ANGLE_DEG, cluster_planes
from .polyhedron import FLAT_SHARE, MAX_PLANES, Polyhedron, merge_corners, validate
from .report import Report
from .topology import EdgeCensus, UnionFind, signed_volume, surface_area

SOLVE_CUTOFF = 1e-3     # corner solve: directions weaker than this (unit-normal scale) keep the vertex


class RebuildResult(object):
    """solid: Polyhedron or None (face labels are region indices). regions: per
    region its triangle indices. planes: per region (unit normal, offset)."""

    def __init__(self):
        self.solid = None
        self.regions = []
        self.planes = []
        self.merged = False
        self.report = Report("rebuild")

    @property
    def status(self):
        return self.report.status


def solve_corner(planes, start, cutoff=SOLVE_CUTOFF):
    """The point closest to the planes (unit normal, offset), least squares.

    Directions in which the planes give less than `cutoff` of grip (an
    eigenvalue of sum(n n^T) below it) keep the start point's position.
    """
    m = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    rhs = [0.0, 0.0, 0.0]
    for n, d in planes:
        for i in range(3):
            rhs[i] += n[i] * d
            for j in range(3):
                m[i][j] += n[i] * n[j]
    # Solve for the step from the start point, one eigen direction at a time.
    residual = tuple(rhs[i] - sum(m[i][j] * start[j] for j in range(3)) for i in range(3))
    values, vectors = vec.eig_sym3(m)
    point = start
    for value, v in zip(values, vectors):
        if value > cutoff:
            point = vec.add(point, vec.scale(v, vec.dot(v, residual) / value))
    return point


def rebuild(vertices, faces, tolerance=0.001, facets=None):
    """A Polyhedron with one exact flat face per planar region of a closed mesh.

    vertices, faces  a cleaned, closed, outward-facing triangle mesh
    tolerance        modelling tolerance: corners closer than this are merged,
                     faces must stay flat within half of it
    facets           facets from facets.cluster_planes(), found when omitted
    Returns a RebuildResult. See the module docs for the steps.
    """
    res = RebuildResult()
    r = res.report
    if tolerance <= 0.0:
        r.fail("the tolerance must be more than zero, got %.4g" % tolerance)
        return res
    if not faces:
        r.fail("the mesh has no faces")
        return res
    census = EdgeCensus(faces)
    if census.naked or census.non_manifold:
        r.fail("the mesh is not closed (%d open, %d non-manifold edges); clean it first" % (
            len(census.naked), len(census.non_manifold)))
        return res
    if census.flipped:
        r.fail("%d edges where neighbouring faces disagree on which side is out" % len(census.flipped))
        return res
    mesh_volume = signed_volume(vertices, faces)
    if mesh_volume <= 0.0:
        r.fail("the mesh is inside out or flat (volume %.4g)" % mesh_volume)
        return res
    if facets is None:
        facets = cluster_planes(vertices, faces)

    owner = _regions(vertices, faces, facets, census, FLAT_SHARE * tolerance)
    count = max(owner) + 1
    regions = [[] for _ in range(count)]
    for t, g in enumerate(owner):
        regions[g].append(t)
    planes = [_fit(vertices, faces, tris) for tris in regions]
    res.regions, res.planes = regions, planes
    if count < 4:
        r.fail("only %d flat regions found; a closed solid needs at least 4" % count)
        return res
    if count > MAX_PLANES:
        r.warn("%d flat faces: the part looks curved rather than faceted" % count)

    shallow = _shallow_pairs(census, owner, planes)
    if shallow:
        r.warn("%d pairs of neighbouring faces meet at less than %.3g degrees but are not flat "
               "together within the tolerance: the mesh is likely noisier than the tolerance "
               "(try a larger one)" % (shallow, ANGLE_DEG))

    loops, problem = _loops(vertices, faces, owner, regions, planes)
    if problem:
        r.fail(problem)
        return res

    # Which regions meet at each vertex; corners are where three or more do.
    meets = [set() for _ in vertices]
    for t, f in enumerate(faces):
        for i in f:
            meets[i].add(owner[t])
    keep = set(i for i, m in enumerate(meets) if len(m) >= 3)
    for outer, holes in loops:
        for loop in [outer] + holes:
            if sum(1 for i in loop if i in keep) < 3:
                keep.update(loop)  # a loop with too few corners keeps all its vertices

    index = {}
    corners, corner_planes = [], []
    for i in sorted(keep):
        index[i] = len(corners)
        ps = [planes[g] for g in sorted(meets[i])]
        corner_planes.append(ps)
        corners.append(solve_corner(ps, vertices[i]))

    poly_faces, poly_holes, labels = [], [], []
    for g, (outer, holes) in enumerate(loops):
        poly_faces.append([index[i] for i in outer if i in keep])
        poly_holes.append([[index[i] for i in h if i in keep] for h in holes])
        labels.append(g)
    exact = Polyhedron(corners, poly_faces, labels, poly_holes)

    flat_tol = FLAT_SHARE * tolerance

    def place(group):
        ps = []
        for k in group:
            for pl in corner_planes[k]:
                if pl not in ps:
                    ps.append(pl)
        return solve_corner(ps, vec.mean([corners[k] for k in group]))

    merged = merge_corners(exact, tolerance, place)
    merge_problems = None  # None: nothing to merge
    if len(merged.vertices) < len(exact.vertices):
        merge_problems, _ = validate(merged, flat_tol)
    if merge_problems == []:
        solid = merged
        res.merged = True
        r.note("%d corners closer than the tolerance merged" % (len(exact.vertices) - len(merged.vertices)))
    else:
        problems, _ = validate(exact, flat_tol)
        if problems:
            r.fail("the rebuilt solid fails its own checks: " + "; ".join(problems))
            return res
        solid = exact
        if merge_problems:
            r.warn("short edges kept: merging corners within %.3g would break the solid (%s)" % (
                tolerance, merge_problems[0]))
    res.solid = solid

    _grade(res, vertices, faces, owner, mesh_volume, tolerance)
    return res


# ---------------------------------------------------------------- regions

def _regions(vertices, faces, facets, census, flat_tol):
    """Region index per triangle: facets split into edge-connected pieces, then
    neighbouring pieces joined while together they stay flat within flat_tol.

    The joining is greedy, flattest pair first, with the plane refitted after
    every join, so the fragments a noisy face falls into grow back into one.
    """
    facet_of = [-1] * len(faces)
    for k, fc in enumerate(facets):
        for t in fc.faces:
            facet_of[t] = k
    next_id = len(facets)
    for t in range(len(faces)):
        if facet_of[t] < 0:  # zero-area triangles are left out of facets
            facet_of[t] = next_id
            next_id += 1
    pairs = []
    for uses in census.uses.values():
        pairs.append((uses[0][0], uses[1][0]))
    uf = UnionFind(len(faces))
    for t0, t1 in pairs:
        if facet_of[t0] == facet_of[t1]:
            uf.union(t0, t1)

    corners = {}     # region root -> set of vertex indices
    for t, f in enumerate(faces):
        corners.setdefault(uf.find(t), set()).update(f)
    near = {}        # region root -> set of neighbouring region roots
    for t0, t1 in pairs:
        a, b = uf.find(t0), uf.find(t1)
        if a != b:
            near.setdefault(a, set()).add(b)
            near.setdefault(b, set()).add(a)
    version = dict((root, 0) for root in corners)

    def flatness(a, b):
        points = [vertices[i] for i in corners[a] | corners[b]]
        _, n, _, _ = vec.plane_fit(points)
        ds = [vec.dot(n, p) for p in points]
        mid = sum(ds) / len(ds)
        return max(abs(x - mid) for x in ds)

    heap = []
    for a in near:
        for b in near[a]:
            if a < b:
                cost = flatness(a, b)
                if cost <= flat_tol:
                    heapq.heappush(heap, (cost, a, b, 0, 0))
    while heap:
        cost, a, b, va, vb = heapq.heappop(heap)
        if version.get(a) != va or version.get(b) != vb:
            continue  # one of them has changed since this pair was measured
        uf.union(a, b)
        root = uf.find(a)
        gone = b if root == a else a
        corners[root] = corners[a] | corners[b]
        near[root] = (near.get(a, set()) | near.get(b, set())) - set([a, b])
        for other in near[root]:
            near[other].discard(gone)
            near[other].add(root)
        del corners[gone]
        near.pop(gone, None)
        del version[gone]
        version[root] += 1
        for other in near[root]:
            c = flatness(root, other)
            if c <= flat_tol:
                lo, hi = (root, other) if root < other else (other, root)
                heapq.heappush(heap, (c, lo, hi, version[lo], version[hi]))

    ids = {}
    owner = []
    for t in range(len(faces)):
        root = uf.find(t)
        if root not in ids:
            ids[root] = len(ids)
        owner.append(ids[root])
    return owner


def _shallow_pairs(census, owner, planes):
    """Neighbouring regions whose planes are less than ANGLE_DEG apart."""
    cos_limit = math.cos(math.radians(ANGLE_DEG))
    seen = set()
    for uses in census.uses.values():
        a, b = owner[uses[0][0]], owner[uses[1][0]]
        if a != b:
            key = (a, b) if a < b else (b, a)
            if key not in seen and vec.dot(planes[a][0], planes[b][0]) > cos_limit:
                seen.add(key)
    return len(seen)


def _fit(vertices, faces, tris):
    """(unit normal, offset) of the plane through a region's vertices, facing
    the way its triangles do. Each vertex weighs as much as the area of its
    triangles in the region, so a sliver joined to a big face barely tilts it."""
    wound = (0.0, 0.0, 0.0)
    weight = {}
    for t in tris:
        a, b, c = [vertices[i] for i in faces[t]]
        cr = vec.cross(vec.sub(b, a), vec.sub(c, a))
        wound = vec.add(wound, cr)
        share = vec.length(cr) / 6.0  # a third of the triangle's area to each corner
        for i in faces[t]:
            weight[i] = weight.get(i, 0.0) + share
    ids = sorted(weight)
    points = [vertices[i] for i in ids]
    weights = [weight[i] for i in ids]
    normal = vec.unit(wound)
    if len(points) >= 3 and sum(weights) > 0.0:
        _, fitted, _, flatness = vec.plane_fit(points, weights)
        if flatness < 0.5 or normal is None:
            normal = fitted
    else:
        weights = [1.0] * len(points)
    if normal is None:
        normal = (0.0, 0.0, 1.0)
    if vec.dot(normal, wound) < 0.0:
        normal = vec.scale(normal, -1.0)
    total = sum(weights)
    offset = sum(w * vec.dot(normal, p) for p, w in zip(points, weights)) / total
    return normal, offset


# ---------------------------------------------------------------- loops

def _loops(vertices, faces, owner, regions, planes):
    """Per region (outer loop, holes) as mesh vertex index lists, or a problem."""
    out = []
    for g, tris in enumerate(regions):
        edges = set()
        for t in tris:
            f = faces[t]
            for k in range(3):
                edges.add((f[k], f[(k + 1) % 3]))
        nexts = {}
        for a, b in edges:
            if (b, a) not in edges:
                nexts.setdefault(a, []).append(b)
        n = planes[g][0]
        u, v = vec.perpendicular_pair(n)

        def flat(i):
            p = vertices[i]
            return (vec.dot(p, u), vec.dot(p, v))

        loops = []
        while nexts:
            start = min(nexts)
            loop = [start]
            prev, current = None, start
            while True:
                options = nexts.get(current)
                if not options:
                    return None, "a flat region's boundary does not close"
                k = 0 if len(options) == 1 or prev is None else _turn(flat, prev, current, options)
                nxt = options.pop(k)
                if not options:
                    del nexts[current]
                if nxt == start:
                    break
                loop.append(nxt)
                prev, current = current, nxt
                if len(loop) > len(edges):
                    return None, "a flat region's boundary does not close"
            loops.append(loop)
        outer, holes = [], []
        for loop in loops:
            area = 0.0
            pts = [flat(i) for i in loop]
            for k in range(len(pts)):
                x0, y0 = pts[k]
                x1, y1 = pts[(k + 1) % len(pts)]
                area += x0 * y1 - x1 * y0
            (outer if area > 0.0 else holes).append(loop)
        if len(outer) != 1:
            return None, "a flat region has %d outer boundaries instead of one" % len(outer)
        out.append((outer[0], holes))
    return out, None


def _turn(flat, prev, current, options):
    """Where a region touches itself at a corner: the next boundary edge that
    keeps to the same side, the first one met turning clockwise from where
    the boundary came from."""
    c = flat(current)
    p = flat(prev)
    back = math.atan2(p[1] - c[1], p[0] - c[0])
    best, best_k = None, 0
    for k, o in enumerate(options):
        q = flat(o)
        a = (back - math.atan2(q[1] - c[1], q[0] - c[0])) % (2.0 * math.pi)
        if a <= 0.0:
            a += 2.0 * math.pi
        if best is None or a < best:
            best, best_k = a, k
    return best_k


# ---------------------------------------------------------------- grading

def _grade(res, vertices, faces, owner, mesh_volume, tolerance):
    r = res.report
    solid = res.solid
    off = 0.0
    for t, f in enumerate(faces):
        n, d = res.planes[owner[t]]
        for i in f:
            off = max(off, abs(vec.dot(n, vertices[i]) - d))
    if off > tolerance:
        r.warn("the mesh is up to %.3g off its flat faces (more than the tolerance)" % off)
    shortest = min(solid.edge_lengths())
    if shortest < tolerance:
        r.warn("shortest edge %.3g is below the tolerance %.3g" % (shortest, tolerance))
    volume = solid.volume()
    allowed = surface_area(vertices, faces) * max(off, tolerance)
    if abs(volume - mesh_volume) > allowed:
        r.messy("the solid's volume %.6g differs from the mesh's %.6g by more than its faces "
                "can explain" % (volume, mesh_volume))
    holes = sum(len(h) for h in solid.holes)
    r.note("%d flat faces%s from %d triangles" % (
        len(solid.faces), " (%d with holes)" % sum(1 for h in solid.holes if h) if holes else "", len(faces)))
    r.note("volume %.6g (mesh %.6g)" % (volume, mesh_volume))
