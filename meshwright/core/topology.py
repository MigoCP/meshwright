"""Mesh topology: welding, edge census, pieces, orientation and a graded health check.

Meshes are plain lists: vertices as (x, y, z) tuples and faces as (i, j, k)
index triples. Quads and n-gons are split into triangles by the adapters
before they reach this module.

The headline function is health(): it cleans what can be cleaned safely
(coincident vertices, zero-area and duplicate faces, loose junk, faces pointing
the wrong way) and grades what cannot (holes, extra pieces, non-manifold edges),
returning the cleaned mesh together with a Report.
"""

from __future__ import absolute_import, division

import math

from . import vec
from .report import Report

# Grading thresholds. Sizes are relative so results do not depend on units.
JUNK_RATIO = 0.01          # a piece below this share of the largest piece's area is junk
SMALL_HOLE_AREA = 0.02     # a hole whose filling patch is below this share of the area is small
FEW_NON_MANIFOLD = 0.01    # non-manifold edges up to this share of all edges are "a few"
DEGENERATE_AREA = 1e-12    # face area below this x size^2 counts as zero


# ---------------------------------------------------------------- union-find

class UnionFind(object):
    """Disjoint sets over 0..n-1 with path halving. Roots are the smallest index."""

    def __init__(self, n):
        self.parent = list(range(n))

    def find(self, i):
        parent = self.parent
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return False
        if ra < rb:
            self.parent[rb] = ra
        else:
            self.parent[ra] = rb
        return True

    def groups(self):
        """Lists of members, one per set, ordered by their smallest member."""
        out = {}
        for i in range(len(self.parent)):
            out.setdefault(self.find(i), []).append(i)
        return [out[k] for k in sorted(out)]


# ---------------------------------------------------------------- measures

def bounding_size(vertices):
    """Diagonal of the axis-aligned bounding box (0 for no vertices)."""
    if not vertices:
        return 0.0
    lo = [min(v[k] for v in vertices) for k in range(3)]
    hi = [max(v[k] for v in vertices) for k in range(3)]
    return vec.length(vec.sub(hi, lo))


def face_area(vertices, face):
    a, b, c = vertices[face[0]], vertices[face[1]], vertices[face[2]]
    return 0.5 * vec.length(vec.cross(vec.sub(b, a), vec.sub(c, a)))


def face_normal(vertices, face):
    """Unit normal following the face's winding, or None for a zero-area face."""
    a, b, c = vertices[face[0]], vertices[face[1]], vertices[face[2]]
    return vec.unit(vec.cross(vec.sub(b, a), vec.sub(c, a)))


def surface_area(vertices, faces):
    return sum(face_area(vertices, f) for f in faces)


def signed_volume(vertices, faces):
    """Volume enclosed by a closed mesh; positive when faces point outward.

    Each face adds the signed volume of the tetrahedron it forms with a
    reference point. Using the mesh centre as that point (instead of the
    origin) keeps the sum accurate for meshes far from the origin.
    """
    if not faces:
        return 0.0
    ref = vec.mean([vertices[i] for f in faces for i in f])
    total = 0.0
    for f in faces:
        a = vec.sub(vertices[f[0]], ref)
        b = vec.sub(vertices[f[1]], ref)
        c = vec.sub(vertices[f[2]], ref)
        total += vec.dot(a, vec.cross(b, c))
    return total / 6.0


# ---------------------------------------------------------------- welding

def weld(vertices, faces, tolerance=0.0):
    """Merge vertices closer than tolerance (exactly equal ones when 0).

    Merging is transitive (a chain of close vertices becomes one). Each
    merged group keeps the position of its lowest-index vertex, so welding
    never moves geometry by more than the chain length. Unused vertices are
    dropped. Returns (vertices, faces, remap) where remap[old] = new index.
    """
    n = len(vertices)
    uf = UnionFind(n)
    if tolerance <= 0.0:
        first = {}
        for i, p in enumerate(vertices):
            key = (p[0], p[1], p[2])
            if key in first:
                uf.union(first[key], i)
            else:
                first[key] = i
    else:
        cell = float(tolerance)
        tol_sq = cell * cell
        grid = {}
        for i, p in enumerate(vertices):
            cx = int(math.floor(p[0] / cell))
            cy = int(math.floor(p[1] / cell))
            cz = int(math.floor(p[2] / cell))
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        for j in grid.get((cx + dx, cy + dy, cz + dz), ()):
                            if vec.length_sq(vec.sub(p, vertices[j])) <= tol_sq:
                                uf.union(i, j)
            grid.setdefault((cx, cy, cz), []).append(i)
    # New indices go to the roots of used vertices, in index order.
    used = set(i for f in faces for i in f)
    remap = [-1] * n
    out_vertices = []
    root_new = {}
    for i in sorted(used):
        root = uf.find(i)
        if root not in root_new:
            root_new[root] = len(out_vertices)
            out_vertices.append(tuple(float(c) for c in vertices[root]))
    for i in range(n):
        remap[i] = root_new.get(uf.find(i), -1)
    out_faces = [tuple(remap[i] for i in f) for f in faces]
    return out_vertices, out_faces, remap


# ---------------------------------------------------------------- edges and pieces

class EdgeCensus(object):
    """How every undirected edge is used by the faces.

    uses[(a, b)] (a < b) lists (face index, True when the face runs a to b).
    naked edges are used once, non_manifold more than twice, and flipped
    edges twice in the same direction (two neighbours that disagree on
    which side is out).
    """

    def __init__(self, faces):
        uses = {}
        for fi, f in enumerate(faces):
            for k in range(3):
                a, b = f[k], f[(k + 1) % 3]
                if a == b:
                    continue
                key = (a, b) if a < b else (b, a)
                uses.setdefault(key, []).append((fi, a < b))
        self.uses = uses
        self.naked = sorted(k for k, u in uses.items() if len(u) == 1)
        self.non_manifold = sorted(k for k, u in uses.items() if len(u) > 2)
        self.flipped = sorted(k for k, u in uses.items() if len(u) == 2 and u[0][1] == u[1][1])

    @property
    def edge_count(self):
        return len(self.uses)

    @property
    def closed(self):
        return not self.naked and not self.non_manifold


def pieces(faces, census=None):
    """Groups of face indices connected through shared EDGES.

    Faces that only touch at a vertex belong to different pieces: two solids
    meeting at a single point are still two solids.
    """
    if census is None:
        census = EdgeCensus(faces)
    uf = UnionFind(len(faces))
    for u in census.uses.values():
        for k in range(1, len(u)):
            uf.union(u[0][0], u[k][0])
    return uf.groups()


def orient(vertices, faces, census=None):
    """Make face windings agree across shared edges, then point closed pieces outward.

    Walks each piece across its two-face edges, reversing any neighbour that
    runs the shared edge in the same direction. A piece that cannot be made
    consistent (a Moebius-like surface) is counted as a conflict. Afterwards a
    closed piece with negative volume is reversed as a whole; an open piece
    keeps the winding most of its faces already had.

    Returns (faces, changed, conflicts) where changed is the number of faces
    whose winding differs from the input.
    """
    if census is None:
        census = EdgeCensus(faces)
    faces = [tuple(f) for f in faces]
    flip = [None] * len(faces)  # None: not visited; else bool "reverse this face"
    conflicts = 0
    neighbours = [[] for _ in faces]
    for u in census.uses.values():
        if len(u) == 2:
            (f0, d0), (f1, d1) = u
            neighbours[f0].append((f1, d0 == d1))
            neighbours[f1].append((f0, d0 == d1))
    groups = pieces(faces, census)
    for group in groups:
        seed = group[0]
        flip[seed] = False
        stack = [seed]
        while stack:
            f = stack.pop()
            for g, same_direction in neighbours[f]:
                # Neighbours must run the shared edge in opposite directions.
                want = flip[f] if not same_direction else (not flip[f])
                if flip[g] is None:
                    flip[g] = want
                    stack.append(g)
                elif flip[g] != want:
                    conflicts += 1
        # Choose the overall side of this piece.
        piece_faces = [_wound(faces[i], flip[i]) for i in group]
        piece_census = EdgeCensus(piece_faces)
        if piece_census.closed:
            reverse_all = signed_volume(vertices, piece_faces) < 0.0
        else:
            reversed_count = sum(1 for i in group if flip[i])
            reverse_all = reversed_count * 2 > len(group)
        if reverse_all:
            for i in group:
                flip[i] = not flip[i]
    out = [_wound(faces[i], flip[i]) for i in range(len(faces))]
    changed = sum(1 for i in range(len(faces)) if flip[i])
    return out, changed, conflicts // 2


def _wound(face, reverse):
    return (face[0], face[2], face[1]) if reverse else face


def boundary_loops(faces, census=None):
    """Rims of the holes: chains of naked edges as vertex index lists.

    Each naked edge is walked in its face's direction, so loops come out
    consistently oriented. A loop that cannot be closed (a rim touching
    itself at a non-manifold vertex) is returned as an open chain.
    """
    if census is None:
        census = EdgeCensus(faces)
    nexts = {}
    for key in census.naked:
        fi, forward = census.uses[key][0]
        a, b = key if forward else (key[1], key[0])
        nexts.setdefault(a, []).append(b)
    loops = []
    while nexts:
        start = min(nexts)
        chain = [start]
        current = start
        while current in nexts:
            options = nexts[current]
            nxt = options.pop()
            if not options:
                del nexts[current]
            if nxt == start:
                break
            chain.append(nxt)
            current = nxt
        loops.append(chain)
    return loops


# ---------------------------------------------------------------- health check

class Health(object):
    """Result of health(): the cleaned mesh, its Report and what was found."""

    def __init__(self):
        self.vertices = []
        self.faces = []
        self.report = Report("mesh health")
        self.stats = {}

    @property
    def status(self):
        return self.report.status

    @property
    def closed(self):
        return bool(self.faces) and self.stats.get("naked_edges", 1) == 0 \
            and self.stats.get("non_manifold_edges", 1) == 0


def health(vertices, faces, tolerance=0.0, junk_ratio=JUNK_RATIO):
    """Check and clean a triangle mesh. Never raises on bad geometry.

    Cleans: vertices closer than tolerance (exact duplicates when 0),
    zero-area faces, duplicate faces, loose pieces smaller than junk_ratio of
    the largest piece, and face windings (made consistent, closed pieces
    pointing out).

    Grades:
      0  closed, one piece, nothing found worth mentioning
      1  warnings: something was cleaned, small holes, a few non-manifold edges
      2  messy: several real pieces, big holes, many non-manifold edges, or
         windings that cannot be made consistent
      3  failed: no faces, invalid indices or non-finite coordinates
    """
    h = Health()
    r = h.report
    st = h.stats
    st["faces_in"] = len(faces)
    st["vertices_in"] = len(vertices)

    # Things that make the input unusable.
    if not faces:
        r.fail("the mesh has no faces")
        return h
    n = len(vertices)
    for f in faces:
        if len(f) != 3 or any((i < 0 or i >= n) for i in f):
            r.fail("faces refer to vertices that do not exist")
            return h
    for p in vertices:
        if any(not _finite(c) for c in p):
            r.fail("some vertex coordinates are not finite numbers")
            return h

    # Weld.
    verts, tris, _ = weld(vertices, faces, tolerance)
    st["welded_vertices"] = len(set(i for f in faces for i in f)) - len(verts)
    size = bounding_size(verts)
    st["size"] = size

    # Zero-area and duplicate faces.
    min_area = DEGENERATE_AREA * size * size
    kept, degenerate = [], 0
    for f in tris:
        if f[0] == f[1] or f[1] == f[2] or f[0] == f[2] or face_area(verts, f) <= min_area:
            degenerate += 1
        else:
            kept.append(f)
    seen, unique, duplicates = set(), [], 0
    for f in kept:
        key = tuple(sorted(f))
        if key in seen:
            duplicates += 1
        else:
            seen.add(key)
            unique.append(f)
    tris = unique
    st["degenerate_faces"] = degenerate
    st["duplicate_faces"] = duplicates
    if not tris:
        r.fail("every face has zero area")
        return h
    if degenerate:
        share = degenerate / float(st["faces_in"])
        text = "removed %d zero-area face%s" % (degenerate, _s(degenerate))
        if share <= 0.01:
            r.note(text)
        else:
            r.warn(text + " (%.0f%% of the faces)" % (100 * share))
    if duplicates:
        r.warn("removed %d duplicate face%s" % (duplicates, _s(duplicates)))

    # Pieces and junk.
    census = EdgeCensus(tris)
    groups = pieces(tris, census)
    areas = [sum(face_area(verts, tris[i]) for i in g) for g in groups]
    largest = max(areas)
    keep_groups, dropped, dropped_area = [], 0, 0.0
    for g, a in zip(groups, areas):
        if a < junk_ratio * largest:
            dropped += 1
            dropped_area += a
        else:
            keep_groups.append(g)
    st["pieces_in"] = len(groups)
    st["junk_dropped"] = dropped
    if dropped:
        r.warn("dropped %d loose bit%s (%.2g%% of the area)" % (
            dropped, _s(dropped), 100.0 * dropped_area / sum(areas)))
    tris = [tris[i] for g in keep_groups for i in g]
    st["pieces"] = len(keep_groups)
    if len(keep_groups) > 1:
        r.messy("%d separate pieces" % len(keep_groups))

    # Orientation.
    census = EdgeCensus(tris)
    st["flipped_edges_in"] = len(census.flipped)
    tris, changed, conflicts = orient(verts, tris, census)
    st["faces_reoriented"] = changed
    st["orientation_conflicts"] = conflicts
    if conflicts:
        r.messy("face directions cannot be made consistent (%d conflicting edge%s)" % (
            conflicts, _s(conflicts)))
    if changed:
        if changed == len(tris):
            r.warn("all faces pointed inward; reversed them")
        else:
            r.warn("reversed %d face%s that disagreed with their neighbours" % (changed, _s(changed)))

    # Holes and non-manifold edges, on the cleaned faces.
    census = EdgeCensus(tris)
    st["edges"] = census.edge_count
    st["naked_edges"] = len(census.naked)
    st["non_manifold_edges"] = len(census.non_manifold)
    loops = boundary_loops(tris, census)
    # A hole's size is the area of the patch that would fill it, compared with
    # the whole surface: one missing sliver is small, a missing end cap is not.
    patches = [vec.length(vec.area_vector([verts[i] for i in loop])) for loop in loops]
    area = surface_area(verts, tris)
    st["holes"] = len(loops)
    st["largest_hole_share"] = max(patches) / area if patches else 0.0
    if loops:
        if st["largest_hole_share"] > SMALL_HOLE_AREA:
            r.messy("%d hole%s, the largest about %.0f%% of the surface" % (
                len(loops), _s(len(loops)), 100.0 * st["largest_hole_share"]))
        else:
            r.warn("%d small hole%s" % (len(loops), _s(len(loops))))
    if census.non_manifold:
        share = len(census.non_manifold) / float(census.edge_count)
        text = "%d non-manifold edge%s (used by more than two faces)" % (
            len(census.non_manifold), _s(len(census.non_manifold)))
        if share <= FEW_NON_MANIFOLD:
            r.warn(text)
        else:
            r.messy(text)

    # Output: compacted vertices and the summary numbers.
    h.vertices, h.faces = _compact(verts, tris)
    st["faces"] = len(h.faces)
    st["area"] = surface_area(h.vertices, h.faces)
    st["volume"] = signed_volume(h.vertices, h.faces) if census.closed else None
    return h


def _compact(vertices, faces):
    used = sorted(set(i for f in faces for i in f))
    new = dict((old, k) for k, old in enumerate(used))
    return [vertices[i] for i in used], [tuple(new[i] for i in f) for f in faces]


def _finite(x):
    return not (math.isnan(x) or math.isinf(x))


def _s(count):
    return "" if count == 1 else "s"
