"""Procedural test geometry.

Every test mesh is generated in code: no files, no scanned data. A mesh here
is a TriMesh (vertex list + triangle list), the same plain format the core
works on. Damage helpers (noise, holes, flipped faces, junk...) return new
meshes so one clean shape can be broken in many controlled ways.
"""

from __future__ import division

import math

import numpy as np
from scipy.spatial import ConvexHull
from scipy.spatial.transform import Rotation


class TriMesh(object):
    """Vertices (n x 3 float array) and triangles (m x 3 int array)."""

    def __init__(self, vertices, faces):
        self.vertices = np.asarray(vertices, dtype=float).reshape(-1, 3)
        self.faces = np.asarray(faces, dtype=int).reshape(-1, 3)

    def copy(self):
        return TriMesh(self.vertices.copy(), self.faces.copy())

    def as_lists(self):
        """Plain Python lists of tuples, the format the core takes."""
        verts = [tuple(float(c) for c in v) for v in self.vertices]
        faces = [tuple(int(i) for i in f) for f in self.faces]
        return verts, faces

    def volume(self):
        """Signed volume (positive when faces point outward)."""
        v = self.vertices[self.faces]
        return float(np.einsum("ij,ij->i", v[:, 0], np.cross(v[:, 1], v[:, 2])).sum() / 6.0)

    def area(self):
        v = self.vertices[self.faces]
        return float(np.linalg.norm(np.cross(v[:, 1] - v[:, 0], v[:, 2] - v[:, 0]), axis=1).sum() / 2.0)

    def edges(self):
        """Undirected edges as a set of sorted vertex index pairs."""
        out = set()
        for a, b, c in self.faces:
            for i, j in ((a, b), (b, c), (c, a)):
                out.add((min(i, j), max(i, j)))
        return out

    def euler_characteristic(self):
        """V - E + F over the vertices actually used. 2 for a closed sphere-like solid."""
        used = len(np.unique(self.faces))
        return used - len(self.edges()) + len(self.faces)


# ---------------------------------------------------------------- clean shapes

def _outward(vertices, faces):
    """Orient every face of a convex mesh away from its centre."""
    c = vertices.mean(axis=0)
    fixed = []
    for f in faces:
        a, b, d = vertices[f[0]], vertices[f[1]], vertices[f[2]]
        n = np.cross(b - a, d - a)
        if np.dot(n, a - c) < 0:
            f = (f[0], f[2], f[1])
        fixed.append(tuple(f))
    return np.array(fixed, dtype=int)


def convex_hull(points):
    """Closed outward triangle mesh of the convex hull of points."""
    pts = np.asarray(points, dtype=float)
    hull = ConvexHull(pts)
    used = np.unique(hull.simplices)
    remap = -np.ones(len(pts), dtype=int)
    remap[used] = np.arange(len(used))
    verts = pts[used]
    faces = remap[hull.simplices]
    return TriMesh(verts, _outward(verts, faces))


def loop(z, sides, radius_x, radius_y=None, phase=0.0):
    """Points of a regular polygon loop at height z (an ellipse when radii differ)."""
    if radius_y is None:
        radius_y = radius_x
    out = []
    for k in range(sides):
        a = phase + 2.0 * math.pi * k / sides
        out.append((radius_x * math.cos(a), radius_y * math.sin(a), z))
    return out


def stacked_loops(levels):
    """Convex hull of several loops. levels: list of loop() argument tuples.

    A point (a loop of radius 0) is allowed as an end, which gives a pointed tip.
    """
    pts = []
    for args in levels:
        z, sides, rx = args[0], args[1], args[2]
        ry = args[3] if len(args) > 3 else rx
        phase = args[4] if len(args) > 4 else 0.0
        if rx == 0 and ry == 0:
            pts.append((0.0, 0.0, z))
        else:
            pts.extend(loop(z, sides, rx, ry, phase))
    return convex_hull(pts)


def faceted_part(sides=16, width=1.0):
    """A convex faceted part: flat top, widest band near the top, pointed bottom.

    Generic stand-in for any faceted convex solid with a strong flat end and an
    obvious natural axis (world Z here, before any random rotation).
    """
    r = 0.5 * width
    return stacked_loops([
        (0.16 * width, sides, 0.55 * r, 0.55 * r, math.pi / sides),  # flat top
        (0.02 * width, sides, r),                                     # widest band, top edge
        (0.0, sides, r),                                              # widest band, bottom edge
        (-0.43 * width, 1, 0.0),                                      # tip
    ])


def pyramid(base=1.0, height=1.0, apex=(0.2, 0.1)):
    """Square-based pyramid with its apex off-centre, so its four sides all slope differently."""
    h = 0.5 * base
    pts = [(-h, -h, 0.0), (h, -h, 0.0), (h, h, 0.0), (-h, h, 0.0), (apex[0], apex[1], height)]
    return convex_hull(pts)


def ball(sides=24, bands=12, radius=1.0):
    """A faceted ball: loops of latitude from pole to pole."""
    levels = [(-radius, 1, 0.0)]
    for k in range(1, bands):
        a = math.pi * k / bands
        levels.append((-radius * math.cos(a), sides, radius * math.sin(a)))
    levels.append((radius, 1, 0.0))
    return stacked_loops(levels)


def cylinder(sides=48, radius=1.0, height=2.0):
    return stacked_loops([(0.0, sides, radius), (height, sides, radius)])


def lofted(levels, sides):
    """Closed solid through loops of equal side count, bottom to top, joined by
    quads and capped by fans: unlike stacked_loops, hollows are kept.
    levels: (z, radius_x, radius_y) bottom to top."""
    verts, faces = [], []
    for z, rx, ry in levels:
        verts.extend(loop(z, sides, rx, ry))
    for k in range(len(levels) - 1):
        a, b = k * sides, (k + 1) * sides
        for i in range(sides):
            j = (i + 1) % sides
            faces.append((a + i, a + j, b + j))
            faces.append((a + i, b + j, b + i))
    bottom, top = len(verts), len(verts) + 1
    verts.append((0.0, 0.0, levels[0][0]))
    verts.append((0.0, 0.0, levels[-1][0]))
    last = (len(levels) - 1) * sides
    for i in range(sides):
        j = (i + 1) % sides
        faces.append((bottom, j, i))
        faces.append((top, last + i, last + j))
    return TriMesh(verts, faces)


def hourglass(sides=32, radius=1.0, waist=0.4, height=2.0):
    """Wide at both ends, narrow in the middle: a profile that goes in and out."""
    return lofted([(0.0, radius, radius), (0.5 * height, waist, waist), (height, radius, radius)], sides)


def wedge(length=2.0, width=1.0, depth=1.0):
    """A flat top tapering down to a straight ridge."""
    a, b = 0.5 * length, 0.5 * width
    return convex_hull([(-a, -b, 0.0), (a, -b, 0.0), (a, b, 0.0), (-a, b, 0.0),
                        (-a, 0.0, -depth), (a, 0.0, -depth)])


def oval_part(sides=48):
    """Like faceted_part, but with an oval widest band."""
    return stacked_loops([(0.16, sides // 2, 0.3, 0.2), (0.02, sides, 0.5, 0.35),
                          (0.0, sides, 0.5, 0.35), (-0.43, 1, 0.0)])


def box(size=(1.0, 1.0, 1.0)):
    sx, sy, sz = size
    pts = [(x, y, z) for x in (0, sx) for y in (0, sy) for z in (0, sz)]
    return convex_hull(pts)


def _ear_clip(poly):
    """Triangulate a simple counter-clockwise 2D polygon (indices into poly)."""
    idx = list(range(len(poly)))
    tris = []

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    guard = 0
    while len(idx) > 3 and guard < 10000:
        guard += 1
        n = len(idx)
        for k in range(n):
            i0, i1, i2 = idx[(k - 1) % n], idx[k], idx[(k + 1) % n]
            a, b, c = poly[i0], poly[i1], poly[i2]
            if cross(a, b, c) <= 0:
                continue  # reflex corner
            inside = False
            for j in idx:
                if j in (i0, i1, i2):
                    continue
                p = poly[j]
                if cross(a, b, p) >= 0 and cross(b, c, p) >= 0 and cross(c, a, p) >= 0:
                    inside = True
                    break
            if not inside:
                tris.append((i0, i1, i2))
                idx.pop(k)
                break
    tris.append(tuple(idx))
    return tris


def prism(outline, height=1.0):
    """Closed prism over a simple counter-clockwise 2D outline (may be non-convex)."""
    n = len(outline)
    verts = [(x, y, 0.0) for x, y in outline] + [(x, y, height) for x, y in outline]
    faces = []
    for a, b, c in _ear_clip(outline):
        faces.append((a, c, b))              # bottom, facing down
        faces.append((a + n, b + n, c + n))  # top, facing up
    for i in range(n):
        j = (i + 1) % n
        faces.append((i, j, j + n))
        faces.append((i, j + n, i + n))
    return TriMesh(verts, faces)


def l_block(size=2.0, notch=1.0, height=1.0):
    """Non-convex L-shaped prism."""
    s, k = size, notch
    return prism([(0, 0), (s, 0), (s, s - k), (s - k, s - k), (s - k, s), (0, s)], height)


def prism_with_holes(outline, holes, height=1.0):
    """Closed prism over a 2D outline with holes running right through it."""
    from meshwright.core import poly2d

    def ccw(loop, want_ccw):
        a = poly2d.signed_area(loop)
        return list(loop) if (a > 0) == want_ccw else list(reversed(loop))

    loops = [ccw(outline, True)] + [ccw(h, False) for h in holes]
    flat = [p for loop in loops for p in loop]
    n = len(flat)
    verts = [(x, y, 0.0) for x, y in flat] + [(x, y, height) for x, y in flat]
    faces = []
    for a, b, c in poly2d.triangulate(loops[0], loops[1:]):
        faces.append((a, c, b))              # bottom, facing down
        faces.append((a + n, b + n, c + n))  # top, facing up
    start = 0
    for loop in loops:  # walls: the outline's face out, the holes' face into the hole
        m = len(loop)
        for k in range(m):
            i, j = start + k, start + (k + 1) % m
            faces.append((i, j, j + n))
            faces.append((i, j + n, i + n))
        start += m
    return TriMesh(verts, faces)


def frame(size=2.0, hole=1.0, height=0.5):
    """Square plate with a square hole through the middle."""
    a, b = 0.5 * size, 0.5 * hole
    return prism_with_holes([(-a, -a), (a, -a), (a, a), (-a, a)], [[(-b, -b), (b, -b), (b, b), (-b, b)]], height)


def plate_with_round_hole(size=4.0, radius=1.0, sides=24, height=0.5):
    a = 0.5 * size
    hole = [(radius * math.cos(2 * math.pi * k / sides), radius * math.sin(2 * math.pi * k / sides))
            for k in range(sides)]
    return prism_with_holes([(-a, -a), (a, -a), (a, a), (-a, a)], [hole], height)


def u_block(width=3.0, depth=2.0, slot=1.0, height=1.0):
    """U-shaped prism: the end faces of its two arms are separate faces on one plane."""
    w, d, k = width, depth, slot
    a = 0.5 * (w - k)
    return prism([(0, 0), (w, 0), (w, d), (w - a, d), (w - a, d - k), (a, d - k), (a, d), (0, d)], height)


def subdivided(mesh, times=1):
    """Split every triangle into four (edge midpoints), keeping the shape exactly.

    Turns a coarse mesh into a dense one like a scan, so noise then roughens
    flat regions the way it does on real data.
    """
    verts = [tuple(v) for v in mesh.vertices]
    faces = [tuple(int(i) for i in f) for f in mesh.faces]
    for _ in range(times):
        midpoint = {}

        def mid(a, b):
            key = (min(a, b), max(a, b))
            if key not in midpoint:
                pa, pb = verts[a], verts[b]
                verts.append(tuple((pa[k] + pb[k]) / 2.0 for k in range(3)))
                midpoint[key] = len(verts) - 1
            return midpoint[key]

        new_faces = []
        for a, b, c in faces:
            ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
            new_faces += [(a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca)]
        faces = new_faces
    return TriMesh(verts, faces)


# ---------------------------------------------------------------- transforms

def rng(seed=0):
    return np.random.default_rng(seed)


def random_rotation(generator):
    return Rotation.random(random_state=int(generator.integers(0, 2 ** 31 - 1))).as_matrix()


def moved(mesh, rotation=None, translation=(0.0, 0.0, 0.0)):
    out = mesh.copy()
    if rotation is not None:
        out.vertices = out.vertices.dot(np.asarray(rotation).T)
    out.vertices = out.vertices + np.asarray(translation, dtype=float)
    return out


def randomly_placed(mesh, generator, spread=50.0):
    """Mesh under a random rotation and translation. Returns (mesh, rotation)."""
    rot = random_rotation(generator)
    t = generator.uniform(-spread, spread, 3)
    return moved(mesh, rot, t), rot


# ---------------------------------------------------------------- damage

def as_float32(mesh):
    """Round vertex coordinates through float32, like many exchange formats do."""
    out = mesh.copy()
    out.vertices = out.vertices.astype(np.float32).astype(float)
    return out


def noisy(mesh, generator, sigma):
    out = mesh.copy()
    out.vertices = out.vertices + generator.normal(0.0, sigma, out.vertices.shape)
    return out


def with_holes(mesh, generator, count=1):
    """Remove count random faces, leaving open edges."""
    out = mesh.copy()
    drop = generator.choice(len(out.faces), size=count, replace=False)
    out.faces = np.delete(out.faces, drop, axis=0)
    return out


def with_flipped_faces(mesh, generator, count=1):
    out = mesh.copy()
    pick = generator.choice(len(out.faces), size=count, replace=False)
    out.faces[pick] = out.faces[pick][:, ::-1]
    return out


def with_duplicate_faces(mesh, generator, count=1):
    out = mesh.copy()
    pick = generator.choice(len(out.faces), size=count, replace=False)
    out.faces = np.vstack([out.faces, out.faces[pick]])
    return out


def combined(*meshes):
    """Several meshes as one (separate pieces)."""
    verts, faces, offset = [], [], 0
    for m in meshes:
        verts.append(m.vertices)
        faces.append(m.faces + offset)
        offset += len(m.vertices)
    return TriMesh(np.vstack(verts), np.vstack(faces))


def with_junk(mesh, generator, count=1, size_ratio=0.01):
    """Add tiny loose tetrahedra near the mesh (scan crumbs)."""
    extent = float(np.ptp(mesh.vertices, axis=0).max())
    pieces = [mesh]
    for _ in range(count):
        centre = mesh.vertices.mean(axis=0) + generator.uniform(-1.0, 1.0, 3) * extent
        tet = convex_hull(centre + generator.normal(0.0, size_ratio * extent, (4, 3)))
        pieces.append(tet)
    return combined(*pieces)


def unwelded(mesh):
    """Every triangle gets its own three vertices (like an STL import)."""
    v = mesh.vertices[mesh.faces].reshape(-1, 3)
    return TriMesh(v, np.arange(len(v)).reshape(-1, 3))
