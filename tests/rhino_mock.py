"""A small stand-in for the parts of RhinoCommon that solidwright's adapters use.

Rhino is not available in CI, so tests install this module as `Rhino`,
`Rhino.Geometry` and `Rhino.Geometry.Intersect`. It mimics names and behaviour
closely enough for the thin adapters in solidwright.rhino, and nothing more.

Covered: Point3d, Point3f, Vector3d, Plane, Transform, Polyline, Line, PolylineCurve,
Mesh (Vertices, Faces, TopologyVertices), Intersect.Intersection.MeshPlane and
a polygon-soup Brep: just enough of Brep (planar faces, join, from mesh, merge
coplanar faces, solid tests) to exercise solidwright.rhino.build. Set
Brep.fail_planar = True to make CreatePlanarBreps fail, as it can in Rhino.
Booleans, BoundingBox, VolumeMassProperties and Intersection.BrepBrep work
for Breps made of axis-aligned boxes (Brep.CreateFromBox), on a grid of cells.
"""

from __future__ import division

import sys
import types

import numpy as np


class Point3d(object):
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.X, self.Y, self.Z = float(x), float(y), float(z)

    def _a(self):
        return np.array([self.X, self.Y, self.Z])

    def __add__(self, other):
        return Point3d(*(self._a() + other._a()))

    def __sub__(self, other):
        r = self._a() - other._a()
        return Vector3d(*r) if isinstance(other, Point3d) else Point3d(*r)

    def __eq__(self, other):
        return isinstance(other, Point3d) and (self.X, self.Y, self.Z) == (other.X, other.Y, other.Z)

    def __hash__(self):
        return hash((self.X, self.Y, self.Z))

    def DistanceTo(self, other):
        return float(np.linalg.norm(self._a() - other._a()))

    def Transform(self, xform):
        self.X, self.Y, self.Z = xform._apply_point(self._a())

    def __repr__(self):
        return "Point3d(%g, %g, %g)" % (self.X, self.Y, self.Z)


class Point3f(Point3d):
    """Single precision point, as stored in Mesh.Vertices."""

    def __init__(self, x=0.0, y=0.0, z=0.0):
        x, y, z = np.float32(x), np.float32(y), np.float32(z)
        Point3d.__init__(self, float(x), float(y), float(z))


class Vector3d(object):
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.X, self.Y, self.Z = float(x), float(y), float(z)

    def _a(self):
        return np.array([self.X, self.Y, self.Z])

    @property
    def Length(self):
        return float(np.linalg.norm(self._a()))

    def Unitize(self):
        n = self.Length
        if n == 0.0:
            return False
        self.X, self.Y, self.Z = self._a() / n
        return True

    def __add__(self, other):
        return Vector3d(*(self._a() + other._a()))

    def __sub__(self, other):
        return Vector3d(*(self._a() - other._a()))

    def __mul__(self, other):
        if isinstance(other, Vector3d):
            return float(np.dot(self._a(), other._a()))
        return Vector3d(*(self._a() * other))

    __rmul__ = __mul__

    @staticmethod
    def CrossProduct(a, b):
        return Vector3d(*np.cross(a._a(), b._a()))

    def __repr__(self):
        return "Vector3d(%g, %g, %g)" % (self.X, self.Y, self.Z)


# RhinoCommon's static members. Under pythonnet they also show on instances
# (a Vector3d "has" a ZAxis), so the mock has them too.
Vector3d.Zero = Vector3d(0.0, 0.0, 0.0)
Vector3d.XAxis = Vector3d(1.0, 0.0, 0.0)
Vector3d.YAxis = Vector3d(0.0, 1.0, 0.0)
Vector3d.ZAxis = Vector3d(0.0, 0.0, 1.0)
Point3d.Origin = Point3d(0.0, 0.0, 0.0)


def _unit(a):
    n = np.linalg.norm(a)
    return a / n if n > 0 else a


class Plane(object):
    def __init__(self, origin, normal_or_x, y_axis=None):
        o = np.array([origin.X, origin.Y, origin.Z])
        if y_axis is None:
            z = _unit(np.array([normal_or_x.X, normal_or_x.Y, normal_or_x.Z]))
            seed = np.eye(3)[int(np.argmin(np.abs(z)))]
            x = _unit(np.cross(seed, z))
        else:
            x = _unit(np.array([normal_or_x.X, normal_or_x.Y, normal_or_x.Z]))
            y = np.array([y_axis.X, y_axis.Y, y_axis.Z])
            y = _unit(y - np.dot(y, x) * x)
            z = np.cross(x, y)
        y = np.cross(z, x)
        self._o, self._x, self._y, self._z = o, x, y, z

    class _WorldXY(object):
        """A static property in RhinoCommon: a fresh plane on every read."""

        def __get__(self, obj, cls):
            return Plane(Point3d(0, 0, 0), Vector3d(1, 0, 0), Vector3d(0, 1, 0))

    WorldXY = _WorldXY()

    @property
    def Origin(self):
        return Point3d(*self._o)

    @property
    def XAxis(self):
        return Vector3d(*self._x)

    @property
    def YAxis(self):
        return Vector3d(*self._y)

    @property
    def ZAxis(self):
        return Vector3d(*self._z)

    @property
    def Normal(self):
        return self.ZAxis

    def DistanceTo(self, point):
        return float(np.dot(point._a() - self._o, self._z))

    def PointAt(self, u, v, w=0.0):
        return Point3d(*(self._o + u * self._x + v * self._y + w * self._z))


def _box_of(points, plane_or_flag):
    pts = np.asarray(points, dtype=float).reshape(-1, 3)
    if isinstance(plane_or_flag, Plane):
        pl = plane_or_flag
        d = pts - pl._o
        pts = np.stack([d.dot(pl._x), d.dot(pl._y), d.dot(pl._z)], axis=1)
    return BoundingBox(Point3d(*pts.min(axis=0)), Point3d(*pts.max(axis=0)))


class Transform(object):
    def __init__(self, matrix=None):
        self.M = np.eye(4) if matrix is None else np.array(matrix, dtype=float)

    @staticmethod
    def Identity():
        return Transform()

    @staticmethod
    def Translation(v):
        m = np.eye(4)
        m[:3, 3] = [v.X, v.Y, v.Z]
        return Transform(m)

    @staticmethod
    def Scale(plane, sx, sy, sz):
        """Scale by sx, sy, sz along the plane's axes, about its origin."""
        f = np.eye(4)
        f[:3, 0], f[:3, 1], f[:3, 2], f[:3, 3] = plane._x, plane._y, plane._z, plane._o
        return Transform(f.dot(np.diag([sx, sy, sz, 1.0])).dot(np.linalg.inv(f)))

    @staticmethod
    def PlaneToPlane(a, b):
        """Maps plane a's frame onto plane b's frame."""
        fa = np.eye(4)
        fa[:3, 0], fa[:3, 1], fa[:3, 2], fa[:3, 3] = a._x, a._y, a._z, a._o
        fb = np.eye(4)
        fb[:3, 0], fb[:3, 1], fb[:3, 2], fb[:3, 3] = b._x, b._y, b._z, b._o
        return Transform(fb.dot(np.linalg.inv(fa)))

    def __mul__(self, other):
        return Transform(self.M.dot(other.M))

    def TryGetInverse(self):
        try:
            return True, Transform(np.linalg.inv(self.M))
        except np.linalg.LinAlgError:
            return False, Transform()

    def _apply_point(self, p):
        return self.M[:3, :3].dot(p) + self.M[:3, 3]


class Polyline(list):
    """List of Point3d. Closed when the last point equals the first."""

    def Add(self, p):
        self.append(p)

    @property
    def Count(self):
        return len(self)

    @property
    def IsClosed(self):
        return len(self) > 2 and self[0] == self[-1]


class _MeshFace(object):
    def __init__(self, a, b, c, d=None):
        self.A, self.B, self.C = a, b, c
        self.D = c if d is None else d

    @property
    def IsQuad(self):
        return self.C != self.D

    @property
    def IsTriangle(self):
        return not self.IsQuad


class _Vertices(list):
    """Like Rhino 8: positions kept in double precision (Point3dAt) unless
    UseDoublePrecisionVertices is switched off; items read as Point3f."""

    def __init__(self):
        list.__init__(self)
        self.UseDoublePrecisionVertices = True
        self._exact = []

    @property
    def Count(self):
        return len(self)

    def Add(self, x, y=None, z=None):
        if y is None:
            x, y, z = x.X, x.Y, x.Z
        self.append(Point3f(x, y, z))
        self._exact.append(Point3d(x, y, z) if self.UseDoublePrecisionVertices else Point3f(x, y, z))
        return len(self) - 1

    def Point3dAt(self, i):
        p = self._exact[i]
        return Point3d(p.X, p.Y, p.Z)


class _Faces(list):
    @property
    def Count(self):
        return len(self)

    def AddFace(self, a, b, c, d=None):
        self.append(_MeshFace(a, b, c, d))
        return len(self) - 1


class _TopologyVertices(object):
    """Welds vertices with identical (float32) coordinates, like Rhino does."""

    def __init__(self, mesh):
        self._mesh = mesh
        self._build()

    def _build(self):
        key_to_top = {}
        self._top_of = []
        self._members = []
        self._points = []
        for i, p in enumerate(self._mesh.Vertices):
            key = (p.X, p.Y, p.Z)
            t = key_to_top.get(key)
            if t is None:
                t = len(self._members)
                key_to_top[key] = t
                self._members.append([])
                self._points.append(Point3f(p.X, p.Y, p.Z))
            self._top_of.append(t)
            self._members[t].append(i)

    @property
    def Count(self):
        return len(self._members)

    def __getitem__(self, i):
        return self._points[i]

    def TopologyVertexIndex(self, vertex_index):
        return self._top_of[vertex_index]

    def MeshVertexIndices(self, topology_index):
        return list(self._members[topology_index])


class _Normals(list):
    def ComputeNormals(self):
        return True


class Mesh(object):
    def __init__(self):
        self.Vertices = _Vertices()
        self.Faces = _Faces()
        self.Normals = _Normals()

    def Compact(self):
        return True

    @property
    def TopologyVertices(self):
        return _TopologyVertices(self)

    @property
    def IsValid(self):
        n = self.Vertices.Count
        return all(0 <= i < n for f in self.Faces for i in (f.A, f.B, f.C, f.D))

    def _triangles(self):
        for f in self.Faces:
            ids = [f.A, f.B, f.C] + ([f.D] if f.IsQuad else [])
            pts = [self.Vertices.Point3dAt(i)._a() for i in ids]
            for k in range(1, len(pts) - 1):
                yield pts[0], pts[k], pts[k + 1]

    def GetBoundingBox(self, accurate_or_plane):
        return _box_of([self.Vertices.Point3dAt(i)._a() for i in range(self.Vertices.Count)], accurate_or_plane)

    def Duplicate(self):
        m = Mesh()
        m.Vertices.UseDoublePrecisionVertices = self.Vertices.UseDoublePrecisionVertices
        for i in range(self.Vertices.Count):
            m.Vertices.Add(self.Vertices.Point3dAt(i))
        for f in self.Faces:
            m.Faces.AddFace(f.A, f.B, f.C, f.D if f.IsQuad else None)
        return m

    def Transform(self, xform):
        exact = [self.Vertices.Point3dAt(i) for i in range(self.Vertices.Count)]
        del self.Vertices[:]
        del self.Vertices._exact[:]
        for p in exact:
            self.Vertices.Add(*xform._apply_point(p._a()))
        return True

    @property
    def IsClosed(self):
        top = self.TopologyVertices
        counts = {}
        for f in self.Faces:
            ids = [f.A, f.B, f.C] + ([f.D] if f.IsQuad else [])
            ids = [top.TopologyVertexIndex(i) for i in ids]
            for k in range(len(ids)):
                e = tuple(sorted((ids[k], ids[(k + 1) % len(ids)])))
                counts[e] = counts.get(e, 0) + 1
        return bool(counts) and all(c == 2 for c in counts.values())

    @staticmethod
    def CreateFromBrep(brep, parameters):
        m = Mesh()
        for poly in brep._polygons():
            ids = [m.Vertices.Add(*p) for p in poly]
            for k in range(1, len(ids) - 1):
                m.Faces.AddFace(ids[0], ids[k], ids[k + 1])
        return [m]

    @staticmethod
    def from_trimesh(trimesh):
        """Test helper (not RhinoCommon): build a mock Mesh from tests.shapes.TriMesh."""
        m = Mesh()
        for v in trimesh.vertices:
            m.Vertices.Add(*v)
        for f in trimesh.faces:
            m.Faces.AddFace(*[int(i) for i in f])
        return m


class Line(object):
    def __init__(self, a, b):
        self.From, self.To = a, b

    @property
    def Length(self):
        return self.From.DistanceTo(self.To)


class PolylineCurve(object):
    """Parameter t runs 0, 1, 2... over the points, like Rhino's polyline curves."""

    def __init__(self, points):
        self.points = [Point3d(p.X, p.Y, p.Z) for p in points]

    def PointAt(self, t):
        k = min(int(np.floor(t)), len(self.points) - 2)
        u = t - k
        a, b = self.points[k]._a(), self.points[k + 1]._a()
        return Point3d(*(a + u * (b - a)))

    def DivideByCount(self, count, include_ends):
        """Parameters at equal lengths (closed curves: count points, the end not repeated)."""
        pts = np.array([p._a() for p in self.points])
        seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
        cum = np.concatenate([[0.0], np.cumsum(seg)])
        n = count if self.IsClosed else count + 1
        out = []
        for i in range(n):
            s = cum[-1] * i / count
            k = min(int(np.searchsorted(cum, s, side="right")) - 1, len(seg) - 1)
            out.append(k + ((s - cum[k]) / seg[k] if seg[k] > 0 else 0.0))
        return out

    @property
    def IsClosed(self):
        return len(self.points) > 3 and self.points[0] == self.points[-1]


class Circle(object):
    def __init__(self, plane, radius):
        self.Plane, self.Radius = plane, float(radius)

    @property
    def Center(self):
        return self.Plane.Origin


class MeshingParameters(object):
    FastRenderMesh = "fast"


class SpaceMorph(object):
    """Calls MorphPoint on every Brep corner or mesh vertex (the mock Brep is
    made of flat polygons, so moving its corners is the whole morph)."""

    def __init__(self):
        self.PreserveStructure = False
        self.QuickPreview = False
        self.Tolerance = 0.0

    def Morph(self, geometry):
        if isinstance(geometry, Mesh):
            pts = [self.MorphPoint(geometry.Vertices.Point3dAt(i)) for i in range(geometry.Vertices.Count)]
            del geometry.Vertices[:]
            del geometry.Vertices._exact[:]
            for p in pts:
                geometry.Vertices.Add(p.X, p.Y, p.Z)
            return True
        for poly in geometry._polygons():
            for k in range(len(poly)):
                q = self.MorphPoint(Point3d(*poly[k]))
                poly[k][:] = (q.X, q.Y, q.Z)
        return True


class BrepSolidOrientation(object):
    None_ = 0
    Outward = 1
    Inward = -1
    Unknown = 2


class _BrepFaces(list):
    @property
    def Count(self):
        return len(self)


class Brep(object):
    """Polygon faces (lists of numpy points), each face a list of polygons
    (one for a planar face, several triangles while coplanar faces are not merged)."""

    fail_planar = False

    def __init__(self, faces):
        self.Faces = _BrepFaces(faces)
        self.tolerance = 1e-6

    def _polygons(self):
        return [poly for face in self.Faces for poly in face]

    @staticmethod
    def CreatePlanarBreps(curves, tolerance):
        """One planar face from a closed curve, or from a list of them (outer loop first, then holes)."""
        if not isinstance(curves, (list, tuple)):
            curves = [curves]
        if Brep.fail_planar or not all(c.IsClosed for c in curves):
            return None
        loops = [[p._a() for p in c.points[:-1]] for c in curves]
        normal = np.zeros(3)
        for i in range(len(loops[0])):
            normal += np.cross(loops[0][i], loops[0][(i + 1) % len(loops[0])])
        if np.linalg.norm(normal) <= tolerance * tolerance:
            return None
        return [Brep([loops])]

    @staticmethod
    def JoinBreps(breps, tolerance):
        joined = Brep([face for b in breps for face in b.Faces])
        joined.tolerance = tolerance
        return [joined]

    @staticmethod
    def CreateFromMesh(mesh, trimmed_triangles):
        faces = []
        for f in mesh.Faces:
            ids = [f.A, f.B, f.C] + ([f.D] if f.IsQuad else [])
            faces.append([[np.array([mesh.Vertices[i].X, mesh.Vertices[i].Y, mesh.Vertices[i].Z])
                           for i in ids]])
        return Brep(faces)

    def _edge_counts(self):
        def key(p):
            return tuple(np.round(p / max(self.tolerance, 1e-12)).astype(int))
        counts = {}
        for poly in self._polygons():
            for i in range(len(poly)):
                e = (key(poly[i]), key(poly[(i + 1) % len(poly)]))
                counts[e] = counts.get(e, 0) + 1
        return counts

    @property
    def IsSolid(self):
        counts = self._edge_counts()
        return bool(counts) and all(c == 1 and counts.get((b, a), 0) == 1 for (a, b), c in counts.items())

    @property
    def IsValid(self):
        return len(self.Faces) > 0

    def _volume(self):
        total = 0.0
        for poly in self._polygons():
            for i in range(1, len(poly) - 1):
                total += float(np.dot(poly[0], np.cross(poly[i], poly[i + 1])))
        return total / 6.0

    @property
    def SolidOrientation(self):
        if not self.IsSolid:
            return BrepSolidOrientation.None_
        return BrepSolidOrientation.Outward if self._volume() > 0 else BrepSolidOrientation.Inward

    def Flip(self):
        for face in self.Faces:
            for poly in face:
                poly.reverse()

    # ---- bounding box, copies and moves

    def _points(self):
        return np.array([p for poly in self._polygons() for p in poly])

    def GetBoundingBox(self, accurate_or_plane):
        return _box_of(self._points(), accurate_or_plane)

    def Duplicate(self):
        return self.DuplicateBrep()

    def MakeDeformable(self):
        return True

    def Transform(self, xform):
        for poly in self._polygons():
            for k in range(len(poly)):
                poly[k][:] = xform._apply_point(poly[k])
        return True

    def DuplicateBrep(self):
        copy = Brep([[[p.copy() for p in poly] for poly in face] for face in self.Faces])
        copy.tolerance = self.tolerance
        return copy

    def Translate(self, v):
        for poly in self._polygons():
            for p in poly:
                p += v._a()
        return True

    @staticmethod
    def CreateFromBox(box):
        """A box as six outward quads."""
        (x0, y0, z0), (x1, y1, z1) = box.Min._a(), box.Max._a()
        return _cells_to_breps(*_grid([x0, x1], [y0, y1], [z0, z1], lambda c: True))[0]

    # ---- Booleans, for Breps whose faces are all axis-aligned rectangles

    @staticmethod
    def CreateBooleanUnion(breps, tolerance):
        breps = list(breps)
        if not all(b._boxy() for b in breps):
            return None
        xs, ys, zs = _grid_of(breps)
        return _cells_to_breps(xs, ys, zs, [c for c in _all_cells(xs, ys, zs)
                                             if any(b._inside(_centre(xs, ys, zs, c)) for b in breps)])

    @staticmethod
    def CreateBooleanIntersection(a, b, tolerance):
        if not (a._boxy() and b._boxy()):
            return None
        xs, ys, zs = _grid_of([a, b])
        cells = [c for c in _all_cells(xs, ys, zs)
                 if a._inside(_centre(xs, ys, zs, c)) and b._inside(_centre(xs, ys, zs, c))]
        return _cells_to_breps(xs, ys, zs, cells) if cells else []

    @property
    def Vertices(self):
        return _BrepVertices(_BrepVertex(p) for p in self._points())

    def IsPointInside(self, point, tolerance, strictly_in):
        return self._inside(point._a())

    def _boxy(self):
        """Every face polygon is a rectangle square to an axis."""
        for poly in self._polygons():
            spread = np.ptp(np.array(poly), axis=0)
            if len(poly) != 4 or np.sum(spread < 1e-12) != 1:
                return False
        return True

    def _inside(self, p):
        """Ray along +X: count the X-facing rectangles it passes through."""
        hits = 0
        for poly in self._polygons():
            a = np.array(poly)
            lo, hi = a.min(axis=0), a.max(axis=0)
            if hi[0] - lo[0] < 1e-12 and lo[0] > p[0] and lo[1] < p[1] < hi[1] and lo[2] < p[2] < hi[2]:
                hits += 1
        return hits % 2 == 1

    def MergeCoplanarFaces(self, tolerance):
        """Join faces that share an edge and lie on the same plane (normal and
        offset rounded to the tolerance)."""
        def key_of(p):
            return tuple(np.round(p / max(self.tolerance, 1e-12)).astype(int))
        plane_keys, edge_faces = [], {}
        for fi, face in enumerate(self.Faces):
            poly = face[0]
            n = np.zeros(3)
            for i in range(len(poly)):
                n += np.cross(poly[i], poly[(i + 1) % len(poly)])
            n = n / np.linalg.norm(n)
            plane_keys.append((tuple(np.round(n / 1e-6).astype(int)), int(round(np.dot(n, poly[0]) / tolerance))))
            for loop in face:
                for i in range(len(loop)):
                    e = tuple(sorted((key_of(loop[i]), key_of(loop[(i + 1) % len(loop)]))))
                    edge_faces.setdefault(e, []).append(fi)
        parent = list(range(len(self.Faces)))

        def find(i):
            while parent[i] != i:
                i = parent[i]
            return i
        for fs in edge_faces.values():
            for a in fs:
                for b in fs:
                    if a != b and plane_keys[a] == plane_keys[b]:
                        parent[find(a)] = find(b)
        groups, order = {}, []
        for fi, face in enumerate(self.Faces):
            r = find(fi)
            if r not in groups:
                groups[r] = []
                order.append(r)
            groups[r].extend(face)
        merged = len(order) < len(self.Faces)
        self.Faces = _BrepFaces(groups[r] for r in order)
        return merged


class _BrepVertex(object):
    def __init__(self, p):
        self.Location = Point3d(*p)


class _BrepVertices(list):
    @property
    def Count(self):
        return len(self)


class BoundingBox(object):
    def __init__(self, a, b):
        self.Min, self.Max = a, b


class VolumeMassProperties(object):
    def __init__(self, volume):
        self.Volume = volume

    @staticmethod
    def Compute(geometry):
        if isinstance(geometry, Mesh):
            if not geometry.IsClosed:
                return None
            return VolumeMassProperties(sum(float(np.dot(a, np.cross(b, c))) for a, b, c in
                                            geometry._triangles()) / 6.0)
        return VolumeMassProperties(geometry._volume()) if geometry.IsSolid else None


class AreaMassProperties(object):
    def __init__(self, area):
        self.Area = area

    @staticmethod
    def Compute(geometry):
        if isinstance(geometry, Mesh):
            tris = geometry._triangles()
        else:
            tris = ((poly[0], poly[i], poly[i + 1]) for poly in geometry._polygons()
                    for i in range(1, len(poly) - 1))
        return AreaMassProperties(sum(0.5 * float(np.linalg.norm(np.cross(b - a, c - a))) for a, b, c in tris))


def _grid_of(breps):
    pts = np.concatenate([b._points() for b in breps])
    return [sorted(set(np.round(pts[:, k], 9))) for k in range(3)]


def _all_cells(xs, ys, zs):
    return [(i, j, k) for i in range(len(xs) - 1) for j in range(len(ys) - 1) for k in range(len(zs) - 1)]


def _centre(xs, ys, zs, c):
    i, j, k = c
    return ((xs[i] + xs[i + 1]) / 2, (ys[j] + ys[j + 1]) / 2, (zs[k] + zs[k + 1]) / 2)


def _grid(xs, ys, zs, keep):
    return xs, ys, zs, [c for c in _all_cells(xs, ys, zs) if keep(c)]


_STEPS = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]


def _boundary_corners(cells):
    """Grid corners (index triples) on the faces between the cells and the rest."""
    out = set()
    for c in cells:
        for axis in range(3):
            for side in (1, -1):
                n = list(c)
                n[axis] += side
                if tuple(n) in cells:
                    continue
                u, v = (axis + 1) % 3, (axis + 2) % 3
                for du in (0, 1):
                    for dv in (0, 1):
                        k = list(c)
                        k[axis] += 1 if side > 0 else 0
                        k[u] += du
                        k[v] += dv
                        out.add(tuple(k))
    return out


def _cells_to_breps(xs, ys, zs, cells):
    """Face-connected groups of grid cells, each as a Brep of outward cell-face quads."""
    cells = set(cells)
    out = []
    while cells:
        todo = [cells.pop()]
        group = set(todo)
        while todo:
            i, j, k = todo.pop()
            for d in _STEPS:
                n = (i + d[0], j + d[1], k + d[2])
                if n in cells:
                    cells.discard(n)
                    group.add(n)
                    todo.append(n)
        faces = []
        g = (xs, ys, zs)
        for c in group:
            for axis in range(3):
                for side in (1, -1):
                    n = list(c)
                    n[axis] += side
                    if tuple(n) in group:
                        continue
                    u, v = (axis + 1) % 3, (axis + 2) % 3
                    corner = [None] * 3
                    quad = []
                    for du, dv in ((0, 0), (1, 0), (1, 1), (0, 1)):
                        corner[axis] = g[axis][c[axis] + (1 if side > 0 else 0)]
                        corner[u] = g[u][c[u] + du]
                        corner[v] = g[v][c[v] + dv]
                        quad.append(np.array(corner, dtype=float))
                    if side < 0:
                        quad.reverse()
                    faces.append([quad])
        out.append(Brep(faces))
    return out


class _Intersection(object):
    @staticmethod
    def BrepBrep(a, b, tolerance):
        """(True, curves, points). Curves exist where the two boundaries meet
        (a piece fully inside the other gives none, as in Rhino)."""
        xs, ys, zs = _grid_of([a, b])
        ina = set(c for c in _all_cells(xs, ys, zs) if a._inside(_centre(xs, ys, zs, c)))
        inb = set(c for c in _all_cells(xs, ys, zs) if b._inside(_centre(xs, ys, zs, c)))
        meet = _boundary_corners(ina) & _boundary_corners(inb)
        return True, ([Line(Point3d(), Point3d(1, 0, 0))] if meet else []), []

    @staticmethod
    def MeshPlane(mesh, plane):
        """Closed or open section polylines of a mesh by a plane, or None.

        Each triangle crossing the plane contributes one segment whose end
        points are keyed by the (welded) edge they lie on, so segments chain
        exactly. Vertices lying on the plane count as above it, a standard
        trick that avoids degenerate cases.
        """
        top = mesh.TopologyVertices
        pts = np.array([[p.X, p.Y, p.Z] for p in top._points])
        d = (pts - plane._o).dot(plane._z)
        side = d >= 0.0
        links = {}
        cross_point = {}

        def edge_point(i, j):
            key = (min(i, j), max(i, j))
            if key not in cross_point:
                t = d[i] / (d[i] - d[j])
                cross_point[key] = Point3d(*(pts[i] + t * (pts[j] - pts[i])))
            return key

        for f in mesh.Faces:
            corners = [top.TopologyVertexIndex(f.A), top.TopologyVertexIndex(f.B),
                       top.TopologyVertexIndex(f.C)]
            if f.IsQuad:
                corners.append(top.TopologyVertexIndex(f.D))
            ends = []
            for k in range(len(corners)):
                i, j = corners[k], corners[(k + 1) % len(corners)]
                if side[i] != side[j]:
                    ends.append(edge_point(i, j))
            for k in range(0, len(ends) - 1, 2):
                a, b = ends[k], ends[k + 1]
                links.setdefault(a, []).append(b)
                links.setdefault(b, []).append(a)
        if not links:
            return None
        seen = set()
        result = []
        # Start open chains at their ends first, then pick up closed loops.
        starts = [k for k, v in links.items() if len(v) == 1] + list(links.keys())
        for start in starts:
            if start in seen:
                continue
            chain = [start]
            seen.add(start)
            current = start
            while True:
                nxt = [n for n in links[current] if n not in seen]
                if not nxt:
                    break
                current = nxt[0]
                seen.add(current)
                chain.append(current)
            line = Polyline(cross_point[k] for k in chain)
            if len(chain) > 2 and chain[0] in links[chain[-1]]:
                line.append(cross_point[chain[0]])
            result.append(line)
        return result


class _GenericList(object):
    """System.Collections.Generic.List: List[T]() gives a list with .Add."""

    class _List(list):
        def Add(self, item):
            self.append(item)

    def __getitem__(self, item_type):
        return self._List


def install():
    """Register the mock as the Rhino, Rhino.Geometry and Rhino.Geometry.Intersect modules."""
    if "Rhino" in sys.modules and not getattr(sys.modules["Rhino"], "IS_MOCK", False):
        return  # real Rhino present: leave it alone
    rhino = types.ModuleType("Rhino")
    rhino.IS_MOCK = True
    geometry = types.ModuleType("Rhino.Geometry")
    intersect = types.ModuleType("Rhino.Geometry.Intersect")
    for cls in (Point3d, Point3f, Vector3d, Plane, Transform, Polyline, Line, PolylineCurve, Mesh,
                Brep, BrepSolidOrientation, BoundingBox, VolumeMassProperties, AreaMassProperties, Circle,
                MeshingParameters, SpaceMorph):
        setattr(geometry, cls.__name__, cls)
    intersect.Intersection = _Intersection
    geometry.Intersect = intersect
    rhino.Geometry = geometry
    if "System" not in sys.modules:
        system = types.ModuleType("System")
        system.Double = float
        collections = types.ModuleType("System.Collections")
        generic = types.ModuleType("System.Collections.Generic")
        generic.List = _GenericList()
        system.Collections = collections
        collections.Generic = generic
        sys.modules["System"] = system
        sys.modules["System.Collections"] = collections
        sys.modules["System.Collections.Generic"] = generic
    sys.modules["Rhino"] = rhino
    sys.modules["Rhino.Geometry"] = geometry
    sys.modules["Rhino.Geometry.Intersect"] = intersect
