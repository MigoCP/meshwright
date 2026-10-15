"""A small stand-in for the parts of RhinoCommon that meshwright's adapters use.

Rhino is not available in CI, so tests install this module as `Rhino`,
`Rhino.Geometry` and `Rhino.Geometry.Intersect`. It mimics names and behaviour
closely enough for the thin adapters in meshwright.rhino, and nothing more.

Covered: Point3d, Point3f, Vector3d, Plane, Transform, Polyline, PolylineCurve,
Mesh (Vertices, Faces, TopologyVertices), Intersect.Intersection.MeshPlane and
a polygon-soup Brep: just enough of Brep (planar faces, join, from mesh, merge
coplanar faces, solid tests) to exercise meshwright.rhino.build. Set
Brep.fail_planar = True to make CreatePlanarBreps fail, as it can in Rhino.
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

    @staticmethod
    def WorldXY():
        return Plane(Point3d(0, 0, 0), Vector3d(1, 0, 0), Vector3d(0, 1, 0))

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

    @staticmethod
    def from_trimesh(trimesh):
        """Test helper (not RhinoCommon): build a mock Mesh from tests.shapes.TriMesh."""
        m = Mesh()
        for v in trimesh.vertices:
            m.Vertices.Add(*v)
        for f in trimesh.faces:
            m.Faces.AddFace(*[int(i) for i in f])
        return m


class PolylineCurve(object):
    def __init__(self, points):
        self.points = [Point3d(p.X, p.Y, p.Z) for p in points]

    @property
    def IsClosed(self):
        return len(self.points) > 3 and self.points[0] == self.points[-1]


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
    def CreatePlanarBreps(curve, tolerance):
        if Brep.fail_planar or not curve.IsClosed:
            return None
        poly = [p._a() for p in curve.points[:-1]]
        normal = np.zeros(3)
        for i in range(len(poly)):
            normal += np.cross(poly[i], poly[(i + 1) % len(poly)])
        if np.linalg.norm(normal) <= tolerance * tolerance:
            return None
        return [Brep([[poly]])]

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

    def MergeCoplanarFaces(self, tolerance):
        """Group faces by plane (normal and offset rounded to the tolerance)."""
        groups = {}
        order = []
        for face in self.Faces:
            poly = face[0]
            n = np.zeros(3)
            for i in range(len(poly)):
                n += np.cross(poly[i], poly[(i + 1) % len(poly)])
            n = n / np.linalg.norm(n)
            key = (tuple(np.round(n / 1e-6).astype(int)), int(round(np.dot(n, poly[0]) / tolerance)))
            if key not in groups:
                groups[key] = []
                order.append(key)
            groups[key].extend(face)
        merged = len(order) < len(self.Faces)
        self.Faces = _BrepFaces(groups[k] for k in order)
        return merged


class _Intersection(object):
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


def install():
    """Register the mock as the Rhino, Rhino.Geometry and Rhino.Geometry.Intersect modules."""
    if "Rhino" in sys.modules and not getattr(sys.modules["Rhino"], "IS_MOCK", False):
        return  # real Rhino present: leave it alone
    rhino = types.ModuleType("Rhino")
    rhino.IS_MOCK = True
    geometry = types.ModuleType("Rhino.Geometry")
    intersect = types.ModuleType("Rhino.Geometry.Intersect")
    for cls in (Point3d, Point3f, Vector3d, Plane, Transform, Polyline, PolylineCurve, Mesh,
                Brep, BrepSolidOrientation):
        setattr(geometry, cls.__name__, cls)
    intersect.Intersection = _Intersection
    geometry.Intersect = intersect
    rhino.Geometry = geometry
    sys.modules["Rhino"] = rhino
    sys.modules["Rhino.Geometry"] = geometry
    sys.modules["Rhino.Geometry.Intersect"] = intersect
