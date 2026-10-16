"""Apply meshwright.core.deform point maps to Rhino geometry.

The point maps work in a local frame (z along the axis); a Plane gives that
frame. Each piece of geometry is copied and morphed through a SpaceMorph
that converts every point to the local frame, maps it and converts it back.
A morph that fails or leaves an invalid result is reported for that input,
never passed through silently.
"""

from __future__ import absolute_import, division

import Rhino.Geometry as rg

from ..core.report import Report


class _Frame(object):
    def __init__(self, plane):
        self.o = (plane.Origin.X, plane.Origin.Y, plane.Origin.Z)
        self.x = (plane.XAxis.X, plane.XAxis.Y, plane.XAxis.Z)
        self.y = (plane.YAxis.X, plane.YAxis.Y, plane.YAxis.Z)
        self.z = (plane.ZAxis.X, plane.ZAxis.Y, plane.ZAxis.Z)

    def local(self, p):
        d = (p[0] - self.o[0], p[1] - self.o[1], p[2] - self.o[2])
        return (_dot(d, self.x), _dot(d, self.y), _dot(d, self.z))

    def world(self, q):
        o, x, y, z = self.o, self.x, self.y, self.z
        return tuple(o[k] + x[k] * q[0] + y[k] * q[1] + z[k] * q[2] for k in range(3))


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


class _LocalMorph(rg.SpaceMorph):
    """Set .frame and .fn after creating it (a .NET subclass takes no constructor arguments)."""

    def MorphPoint(self, p):
        q = self.frame.world(self.fn(self.frame.local((p.X, p.Y, p.Z))))
        return rg.Point3d(q[0], q[1], q[2])


def local_points(items, plane):
    """Points of the geometry in the plane's frame: mesh vertices, or the
    vertices of a quick render mesh of each Brep."""
    frame = _Frame(plane)
    out = []
    for g in items:
        for m in _meshes(g):
            for i in range(m.Vertices.Count):
                p = m.Vertices[i]
                out.append(frame.local((p.X, p.Y, p.Z)))
    return out


def outline(curve, plane, count=720):
    """A closed curve as (x, y) points in the plane's frame (its height along the axis dropped)."""
    frame = _Frame(plane)
    ts = curve.DivideByCount(count, True) or []
    out = []
    for t in ts:
        p = curve.PointAt(t)
        q = frame.local((p.X, p.Y, p.Z))
        out.append((q[0], q[1]))
    return out


def _meshes(g):
    if isinstance(g, rg.Mesh):
        return [g]
    if isinstance(g, rg.Brep):
        return list(rg.Mesh.CreateFromBrep(g, rg.MeshingParameters.FastRenderMesh) or [])
    return []


def extent(items, plane):
    """(low, high) of all the geometry along the plane's Z, in the plane's frame, or None."""
    lows, highs = [], []
    for g in items:
        if g is None:
            continue
        box = g.GetBoundingBox(plane)
        lows.append(box.Min.Z)
        highs.append(box.Max.Z)
    if not lows:
        return None
    return min(lows), max(highs)


def morph(items, point_map, plane, tolerance, affine=None):
    """(copies, report): every input morphed by point_map in the plane's frame.

    affine, when given, is a Transform that does exactly what the map does on
    this geometry (a compress starting below it is a plain scale); it is tried
    first because it keeps the geometry exact, with the morph as fallback.
    """
    report = Report("deform")
    shape = _LocalMorph()
    shape.frame = _Frame(plane)
    shape.fn = point_map
    shape.PreserveStructure = False
    shape.QuickPreview = False
    shape.Tolerance = tolerance
    out = []
    for i, g in enumerate(items):
        if g is None:
            out.append(None)
            continue
        was_solid = _solid(g)
        moved = None
        if affine is not None:
            moved = _transformed(g, affine)
        if moved is None:
            moved = g.Duplicate()
            try:
                ok = shape.Morph(moved)
            except Exception:
                ok = False
            if not ok:
                report.messy("input %d could not be morphed: passed through unchanged" % i)
                out.append(g.Duplicate())
                continue
        if not moved.IsValid:
            report.messy("input %d is not valid after the morph" % i)
        elif was_solid and not _solid(moved):
            report.messy("input %d was closed and is open after the morph" % i)
        out.append(moved)
    return out, report


def _transformed(g, xform):
    d = g.Duplicate()
    if isinstance(d, rg.Brep) and not d.MakeDeformable():
        return None
    if not d.Transform(xform) or not d.IsValid:
        return None
    return d


def _solid(g):
    if isinstance(g, rg.Brep):
        return bool(g.IsSolid)
    if isinstance(g, rg.Mesh):
        return bool(g.IsClosed)
    return False
