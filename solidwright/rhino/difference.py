"""Rhino Brep calls for solidwright.core.difference_plan.

BrepOps answers difference_plan's questions with RhinoCommon: Boolean
difference, volume, bounding box, intersection curves, point inclusion,
points spread over the faces, a scaled copy and a sealed hollow. The solid tests, volume,
box and intersection curves are union's (solidwright.rhino.union).
Every call that Rhino can fail returns None (or False) instead of raising,
so a failure is just one more refusing cutter.
"""

from __future__ import absolute_import, division

import Rhino.Geometry as rg

from ..core import difference_plan
from . import union

GRIDS = (3, 7)  # samples per face: an n x n grid over its parameters, the finer one for small trimmed faces


class BrepOps(union.BrepOps, difference_plan.Ops):
    def difference(self, solids, cutters, tolerance):
        try:
            out = rg.Brep.CreateBooleanDifference(list(solids), list(cutters), tolerance, True)
        except Exception:
            return None
        if out is None:
            return None
        out = list(out)
        for b in out:
            if b is not None and b.IsSolid and b.SolidOrientation == rg.BrepSolidOrientation.Inward:
                b.Flip()
        return out

    def inside(self, brep, point):
        try:
            return bool(brep.IsPointInside(rg.Point3d(point[0], point[1], point[2]), self.tolerance, True))
        except Exception:
            return False

    def samples(self, brep):
        found = []
        for face in brep.Faces:
            for n in GRIDS:
                got = _face_samples(face, n)
                if got:
                    found.extend(got)
                    break
        return found

    def scaled(self, brep, centre, factors):
        plane = rg.Plane(rg.Point3d(centre[0], centre[1], centre[2]), rg.Vector3d.XAxis, rg.Vector3d.YAxis)
        copy = brep.DuplicateBrep()
        try:
            if not copy.Transform(rg.Transform.Scale(plane, factors[0], factors[1], factors[2])):
                return None
        except Exception:
            return None
        return copy

    def hollowed(self, brep, cutter):
        """Rhino's difference returns nothing for a cutter wholly inside the
        solid; the solid with the cutter turned inside out as a second shell is
        that difference."""
        shell = brep.DuplicateBrep()
        hole = cutter.DuplicateBrep()
        try:
            hole.Flip()
            shell.Append(hole)
        except Exception:
            return None
        return shell


def _face_samples(face, n):
    """(point, outward unit normal) on an n x n grid of the face's parameters,
    the points inside its trims only. BrepFace.NormalAt already points out of
    the solid on a reversed face (measured in Rhino 8), so it is not flipped."""
    found = []
    try:
        du, dv = face.Domain(0), face.Domain(1)
    except Exception:
        return found
    for i in range(n):
        for j in range(n):
            u = du.ParameterAt((i + 0.5) / n)
            v = dv.ParameterAt((j + 0.5) / n)
            try:
                if face.IsPointOnFace(u, v) != rg.PointFaceRelation.Interior:
                    continue
                p = face.PointAt(u, v)
                nv = face.NormalAt(u, v)
            except Exception:
                continue
            size = (nv.X * nv.X + nv.Y * nv.Y + nv.Z * nv.Z) ** 0.5
            if size > 0.0:
                s = 1.0 / size
                found.append(((p.X, p.Y, p.Z), (nv.X * s, nv.Y * s, nv.Z * s)))
    return found
