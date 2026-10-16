"""Rhino Brep calls for solidwright.core.union_plan.

BrepOps answers union_plan's questions with RhinoCommon: Boolean union and
intersection, volume, bounding box, intersection curves, point inclusion and
a moved copy.
Every call that Rhino can fail returns None (or 0 for a volume) instead of
raising, so a failure is just one more refusing piece.
"""

from __future__ import absolute_import, division

import Rhino.Geometry as rg
from Rhino.Geometry.Intersect import Intersection

from ..core import union_plan


class BrepOps(union_plan.Ops):
    def __init__(self, tolerance):
        self.tolerance = tolerance

    def problem(self, brep):
        if not brep.IsValid:
            return "not valid"
        if not brep.IsSolid:
            return "open"
        return None

    def volume(self, brep):
        try:
            props = rg.VolumeMassProperties.Compute(brep)
        except Exception:
            return 0.0
        return 0.0 if props is None else abs(props.Volume)

    def box(self, brep):
        b = brep.GetBoundingBox(True)
        return (b.Min.X, b.Min.Y, b.Min.Z), (b.Max.X, b.Max.Y, b.Max.Z)

    def union(self, a, b):
        try:
            out = rg.Brep.CreateBooleanUnion([a, b], self.tolerance)
        except Exception:
            return None
        return None if out is None else list(out)

    def intersection_volume(self, a, b):
        try:
            out = rg.Brep.CreateBooleanIntersection(a, b, self.tolerance)
        except Exception:
            return None
        if out is None or len(out) == 0:
            return None
        return sum(self.volume(c) for c in out)

    def touch(self, a, b):
        """True when intersection curves exist. Touching at a point does not
        count: two solids meeting at a point cannot become one closed solid."""
        try:
            found = Intersection.BrepBrep(a, b, self.tolerance)
        except Exception:
            return False
        ok, curves = found[0], found[1]
        return bool(ok) and curves is not None and len(curves) > 0

    def encloses(self, outer, inner):
        """No intersection curves, and a point of inner is inside outer."""
        if self.touch(outer, inner):
            return False
        if inner.Vertices.Count > 0:
            point = inner.Vertices[0].Location
        else:
            face = inner.Faces[0]  # no vertices: an untrimmed closed face, so its middle is on it
            point = face.PointAt(face.Domain(0).Mid, face.Domain(1).Mid)
        try:
            return bool(outer.IsPointInside(point, self.tolerance, False))
        except Exception:
            return False

    def moved(self, brep, vector):
        copy = brep.DuplicateBrep()
        if not copy.Translate(rg.Vector3d(vector[0], vector[1], vector[2])):
            return None
        return copy
