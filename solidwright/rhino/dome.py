"""Rounded caps on flat outlines, and parts cut off and capped with them.

dome() lofts an outline through shrinking, rising copies of itself (a
quarter ellipse in profile) to an apex above its centroid, so the cap meets
the outline's plane square on. cap_part() cuts a part with a plane, keeps
the side below it, caps every section with a dome and joins the lot into
one closed Brep, saying what went wrong when it does not close.
"""

from __future__ import absolute_import, division

import Rhino.Geometry as rg
from Rhino.Geometry.Intersect import Intersection

from ..core import deform
from ..core.report import Report
from .sweep import curve_plane


def dome(outline, ratio, steps, tol, up=None):
    """(list of Breps or None, height, reason) for a dome on a closed planar outline.

    up picks the side: the dome rises on the side of the outline's plane it
    points to (default: the side of world +Z; the plane normal for a
    vertical outline).
    """
    if not outline.IsClosed:
        return None, 0.0, "the outline is not closed"
    plane = curve_plane(outline, tol)
    if plane is None:
        return None, 0.0, "the outline is not planar"
    props = rg.AreaMassProperties.Compute(outline)
    if props is None or props.Area <= 0.0:
        return None, 0.0, "the outline has no area"
    n = plane.ZAxis
    hint = up if up is not None else rg.Vector3d.ZAxis
    if n * hint < 0.0:
        n = -n
    base = rg.Plane(props.Centroid, n)
    height = deform.dome_height(props.Area, ratio)
    profiles = [outline]
    for scale, z in deform.dome_levels(height, steps):
        c = outline.DuplicateCurve()
        c.Transform(rg.Transform.Translation(n * z) * rg.Transform.Scale(base, scale, scale, 1.0))
        profiles.append(c)
    apex = props.Centroid + n * height
    lofts = rg.Brep.CreateFromLoft(profiles, rg.Point3d.Unset, apex, rg.LoftType.Normal, False)
    if not lofts:
        return None, height, "the loft failed"
    return list(lofts), height, None


def sections(brep, plane, tol, report):
    """Closed outlines where the plane cuts the Brep, holes refused."""
    rc, curves, _ = Intersection.BrepPlane(brep, plane, tol)
    if not rc or not curves:
        return []
    raw = list(curves)
    loops = rg.Curve.JoinCurves(raw, 5.0 * tol) or raw
    closed = []
    for c in loops:
        if not c.IsClosed:
            c.MakeClosed(10.0 * tol)  # heal hairline gaps at seams
        if c.IsClosed:
            closed.append(rg.Curve.ProjectToPlane(c, plane) or c)
        else:
            report.messy("a section stays open (%.4g long): left without a dome" % c.GetLength())
    for i in range(len(closed)):
        for j in range(len(closed)):
            if i != j and rg.Curve.PlanarClosedCurveRelationship(closed[i], closed[j], plane, tol) == \
                    rg.RegionContainment.AInsideB:
                report.fail("the cut has a hole (the part is hollow there): a dome would close over it")
                return []
    if len(closed) > 1:  # overlapping outlines (overlapping pieces in the input) become one
        merged = rg.Curve.CreateBooleanUnion(closed, tol)
        if merged and len(merged) > 0:
            closed = list(merged)
    out = []
    for c in closed:
        props = rg.AreaMassProperties.Compute(c)
        if props is None or props.Area < 100.0 * tol * tol:
            report.note("skipped a section too small to cap")
            continue
        out.append(c)
    return out


def cap_part(brep, plane, ratio, steps, tol):
    """(breps, domes, outlines, report): the part below the plane, capped with domes and joined."""
    report = Report("dome cap")
    box = brep.GetBoundingBox(True)
    size = 2.0 * box.Diagonal.Length + abs(plane.DistanceTo(box.Center))
    cutter = rg.PlaneSurface(plane, rg.Interval(-size, size), rg.Interval(-size, size)).ToBrep()
    outlines = sections(brep, plane, tol, report)
    if report.failed:
        return [], [], outlines, report
    if not outlines:
        report.fail("the plane does not cut the part")
        return [], [], [], report
    pieces = brep.Split(cutter, tol)
    if not pieces:
        report.fail("Rhino could not split the part with the plane")
        return [], [], outlines, report
    lows = [p for p in pieces if plane.DistanceTo(p.GetBoundingBox(True).Center) < 0.0]
    domes = []
    for k, c in enumerate(outlines):
        made, height, why = dome(c, ratio, steps, tol, up=plane.ZAxis)
        if made is None:
            report.messy("section %d: %s" % (k, why))
            continue
        domes.extend(made)
        report.note("section %d: dome %.4g high" % (k, height))
    group = lows + domes
    joined = rg.Brep.JoinBreps(group, tol)
    ok = bool(joined) and all(j.IsSolid for j in joined)
    if not ok:
        wider = rg.Brep.JoinBreps(group, 3.0 * tol)
        if wider and all(j.IsSolid for j in wider):
            joined, ok = wider, True
            report.warn("the domes joined only at 3 x the tolerance")
    if not joined:
        report.messy("the pieces did not join: output separately")
        return group, domes, outlines, report
    joined = list(joined)
    naked, naked_len = 0, 0.0
    for j in joined:
        if j.IsSolid and j.SolidOrientation == rg.BrepSolidOrientation.Inward:
            j.Flip()
        edges = j.DuplicateNakedEdgeCurves(True, True) or []
        naked += len(edges)
        naked_len += sum(e.GetLength() for e in edges)
    if not ok:
        report.messy("not closed: %d naked edge%s, %.4g long in all" % (
            naked, "" if naked == 1 else "s", naked_len))
    report.note("%d piece%s joined into %d Brep%s" % (
        len(group), "" if len(group) == 1 else "s", len(joined), "" if len(joined) == 1 else "s"))
    return joined, domes, outlines, report
