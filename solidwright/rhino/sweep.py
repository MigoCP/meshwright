"""A one-rail sweep that tries hard and says how it went.

robust_sweep() runs Rhino's Sweep1 the careful way and falls back step by step:

1. Cleans the rail and sections: nested polycurves, segments shorter than
   the tolerance, ends that almost touch get closed (and reported).
2. Finds where each section sits on the rail (the rail through the section's
   plane, then curve intersections, then closest points) and hands those
   parameters to the solver instead of letting it guess.
3. Sorts the sections along the rail. On a closed rail the seam moves to the
   first section and the sweep is closed.
4. Matches the direction (and seam, for closed sections) of every section to
   the first one, compared in rotation-minimizing frames along the rail, so
   rail twist does not confuse it.
5. Tries solvers in order and keeps the first clean result (one valid Brep,
   closed when it should be): Sweep1 (with the other miter types when the
   rail has kinks), Sweep1 with sections refit, Sweep1 on a rebuilt rail, its
   own frames + loft (kinks mitered, closed-rail twist spread evenly), Sweep1
   at 10 x tolerance, CreateFromSweepSegmented, CreateFromSweep. On a tie the
   earlier method wins: the earlier ones stay closer to what was asked.

Every result is joined, kinky faces split, capped when asked, and flipped
outward when it is a solid.
"""

from __future__ import absolute_import, division

import math

import System
from System.Collections.Generic import List

import Rhino.Geometry as rg
from Rhino.Geometry.Intersect import Intersection

from ..core import frames as fr
from ..core import vec
from ..core.report import Report

MITERS = (0, 1, 2)  # 0 none, 1 trimmed, 2 untrimmed ("rotated" in Grasshopper)


class SweepResult(object):
    """breps (list), frames (Planes where each section attaches), method (the solver that won)."""

    def __init__(self, report):
        self.report = report
        self.breps = []
        self.frames = []
        self.method = None

    @property
    def status(self):
        return self.report.status


def _v3(v):
    return (v.X, v.Y, v.Z)


def _curve_list(curves):
    out = List[rg.Curve]()
    for c in curves:
        out.Add(c)
    return out


def _double_list(values):
    out = List[System.Double]()
    for v in values:
        out.Add(float(v))
    return out


def _brep_list(breps):
    out = List[rg.Brep]()
    for b in breps:
        out.Add(b)
    return out


class RhinoRail(object):
    """A Rhino curve as the rail solidwright.core.frames expects."""

    def __init__(self, curve):
        self.curve = curve
        d = curve.Domain
        self.domain = (d.T0, d.T1)
        self.closed = curve.IsClosed
        self.length = curve.GetLength()

    def point(self, t):
        return _v3(self.curve.PointAt(t))

    def tangent(self, t):
        return vec.unit(_v3(self.curve.TangentAt(t))) or (0.0, 0.0, 0.0)

    def side_tangents(self, t):
        return side_tangents(self.curve, t)

    def start_x(self):
        try:
            rc, plane = self.curve.PerpendicularFrameAt(self.domain[0])
        except Exception:
            return None
        return _v3(plane.XAxis) if rc else None


def to_plane(st):
    """A station as a Plane, Z along the rail tangent."""
    o, x, y = st.plane_axes()
    return rg.Plane(rg.Point3d(*o), rg.Vector3d(*x), rg.Vector3d(*y))


# ---------------------------------------------------------------- curve prep and placement

def clean_curve(crv, tol, name, report):
    c = crv.DuplicateCurve()
    try:
        c.RemoveNesting()  # polycurves only, harmless otherwise
    except Exception:
        pass
    try:
        if c.RemoveShortSegments(tol):
            report.note("%s: removed segments shorter than the tolerance" % name)
    except Exception:
        pass
    if not c.IsClosed:
        gap = c.PointAtStart.DistanceTo(c.PointAtEnd)
        if 0.0 < gap <= 10.0 * tol and c.GetLength() > 100.0 * tol:
            if c.MakeClosed(10.0 * tol):
                report.warn("%s: its ends were %.4g apart, so it was closed" % (name, gap))
    return c


def span_snap(crv, t, tol):
    """t moved onto an end of its own span when it lies within tol of it, so a
    new seam never leaves a sliver span behind."""
    try:
        p = crv.PointAt(t)
        for i in range(crv.SpanCount):
            dom = crv.SpanDomain(i)
            if dom.T0 <= t <= dom.T1:
                for u in (dom.T0, dom.T1):
                    if crv.PointAt(u).DistanceTo(p) <= tol:
                        return u
                break
    except Exception:
        pass
    return t


def reseam(crv, pt, tol):
    """Copy of a closed curve with its seam at the point closest to pt, or None."""
    for c in (crv.DuplicateCurve(), crv.ToNurbsCurve()):
        if c is None:
            continue
        rc, t = c.ClosestPoint(pt)
        if rc and c.ChangeClosedCurveSeam(span_snap(c, t, tol)):
            return c
    return None


def side_tangents(crv, t):
    """Unit tangents just below and just above parameter t."""
    try:
        below = crv.DerivativeAt(t, 1, rg.CurveEvaluationSide.Below)[1]
        above = crv.DerivativeAt(t, 1, rg.CurveEvaluationSide.Above)[1]
        tb, ta = vec.unit(_v3(below)), vec.unit(_v3(above))
        if tb is not None and ta is not None:
            return tb, ta
    except Exception:
        pass
    d = crv.Domain
    eps = 1e-6 * d.Length
    tb = vec.unit(_v3(crv.TangentAt(max(d.T0, t - eps)))) or (0.0, 0.0, 0.0)
    ta = vec.unit(_v3(crv.TangentAt(min(d.T1, t + eps)))) or (0.0, 0.0, 0.0)
    return tb, ta


def find_kinks(crv, ang):
    """Inner parameters where the tangent turns more than the angle tolerance."""
    d = crv.Domain
    cos_tol = math.cos(ang)
    kinks = []
    t0 = d.T0
    for _ in range(100000):
        rc, t = crv.GetNextDiscontinuity(rg.Continuity.C1_continuous, t0, d.T1)
        if not rc or t <= t0 or t >= d.T1:
            break
        tb, ta = side_tangents(crv, t)
        if vec.dot(tb, ta) < cos_tol:
            kinks.append(t)
        t0 = t
    return kinks


def seam_is_kink(crv, ang):
    ts = vec.unit(_v3(crv.TangentAtStart)) or (0.0, 0.0, 0.0)
    te = vec.unit(_v3(crv.TangentAtEnd)) or (0.0, 0.0, 0.0)
    return vec.dot(ts, te) < math.cos(ang)


def curve_plane(crv, tol):
    """Plane of a planar curve that is not a line, else None. TryGetPlane has its
    out parameter first: IronPython takes TryGetPlane(tol), Rhino 8's Python 3
    needs a placeholder plane, so both call styles are tried."""
    try:
        if crv.IsLinear(tol) or not crv.IsPlanar(tol):
            return None
    except Exception:
        return None
    for args in ((tol,), (rg.Plane.WorldXY, tol)):
        try:
            res = crv.TryGetPlane(*args)
            if res[0]:
                return res[1]
        except Exception:
            continue
    return None


def attach_param(rail, sec, tol):
    """Rail parameter a section belongs to:
    1. planar section: where the rail pierces its plane, nearest to the curve;
    2. otherwise: rail and section intersections;
    3. otherwise: the closest points between the two curves."""
    size = max(sec.GetBoundingBox(True).Diagonal.Length, tol)
    cands = []
    plane = curve_plane(sec, tol)
    if plane is not None:
        events = Intersection.CurvePlane(rail, plane, tol)
        if events is not None:
            for e in events:
                if e.IsPoint:
                    cands.append(e.ParameterA)
    if not cands:
        events = Intersection.CurveCurve(rail, sec, tol, tol)
        if events is not None:
            for e in events:
                cands.append(e.ParameterA)
    best = None
    for t in cands:
        p = rail.PointAt(t)
        rc, ts = sec.ClosestPoint(p)
        dist = p.DistanceTo(sec.PointAt(ts))
        if best is None or dist < best[0]:
            best = (dist, t)
    if best is not None and best[0] <= size:
        return best[1]
    rc, pa, pb = rail.ClosestPoints(sec)
    if rc:
        rc, t = rail.ClosestPoint(pa)
        if rc:
            return t
    rc, t = rail.ClosestPoint(sec.PointAtStart)
    return t


def place_sections(rail, secs, tol):
    """Sort sections along the rail; on a closed rail move its seam to the
    first section. Returns (rail, secs, params)."""
    params = [attach_param(rail, s, tol) for s in secs]
    if rail.IsClosed:
        # A section at the seam counts as the start, not the end: decided by distance on the rail.
        d = rail.Domain
        start, mid = rail.PointAtStart, d.Mid
        params = [d.T0 if t > mid and rail.PointAt(t).DistanceTo(start) <= tol else t for t in params]
    order = sorted(range(len(secs)), key=lambda i: params[i])
    secs = [secs[i] for i in order]
    params = [params[i] for i in order]
    if rail.IsClosed:
        moved = reseam(rail, rail.PointAt(params[0]), tol)
        if moved is not None:
            rail = moved
            d = rail.Domain
            start, mid = rail.PointAtStart, d.Mid
            rest = [attach_param(rail, s, tol) for s in secs[1:]]
            rest = [d.T1 if t < mid and rail.PointAt(t).DistanceTo(start) <= tol else t for t in rest]
            order = sorted(range(len(rest)), key=lambda i: rest[i])
            secs = [secs[0]] + [secs[1 + i] for i in order]
            params = [d.T0] + [rest[i] for i in order]
    return rail, secs, params


def align_sections(secs, planes, tol, report):
    """Match direction (and seam, when closed) of every section to the first,
    comparing them in their own rail frames so rail twist does not matter."""
    if len(secs) < 2:
        return secs
    ref, f0 = secs[0], planes[0]
    up = f0.ZAxis
    out = [ref]
    for i in range(1, len(secs)):
        s = secs[i].DuplicateCurve()
        probe = s.DuplicateCurve()
        probe.Transform(rg.Transform.PlaneToPlane(planes[i], f0))
        if ref.IsClosed:
            o_ref = ref.ClosedCurveOrientation(up)
            o_s = probe.ClosedCurveOrientation(up)
            und = rg.CurveOrientation.Undefined
            if o_ref != und and o_s != und and o_ref != o_s:
                s.Reverse()
                probe.Reverse()
                report.note("section %d: reversed to match section 0" % i)
            rc, t = probe.ClosestPoint(ref.PointAtStart)
            if rc:
                moved = reseam(s, s.PointAt(t), tol)
                if moved is not None:
                    s = moved
        else:
            a0, a1 = ref.PointAtStart, ref.PointAtEnd
            b0, b1 = probe.PointAtStart, probe.PointAtEnd
            if a0.DistanceTo(b1) + a1.DistanceTo(b0) < a0.DistanceTo(b0) + a1.DistanceTo(b1):
                s.Reverse()
                report.note("section %d: reversed to match section 0" % i)
        out.append(s)
    return out


# ---------------------------------------------------------------- solvers

def sweep_one_rail(rail, secs, params, closed, tol, ang, miter, refit=None):
    sweep = rg.SweepOneRail()
    sweep.SweepTolerance = tol
    sweep.AngleToleranceRadians = ang
    sweep.ClosedSweep = closed
    sweep.MiterType = miter
    if refit is None:
        return sweep.PerformSweep(rail, _curve_list(secs), _double_list(params))
    return sweep.PerformSweepRefit(rail, _curve_list(secs), _double_list(params), refit)


def rebuilt_rail(rail, tol):
    """(smooth evenly parameterized copy within tolerance, point count, deviation), or None."""
    for n in (12, 24, 48, 96, 192):
        try:
            r = rail.Rebuild(n, 3, True)
        except Exception:
            r = None
        if r is None:
            continue
        if rail.IsClosed and not r.IsClosed:
            r.MakeClosed(tol)
        res = rg.Curve.GetDistancesBetweenCurves(rail, r, 0.1 * tol)
        if res[0] and res[1] <= tol:
            return r, n, res[1]
    return None


def compatible(curves, tol, ang, report):
    """Section curves with the same degree, knots and point count, so they blend."""
    try:
        res = rg.NurbsCurve.MakeCompatible(_curve_list(curves), rg.Point3d.Unset, rg.Point3d.Unset,
                                           0, 0, tol, ang)
        if res is not None and len(res) == len(curves):
            res = list(res)
            n, deg = res[0].Points.Count, res[0].Degree
            if all(r.Points.Count == n and r.Degree == deg for r in res):
                return res
    except Exception:
        pass
    report.warn("frames + loft: sections rebuilt to 64 points so they can blend (sharp corners soften)")
    return [c.Rebuild(64, 3, False) for c in curves]


def blend(a, b, f):
    """Point-wise blend of two compatible NURBS curves."""
    if f <= 0.0:
        return a
    if f >= 1.0:
        return b
    c = rg.NurbsCurve(a)
    for i in range(a.Points.Count):
        pa, pb = a.Points[i], b.Points[i]
        la, lb = pa.Location, pb.Location
        loc = rg.Point3d(la.X + (lb.X - la.X) * f, la.Y + (lb.Y - la.Y) * f, la.Z + (lb.Z - la.Z) * f)
        c.Points[i] = rg.ControlPoint(loc, pa.Weight + (pb.Weight - pa.Weight) * f)
    return c


def manual_sweep(rail, secs, stations, attach, tol, ang, report):
    """Place the sections on the rotation-minimizing frames and loft. Kinks are
    mitered by projecting onto the bisector plane, so the pieces meet exactly."""
    world = rg.Plane.WorldXY
    closed = rail.IsClosed
    locs = []
    for s, st in zip(secs, attach):
        c = s.ToNurbsCurve()
        c.Transform(rg.Transform.PlaneToPlane(to_plane(st), world))
        locs.append(c)
    if len(locs) > 1:
        locs = compatible(locs, tol, ang, report)

    extent = 0.0
    for c in locs:
        for i in range(c.Points.Count):
            q = c.Points[i].Location
            extent = max(extent, math.sqrt(q.X * q.X + q.Y * q.Y + q.Z * q.Z))

    pos = [st.s for st in attach]
    total = stations[-1].s
    runs, loft_closed = fr.split_runs(stations, closed)
    max_turn = math.radians(4.0)
    max_step = total / 24.0 if len(locs) > 1 else float("inf")

    pieces = []
    for run in runs:
        straight = all(vec.dot(st.T, run[0].T) > 1.0 - 1e-9 for st in run)
        kept = fr.thin_run(run, extent, max_turn, max_step)
        if loft_closed:
            kept = kept[:-1]  # the last station repeats the first
        crvs = []
        for st in kept:
            k0, k1, f = fr.blend_at(pos, total, closed, st.s)
            c = rg.NurbsCurve(blend(locs[k0], locs[k1], f))
            c.Transform(rg.Transform.PlaneToPlane(world, to_plane(st)))
            if st.kink and st.M is not None:
                if vec.dot(st.T, st.M) > 0.1:
                    miter = rg.Plane(rg.Point3d(*st.p), rg.Vector3d(*st.M))
                    c.Transform(rg.Transform.ProjectAlong(miter, rg.Vector3d(*st.T)))
                else:
                    report.warn("frames + loft: hairpin kink near t = %.4g left unmitered" % st.t)
            crvs.append(c)
        if len(crvs) < 2:
            continue
        loft_type = rg.LoftType.Straight if straight else rg.LoftType.Normal
        lofts = rg.Brep.CreateFromLoft(_curve_list(crvs), rg.Point3d.Unset, rg.Point3d.Unset,
                                       loft_type, loft_closed)
        if lofts is None or len(lofts) == 0:
            raise ValueError("the loft failed on a rail piece")
        pieces.extend(list(lofts))
    return pieces


# ---------------------------------------------------------------- grading and main

def assess(result, tol, ang, closed_secs, closed_rail, cap_ends):
    """Join, tidy and grade a solver result: (score, breps), 0 is clean, or None.

    1000 per invalid Brep, 10 per extra piece, 5 when it should be a closed
    solid and is not."""
    if result is None:
        return None
    breps = [b for b in result if b is not None]
    if not breps:
        return None
    if len(breps) > 1:
        joined = rg.Brep.JoinBreps(_brep_list(breps), tol)
        if joined is not None and len(joined) > 0:
            breps = list(joined)
    out = []
    for b in breps:
        if not b.IsValid:
            try:
                b.Repair(tol)
            except Exception:
                pass
        if not b.IsSolid:
            try:
                b.JoinNakedEdges(tol)  # both ends of a sweep meeting at a seam miter
            except Exception:
                pass
        try:
            dup = b.DuplicateBrep()
            if dup.Faces.SplitKinkyFaces(ang, True) and dup.IsValid:
                b = dup
        except Exception:
            pass
        if cap_ends and closed_secs and not closed_rail and not b.IsSolid:
            capped = b.CapPlanarHoles(tol)
            if capped is not None:
                b = capped
        if b.IsSolid and b.SolidOrientation == rg.BrepSolidOrientation.Inward:
            b.Flip()
        out.append(b)
    score = 1000 * len([b for b in out if not b.IsValid]) + 10 * (len(out) - 1)
    if closed_secs and (closed_rail or cap_ends) and not all(b.IsSolid for b in out):
        score += 5
    return score, out


def robust_sweep(rail_curve, section_curves, miter=1, cap=False, tol=0.001, ang=math.radians(1.0)):
    """SweepResult for one rail and its sections (any order)."""
    report = Report("robust sweep")
    res = SweepResult(report)
    if miter not in MITERS:
        report.fail("the miter must be 0 (none), 1 (trimmed) or 2 (untrimmed), got %r" % (miter,))
        return res
    if rail_curve is None:
        report.fail("no rail curve")
        return res
    if not section_curves:
        report.fail("no section curves")
        return res
    report.note("tolerance %.4g, angle tolerance %.2f degrees" % (tol, math.degrees(ang)))
    rail = clean_curve(rail_curve, tol, "rail", report)
    secs = [clean_curve(s, tol, "section %d" % i, report) for i, s in enumerate(section_curves)]
    if len(set(s.IsClosed for s in secs)) > 1:
        report.fail("the sections must be all closed or all open")
        return res
    closed_secs = secs[0].IsClosed

    rail, secs, params = place_sections(rail, secs, tol)
    closed_rail = rail.IsClosed
    # Measured along the rail, so a rail that passes near itself is fine.
    gaps = [rail.GetLength(rg.Interval(a, b)) for a, b in zip(params[:-1], params[1:])]
    if closed_rail and len(params) > 1:
        gaps.append(rail.GetLength() - rail.GetLength(rg.Interval(params[0], params[-1])))
    if any(g <= tol for g in gaps):
        report.fail("two sections sit on the same spot of the rail")
        return res

    kinks = find_kinks(rail, ang)
    seam_kink = closed_rail and seam_is_kink(rail, ang)
    frames = fr.build_stations(RhinoRail(rail), params, kinks, seam_kink)
    res.frames = [to_plane(st) for st in frames.attach]
    secs = align_sections(secs, res.frames, tol, report)
    if frames.twist:
        report.note("rail: spread %.3f degrees of closed-rail twist evenly" % math.degrees(frames.twist))
    report.note("rail: %s, %d kink%s%s" % ("closed" if closed_rail else "open", len(kinks),
                                         "" if len(kinks) == 1 else "s", ", kink at the seam" if seam_kink else ""))
    if len(kinks) > 200:
        report.warn("the rail has %d kinks: if it should be smooth, Interpolate or Fit it first" % len(kinks))
    report.note("sections: %d %s, at rail t = %s" % (len(secs), "closed" if closed_secs else "open",
                                                    ", ".join("%.5g" % t for t in params)))

    def from_rebuilt_rail():
        got = rebuilt_rail(rail, tol)
        if got is None:
            raise ValueError("could not rebuild the rail within the tolerance")
        r, n, dev = got
        r, s2, p2 = place_sections(r, secs, tol)
        report.note("rail rebuilt with %d points, deviation %.2g" % (n, dev))
        return sweep_one_rail(r, s2, p2, r.IsClosed, tol, ang, miter)

    attempts = []
    miters = [miter]
    if kinks or seam_kink:  # no miter at all only when asked for
        miters += [m for m in (1, 2) if m != miter]
    for m in miters:
        attempts.append(("Sweep1, miter %d" % m,
                         lambda m=m: sweep_one_rail(rail, secs, params, closed_rail, tol, ang, m)))
    attempts.append(("Sweep1, sections refit",
                     lambda: sweep_one_rail(rail, secs, params, closed_rail, tol, ang, miter, refit=tol)))
    if not kinks and not seam_kink:
        attempts.append(("Sweep1, rebuilt rail", from_rebuilt_rail))
    attempts.append(("frames + loft",
                     lambda: manual_sweep(rail, secs, frames.stations, frames.attach, tol, ang, report)))
    attempts.append(("Sweep1, 10 x tolerance",
                     lambda: sweep_one_rail(rail, secs, params, closed_rail, 10.0 * tol, ang, miter)))
    attempts.append(("CreateFromSweepSegmented",
                     lambda: rg.Brep.CreateFromSweepSegmented(rail, _curve_list(secs), closed_rail, tol)))
    attempts.append(("CreateFromSweep",
                     lambda: rg.Brep.CreateFromSweep(rail, _curve_list(secs), closed_rail, tol)))

    best = None
    tried = []
    for name, fn in attempts:
        try:
            result = fn()
        except Exception as e:
            report.note("%s: %s" % (name, e))
            tried.append(name)
            continue
        graded = assess(result, tol, ang, closed_secs, closed_rail, cap)
        if graded is None:
            report.note("%s: no result" % name)
            tried.append(name)
            continue
        score, breps = graded
        report.note("%s: %d Brep%s, score %d" % (name, len(breps), "" if len(breps) == 1 else "s", score))
        if best is None or score < best[0]:
            best = (score, breps, name)
        if score == 0:
            break
        tried.append(name)

    if best is None:
        report.fail("every sweep method failed")
        return res
    score, breps, name = best
    res.breps, res.method = breps, name
    if score >= 1000:
        report.fail("only invalid Breps came out (best: %s)" % name)
    elif score > 0:
        report.messy("best effort with %s: the result is split or open" % name)
    elif name.startswith("Sweep1, 10 x"):
        report.messy("swept only at 10 x the tolerance")
    elif tried:
        report.warn("swept with %s (%s did not give a clean result)" % (name, ", ".join(tried)))
    else:
        report.note("swept with %s" % name)
    return res
