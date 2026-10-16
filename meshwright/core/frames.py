"""Rotation-minimizing frames along a rail, and the pieces of a frames + loft sweep.

A rotation-minimizing frame (RMF) turns with the rail and never spins about
it, so a section carried along it keeps its orientation the way a hand
would, with no twist the rail did not ask for. Frames are marched with the
double reflection method (Wang, Juettler, Zheng and Liu, 2008), which is
exact for circles and helices and very close elsewhere.

The rail is any object with these members (meshwright.rhino.sweep wraps a
Rhino curve; tests use exact circles and helices):

    domain               (t0, t1)
    closed               True for a closed rail
    length               its length
    point(t)             (x, y, z)
    tangent(t)           unit tangent
    side_tangents(t)     (unit tangent just below t, just above t)
    start_x()            a preferred x axis at the start, or None

Stations are sampled densely where the rail bends. A kink gives two stations
at the same point, one per side ("in" then "out"). On a closed rail the
frames come back around turned by some angle (the rail's total twist); that
angle is spread evenly along the rail so the seam closes without a jump.
"""

from __future__ import absolute_import, division

import math

from . import vec

MAX_TURN = math.radians(1.0)  # one station per degree of rail turning (or more)
PIECES = 240  # and one per 1/240 of its length
BASE = 2000  # parameter steps examined


class Station(object):
    """One frame along the rail.

    t, p    rail parameter and point
    T, X    unit tangent and frame x axis (the frame y is T x X)
    s       arc length from the start (measured along the stations)
    kink    None, "in" (arriving at a kink) or "out" (leaving it)
    q       indices of the sections attached here
    M       miter plane normal at a kink (the bisector of both tangents), or None
    """

    def __init__(self, t, p, T, kink=None, q=None):
        self.t, self.p, self.T = t, p, T
        self.kink = kink
        self.q = q or []
        self.X = None
        self.s = 0.0
        self.M = None

    def plane_axes(self):
        """(origin, x, y) of the section plane: z along the tangent."""
        return self.p, self.X, vec.cross(self.T, self.X)


class Frames(object):
    """stations; attach[i] is the station of section i; twist is the closed-rail
    angle spread along the rail (radians, 0 for an open rail)."""

    def __init__(self, stations, attach, twist):
        self.stations = stations
        self.attach = attach
        self.twist = twist


def any_perp(t):
    helper = (1.0, 0.0, 0.0) if abs(t[0]) < 0.9 else (0.0, 1.0, 0.0)
    return vec.unit(vec.cross(t, helper)) or (0.0, 0.0, 0.0)


def rotate_min(v, ta, tb):
    """v turned by the smallest rotation that takes unit vector ta onto tb."""
    axis = vec.cross(ta, tb)
    s = vec.length(axis)
    if s < 1e-12:
        return v
    return vec.rotate(v, axis, math.atan2(s, vec.dot(ta, tb)))


def sample_params(rail, max_turn=MAX_TURN, pieces=PIECES, base=BASE):
    """Rail parameters, denser where the rail bends."""
    t0, t1 = rail.domain
    step = rail.length / float(pieces)
    ts = [t0]
    last_tan, last_pt = rail.tangent(t0), rail.point(t0)
    turn = run = 0.0
    for i in range(1, base + 1):
        t = t0 + (t1 - t0) * i / float(base)
        tan, pt = rail.tangent(t), rail.point(t)
        turn += math.acos(max(-1.0, min(1.0, vec.dot(last_tan, tan))))
        run += vec.distance(pt, last_pt)
        last_tan, last_pt = tan, pt
        if turn >= max_turn or run >= step or i == base:
            ts.append(t)
            turn = run = 0.0
    return ts


def build_stations(rail, params, kinks=(), seam_kink=False, samples=None):
    """Frames along the rail with a station at every section parameter (params),
    two at every kink, and the samples in between (default: sample_params)."""
    t0, t1 = rail.domain
    items = [(t, 0, -1) for t in (samples if samples is not None else sample_params(rail))]
    items += [(t, 1, i) for i, t in enumerate(params)]
    items += [(t, 2, -1) for t in kinks]
    items += [(t0, 3, -1), (t1, 3, -1)]
    items.sort(key=lambda it: it[0])

    # Merge items at (almost) the same parameter; an exact one (section, kink, end) sets t.
    eps = 1e-7 * (t1 - t0)
    groups = []
    for t, kind, idx in items:
        if groups and abs(t - groups[-1][0]) <= eps:
            g = groups[-1]
        else:
            g = [t, [], False, False]  # t, sections, kink, exact
            groups.append(g)
        if kind != 0 and not g[3]:
            g[0], g[3] = t, True
        if kind == 1:
            g[1].append(idx)
        elif kind == 2:
            g[2] = True

    stations = []
    for t, q, kink, _ in groups:
        p = rail.point(t)
        if kink:
            tb, ta = rail.side_tangents(t)
            stations.append(Station(t, p, tb, "in"))
            stations.append(Station(t, p, ta, "out", q))
        else:
            tan = rail.tangent(t)
            if vec.length(tan) < 0.5 and stations:
                tan = stations[-1].T
            stations.append(Station(t, p, tan, None, q))

    # Starting x axis: the rail's own perpendicular frame when it has one.
    first = stations[0]
    x0 = rail.start_x() or any_perp(first.T)
    x0 = vec.unit(vec.sub(x0, vec.scale(first.T, vec.dot(x0, first.T))))
    first.X = x0 if x0 is not None else any_perp(first.T)

    tiny2 = (1e-9 * max(rail.length, 1.0)) ** 2
    for i in range(1, len(stations)):
        a, b = stations[i - 1], stations[i]
        v1 = vec.sub(b.p, a.p)
        c1 = vec.dot(v1, v1)
        if c1 <= tiny2:  # same point: the tangent jumps at a kink
            r = rotate_min(a.X, a.T, b.T)
        else:  # double reflection
            rl = vec.sub(a.X, vec.scale(v1, 2.0 * vec.dot(v1, a.X) / c1))
            tl = vec.sub(a.T, vec.scale(v1, 2.0 * vec.dot(v1, a.T) / c1))
            v2 = vec.sub(b.T, tl)
            c2 = vec.dot(v2, v2)
            r = rl if c2 < 1e-24 else vec.sub(rl, vec.scale(v2, 2.0 * vec.dot(v2, rl) / c2))
        r = vec.unit(vec.sub(r, vec.scale(b.T, vec.dot(r, b.T))))
        b.X = r if r is not None else any_perp(b.T)

    acc = 0.0
    for i in range(1, len(stations)):
        acc += vec.distance(stations[i].p, stations[i - 1].p)
        stations[i].s = acc

    twist = 0.0
    if rail.closed:
        if seam_kink:
            stations[0].kink = "out"
            stations[-1].kink = "in"
        a, b = stations[-1], stations[0]
        xe = rotate_min(a.X, a.T, b.T)
        phi = math.atan2(vec.dot(vec.cross(xe, b.X), b.T), vec.dot(xe, b.X))
        if acc > 0.0 and abs(phi) > 1e-12:
            for st in stations:
                st.X = vec.unit(vec.rotate(st.X, st.T, phi * st.s / acc))
            twist = phi

    n = len(stations)
    for i, st in enumerate(stations):
        if st.kink == "in":
            nxt = stations[(i + 1) % n]
            m = vec.unit(vec.add(st.T, nxt.T))
            st.M = m
            nxt.M = m

    attach = [None] * len(params)
    for st in stations:
        for q in st.q:
            attach[q] = st
    for i in range(len(attach)):
        if attach[i] is None:
            attach[i] = min(stations, key=lambda st: abs(st.t - params[i]))
    return Frames(stations, attach, twist)


def blend_at(pos, total, closed, sv):
    """Sections to blend at arc length sv, given the sections' arc lengths pos
    (sorted) and the rail length: (k0, k1, f), f of the way from k0 to k1."""
    n = len(pos)
    if n == 1:
        return 0, 0, 0.0
    if closed:
        for k in range(n):
            a = pos[k]
            b = pos[k + 1] if k + 1 < n else total + pos[0]
            for x in (sv, sv + total):
                if b > a and a - 1e-9 <= x <= b + 1e-9:
                    return k, (k + 1) % n, min(1.0, max(0.0, (x - a) / (b - a)))
        return 0, 0, 0.0
    if sv <= pos[0]:
        return 0, 0, 0.0
    for k in range(n - 1):
        if sv <= pos[k + 1]:
            span = pos[k + 1] - pos[k]
            return k, k + 1, (sv - pos[k]) / span if span > 0.0 else 1.0
    return n - 1, n - 1, 0.0


def split_runs(stations, closed):
    """Cut the stations at kinks: (runs, loft_closed). A smooth closed rail
    gives one run whose last station repeats the first."""
    runs = [[]]
    for st in stations:
        if st.kink == "out" and runs[-1]:
            runs.append([])
        runs[-1].append(st)
    if closed and stations[0].kink is None:
        if len(runs) == 1:
            return runs, True
        # The last station repeats the first; keep the first, it carries section 0.
        return [runs[-1][:-1] + runs[0]] + runs[1:-1], False
    return runs, False


def thin_run(run, extent, max_turn, max_step):
    """The stations a loft needs: both ends, every section, one per max_turn of
    tangent turning (and per max_step when blending sections). Stations right
    next to a mitered end are dropped so the loft cannot fold there."""
    n = len(run)
    if n <= 2:
        return list(run)
    dist = [0.0]
    for i in range(1, n):
        dist.append(dist[-1] + vec.distance(run[i].p, run[i - 1].p))

    def clearance(st):
        if not st.kink or st.M is None:
            return 0.0
        c = max(vec.dot(st.T, st.M), 1e-3)
        return 1.5 * extent * math.sqrt(max(0.0, 1.0 - c * c)) / c

    clear_start, clear_end = clearance(run[0]), clearance(run[-1])
    keep = [0]
    for i in range(1, n - 1):
        st = run[i]
        if st.q:
            keep.append(i)
            continue
        if dist[i] < clear_start or dist[-1] - dist[i] < clear_end:
            continue
        last = run[keep[-1]]
        turn = math.acos(max(-1.0, min(1.0, vec.dot(last.T, st.T))))
        if turn >= max_turn or dist[i] - dist[keep[-1]] >= max_step:
            keep.append(i)
    # No sliver span before the end: drop the last kept station if it sits
    # right next to it (sections always stay).
    if len(keep) > 1 and not run[keep[-1]].q:
        last, end = run[keep[-1]], run[n - 1]
        turn = math.acos(max(-1.0, min(1.0, vec.dot(last.T, end.T))))
        if turn < 0.5 * max_turn and dist[-1] - dist[keep[-1]] < 0.5 * max_step:
            keep.pop()
    keep.append(n - 1)
    return [run[i] for i in keep]
