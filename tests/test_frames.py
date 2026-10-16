"""Rotation-minimizing frames on exact rails: circle, helix, a closed twisted loop, a kinked polyline."""

import math

import pytest

from meshwright.core import frames as fr
from meshwright.core import vec


class Analytic(object):
    """A rail from a point function and its derivative, on [0, 1]."""

    def __init__(self, point, deriv, closed, length):
        self._p, self._d = point, deriv
        self.closed = closed
        self.length = length
        self.domain = (0.0, 1.0)

    def point(self, t):
        return self._p(t)

    def tangent(self, t):
        return vec.unit(self._d(t))

    def side_tangents(self, t):
        return self.tangent(t), self.tangent(t)

    def start_x(self):
        return None


def circle(r=2.0):
    w = 2.0 * math.pi
    return Analytic(lambda t: (r * math.cos(w * t), r * math.sin(w * t), 0.0),
                    lambda t: (-r * w * math.sin(w * t), r * w * math.cos(w * t), 0.0), True, w * r)


def helix(r=1.0, pitch=0.5, turns=2.0):
    w = 2.0 * math.pi * turns
    h = pitch * turns
    return Analytic(lambda t: (r * math.cos(w * t), r * math.sin(w * t), h * t),
                    lambda t: (-r * w * math.sin(w * t), r * w * math.cos(w * t), h), False,
                    math.hypot(w * r, h))


def trefoil():
    """A trefoil knot: closed and twisted, so frames come back turned (about 2.2 radians)."""
    w = 2.0 * math.pi
    return Analytic(lambda t: (math.sin(w * t) + 2 * math.sin(2 * w * t), math.cos(w * t) - 2 * math.cos(2 * w * t),
                               -math.sin(3 * w * t)),
                    lambda t: (w * math.cos(w * t) + 4 * w * math.cos(2 * w * t),
                               -w * math.sin(w * t) + 4 * w * math.sin(2 * w * t), -3 * w * math.cos(3 * w * t)),
                    True, 30.0)


class Polyline(object):
    """Parameter k..k+1 runs along segment k; kinks at the inner corners."""

    def __init__(self, points):
        self.pts = points
        self.closed = False
        self.domain = (0.0, float(len(points) - 1))
        self.length = sum(vec.distance(a, b) for a, b in zip(points, points[1:]))

    def _seg(self, t, side=0.0):
        k = int(math.floor(t + side))
        return max(0, min(len(self.pts) - 2, k))

    def point(self, t):
        k = self._seg(t)
        return vec.lerp(self.pts[k], self.pts[k + 1], t - k)

    def tangent(self, t):
        k = self._seg(t)
        return vec.unit(vec.sub(self.pts[k + 1], self.pts[k]))

    def side_tangents(self, t):
        kb, ka = self._seg(t, -1e-9), self._seg(t, 1e-9)
        return (vec.unit(vec.sub(self.pts[kb + 1], self.pts[kb])),
                vec.unit(vec.sub(self.pts[ka + 1], self.pts[ka])))

    def start_x(self):
        return (0.0, 0.0, 1.0)


def check_frames(stations):
    for st in stations:
        assert abs(vec.length(st.X) - 1.0) < 1e-9 and abs(vec.dot(st.X, st.T)) < 1e-9


def spin(stations):
    """Largest turn about the tangent between neighbours, beyond the minimal rotation."""
    worst = 0.0
    for a, b in zip(stations, stations[1:]):
        r = fr.rotate_min(a.X, a.T, b.T)
        worst = max(worst, abs(math.atan2(vec.dot(vec.cross(r, b.X), b.T), vec.dot(r, b.X))))
    return worst


def test_circle_frames_do_not_spin_and_close():
    f = fr.build_stations(circle(), [0.0, 0.5])
    check_frames(f.stations)
    assert abs(f.twist) < 1e-6  # a flat closed curve has no twist to spread
    assert spin(f.stations) < 1e-6
    assert f.attach[0].t == 0.0 and f.attach[1].t == 0.5


def test_helix_frames_are_rotation_minimizing():
    f = fr.build_stations(helix(), [0.0, 0.3, 1.0])
    check_frames(f.stations)
    assert spin(f.stations) < 1e-6
    assert f.stations[-1].s == pytest.approx(helix().length, rel=1e-4)


def test_closed_knot_spreads_its_twist():
    rail = trefoil()
    f = fr.build_stations(rail, [0.0])
    check_frames(f.stations)
    assert abs(f.twist) > 1.0
    # After spreading, the frame at the end turns back onto the frame at the start.
    a, b = f.stations[-1], f.stations[0]
    back = fr.rotate_min(a.X, a.T, b.T)
    assert vec.distance(back, b.X) < 1e-6
    # ...and the spin per station is the twist shared out evenly.
    assert spin(f.stations) < abs(f.twist) / 50.0


def test_kinks_give_two_stations_a_miter_and_two_runs():
    rail = Polyline([(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (2.0, 2.0, 0.0)])
    f = fr.build_stations(rail, [0.0, 1.5], kinks=[1.0])
    kinked = [st for st in f.stations if st.t == 1.0]
    assert [st.kink for st in kinked] == ["in", "out"]
    m = kinked[0].M
    assert vec.distance(m, vec.unit((1.0, 1.0, 0.0))) < 1e-12 and kinked[1].M == m
    assert kinked[1].X == pytest.approx((0.0, 0.0, 1.0))  # turned with the rail, not spun
    runs, loft_closed = fr.split_runs(f.stations, False)
    assert len(runs) == 2 and not loft_closed
    assert runs[1][0].kink == "out"


def test_closed_smooth_rail_is_one_closed_loft():
    f = fr.build_stations(circle(), [0.0])
    runs, loft_closed = fr.split_runs(f.stations, True)
    assert len(runs) == 1 and loft_closed


def test_thin_run_keeps_ends_sections_and_turns():
    f = fr.build_stations(circle(), [0.0, 0.25])
    run = f.stations
    kept = fr.thin_run(run, 0.1, math.radians(10.0), float("inf"))
    assert kept[0] is run[0] and kept[-1] is run[-1]
    assert f.attach[1] in kept
    assert 34 <= len(kept) <= 40  # about one per 10 degrees around the circle
    straight = fr.build_stations(Polyline([(0, 0, 0), (5, 0, 0)]), [0.0])
    assert len(fr.thin_run(straight.stations, 0.1, math.radians(4.0), float("inf"))) == 2


def test_blend_at():
    assert fr.blend_at([0.0], 10.0, False, 4.0) == (0, 0, 0.0)
    assert fr.blend_at([0.0, 10.0], 10.0, False, 2.5) == (0, 1, 0.25)
    assert fr.blend_at([2.0, 6.0], 10.0, False, 1.0) == (0, 0, 0.0)
    assert fr.blend_at([2.0, 6.0], 10.0, False, 8.0) == (1, 1, 0.0)
    # Closed: past the last section it blends back to the first across the seam.
    k0, k1, f = fr.blend_at([0.0, 4.0], 10.0, True, 7.0)
    assert (k0, k1) == (1, 0) and abs(f - 0.5) < 1e-12
