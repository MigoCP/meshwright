import math

from meshwright.core import deform


def close(a, b, eps=1e-12):
    return all(abs(x - y) < eps for x, y in zip(a, b))


def test_compress_moves_only_what_is_above_the_start():
    m, report = deform.compress(1.0, 0.5)
    assert report.status == 0
    assert m((2.0, 3.0, 0.5)) == (2.0, 3.0, 0.5)
    assert close(m((2.0, 3.0, 3.0)), (2.0, 3.0, 2.0))
    assert deform.compress(0.0, 0.0)[0] is None and deform.compress(0.0, -1.0)[1].failed


def test_curl_eases_in_and_keeps_lengths_on_the_hinge_circle():
    m = deform.Curl(start=0.0, length=2.0, radius=1.0, angle=math.radians(90))
    # Below the hinge and on the axis: no change.
    assert m((1.0, 0.0, -0.5)) == (1.0, 0.0, -0.5) and m((0.0, 0.0, 1.0)) == (0.0, 0.0, 1.0)
    # On the hinge circle at the top: rolled 90 degrees in towards the axis, height kept as a length.
    assert close(m((1.0, 0.0, 2.0)), (1.0 - 2.0, 0.0, 0.0))
    # Just above the hinge the turn is tiny (no crease).
    x, y, z = m((1.0, 0.0, 0.01))
    assert abs(x - 1.0) < 1e-6 and abs(z - 0.01) < 1e-9
    # Direction around the axis is kept.
    x, y, z = m((0.0, 2.0, 1.0))
    assert abs(x) < 1e-12 and y < 2.0


def test_curl_finds_the_hinge_radius_and_fails_cleanly():
    pts = [(1.0, 0.0, 1.0), (0.0, 3.0, 2.0), (5.0, 0.0, -1.0)]  # the last is below the start
    m, report = deform.curl(0.0, 2.0, math.radians(30), points=pts)
    assert m.radius == 2.0 and "hinge radius 2" in str(report)
    assert deform.curl(1.0, 1.0, 0.3, radius=1.0)[1].failed
    assert deform.curl(0.0, 1.0, 0.3, points=[(1.0, 0.0, -1.0)])[1].failed
    assert deform.curl(0.0, 1.0, 0.3, radius=0.0)[1].failed


def square(h):
    return [(-h, -h), (h, -h), (h, h), (-h, h)]


def circle(r, n=64):
    return [(r * math.cos(2 * math.pi * k / n), r * math.sin(2 * math.pi * k / n)) for k in range(n)]


def test_radial_fit_lands_one_outline_on_another():
    m, report = deform.radial_fit(square(1.0), circle(2.0, 720), start=0.0, fade=1.0)
    assert report.status == 0
    # A point of the square's side lands on the circle; at and below the start fully.
    x, y, z = m((1.0, 0.5, -0.3))
    assert abs(math.hypot(x, y) - 2.0) < 1e-3 and z == -0.3
    # Right at a sharp corner the bins (1 degree wide) round it a little: about 1% here.
    x, y, _ = m((1.0, 1.0, 0.0))
    assert abs(math.hypot(x, y) - 2.0) < 0.02
    # Halfway up the fade it is halfway (smoothstep at 0.5), and past the fade nothing moves.
    x, y, _ = m((1.0, 0.0, 0.5))
    assert abs(x - 1.5) < 1e-3
    assert close(m((1.0, 0.0, 1.5)), (1.0, 0.0, 1.5))


def test_radial_fit_reports_outlines_off_the_axis():
    off = [(x + 5.0, y) for x, y in square(1.0)]
    m, report = deform.radial_fit(off, circle(2.0), 0.0, 1.0)
    assert report.status == 2 and "miss an outline" in str(report)
    far = [(x, y + 9.0) for x, y in circle(1.0)]  # no ray from the axis meets both
    assert deform.radial_fit(off, far, 0.0, 1.0)[0] is None
    assert deform.radial_fit([(0, 0), (1, 0)], circle(1.0), 0.0, 1.0)[1].failed


def test_ray_hit_and_dome_levels():
    assert abs(deform.ray_hit(square(1.0), (1.0, 0.0)) - 1.0) < 1e-12
    assert abs(deform.ray_hit(square(1.0), (math.sqrt(0.5), math.sqrt(0.5))) - math.sqrt(2.0)) < 1e-12
    levels = deform.dome_levels(2.0, steps=2)
    assert close(levels[0], (math.cos(math.pi / 6), 1.0)) and close(levels[1], (0.5, 2.0 * math.sin(math.pi / 3)))
