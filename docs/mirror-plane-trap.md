# The mirror-plane trap

`NaturalAxis` finds the direction a part "obviously" sits along, from the
mesh alone. Its main test is whether a direction looks *level*: look down
the direction, take the points on the outline (the silhouette), fit a plane
through them, and check that the plane is square to the direction. Looking
down the axis of a turned or faceted part, the outline is its widest band,
and that band is flat and square to the axis.

The trap: **for a symmetric part, the view square to any mirror plane also
looks level.** The silhouette seen square to a mirror plane is symmetric
about that plane, so the plane fitted through it is the mirror plane itself,
exactly square to the view. A part with 16 sides has 16 such side views, and
each one passes the level test as well as the true axis does.

## Numbers

A 16-sided part with a flat top, its widest band just below the top and a
long point 2 units below that (`tall_tip_part` in `tests/test_orient.py`):

| View | Tilt from level | Level? | Outline area | Strong flat end |
|---|---:|---|---:|---|
| down the true axis | 0.000 deg | yes | 0.765 | yes |
| square to a mirror plane | 0.000 deg | yes | **1.128** | no |
| 3 deg off a mirror plane | 1.496 deg | yes (limit 1.5) | 1.127 | no |
| tilted 20 deg from the axis | no plane | no | 0.788 | no |

Sweeping the side view around the axis in 1 degree steps, **168 of 360
side views count as level**: the 16 mirror views, plus everything close
enough to one that the tilt stays under the 1.5 degree limit.

Worse, the side view shows 47% *more* outline than the true axis, because
the long point makes the part tall. A rule like "among level directions,
pick the largest outline" picks a side view here, 90 degrees off.

## What decides it

1. **A strong flat end wins.** A flat facet square to the direction that
   fills most of that end of the outline is the best evidence there is. The
   top view has one; no side view does. This picks the true axis for the
   part above.
2. **Otherwise the largest outline wins**, which is right for parts with no
   flat end (a lens, a double cone).
3. **A near tie is reported, not guessed.** A tall square bar has two side
   axes 90 degrees apart that fit equally well; `NaturalAxis` returns status
   2 (messy) with "two directions 90 degrees apart fit about equally well"
   instead of picking one at random.

The test `test_mirror_plane_view_looks_level_but_is_not_chosen` pins all of
this down: the side view is level, it shows more outline, and the flat end
is what decides. The same part placed at 8 random orientations, stored in
single precision, comes back within 0.01 degrees of the true axis.
