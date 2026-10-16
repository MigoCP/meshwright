# Why interpolating between sections cuts into the part

`ClearanceSolid` describes a part by cross sections square to its axis and
joins neighbouring sections with flat faces. Between two levels, every
section edge slides linearly with height. That is exact across a flat side
facet: in cross section a facet is a line whose offset changes linearly with
height. It stops being exact wherever the side **bends between two levels**.

A convex profile bulges outward past the straight line (the chord) between
any two of its points. So when the side bends between two levels, the
chord, and the solid built on it, lies *inside* the part there. It is not a
small rounding error: the solid cuts into the part by the full bulge.

## Numbers

A 16-sided part: flat top (radius 0.3, height 0.2), widest band (radius 0.5,
from height 0 to 0.02), a bend in the lower taper (radius 0.35 at height
-0.2) and a point at -0.6. The chord from the band to the point passes
height -0.2 at radius 0.333, while the part reaches 0.35 there.

Smallest sideways clearance measured at every part vertex (negative: the
part pokes out of the solid):

| Levels used | Clearance 0 | Clearance 0.03 |
|---|---:|---:|
| band and ends only | **-0.0175** (4 levels, 48 vertices out) | 0.0125, 42% of it (4 levels) |
| + bends found in the area profile | -0.0022 (5 levels, 32 vertices out) | 0.0278 (5 levels) |
| + refinement at the worst vertices | -0.0002, within tolerance (6 levels) | 0.0278 (5 levels: none needed) |

A ball (24 sides, 12 bands) with clearance 0.03, where the whole side bends:

| Levels used | Levels | Smallest clearance |
|---|---:|---:|
| band and ends only | 4 | -0.209 |
| + bends | 14 | 0.0260 |
| + refinement | 16 | 0.0283 |

## What the solid does about it

1. **Bends.** The square root of the section area is a straight line along
   a run of similar sections, so where it stops being straight, the profile
   bends. A level goes there.
2. **Refinement.** The area averages all sides together, so a bend on one
   side can hide. After the solid is built, every part vertex is checked;
   any that comes out closer than 90% of the clearance asks for a measured
   level at its own height (the worst one per span), and the check repeats.
3. **Grading.** If a vertex is still outside after refinement, the result
   fails ("the part pokes out of the solid at N vertices"). If the clearance
   drops below half of what was asked, it is a warning that says how far.

The same lesson applies to any model built from sections or tiers: it needs
a level at every height where any side's profile bends, and the only
reliable way to find them all is to check the result against the part.
