# Verified union: "it ran" is not "it worked"

A Boolean union of many solids has two ways to go wrong. The loud one: one
pair with flush (coincident) faces is enough for a union over everything to
fail. The quiet one: the union returns a valid closed solid that is really
one input with the other silently dropped. The second is worse, because
nothing looks wrong until the part is missing.

`VerifiedUnion` (ported from an earlier C# Rhino command, BetterUnion)
unites one piece at a time, largest first, and checks every step before
keeping it.

## What a union must pass

One valid closed solid comes back, and:

1. **Its bounding box still covers both inputs' boxes.** The cheap tell of a
   dropped input.
2. **Its volume adds up**: at least as big as either input, no bigger than
   both together (within 1% of the smaller one).
3. **Each input adds volume.** An input that adds almost nothing is either
   inside the other (fine) or was dropped. A dropped piece can sit inside
   the other's bounding box, so check 1 misses it. The volume of the common
   part decides: union = a + b - common.

When the union refuses, the piece is set aside and retried once everything
else is in (it often unites fine once its neighbours have merged). Only
then is it nudged by 10 x the tolerance along +Z, -Z, +X, -X, +Y, -Y in
turn, so flush faces stop being coincident. A piece that still refuses is
left out and named, never dropped quietly. Results that grew into each
other are merged pairwise at the end, through the same checks.

## What Rhino actually returns

Measured with `Brep.CreateBooleanUnion` on two boxes in Rhino 8:

| Case | Result |
|---|---|
| overlapping | 1 Brep, volume 12 (correct) |
| sharing a face | 1 Brep (correct) |
| apart | 2 Breps, the inputs |
| one inside the other | 1 Brep, the outer one (correct) |
| touching along an edge only | 1 Brep that is **not a closed solid** |

Two of these shaped the checks:

- **Touching along an edge** gives a Brep that is not a closed solid. It is
  rejected, and the piece ends up reported as one that could not be united.
- **One inside the other**: the union is right, but the inner piece adds
  no volume, so the common part has to confirm it. `CreateBooleanIntersection`
  returns *nothing* for a piece fully inside another. The original command
  therefore rejected every enclosed piece. Now, when there is no common part
  to compute, the piece must lie inside the other: no intersection curves,
  and a point of it inside the other solid. Tested with a ball and a torus
  inside a box.

## Cost

Checking every step costs time. 51 parts (a plate, 25 pins standing on it
and 25 caps on the pins, every contact a flush face):

| | Result | Time |
|---|---|---:|
| one plain union | 1 solid, volume 104.6327 | 0.1 s |
| `VerifiedUnion` | 1 solid, volume 104.6327 | 3.2 s |

About 30 times slower, for a result that is known to contain every piece.
In this test the plain union happened to work; the point is that nothing
tells you when it does not.

## Testing without Rhino

The scheduler in `solidwright.core.union_plan` never calls Rhino. It takes the
geometry calls (union, volume, bounding box, intersection, touching, moving)
as an object, so the tests use stand-in solids made of unit cubes, where
union and volume are exact, and play Rhino's failures on purpose: a union
that refuses until a neighbour is in, one that refuses until the piece is
nudged, one that drops an input inside or outside the other's box, one with
a volume that does not add up, and one that returns an open solid.
