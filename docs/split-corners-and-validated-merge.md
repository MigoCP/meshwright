# Split corners and the validated merge

`OffsetSolid` pushes every facet of a part out by a clearance. With an axis,
a facet moves by the clearance times the sideways share of its normal: a
side wall moves the full clearance, a sloped facet less, a flat end not at
all. Each facet keeps its direction, and the solid is rebuilt from the moved
planes.

Where three planes meet at a corner, moving them still leaves one corner.
Where **four or more** meet, moving them by different amounts means they no
longer pass through one point: the corner **splits** into a cluster of
corners a hair apart, joined by edges far shorter than anything you would
model. The result is a valid solid, but those edges are trouble downstream
(Booleans, fillets, meshing).

## Numbers

A faceted part with 65 facets, many of its corners shared by four or more:

| Clearance (sideways) | Corners after the offset | Shortest edge | After the merge | Shortest edge |
|---|---:|---:|---:|---:|
| 0.05 | 97 | 0.005 | 49 | 0.020 |
| 0.2 | 97 | 0.018 | 33 | 0.185 |

At clearance 0.2, the merge also removed 16 faces that had shrunk to
slivers (65 faces down to 49), and the result still passes every check
below.

## Why the merge has to be checked

Merging means moving corners, and moving a corner can bend the faces around
it out of flat, flip an edge, or pull a face in closer than the clearance.
So a merge is only kept when the merged solid passes `validate()`:

- it is still closed (every edge paired with one running the other way),
- every face is flat within half the tolerance (its plane is refitted
  through its new corners),
- every face still keeps at least 90% of the clearance from the part.

`offset_solid()` tries merge distances from large to small (a quarter of
the clearance, a tenth, then the tolerance) and keeps the **biggest one that
passes**. Merging is transitive (union-find), so a chain of close corners
becomes one corner even when its ends are further apart than the distance.

An oval part with 169 facets and clearance 0.1 in every direction shows the
check at work:

| Merge within | Corners | Result |
|---|---:|---|
| 0.025 | 301 to 147 | refused: a face 0.0064 off flat (allowed 0.0005) |
| 0.01 | 301 to 167 | refused: a face 0.00075 off flat |
| 0.001 (the tolerance) | 301 to 289 | kept, shortest edge 0.00103 |

At clearance 0.4 the biggest merge also breaks the solid itself ("8 edges
run the same way twice") and is refused the same way. When even the merge
at the tolerance fails, the unmerged solid is returned with a warning that
says why, instead of a merged solid that is quietly wrong.

`MeshToSolid` meets the same problem from the other side: noise splits the
corners of a rebuilt part, and its corners are merged within the tolerance
only when the result validates.
