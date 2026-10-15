# meshwright

Turn rough meshes into clean, verified solids for Rhino and Grasshopper.

Work in progress. The full README (what it does, a before/after comparison,
install steps for GHPython and Rhino 8, supported Rhino versions) arrives with
the first release.

## Principles

- **Every result is graded.** Tools return a report with a status (0 good,
  1 warning, 2 messy, 3 failed) and readable reasons, instead of raising or
  quietly handing back bad geometry.
- **Pure-Python core.** `meshwright.core` has no Rhino and no numpy imports, so
  it runs in GHPython (IronPython 2.7), Rhino 8 CPython and plain CPython.
  Thin adapters in `meshwright.rhino` convert to and from RhinoCommon.
- **Tested without Rhino.** Tests use procedural geometry and a small mock of
  the RhinoCommon pieces the adapters touch.

## Planar mesh to solid, compared with Rhino

`MeshToSolid` turns a closed mesh of a part made of flat faces into a Brep
with one exact flat face per planar region (parts that are not convex and
faces with holes included). Rhino's own route is `ToNURBS` followed by
`MergeAllCoplanarFaces`. Both were run in Rhino 8 on the same meshes, each
tilted 37 degrees and moved off the origin, with a tolerance of 0.001 and
Rhino's default 1 degree angle tolerance. The test parts are 1 to 4 units
across; the meshes are dense (every face split into many triangles), with
random noise added to the vertices or stored in single precision far from
the origin, as an STL export there would be.

*Corners off* is how far the result's corners are from the true solid (the
largest distance, in model units).

| Part | Mesh | Right faces | meshwright faces | corners off | time | Rhino faces | corners off | time |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| box | noise 1e-4 | 6 | 6 | 7.5e-5 | 0.01 s | 6 | 1.6e-4 | 0.06 s |
| L block | noise 1e-4 | 8 | 8 | 1.3e-4 | 0.02 s | 8 | 1.8e-4 | 0.13 s |
| U block | noise 1e-4 | 10 | 10 | 1.0e-4 | 0.03 s | 10 | 2.4e-4 | 0.22 s |
| frame (square hole) | noise 1e-4 | 10 | 10 | 1.3e-4 | 0.04 s | 10 | 1.7e-4 | 0.29 s |
| faceted part | noise 1e-4 | 65 | 65 | 1.5e-4 | 0.16 s | **221** | 3.0e-4 | 3.0 s |
| plate, round hole | noise 1e-4 | 30 | 30 | 8.7e-5 | 0.16 s | **74** | 2.9e-4 | 4.3 s |
| box | single precision, far | 6 | 6 | 9.9e-6 | 0.01 s | 6 | 6.1e-5 | 0.06 s |
| L block | single precision, far | 8 | 8 | 2.5e-5 | 0.02 s | 8 | 3.6e-5 | 0.13 s |
| U block | single precision, far | 10 | 10 | 2.5e-5 | 0.02 s | 10 | 5.0e-5 | 0.21 s |
| frame (square hole) | single precision, far | 10 | 10 | 1.3e-5 | 0.03 s | 10 | 6.1e-5 | 0.27 s |
| faceted part | single precision, far | 65 | 65 | 3.0e-5 | 0.10 s | **91** | 5.9e-5 | 2.6 s |
| plate, round hole | single precision, far | 30 | 30 | 2.9e-5 | 0.15 s | **35** | 6.2e-5 | 4.0 s |

- **Clean meshes:** both give the same faces and edges, and both are exact
  (corners within 2e-14).
- **Light noise (1e-5):** both find the right faces. meshwright's corners
  are 1.1 to 2.5 times closer, because each corner is solved from the
  refitted planes that meet there instead of being taken from one noisy
  vertex.
- **Heavier noise (1e-4) and single precision far from the origin:**
  Rhino leaves faces split into pieces on the parts with many or narrow
  faces. meshwright keeps them whole by joining neighbouring pieces while
  they stay flat within half the tolerance. Its corners are closer in every
  case, up to 6 times.
- **Validity:** every result from both is a valid closed solid with no open
  edges. The largest edge tolerance in a meshwright result was 7.9e-5 (the
  noisy faceted part); Rhino's were 0.
- **A Rhino-meshed CAD part** (a box with a cylinder hole, single precision
  render mesh): both give 46 faces; corners off 1.6e-7 (meshwright) and
  3.4e-7 (Rhino).
- **Speed:** times are for Rhino 8's CPython. In IronPython (GHPython)
  meshwright took about the same, after about 0.3 s of start-up on its
  first run. The real `ToNURBS` and `MergeAllCoplanarFaces` commands gave
  the same results as their RhinoCommon calls, slightly slower.
- **Noise beyond the tolerance:** when the noise is too large for faces to
  be flat within half the tolerance, meshwright cannot join them either. It
  then reports that the mesh looks noisier than the tolerance instead of
  grading the result good.

## Running the tests

```
python -m pip install pytest numpy scipy
python -m pytest
```
