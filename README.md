# solidwright

[![tests](https://github.com/MigoCP/solidwright/actions/workflows/tests.yml/badge.svg)](https://github.com/MigoCP/solidwright/actions/workflows/tests.yml)
[![version](https://img.shields.io/badge/version-0.2.0-blue)](https://github.com/MigoCP/solidwright/tags)
![Rhino 8](https://img.shields.io/badge/Rhino-8-informational)
![Python](https://img.shields.io/badge/python-IronPython%202.7%20%7C%203.9%2B-informational)
[![license](https://img.shields.io/badge/license-MIT-blue)](LICENSE)

Verified solid modelling for Rhino and Grasshopper: mesh to solid, checked
Booleans, robust sweeps. Every result graded.

![A noisy mesh, Rhino's ToNURBS + MergeAllCoplanarFaces result and solidwright's MeshToSolid result turning side by side, one colour per Brep face](docs/demo.gif)

*A noisy mesh of a plate with a round hole, rebuilt by Rhino (74 faces, its big
faces left split into pieces) and by solidwright's `MeshToSolid` (30 faces, one
per flat region). Each Brep face has its own colour. Captured in Rhino 8.*

solidwright is a set of Grasshopper components and the Python library behind
them. It checks meshes, finds how a part sits, rebuilds faceted meshes as
exact solids, builds clearance solids around parts, and wraps Rhino
operations that fail quietly (Boolean union and difference, Sweep1) in
checks that say when they did. Every component reports how it went.

## Components

**Meshes and parts**

| Component | What it does | More |
|---|---|---|
| **MeshHealth** | Checks a mesh and cleans what is safe to clean (welds vertices, removes zero-area and duplicate faces, drops loose junk, fixes face directions), then grades what is left: holes, extra pieces, non-manifold edges. | |
| **NaturalAxis** | The direction a part obviously sits along (the axis of a turned or faceted body, the normal of a plate), at any orientation, or a check of an axis you give it. | [mirror-plane trap](docs/mirror-plane-trap.md) |
| **ProfileFit** | A part's side outline and where it bends, as a three-piece broken line. | |
| **SizeGroups** | Sorts parts into size groups, by reference parts, by a count, or by gaps in size. | |

**Solids**

| Component | What it does | More |
|---|---|---|
| **MeshToSolid** | A closed Brep with one exact flat face per planar region of a mesh, convex or not, holes included. Where Rhino's `ToNURBS` + `MergeAllCoplanarFaces` leaves big faces split (74 faces for the plate above), it gives the right 30, with corners up to 6 times closer. | [compared with Rhino](docs/mesh-to-solid-vs-rhino.md) |
| **OffsetSolid** | A closed solid around a faceted convex part, every facet pushed out by a clearance (sideways only, or all round). Corners the offset splits are merged only when the result passes its checks. | [split corners](docs/split-corners-and-validated-merge.md) |
| **ClearanceSolid** | A closed solid around any part, curved and scanned ones included, built from cross sections kept a clearance away sideways. | [cutting in between sections](docs/interpolation-cuts-between-kinks.md) |
| **DomeCap** | Rounded caps on flat outlines, or a part cut by a plane and capped. | |

**Booleans, sweeps and shaping**

| Component | What it does | More |
|---|---|---|
| **VerifiedUnion** | Boolean union of many solids, one at a time, every result checked. Catches a piece the union dropped without saying so, and solids touching along an edge that give an open Brep. | [verified union](docs/verified-union.md) |
| **VerifiedDifference** | Boolean difference with many cutters, the ones that break it isolated and rescued, every result checked. Catches a cutter the difference skipped without saying so, and open Breps. | [verified difference](docs/verified-difference.md) |
| **RobustSweep** | A one-rail sweep that places and aligns the sections itself and falls back through seven methods, reporting which one worked; Info lists every attempt. | |
| **Deform** | Compresses, curls or radially fits parts above a height along an axis. | |
| **StackAlongAxis** | Lines parts up end to end along any axis with a fixed gap. | |

![A plate with 25 holes and a ball inside, cut by Rhino's CreateBooleanDifference and by VerifiedDifference, both shown with a section plane moving in to the ball's centre: Rhino's result is solid there, VerifiedDifference's shows the hollow](docs/difference.gif)

*One Boolean difference with 26 cutters: 25 holes and a ball of radius 0.6
inside the plate. Rhino returns a valid closed solid with the ball skipped
and no error (volume 174.867). `VerifiedDifference` catches the skip and cuts
the ball as a sealed hollow (volume 173.962). Captured in Rhino 8.*

Every component has three extra outputs:

| Output | Meaning |
|---|---|
| `Status` | 0 good, 1 warning (usable, worth a look), 2 messy (made from questionable input or by a fallback), 3 failed (no trustworthy result) |
| `OK` | True for status 0 or 1 |
| `Info` | every reason, one per line |

Warnings and errors also show as the component's orange or red balloon.

## Install

solidwright is plain Python with no dependencies. Download or clone this
repository, then make the `solidwright` folder (the one inside the
repository, holding `core`, `rhino` and `gh`) importable:

- **Rhino 8**: copy the `solidwright` folder into
  `%APPDATA%\McNeel\Rhinoceros\8.0\scripts`. Rhino 8 puts that folder on the
  search path of both its Pythons (IronPython 2 and Python 3); checked with
  `RunPythonScript` in each.
- **Anywhere else** (Rhino 7 GHPython, a folder of your choice): add the
  repository folder to the search path at the top of a script:

  ```python
  import sys
  sys.path.append(r"C:\path\to\solidwright-repository")
  ```

Rhino keeps imported modules loaded for the whole session, so restart Rhino
after updating the folder.

## Using a component

Each file in `components/` is one Grasshopper component.

1. Add a script component: in Rhino 8 the **Script** component (Python 3 or
   IronPython 2), in Rhino 7 the **GhPython Script** component.
2. Open its editor and paste the whole file in.
3. Zoom in on the component and set its inputs and outputs to match the
   list at the top of the file: same names, and **List Access** where it says
   "list access". Type hints help (Mesh, Brep, Curve...).

Inputs marked optional can be left unconnected or removed from the
component.

## Using the library from Python

The components are thin wrappers, so the same tools work in any Rhino
Python script. Every tool returns its result with a report:

```python
from solidwright.core import difference_plan
from solidwright.rhino.difference import BrepOps

res = difference_plan.cut(plate, cutters, BrepOps(0.001), 0.001)
print(res.report)            # status and every reason, one per line
if res.report.ok:            # good or warning
    solids = res.solids
else:
    print("cutters not applied:", sorted(res.failed))
```

The geometry checks and planning live in `solidwright.core`, which never
imports Rhino; `solidwright.rhino` holds the RhinoCommon calls.

## Supported versions

| Rhino | Python | Status |
|---|---|---|
| Rhino 8 | IronPython 2.7.12 | tested |
| Rhino 8 | Python 3.9.10 | tested |
| Rhino 7 | IronPython 2.7 (GhPython) | written for it, not tested yet |

"Tested" means every component script was run inside Rhino 8 on real
geometry, from a Python script that feeds it inputs the way Grasshopper
does.

The library itself (`solidwright.core`) is tested in CI on Python 3.9 and
3.12 and is written to run unchanged in IronPython 2.7.

## Principles

- **Every result is graded.** Tools return a report with a status (0 good,
  1 warning, 2 messy, 3 failed) and readable reasons, instead of raising or
  quietly handing back bad geometry.
- **Pure-Python core.** `solidwright.core` has no Rhino and no numpy imports, so
  it runs in GHPython (IronPython 2.7), Rhino 8 CPython and plain CPython.
  Thin adapters in `solidwright.rhino` convert to and from RhinoCommon.
- **Tested without Rhino.** Tests use procedural geometry and a small mock of
  the RhinoCommon pieces the adapters touch.

## Write-ups

Short notes on the problems that took the most thought, with numbers:

- [Planar mesh to solid, compared with Rhino](docs/mesh-to-solid-vs-rhino.md):
  `MeshToSolid` and Rhino's `ToNURBS` + `MergeAllCoplanarFaces` on the same
  noisy and single precision meshes: faces, accuracy and time.
- [The mirror-plane trap](docs/mirror-plane-trap.md): why a side view of a
  symmetric part looks as level as its true axis, and what decides.
- [Why interpolating between sections cuts into the part](docs/interpolation-cuts-between-kinks.md):
  a section-based solid needs a level wherever any side bends.
- [Split corners and the validated merge](docs/split-corners-and-validated-merge.md):
  why pushing faces out splits shared corners, and how the merge is checked.
- [Verified union](docs/verified-union.md): why "it ran" is not "it worked"
  for a Boolean union, and what Rhino actually returns.
- [Verified difference](docs/verified-difference.md): a closed solid is not
  proof of a complete cut; what Rhino returns, and the checks that catch a
  skipped cutter.

## Layout

```
solidwright/
  core/     pure Python geometry: no Rhino, no numpy (IronPython 2.7 and Python 3)
  rhino/    thin adapters between RhinoCommon and the core
  gh/       input and message helpers for Grasshopper components
components/ one Grasshopper component script per tool
tests/      pytest, with procedural test geometry and a small RhinoCommon mock
docs/       the write-ups
```

## Running the tests

```
python -m pip install pytest numpy scipy
python -m pytest
```

## License

MIT. See [LICENSE](LICENSE).
