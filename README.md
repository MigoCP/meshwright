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

## Running the tests

```
python -m pip install pytest numpy scipy
python -m pytest
```
