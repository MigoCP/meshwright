"""Run the Grasshopper component scripts the way Grasshopper does: inputs as globals."""

import io
import os

import Rhino.Geometry as rg

import shapes

COMPONENTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "components")


def run_component(name, **inputs):
    path = os.path.join(COMPONENTS, name + ".py")
    with io.open(path, encoding="utf-8") as f:
        code = compile(f.read(), path, "exec")
    namespace = {"__name__": "__main__", "ghenv": None}
    namespace.update(inputs)
    exec(code, namespace)
    return namespace


def test_mesh_health_component_clean_mesh():
    out = run_component("MeshHealth", M=rg.Mesh.from_trimesh(shapes.faceted_part()))
    assert out["Status"] == 0 and out["OK"] is True
    assert out["Clean"].Faces.Count == len(shapes.faceted_part().faces)


def test_mesh_health_component_defaults_when_inputs_missing():
    # Tol and Junk are not defined at all (inputs removed from the component).
    g = shapes.rng(0)
    out = run_component("MeshHealth", M=rg.Mesh.from_trimesh(shapes.with_junk(shapes.box(), g, 2)))
    assert out["Status"] == 1
    assert "dropped 2 loose bits" in out["Info"]


def test_mesh_health_component_no_mesh():
    out = run_component("MeshHealth", M=None, Tol=None)
    assert out["Status"] == 3 and out["OK"] is False and out["Clean"] is None


def test_natural_axis_component_detects_and_checks():
    g = shapes.rng(1)
    placed, rot = shapes.randomly_placed(shapes.faceted_part(), g)
    truth = rot.dot((0.0, 0.0, 1.0))
    mesh = rg.Mesh.from_trimesh(placed)
    out = run_component("NaturalAxis", M=mesh)
    assert out["Status"] == 0
    axis = (out["Axis"].X, out["Axis"].Y, out["Axis"].Z)
    assert abs(sum(a * b for a, b in zip(axis, truth)) - 1.0) < 1e-6
    assert out["Frame"].ZAxis.X == out["Axis"].X
    # Checking a wrong axis given as a plane: messy, and the given axis is kept.
    plane = rg.Plane(rg.Point3d(0, 0, 0), rg.Vector3d(1.0, 0.0, 0.0))
    out = run_component("NaturalAxis", M=mesh, A=plane)
    assert out["Status"] == 2 and "off the detected one" in out["Info"]
    assert abs(abs(out["Axis"].X) - 1.0) < 1e-9


def test_natural_axis_component_reports_mesh_problems_too():
    g = shapes.rng(2)
    mesh = rg.Mesh.from_trimesh(shapes.with_junk(shapes.faceted_part(), g, 1))
    out = run_component("NaturalAxis", M=mesh)
    assert out["Status"] == 1 and "mesh: dropped 1 loose bit" in out["Info"]
    assert run_component("NaturalAxis", M=None)["Status"] == 3
