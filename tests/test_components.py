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
    assert abs(out["Frame"].ZAxis.X - out["Axis"].X) < 1e-12
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


def test_offset_solid_component():
    mesh = rg.Mesh.from_trimesh(shapes.faceted_part())
    axis = rg.Plane(rg.Point3d(0, 0, 0), rg.Vector3d(0.0, 0.0, 1.0))
    out = run_component("OffsetSolid", M=mesh, Clr=0.05, A=axis, Tol=0.001)
    assert out["Status"] == 0 and out["OK"] is True
    assert out["Solid"].IsSolid
    assert "solid: split corners merged" in out["Info"]


def test_offset_solid_component_defaults_and_failures():
    # Clr, A and Tol missing: zero clearance, document tolerance (0.001 outside Rhino).
    out = run_component("OffsetSolid", M=rg.Mesh.from_trimesh(shapes.box()))
    assert out["Status"] == 0 and out["Solid"].Faces.Count == 6
    out = run_component("OffsetSolid", M=rg.Mesh.from_trimesh(shapes.l_block()), Clr=0.1)
    assert out["Status"] == 2 and "solid: the part is not convex" in out["Info"]
    out = run_component("OffsetSolid", M=None)
    assert out["Status"] == 3 and out["Solid"] is None


def test_mesh_to_solid_component():
    g = shapes.rng(9)
    placed, _ = shapes.randomly_placed(shapes.subdivided(shapes.frame(), 2), g)
    out = run_component("MeshToSolid", M=rg.Mesh.from_trimesh(placed), Tol=0.001)
    assert out["Status"] == 0 and out["OK"] is True
    assert out["Solid"].IsSolid and out["Solid"].Faces.Count == 10
    assert "rebuild: 10 flat faces (2 with holes)" in out["Info"]


def test_mesh_to_solid_component_failures():
    v, f = shapes.box().as_lists()
    out = run_component("MeshToSolid", M=rg.Mesh.from_trimesh(shapes.TriMesh(v, f[:-2])))
    assert out["Status"] == 3 and out["Solid"] is None
    assert run_component("MeshToSolid", M=None)["Status"] == 3


def test_clearance_solid_component():
    g = shapes.rng(11)
    placed, rot = shapes.randomly_placed(shapes.faceted_part(), g)
    mesh = rg.Mesh.from_trimesh(placed)
    up = rot.dot((0.0, 0.0, 1.0))
    out = run_component("ClearanceSolid", M=mesh, Clr=0.03, A=rg.Vector3d(*up), Tol=0.001)
    assert out["Status"] == 0, out["Info"]
    assert out["Solid"].IsSolid
    assert len(out["Sections"]) >= 4 and all(s.IsClosed for s in out["Sections"])
    # Without A the part's natural axis is found first.
    out = run_component("ClearanceSolid", M=mesh, Clr=0.03)
    assert out["Status"] == 0 and "solid: the bottom tapers to a point" in out["Info"]
    assert run_component("ClearanceSolid", M=None)["Status"] == 3


def test_profile_fit_component():
    mesh = rg.Mesh.from_trimesh(shapes.faceted_part())
    out = run_component("ProfileFit", M=mesh, A=rg.Vector3d(0, 0, 1), D=rg.Vector3d(1, 0, 0))
    assert out["Status"] == 0
    assert len(out["Lines"]) == 3 and len(out["Knees"]) == 2
    zs = sorted(k.Z for k in out["Knees"])
    assert abs(zs[0] - 0.0) < 1e-6 and abs(zs[1] - 0.02) < 1e-6
    assert out["Profile"].Count > 8
    assert run_component("ProfileFit", M=None)["Status"] == 3



def _box(x0, y0, z0, x1, y1, z1):
    return rg.Brep.CreateFromBox(rg.BoundingBox(rg.Point3d(x0, y0, z0), rg.Point3d(x1, y1, z1)))


def test_verified_union_component():
    parts = [_box(0, 0, 0, 1, 1, 1), _box(0, 0, 0, 4, 4, 1), _box(3, 3, 1, 4, 4, 2), _box(8, 8, 8, 9, 9, 9)]
    out = run_component("VerifiedUnion", B=parts, Tol=0.001)
    assert out["Status"] == 1 and out["OK"] is True  # one input touches nothing
    assert len(out["Solids"]) == 2 and out["Failed"] == []
    assert "3 of 4 inputs united into 1 solid" in out["Info"]
    assert "solid 0: inputs 0, 1, 2" in out["Info"]


def test_verified_union_component_failures():
    out = run_component("VerifiedUnion")
    assert out["Status"] == 3 and out["Solids"] == []
    opened = _box(0, 0, 0, 1, 1, 1)
    del opened.Faces[0]
    out = run_component("VerifiedUnion", B=[_box(0, 0, 0, 2, 2, 2), opened])
    assert out["Status"] == 3 and "input 1 is open" in out["Info"]
