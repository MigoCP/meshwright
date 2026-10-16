import Rhino.Geometry as rg

from meshwright.core import union_plan
from meshwright.rhino.union import BrepOps


def box(x0, y0, z0, x1, y1, z1):
    return rg.Brep.CreateFromBox(rg.BoundingBox(rg.Point3d(x0, y0, z0), rg.Point3d(x1, y1, z1)))


def test_brep_ops_answers():
    ops = BrepOps(0.001)
    a, b = box(0, 0, 0, 2, 2, 2), box(1, 0, 0, 3, 2, 1)
    assert ops.problem(a) is None and abs(ops.volume(a) - 8.0) < 1e-12
    assert ops.box(b) == ((1.0, 0.0, 0.0), (3.0, 2.0, 1.0))
    out = ops.union(a, b)
    assert len(out) == 1 and abs(ops.volume(out[0]) - 10.0) < 1e-12
    assert abs(ops.intersection_volume(a, b) - 2.0) < 1e-12
    assert ops.touch(a, b) and ops.touch(a, box(2, 0, 0, 3, 1, 1))  # sharing a face counts
    assert not ops.touch(a, box(5, 0, 0, 6, 1, 1))
    assert ops.intersection_volume(a, box(5, 0, 0, 6, 1, 1)) is None
    inner = box(0.5, 0.5, 0.5, 1.0, 1.0, 1.0)
    assert ops.encloses(a, inner) and not ops.encloses(inner, a) and not ops.encloses(a, b)
    moved = ops.moved(b, (0.0, 0.0, 0.5))
    assert ops.box(moved)[0] == (1.0, 0.0, 0.5) and ops.box(b)[0] == (1.0, 0.0, 0.0)  # a copy


def test_failing_rhino_calls_become_refusals(monkeypatch):
    def boom(*args):
        raise RuntimeError("Rhino said no")
    monkeypatch.setattr(rg.Brep, "CreateBooleanUnion", staticmethod(boom))
    ops = BrepOps(0.001)
    assert ops.union(box(0, 0, 0, 1, 1, 1), box(0, 0, 0, 1, 1, 1)) is None


def test_union_all_with_breps():
    parts = [box(0, 0, 0, 4, 4, 1), box(1, 1, 1, 2, 2, 2), box(3, 3, 1, 4, 4, 3), box(9, 9, 9, 10, 10, 10)]
    res = union_plan.union_all(parts, BrepOps(0.001))
    assert res.status == 1, str(res.report)
    assert len(res.islands) == 2 and sorted(res.islands[0].members) == [0, 1, 2]
    assert res.solids[0].IsSolid and abs(BrepOps(0.001).volume(res.solids[0]) - 19.0) < 1e-9
    assert res.solids[1] is parts[3]
