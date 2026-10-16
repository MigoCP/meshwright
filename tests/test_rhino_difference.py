import Rhino.Geometry as rg

from solidwright.core import difference_plan
from solidwright.rhino.difference import BrepOps


def box(x0, y0, z0, x1, y1, z1):
    return rg.Brep.CreateFromBox(rg.BoundingBox(rg.Point3d(x0, y0, z0), rg.Point3d(x1, y1, z1)))


def test_brep_ops_answers():
    ops = BrepOps(0.001)
    a, b = box(0, 0, 0, 4, 4, 2), box(1, 1, -1, 2, 2, 3)
    out = ops.difference([a], [b], 0.001)
    assert len(out) == 1 and abs(ops.volume(out[0]) - 30.0) < 1e-12
    assert ops.inside(a, (0.5, 0.5, 0.5)) and not ops.inside(a, (5.0, 0.5, 0.5))
    found = ops.samples(a)
    assert len(found) == 6 * 9  # a 3 x 3 grid on each face
    for p, n in found:
        assert ops.inside(a, tuple(p[k] - 0.01 * n[k] for k in range(3)))
        assert not ops.inside(a, tuple(p[k] + 0.01 * n[k] for k in range(3)))
    grown = ops.scaled(b, (1.5, 1.5, 1.0), (2.0, 2.0, 1.0))
    assert ops.box(grown) == ((0.5, 0.5, -1.0), (2.5, 2.5, 3.0)) and ops.box(b)[0] == (1.0, 1.0, -1.0)  # a copy


def test_failing_rhino_calls_become_refusals(monkeypatch):
    def boom(*args):
        raise RuntimeError("Rhino said no")
    monkeypatch.setattr(rg.Brep, "CreateBooleanDifference", staticmethod(boom))
    assert BrepOps(0.001).difference([box(0, 0, 0, 2, 2, 2)], [box(1, 1, 1, 3, 3, 3)], 0.001) is None


def test_cut_with_breps():
    # a hole, a pocket, a sealed cavity and a cutter far away
    solid = box(0, 0, 0, 6, 6, 2)
    cutters = [box(1, 1, -1, 2, 2, 3), box(3, 3, 1, 4, 4, 3), box(1, 4, 0.5, 2, 5, 1.5), box(9, 9, 9, 10, 10, 10)]
    res = difference_plan.cut(solid, cutters, BrepOps(0.001))
    assert res.status == 1, str(res.report)  # the cavity is worth a look
    assert sorted(res.applied) == [0, 1, 2] and res.missed == [3] and res.rescued == {2: "as a sealed hollow"}
    assert len(res.solids) == 1 and res.solids[0].IsSolid
    assert abs(BrepOps(0.001).volume(res.solids[0]) - (72.0 - 2.0 - 1.0 - 1.0)) < 1e-9


def test_cut_caught_a_non_manifold_result():
    # The cavity touches the pocket along an edge: the plain result is not a
    # closed solid, and only a grown cavity clears the edge.
    solid = box(0, 0, 0, 6, 6, 2)
    res = difference_plan.cut(solid, [box(3, 3, 1, 4, 4, 3), box(2, 4, 0.5, 3, 5, 1.5)], BrepOps(0.001))
    assert res.status == 1 and res.rescued == {1: "grown per axis"}, str(res.report)
    assert res.solids[0].IsSolid
