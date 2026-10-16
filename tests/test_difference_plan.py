"""difference_plan with stand-in solids: sets of unit cubes, so difference, volume and inclusion are exact.

Rhino's failures are played by FakeOps hooks: a difference that refuses
(always, at one tolerance, unless the job is scaled up, unless the cutter
is grown), one that skips a cutter quietly, one that returns an open piece.
"""

from solidwright.core import difference_plan as dp


class Cells(object):
    def __init__(self, cells, scale=(1.0, 1.0, 1.0), problem=None, name=None):
        self.cells = frozenset(cells)
        self.scale = scale  # what scaled() did to it; the cells stay the same
        self.problem = problem
        self.name = name


def block(x0, y0, z0, x1, y1, z1, name=None):
    return Cells([(x, y, z) for x in range(x0, x1) for y in range(y0, y1) for z in range(z0, z1)], name=name)


def _neighbours(c):
    x, y, z = c
    return [(x + 1, y, z), (x - 1, y, z), (x, y + 1, z), (x, y - 1, z), (x, y, z + 1), (x, y, z - 1)]


def _groups(cells):
    cells = set(cells)
    out = []
    while cells:
        todo = [cells.pop()]
        group = set(todo)
        while todo:
            for n in _neighbours(todo.pop()):
                if n in cells:
                    cells.discard(n)
                    group.add(n)
                    todo.append(n)
        out.append(group)
    return sorted(out, key=lambda g: min(g))


class FakeOps(dp.Ops):
    def __init__(self, refuse=None, hook=None, enclosed=False):
        self.refuse = refuse  # refuse(solids, cutters, tolerance) -> True makes the call fail (None)
        self.hook = hook  # hook(solids, cutters, honest result) -> what to return
        self.enclosed = enclosed  # True: nothing at all for a cutter wholly inside, as in Rhino
        self.calls = []  # names of the cutters of each call

    def problem(self, s):
        return s.problem

    def volume(self, s):
        return float(len(s.cells))

    def box(self, s):
        xs, ys, zs = zip(*s.cells)
        return (min(xs), min(ys), min(zs)), (max(xs) + 1, max(ys) + 1, max(zs) + 1)

    def difference(self, solids, cutters, tolerance):
        self.calls.append([c.name for c in cutters])
        if self.refuse and self.refuse(solids, cutters, tolerance):
            return None
        if self.enclosed and any(not self.touch(s, c) and c.cells <= s.cells for s in solids for c in cutters):
            return []
        left = set()
        for s in solids:
            left |= s.cells
        for c in cutters:
            left -= c.cells
        out = [Cells(g, scale=solids[0].scale) for g in _groups(left)]
        if self.hook:
            return self.hook(solids, cutters, out)
        return out

    def touch(self, a, b):
        """The surfaces meet: they share cells or faces, and neither lies
        deep inside the other (no intersection curves then, as in Rhino)."""
        def near(s):
            return set(n for c in s.cells for n in _neighbours(c) + [c])
        if not near(a) & b.cells:
            return False
        return not (near(a) <= b.cells or near(b) <= a.cells)

    def inside(self, s, p):
        return (int(p[0] // 1), int(p[1] // 1), int(p[2] // 1)) in s.cells

    def samples(self, s):
        found = []
        for c in sorted(s.cells):
            for n in _neighbours(c):
                if n not in s.cells:
                    d = tuple(n[k] - c[k] for k in range(3))
                    found.append((tuple(c[k] + 0.5 + 0.5 * d[k] for k in range(3)), tuple(float(x) for x in d)))
        return found

    def hollowed(self, s, cutter):
        return Cells(s.cells - cutter.cells, s.scale)

    def scaled(self, s, centre, factors):
        return Cells(s.cells, tuple(a * b for a, b in zip(s.scale, factors)), s.problem, s.name)


def _is_one(scale):
    return all(abs(f - 1.0) < 1e-9 for f in scale)


def _volume(res):
    return sum(len(s.cells) for s in res.solids)


def test_clean_cut_in_one_call():
    ops = FakeOps()
    res = dp.cut(block(0, 0, 0, 6, 6, 2), [block(1, 1, 0, 2, 2, 2, "a"), block(4, 4, 0, 5, 5, 2, "b")], ops)
    assert res.status == 0, str(res.report)
    assert _volume(res) == 72 - 4 and len(res.solids) == 1
    assert sorted(res.applied) == [0, 1] and res.failed == {}
    assert ops.calls == [["a", "b"]]
    assert "2 of 2 cutters cut the solid" in str(res.report)
    assert "1 Boolean difference call" in str(res.report)


def test_cut_apart_gives_several_solids():
    res = dp.cut(block(0, 0, 0, 5, 1, 1), [block(2, 0, 0, 3, 1, 1)], FakeOps())
    assert res.status == 0 and len(res.solids) == 2 and _volume(res) == 4


def test_cutters_that_miss():
    # An L: the notch cutter is inside the L's box but touches nothing; the far one is outside the box.
    solid = Cells(set(block(0, 0, 0, 4, 1, 1).cells) | set(block(0, 0, 0, 1, 4, 1).cells))
    notch, far = block(2, 2, 0, 3, 3, 1, "notch"), block(9, 9, 9, 10, 10, 10, "far")
    ops = FakeOps()
    res = dp.cut(solid, [notch, far], ops)
    assert res.status == 0 and sorted(res.missed) == [0, 1] and res.applied == []
    assert res.solids == [solid] and ops.calls == []
    assert "do not reach the solid: cutters 0, 1" in str(res.report)


def test_cavity_from_a_cutter_inside():
    res = dp.cut(block(0, 0, 0, 3, 3, 3), [block(1, 1, 1, 2, 2, 2)], FakeOps())
    assert res.status == 0 and _volume(res) == 26 and res.applied == [0]


def test_cavity_when_the_difference_gives_nothing():
    # Rhino returns nothing for a cutter wholly inside; the hollow rescue makes the cavity.
    cutters = [block(1, 1, 1, 2, 2, 2), block(3, 0, 0, 4, 4, 1)]
    res = dp.cut(block(0, 0, 0, 4, 4, 4), cutters, FakeOps(enclosed=True))
    assert res.status == 1 and res.rescued == {0: "as a sealed hollow"}, str(res.report)
    assert _volume(res) == 64 - 1 - 4 and sorted(res.applied) == [0, 1]
    assert "cutter 0 lies wholly inside the solid: it leaves a sealed hollow" in str(res.report)


def test_bisection_isolates_the_refusing_cutter():
    cutters = [block(k, 0, 0, k + 1, 1, 2, "c%d" % k) for k in range(0, 8, 2)]
    ops = FakeOps(refuse=lambda s, c, t: any(x.name == "c4" for x in c))
    res = dp.cut(block(0, 0, 0, 8, 3, 2), cutters, ops)
    assert res.status == 2, str(res.report)
    assert sorted(res.failed) == [2] and sorted(res.applied) == [0, 1, 3]
    assert _volume(res) == 48 - 6  # cutter 2's cut is missing
    # one pass at 3 tolerances, then the halves: [c0, c2] in bulk, [c4, c6] refuses, then each alone
    assert ops.calls[:5] == [["c0", "c2", "c4", "c6"]] * 3 + [["c0", "c2"], ["c4", "c6"]]
    assert "1 cutter could not be applied, the result is missing its cut: 2" in str(res.report)
    assert "cutter 2, at the tolerance: the difference failed" in str(res.report)
    assert "cutter 2, grown per axis: the difference failed" in str(res.report)


def test_rescue_at_twice_the_tolerance():
    ops = FakeOps(refuse=lambda s, c, t: abs(t - 0.001) < 1e-12)
    res = dp.cut(block(0, 0, 0, 4, 4, 1), [block(1, 1, 0, 2, 2, 1)], ops)
    assert res.status == 0 and res.rescued == {0: "twice the tolerance"}
    assert "cutter 0 cut at twice the tolerance" in str(res.report)


def test_rescue_scaled_up():
    ops = FakeOps(refuse=lambda s, c, t: not all(abs(x.scale[0] - 10.0) < 1e-9 for x in list(s) + list(c)))
    res = dp.cut(block(0, 0, 0, 4, 4, 1), [block(1, 1, 0, 2, 2, 1)], ops)
    assert res.status == 0 and res.rescued == {0: "scaled up 10 times"}
    assert all(_is_one(s.scale) for s in res.solids)  # scaled back
    assert _volume(res) == 15


def test_rescue_by_growing_the_cutter_warns_of_the_over_cut():
    ops = FakeOps(refuse=lambda s, c, t: not 1.0 < c[0].scale[0] < 2.0)  # grown, not scaled up
    res = dp.cut(block(0, 0, 0, 4, 4, 1), [block(1, 1, 0, 3, 3, 1)], ops)
    assert res.status == 1 and res.rescued == {0: "grown per axis"}
    assert "cutter 0 cut only when grown per axis, over-cuts up to 0.00866" in str(res.report)


def test_a_quietly_skipped_cutter_is_caught():
    def skip_b(solids, cutters, out):
        if any(c.name == "b" for c in cutters):
            kept = set()
            for s in solids:
                kept |= s.cells
            for c in cutters:
                if c.name != "b":
                    kept -= c.cells
            return [Cells(g) for g in _groups(kept)]
        return out
    res = dp.cut(block(0, 0, 0, 6, 6, 1), [block(1, 1, 0, 2, 2, 1, "a"), block(4, 4, 0, 5, 5, 1, "b")],
                 FakeOps(hook=skip_b))
    assert res.status == 2 and sorted(res.failed) == [1] and res.applied == [0]
    assert _volume(res) == 35
    assert "cutter 1 was skipped: its material is still inside the result" in str(res.report)


def test_open_or_bloated_results_are_refused():
    def open_piece(solids, cutters, out):
        return [Cells(out[0].cells, problem="open")]
    res = dp.cut(block(0, 0, 0, 4, 4, 1), [block(1, 1, 0, 2, 2, 1)], FakeOps(hook=open_piece))
    assert res.status == 2 and "a piece that is open" in str(res.report)

    def bloated(solids, cutters, out):
        return [Cells(set(out[0].cells) | set(block(0, 0, 1, 4, 4, 2).cells))]
    res = dp.cut(block(0, 0, 0, 4, 4, 1), [block(1, 1, 0, 2, 2, 1)], FakeOps(hook=bloated))
    assert res.status == 2 and "the volume grew" in str(res.report)


def test_result_inside_a_cutter_is_refused():
    # A cut that went only half way down. Its probes would catch it too, so check without them.
    ops = FakeOps()
    solid, cutter = block(0, 0, 0, 4, 4, 2), block(1, 1, 0, 2, 2, 2)
    plan = dp._Plan(ops, 0.001, dp.DifferenceResult(None))
    c = dp._Cutter(0, cutter, ops.volume(cutter), ops.box(cutter))
    half_way = Cells(set(solid.cells) - {(1, 1, 1)})
    assert plan.accept([half_way], [solid], [c]) == "part of the result lies inside cutter 0"
    assert plan.accept(ops.difference([solid], [cutter], 0.001), [solid], [c]) is None


def test_a_cutter_that_swallows_the_solid():
    ops = FakeOps()
    res = dp.cut(block(1, 1, 1, 2, 2, 2), [block(0, 0, 0, 3, 3, 3)], ops)
    assert res.status == 1 and res.solids == [] and ops.calls == []
    assert "cutter 0 swallows the whole solid: nothing is left" in str(res.report)


def test_bad_inputs():
    good = block(0, 0, 0, 4, 4, 1)
    assert dp.cut(None, [], FakeOps()).status == 3
    assert dp.cut(good, [], FakeOps(), tolerance=0.0).status == 3
    res = dp.cut(Cells(good.cells, problem="open"), [block(1, 1, 0, 2, 2, 1)], FakeOps())
    assert res.status == 3 and res.solids == [] and "the solid is open" in str(res.report)
    # an open cutter that reaches the solid fails, one that misses it is just skipped, None is noted
    cutters = [Cells(block(1, 1, 0, 2, 2, 1).cells, problem="open"), Cells({(9, 9, 9)}, problem="open"),
               None, block(3, 3, 0, 4, 4, 1)]
    res = dp.cut(good, cutters, FakeOps())
    assert res.status == 2 and sorted(res.failed) == [0] and res.missed == [1] and res.applied == [3]
    assert "cutter 0, the cutter is open" in str(res.report)
    assert "empty input skipped: 2" in str(res.report)
    assert _volume(res) == 15


def test_no_cutters_leaves_the_solid():
    solid = block(0, 0, 0, 2, 2, 2)
    res = dp.cut(solid, None, FakeOps())
    assert res.status == 0 and res.solids == [solid]
