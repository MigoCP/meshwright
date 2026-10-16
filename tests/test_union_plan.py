"""union_plan with stand-in solids: sets of unit cubes, so union, volume and touching are exact.

Rhino's failures are played by FakeOps hooks: a union that refuses (until a
neighbour is in, or until the piece is nudged), one that drops an input, one
that returns an open solid.
"""

from solidwright.core import union_plan as up


def block(x0, y0, z0, x1, y1, z1, name=None):
    cells = [(x, y, z) for x in range(x0, x1) for y in range(y0, y1) for z in range(z0, z1)]
    return Cells(cells, [name] if name is not None else [])


class Cells(object):
    def __init__(self, cells, names, moved=False, problem=None):
        self.cells = frozenset(cells)
        self.names = frozenset(names)
        self.moved = moved
        self.problem = problem


def _neighbours(c):
    x, y, z = c
    return [(x + 1, y, z), (x - 1, y, z), (x, y + 1, z), (x, y - 1, z), (x, y, z + 1), (x, y, z - 1)]


def _connected(cells):
    cells = set(cells)
    if not cells:
        return False
    todo = [next(iter(cells))]
    seen = set(todo)
    while todo:
        for n in _neighbours(todo.pop()):
            if n in cells and n not in seen:
                seen.add(n)
                todo.append(n)
    return len(seen) == len(cells)


class FakeOps(up.Ops):
    def __init__(self, refuse=None, union_hook=None, volume_hook=None, common=True):
        self.refuse = refuse  # refuse(a, b) -> True makes the union fail (None)
        self.union_hook = union_hook  # union_hook(a, b, honest result) -> what to return
        self.volume_hook = volume_hook
        self.common = common  # False: intersections of a piece fully inside another give nothing, as in Rhino
        self.unions = 0

    def problem(self, s):
        return s.problem

    def volume(self, s):
        if self.volume_hook:
            return self.volume_hook(s)
        return float(len(s.cells))

    def box(self, s):
        xs, ys, zs = zip(*s.cells)
        return (min(xs), min(ys), min(zs)), (max(xs) + 1, max(ys) + 1, max(zs) + 1)

    def union(self, a, b):
        self.unions += 1
        if self.refuse and self.refuse(a, b):
            return None
        both = a.cells | b.cells
        out = [Cells(both, a.names | b.names)] if _connected(both) else [a, b]
        if self.union_hook:
            return self.union_hook(a, b, out)
        return out

    def intersection_volume(self, a, b):
        if not self.common and (a.cells <= b.cells or b.cells <= a.cells):
            return None
        return float(len(a.cells & b.cells))

    def encloses(self, outer, inner):
        return inner.cells <= outer.cells

    def touch(self, a, b):
        return bool(a.cells & b.cells) or any(n in b.cells for c in a.cells for n in _neighbours(c))

    def moved(self, s, vector):
        return Cells(s.cells, s.names, moved=True)


def names(res):
    return [sorted(isl.shape.names) for isl in res.islands]


def notes(res):
    return " | ".join(t for _, t in res.report.reasons)


# ---------------------------------------------------------------- the plain cases

def test_overlapping_chain_becomes_one_solid_largest_first():
    parts = [block(0, 0, 0, 2, 2, 2, "a"), block(1, 0, 0, 6, 2, 2, "big"), block(5, 0, 0, 7, 1, 1, "c")]
    res = up.union_all(parts, FakeOps())
    assert res.status == 0, str(res.report)
    assert len(res.islands) == 1 and res.islands[0].members == [1, 0, 2]  # the biggest seeds it
    assert len(res.solids[0].cells) == 8 + 20 - 4 + 2 - 1
    assert "3 of 3 inputs united into 1 solid" in notes(res)


def test_flush_faces_touch_and_unite():
    parts = [block(0, 0, 0, 2, 2, 2, "a"), block(2, 0, 0, 3, 2, 2, "b")]
    res = up.union_all(parts, FakeOps())
    assert res.status == 0 and names(res) == [["a", "b"]]


def test_parts_that_touch_nothing_stay_as_given():
    a, b, c = block(0, 0, 0, 2, 2, 2, "a"), block(1, 1, 1, 3, 3, 3, "b"), block(10, 0, 0, 11, 1, 1, "c")
    res = up.union_all([a, b, c], FakeOps())
    assert res.status == 1 and "1 input touches nothing, left as they are: 2" in notes(res)
    assert len(res.islands) == 2 and res.solids[1] is c  # the very same object, not a copy


def test_boxes_overlap_but_solids_do_not():
    # An L and a cube sitting in its notch, one cell away: the boxes overlap, the solids never meet.
    ell = Cells(block(0, 0, 0, 4, 1, 1).cells | block(0, 0, 0, 1, 4, 1).cells, ["L"])
    cube = block(2, 2, 0, 3, 3, 1, "c")
    res = up.union_all([ell, cube], FakeOps())
    assert len(res.islands) == 2 and res.failed == []


def test_enclosed_piece_is_accepted_through_the_common_volume():
    outer, inner = block(0, 0, 0, 4, 4, 4, "outer"), block(1, 1, 1, 2, 2, 2, "inner")
    res = up.union_all([inner, outer], FakeOps())
    assert res.status == 0 and names(res) == [["inner", "outer"]]
    assert len(res.solids[0].cells) == 64


def test_enclosed_piece_without_a_common_part():
    # Rhino's Boolean intersection gives nothing for a piece fully inside another.
    outer, inner = block(0, 0, 0, 4, 4, 4, "outer"), block(1, 1, 1, 2, 2, 2, "inner")
    res = up.union_all([inner, outer], FakeOps(common=False))
    assert res.status == 0 and names(res) == [["inner", "outer"]]
    # The same union dropping a piece that is not inside is still caught.
    ell = Cells(block(0, 0, 0, 4, 1, 1).cells | block(0, 0, 0, 1, 4, 1).cells, ["L"])
    notch = block(1, 1, 0, 3, 2, 1, "notch")
    ops = FakeOps(common=False, union_hook=lambda a, b, out: [Cells(a.cells, a.names)])
    ops.intersection_volume = lambda a, b: None
    assert up.union_all([ell, notch], ops).failed == [1]


# ---------------------------------------------------------------- refusals, nudges, islands

def test_stubborn_piece_unites_on_the_second_pass():
    # c only unites once d is in (as when its neighbours must merge first).
    parts = [block(0, 0, 0, 4, 4, 1, "body"), block(1, 1, 1, 3, 3, 2, "c"), block(1, 1, 1, 2, 2, 3, "d")]
    ops = FakeOps(refuse=lambda a, b: "c" in b.names and "d" not in a.names)
    res = up.union_all(parts, ops)
    assert res.status == 0 and len(res.islands) == 1
    assert res.islands[0].members == [0, 2, 1]
    assert "1 input united on the second pass" in notes(res) and res.nudged == {}


def test_flush_piece_needs_a_nudge():
    parts = [block(0, 0, 0, 4, 4, 1, "body"), block(1, 1, 1, 2, 2, 2, "flush")]
    ops = FakeOps(refuse=lambda a, b: "flush" in b.names and not b.moved)
    res = up.union_all(parts, ops, tolerance=0.001)
    assert res.status == 1 and res.nudged == {1: "Z+"}
    assert "moved 0.01 to unite: input 1 (Z+)" in notes(res)


def test_nudge_tries_every_direction_in_order():
    tried = []

    class Picky(FakeOps):
        def moved(self, s, vector):
            tried.append(vector)
            return None if len(tried) < 3 else FakeOps.moved(self, s, vector)
    ops = Picky(refuse=lambda a, b: "p" in b.names and not b.moved)
    res = up.union_all([block(0, 0, 0, 2, 2, 2, "body"), block(2, 0, 0, 3, 1, 1, "p")], ops, nudge=0.5)
    assert res.nudged == {1: "X+"} and tried == [(0, 0, 0.5), (0, 0, -0.5), (0.5, 0, 0)]


def test_a_piece_that_never_unites_is_left_out():
    parts = [block(0, 0, 0, 4, 4, 1, "body"), block(1, 1, 1, 2, 2, 2, "bad"), block(3, 3, 1, 4, 4, 2, "ok")]
    res = up.union_all(parts, FakeOps(refuse=lambda a, b: "bad" in b.names))
    assert res.status == 2 and res.failed == [1]
    assert names(res) == [["body", "ok"]]
    assert "1 input could not be united even after nudging, left out: 1" in notes(res)


def test_islands_that_grew_into_each_other_are_merged():
    # Sorted by volume: body, far, bridge. far meets nothing yet and starts its own
    # island; bridge joins body and reaches far, so the merge pass unites them.
    body, far, bridge = block(0, 0, 0, 5, 5, 1, "body"), block(8, 0, 0, 11, 3, 1, "far"), \
        block(4, 0, 0, 9, 1, 1, "bridge")
    res = up.union_all([body, far, bridge], FakeOps())
    assert res.status == 0 and len(res.islands) == 1
    assert sorted(res.islands[0].members) == [0, 1, 2]


# ---------------------------------------------------------------- results that lie

def test_a_dropped_input_outside_the_box_is_caught():
    parts = [block(0, 0, 0, 4, 4, 1, "body"), block(3, 3, 1, 4, 4, 5, "tall")]
    ops = FakeOps(union_hook=lambda a, b, out: [Cells(a.cells, a.names)] if "tall" in b.names else out)
    res = up.union_all(parts, ops)
    assert res.failed == [1] and res.status == 2


def test_a_dropped_input_inside_the_box_is_caught_by_volume():
    # An L and a piece in its notch, touching it: dropping the piece keeps the box the same.
    ell = Cells(block(0, 0, 0, 4, 1, 1).cells | block(0, 0, 0, 1, 4, 1).cells, ["L"])
    notch = block(1, 1, 0, 3, 2, 1, "notch")
    ops = FakeOps(union_hook=lambda a, b, out: [Cells(a.cells, a.names)] if "notch" in b.names else out)
    res = up.union_all([ell, notch], ops)
    assert res.failed == [1]
    honest = up.union_all([ell, notch], FakeOps())
    assert honest.status == 0


def test_a_volume_that_does_not_add_up_is_rejected():
    parts = [block(0, 0, 0, 2, 2, 2, "a"), block(1, 0, 0, 3, 2, 2, "b")]
    ops = FakeOps(volume_hook=lambda s: 99.0 if len(s.names) > 1 else float(len(s.cells)))
    res = up.union_all(parts, ops)
    assert res.failed == [1]


def test_an_open_result_is_rejected():
    parts = [block(0, 0, 0, 2, 2, 2, "a"), block(1, 0, 0, 3, 2, 2, "b")]
    ops = FakeOps(union_hook=lambda a, b, out: [Cells(o.cells, o.names, problem="open") for o in out])
    assert up.union_all(parts, ops).failed == [1]


def test_adds_up_rules():
    plan = up._Plan([], [], FakeOps(), 0.001, 0.01, None)
    assert plan.adds_up(None, 10.0, None, 5.0, 12.0)
    assert not plan.adds_up(None, 10.0, None, 5.0, 9.0)  # smaller than an input
    assert not plan.adds_up(None, 10.0, None, 5.0, 15.1)  # more than both
    assert not plan.adds_up(None, 0.0, None, 5.0, 5.0)


# ---------------------------------------------------------------- inputs

def test_bad_inputs_stop_before_anything_is_united():
    ops = FakeOps()
    good = block(0, 0, 0, 2, 2, 2, "a")
    opened = Cells(good.cells, ["o"], problem="open")
    res = up.union_all([good, opened, None], ops)
    assert res.status == 3 and ops.unions == 0
    assert "input 1 is open; input 2 is empty" in str(res.report)
    # A volume that cannot be computed would never pass the volume check: say so up front.
    res = up.union_all([good, block(1, 1, 1, 3, 3, 3, "z")], FakeOps(volume_hook=lambda s: 0.0))
    assert res.status == 3 and "input 0 has no volume" in str(res.report)


def test_no_inputs_one_input_and_bad_tolerance():
    assert up.union_all([], FakeOps()).status == 3
    one = block(0, 0, 0, 1, 1, 1, "a")
    res = up.union_all([one], FakeOps())
    assert res.status == 0 and res.solids == [one]
    assert "zero" in str(up.union_all([one, one], FakeOps(), tolerance=0.0).report)
