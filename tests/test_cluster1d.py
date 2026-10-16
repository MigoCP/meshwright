from solidwright.core import cluster1d as c1


def notes(res):
    return " | ".join(t for _, t in res.report.reasons)


def cubes(*sides):
    return [s ** 3 for s in sides]


def test_gaps_over_the_tolerance_start_groups():
    # Sides 1.0, 1.02 (2% apart: same size), 1.5, 1.52, 2.0; given out of order.
    res = c1.group_by_size(cubes(1.5, 1.0, 2.0, 1.02, 1.52), power=1.0 / 3.0)
    assert res.status == 0
    assert res.index == [1, 0, 2, 0, 1] and res.count == 3
    assert abs(res.sizes[2] - 2.0) < 1e-12
    assert res.members() == [[1, 3], [0, 4], [2]]
    assert "3 groups (gaps over 5%)" in notes(res)


def test_count_cuts_at_the_largest_gaps():
    vals = [1.0, 1.1, 1.3, 2.0, 2.05]
    res = c1.group_by_size(vals, count=2)
    assert res.index == [0, 0, 0, 1, 1] and res.status == 1  # every gap but the last is over 5%
    assert "4 distinct sizes found but 2 asked for; some sizes were merged" in notes(res)
    res = c1.group_by_size([1.0, 1.01, 1.02], count=2)
    assert res.count == 2 and "split near-identical inputs" in notes(res)
    res = c1.group_by_size([1.0, 2.0], count=4)
    assert res.count == 2 and "only 2 inputs" in notes(res)


def test_references_on_a_log_scale_with_empty_groups():
    # 1.4 is nearer 2 than 1 on a log scale? log(1.4) = 0.336, log 2 - log 1.4 = 0.357: nearer 1.
    res = c1.group_by_size([1.4, 2.1, 0.9], references=[2.0, 1.0, 8.0])
    assert res.index == [0, 1, 0] and res.count == 3
    assert res.sizes[2] is None and "empty group 2" in notes(res)


def test_bad_values_and_references_are_named_not_dropped_silently():
    res = c1.group_by_size([1.0, None, 0.0, 1.5, float("nan")])
    assert res.index == [0, -1, -1, 1, -1]
    assert "skipped 3 inputs with no size: 1, 2, 4" in notes(res)
    res = c1.group_by_size([1.0, 2.0], references=[1.0, 0.0, 2.0], count=3)
    assert res.count == 2 and "skipped 1 reference with no size" in notes(res)
    assert "the count is ignored" in notes(res)


def test_failures():
    assert c1.group_by_size([]).status == 3
    assert c1.group_by_size([None, 0.0]).status == 3
    assert c1.group_by_size([1.0], tolerance=0.0).status == 3
    assert c1.group_by_size([1.0], references=[0.0]).status == 3
    one = c1.group_by_size([3.0])
    assert one.status == 0 and one.index == [0] and one.sizes == [3.0]
