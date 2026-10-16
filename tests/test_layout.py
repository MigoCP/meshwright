import pytest

from solidwright.core import layout


def test_lowest_first_and_the_first_stays_put():
    res = layout.stack([(5.0, 6.0), (0.0, 2.0), (1.0, 1.5)], 0.3)
    assert res.status == 0 and res.order == [1, 2, 0]
    assert res.shifts == pytest.approx([2.0 + 0.3 + 0.5 + 0.3 - 5.0, 0.0, 2.3 - 1.0])  # tops: 2, 2.8, 4.1
    assert "3 parts stacked, 4.1 long" in str(res.report)


def test_keep_order_and_negative_gap():
    res = layout.stack([(5.0, 6.0), (0.0, 2.0)], -0.5, keep_order=True)
    assert res.order == [0, 1] and res.shifts == [0.0, 5.5]
    assert res.status == 1 and "overlap by 0.5" in str(res.report)


def test_missing_parts_and_bad_inputs():
    res = layout.stack([None, (0.0, 1.0), None, (0.0, 1.0)], 0.0)
    assert res.shifts == [None, 0.0, None, 1.0]
    assert "skipped 2 empty inputs: 0, 2" in str(res.report)
    assert layout.stack([None], 0.1).status == 3
    assert layout.stack([], 0.1).status == 3
    assert layout.stack([(0.0, 1.0)], None).status == 3
