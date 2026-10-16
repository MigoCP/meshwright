import pytest

from solidwright.core.report import FAILED, GOOD, MESSY, WARNING, Report


def test_empty_report_is_good():
    r = Report()
    assert r.status == GOOD
    assert r.ok and not r.failed
    assert r.summary() == "0 good"


def test_status_is_worst_reason():
    r = Report().warn("small hole").messy("two pieces").warn("flipped faces")
    assert r.status == MESSY
    assert not r.ok
    assert r.worst_reasons() == ["two pieces"]


def test_notes_do_not_change_status():
    r = Report().note("dropped 3 crumbs")
    assert r.status == GOOD
    assert r.worst_reasons() == []


def test_failed():
    r = Report("rebuild").warn("noisy").fail("not closed")
    assert r.failed
    assert r.summary() == "rebuild: 3 failed (not closed)"


def test_merge_keeps_reasons_and_prefix():
    inner = Report().messy("no level axis")
    outer = Report().warn("noisy").merge(inner, prefix="axis: ")
    assert outer.status == MESSY
    assert (MESSY, "axis: no level axis") in outer.reasons


def test_rejects_bad_severity():
    with pytest.raises(ValueError):
        Report().add(4, "nope")
    with pytest.raises(ValueError):
        Report().add(-1, "nope")


def test_str_lists_every_reason():
    text = str(Report("check").warn("a").note("b"))
    assert text.splitlines() == ["check: 1 warning (a)", "  [1] a", "  [0] b"]


def test_severity_constants():
    assert (GOOD, WARNING, MESSY, FAILED) == (0, 1, 2, 3)
