"""W7 initial executable acceptance vectors; authored without running tests.

Uses the existing durable decision loader, independent of W7 transform APIs.
"""
from copy import deepcopy
import json

import pytest

import backend


def decision(rrm=90):
    row = {"key": "anchor-A", "rrm": rrm}
    if rrm >= 90:
        row.update(tilt30=True, disp45=False)
    return {"schema": 1, "basis": ["anchor-A"], "active": True, "rows": [row]}


@pytest.mark.parametrize("rrm", [50, 60, 70, 80, 90, 95])
@pytest.mark.parametrize("active", [True, False])
def test_w7_class_decision_roundtrips_including_set_aside(tmp_path, rrm, active):
    value = decision(rrm)
    value["active"] = active
    path = tmp_path / "gaussian-classes.json"
    path.write_text(json.dumps(value))
    assert backend.load_marker_payload(path, "gaussian_classes") == (value, [])


def invalid_decisions():
    value = decision()
    yield {**value, "rows": value["rows"] * 2}
    yield {**value, "basis": ["anchor-B"]}
    yield {**value, "basis": ["anchor-A", "anchor-A"]}
    yield {**value, "rows": []}
    for rrm in [True, 89, 100, "90", 90.0]:
        other = deepcopy(value)
        other["rows"][0]["rrm"] = rrm
        yield other
    for flag in ["tilt30", "disp45"]:
        other = deepcopy(value)
        del other["rows"][0][flag]
        yield other
        for bad in [None, 0, 1, "false"]:
            other = deepcopy(value)
            other["rows"][0][flag] = bad
            yield other
    for rrm in [50, 60, 70, 80]:
        other = decision(rrm)
        other["rows"][0]["tilt30"] = False
        yield other


@pytest.mark.parametrize("value", list(invalid_decisions()))
def test_w7_invalid_class_decision_is_diagnosed_without_rewriting_file(tmp_path, value):
    path = tmp_path / "gaussian-classes.json"
    raw = json.dumps(value)
    path.write_text(raw)
    loaded, diagnostics = backend.load_marker_payload(path, "gaussian_classes")
    assert loaded is None
    assert diagnostics
    assert path.read_text() == raw
