"""Focused regressions for interpreted W7 runtime dependencies."""

import inspect
import sys
import types
import builtins

import backend


def test_best_fit_threshold_bundle_is_self_contained_without_file_or_sibling_import(monkeypatch):
    modules = backend.threshold_module_code(backend.resolve_alaric_dir())
    refe_source = modules[0]["code"]
    original_import = builtins.__import__

    def reject_sibling_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "threshold_rounding":
            raise AssertionError("the interpreted T8 module must not import a sibling attachment")
        return original_import(name, globals, locals, fromlist, level)

    module = types.ModuleType("refe_thresholds")
    module.__dict__["__builtins__"] = dict(vars(builtins), __import__=reject_sibling_import)
    assert "from threshold_rounding import" not in refe_source
    assert "Path(__file__)" not in refe_source
    exec(compile(refe_source, "<refe_thresholds>", "exec"), module.__dict__)
    monkeypatch.setitem(sys.modules, "refe_thresholds", module)
    monkeypatch.delitem(sys.modules, "threshold_rounding", raising=False)

    namespace = {}
    exec(compile(inspect.getsource(backend.thresholds_from_tables), "<threshold-transform>", "exec"), namespace)

    pair_table = (
        "down\tup\tsequence\tovRMSD\tcRMSD\n"
        "1\t2\tAC\t0.10\t0.20\n"
    ).encode()
    base_table = (
        "fragment\tsequence\tnucleotide1_ovRMSD\tnucleotide2_ovRMSD\n"
        "1\tAC\t0.10\t0.20\n"
    ).encode()
    result = namespace["thresholds_from_tables"](
        {"pairs.tsv": pair_table, "bases.tsv": base_table}, {"fit": "ready"}
    )

    assert result["pairs"][0]["ovrmsd"] == "0.11"
    assert result["bases"][0]["nuc1_ovrmsd"] == "0.11"
    assert result["echo"] == {"fit": "ready"}


def test_gaussian_jobs_json_is_self_contained_as_an_interpreted_transformer():
    namespace = {}
    exec(compile(inspect.getsource(backend._stage6_jobs_json), "<jobs-json-transform>", "exec"), namespace)

    assert namespace["_stage6_jobs_json"]({"jobs": [{"key": "anchor-A"}]}) == '[{"key":"anchor-A"}]'


def test_energy_threshold_transformer_bundle_inlines_rounding_without_sibling_import(monkeypatch):
    import stage5 as preparic_stage5

    energy_source = backend.energy_threshold_module_code(backend.resolve_alaric_dir())["code"]
    original_import = builtins.__import__

    def reject_rounding_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name in {"threshold_rounding", ".threshold_rounding"}:
            raise AssertionError("T16 must not import an unattached sibling module")
        return original_import(name, globals, locals, fromlist, level)

    module = types.ModuleType("energy_threshold")
    module.__dict__["__builtins__"] = dict(vars(builtins), __import__=reject_rounding_import)
    exec(compile(energy_source, "<energy_threshold>", "exec"), module.__dict__)
    monkeypatch.setitem(sys.modules, "energy_threshold", module)
    monkeypatch.setitem(sys.modules, "preparic_stage5", preparic_stage5)
    monkeypatch.delitem(sys.modules, "threshold_rounding", raising=False)

    namespace = {}
    exec(compile(inspect.getsource(backend._stage5_energy_thresholds), "<energy-threshold-transform>", "exec"), namespace)
    best_fits = {"echo": {"fit": "current"}}
    energetics = {
        "schema": 1,
        "echo": {"best_fits": best_fits["echo"]},
        "ranges": [{"range": "range-2", "rows": [{
            "key": "frag-1-GU", "bestfit": {"energy": -1.234, "error": None},
        }]}],
    }

    result = namespace["_stage5_energy_thresholds"](energetics, best_fits)
    # -1.234 + 0.05 = -1.184, so the computed threshold is -1.18.
    assert result["values"] == {"frag-1-GU": {"range-2": "-1.18"}}


def _pinned_energy_threshold_module():
    module = types.ModuleType("energy_threshold")
    source = backend.energy_threshold_module_code(backend.resolve_alaric_dir())["code"]
    exec(compile(source, "<energy_threshold>", "exec"), module.__dict__)
    return module


def test_computed_energy_threshold_lies_005_to_006_above_the_energy():
    from fractions import Fraction
    import random

    import stage5 as preparic_stage5

    round_energy = _pinned_energy_threshold_module().energy_threshold
    compute = preparic_stage5.computed_energy_threshold
    # The 1B7F frag-3-UG / range-2 best fit: the plain A6 value (-4.58) lay only
    # 0.0005 above it, less than the 0.001 A coordinate rounding moves it.
    assert compute(-4.580465468932903, round_energy) == "-4.53"
    assert compute(-4.58, round_energy) == "-4.52"  # exactly on a cent: strictly above
    assert compute(-0.051, round_energy) == "0.00"
    assert compute(-0.05, round_energy) == "0.01"
    assert compute(0.0, round_energy) == "0.06"
    assert compute(1.234, round_energy) == "1.29"

    rng = random.Random(1)
    energies = [rng.uniform(-60, 20) for _ in range(5000)]
    energies += [round(value, 2) for value in energies[:1000]]  # values on cent boundaries
    energies += [rng.uniform(-1e-3, 1e-3) for _ in range(200)] + [1e-300, -1e-300]
    for energy in energies:
        threshold = compute(energy, round_energy)
        assert preparic_stage5._ENERGY_VALUE.fullmatch(threshold), threshold
        # Exact arithmetic: a Decimal context would round 0.05 + 1e-300 to 0.05.
        margin = Fraction(threshold) - Fraction(repr(energy))
        assert Fraction(5, 100) < margin <= Fraction(6, 100), (energy, threshold)
