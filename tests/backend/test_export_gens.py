"""Pure generator-planner contract tests using committed normalized fixtures.

These tests intentionally exercise ``plan_gens`` with values shaped like an
E2 ExportSnapshot.  They do not construct a Seamless Context or inspect any
current project mount.
"""

from __future__ import annotations

from copy import deepcopy
import base64
import json
from pathlib import Path
import re

import pytest

import export_gens


_DATA = Path(__file__).parent / "data" / "export-gens"


def _restore_fixture_values(value):
    if isinstance(value, dict):
        if set(value) == {"$bytes_base64"}:
            return base64.b64decode(value["$bytes_base64"], validate=True)
        return {key: _restore_fixture_values(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_restore_fixture_values(item) for item in value]
    return value


@pytest.fixture
def one_b_seven_f_snapshot():
    """The independent §10 input, never copied from a mounted Preparic run."""
    return _restore_fixture_values(json.loads((_DATA / "1b7f.snapshot.json").read_text()))


def _gens(snapshot):
    return export_gens.plan_gens(snapshot)


def _by_path(gens):
    result = {gen["path"]: gen for gen in gens}
    assert len(result) == len(gens), "generator paths must be unique"
    return result


def _action_index(gen):
    actions = gen["actions"]
    names = [item["name"] for item in actions]
    assert len(names) == len(set(names)), f"duplicate local action directory in {gen['path']}"
    return {item["name"]: item for item in actions}


def _inputs(parameters):
    action = parameters["action"]
    if action in {"grow", "energy-filter", "rmsd"}:
        return [parameters[key] for key in ("input", "restrict_input") if key in parameters]
    if action in {"identity", "mask-common-conformer"}:
        return [parameters[key] for key in ("input1", "input2") if key in parameters]
    if action == "apply-mask":
        return [parameters["input"], parameters["mask_input"]]
    return []


def _assert_recipe_dag(gen):
    """Check local producer order and explicit external pool declarations."""
    actions = gen["actions"]
    by_name = _action_index(gen)
    preceding = set()
    external = {row["pool"] for row in gen["external_inputs"]}
    for action in actions:
        for pool in _inputs(action["parameters"]):
            if pool in by_name:
                assert pool in preceding, f"{pool} must precede {action['name']} in {gen['path']}"
            else:
                assert pool in external, f"{pool} is undeclared in {gen['path']}"
        preceding.add(action["name"])
    assert set(gen["final_pools"]) <= set(by_name), gen["path"]
    for pool in external:
        assert pool not in by_name, f"{pool} cannot be both local and external"


def _summary(gens):
    """Stable golden projection of the §10 catalogue contract."""
    def pool_key(name):
        match = re.match(r"frag([1-9][0-9]*)", name)
        return (int(match.group(1)), name)

    return [
        {
            "path": gen["path"],
            "modes": sorted(gen["modes"]),
            "range": gen["range"],
            "direction": gen["direction"],
            "actions": sorted(action["name"] for action in gen["actions"]),
            "final_pools": sorted(gen["final_pools"], key=pool_key),
            "external_inputs": sorted([
                {
                    "pool": row["pool"],
                    "fragment": row["fragment"],
                    "sequence": row["sequence"],
                    "frame": row["frame"],
                    "producers": sorted(row["producers"]),
                }
                for row in gen["external_inputs"]
            ], key=lambda row: row["pool"]),
        }
        for gen in sorted(gens, key=lambda item: item["path"])
    ]


def _named(gens, *, prefix):
    return [gen for gen in gens if gen["path"].startswith(prefix)]


def _actions_of(gen, action_type):
    return [item for item in gen["actions"] if item["parameters"]["action"] == action_type]


def test_1b7f_catalogue_matches_committed_golden(one_b_seven_f_snapshot):
    gens = _gens(one_b_seven_f_snapshot)
    golden = json.loads((_DATA / "1b7f.catalogue.golden.json").read_text())
    assert _summary(gens) == golden


def test_1b7f_exact_physical_gen_count_and_mode_catalogue(one_b_seven_f_snapshot):
    gens = _gens(one_b_seven_f_snapshot)
    by_path = _by_path(gens)
    assert len(gens) == 49

    # §10 retains both anchor-tree and direction alternatives.  The clipped
    # nuc11 outer fragment produces a single-fragment gen, not gen-frag10-11.
    for tree in ("anchor", "anchor-refe"):
        assert f"{tree}/gen-frag10" in by_path
        assert f"{tree}/fwd/fixed/gen-frag3-6" in by_path
        assert f"{tree}/bwd/fixed/gen-frag3-6" in by_path
        assert f"{tree}/fwd/unfixed/range-3/gen-frag3-5" in by_path
        assert f"{tree}/bwd/unfixed/range-3/gen-frag3-5" in by_path
        assert f"{tree}/fwd/unfixed/range-2/gen-frag5-6" in by_path
        assert f"{tree}/bwd/unfixed/range-2/gen-frag5-6" in by_path
    assert not any("gen-frag10-11" in path for path in by_path)
    assert not any("frag11" in pool or "frag11" in action["name"]
                   for gen in gens for pool in gen["final_pools"] for action in gen["actions"])

    # All four dual-range fragments receive pairwise common-conformer gens.
    assert {
        path for path in by_path
        if path.startswith("common-conformer/")
    } == {
        f"common-conformer/gen-frag{n}-range-2-3" for n in range(5, 9)
    }

    # The six recipes whose rendered graphs agree across modes are hoisted.
    for path in (
        "anchor/gen-frag10", "anchor-refe/gen-frag10", "grow/bwd/gen-frag9",
        "merge/fwd/anchor/gen-frag10-merged", "merge/fwd/anchor-refe/gen-frag10-merged",
        "merge/bwd/gen-frag9-merged",
    ):
        assert path in by_path
        assert set(by_path[path]["modes"]) == {"fixed", "unfixed"}


def test_1b7f_anchor_runs_clipping_and_alternatives_follow_r3(one_b_seven_f_snapshot):
    by_path = _by_path(_gens(one_b_seven_f_snapshot))

    # Fixed mode joins the three adjacent selected nucleotides across ranges.
    fixed = by_path["anchor/fwd/fixed/gen-frag3-6"]
    assert fixed["final_pools"] == ["frag3", "frag4", "frag5", "frag6"]
    assert {"anchor/fwd/fixed/gen-frag3-6", "anchor/bwd/fixed/gen-frag3-6"} <= set(by_path)

    # In unfixed mode the range-3 run is K=2 and the range-2 nuc6 run is K=1.
    r3 = by_path["anchor/fwd/unfixed/range-3/gen-frag3-5"]
    r2 = by_path["anchor/fwd/unfixed/range-2/gen-frag5-6"]
    assert r3["final_pools"] == ["frag3", "frag4", "frag5-range3"]
    assert r2["final_pools"] == ["frag5-range2", "frag6-range2"]
    assert by_path["anchor/gen-frag10"]["final_pools"] == ["frag10"]

    # The fwd and bwd trees deliberately collide in local final pool names.
    for tree in ("anchor", "anchor-refe"):
        fwd = by_path[f"{tree}/fwd/fixed/gen-frag3-6"]
        bwd = by_path[f"{tree}/bwd/fixed/gen-frag3-6"]
        assert fwd["final_pools"] == bwd["final_pools"]
        fwd_actions = _action_index(fwd)
        bwd_actions = _action_index(bwd)
        assert set(fwd_actions) == set(bwd_actions)
        # R0 deliberately reuses inner pool names across these alternatives;
        # the first-pass root and local dependency edges change instead.
        assert fwd_actions["frag5"]["parameters"]["direction"] == "fwd"
        assert fwd_actions["frag5"]["parameters"]["input"] == "frag4-anchor-scored"
        assert fwd_actions["frag5"]["parameters"]["restrict_input"] == "frag5-anchor-scored"
        assert fwd_actions["frag4"]["parameters"]["input"] == "frag5"
        assert fwd_actions["frag4"]["parameters"]["restrict_input"] == "frag4-anchor-scored"
        assert bwd_actions["frag4"]["parameters"]["input"] == "frag5-anchor-scored"
        assert bwd_actions["frag4"]["parameters"]["restrict_input"] == "frag4-anchor-scored"
        assert bwd_actions["frag5"]["parameters"]["input"] == "frag4"
        assert bwd_actions["frag5"]["parameters"]["restrict_input"] == "frag5-anchor-scored"
        fwd_names = [action["name"] for action in fwd["actions"]]
        bwd_names = [action["name"] for action in bwd["actions"]]
        assert fwd_names.index("frag5") < fwd_names.index("frag4")
        assert bwd_names.index("frag4") < bwd_names.index("frag5")


def test_r0_r9_pool_names_are_canonical_and_unambiguous(one_b_seven_f_snapshot):
    gens = _gens(one_b_seven_f_snapshot)
    by_path = _by_path(gens)
    all_finals = [pool for gen in gens for pool in gen["final_pools"]]
    assert all_finals
    for gen in gens:
        for final in gen["final_pools"]:
            if "-merged" in final or final.endswith("-common-conf"):
                continue
            assert re.fullmatch(r"frag[1-9][0-9]*(?:-range[1-9][0-9]*)?", final), final
            assert not any(role in final for role in ("-anchor", "-fwd", "-bwd", "-grown", "-unique")), final

    # A run selects one producer per plain pool.  Duplicate names are only
    # intentional across mutually exclusive mode/range alternatives.
    producers = {}
    for gen in gens:
        for pool in gen["final_pools"]:
            producers.setdefault(pool, []).append(gen["path"])
    for pool, paths in producers.items():
        if len(paths) > 1:
            assert all(any(token in path for token in ("fixed", "unfixed", "anchor", "anchor-refe", "grow", "merge"))
                       for path in paths), pool

    unfixed_pools = [pool for gen in gens if "unfixed" in gen["modes"] for pool in gen["final_pools"]]
    assert any(pool.startswith("frag5-range2") for pool in unfixed_pools)
    assert any(pool.startswith("frag5-range3") for pool in unfixed_pools)
    assert any(pool == "frag9" for pool in unfixed_pools)
    assert any(pool == "frag10" for pool in unfixed_pools)
    assert not any("-range-" in pool for pool in unfixed_pools)
    assert all(path.split("/")[-1].startswith("gen-") for path in by_path)


def test_traversal_alternatives_stop_at_clips_and_matrix_boundaries(one_b_seven_f_snapshot):
    by_path = _by_path(_gens(one_b_seven_f_snapshot))

    # Range 3 has no bwd traversal before the clip and stops at false frag9.
    r3_grows = _named(by_path.values(), prefix="grow/fwd/unfixed/range-3/")
    assert {gen["final_pools"][0] for gen in r3_grows} == {
        "frag6-range3", "frag7-range3", "frag8-range3",
    }
    assert not _named(by_path.values(), prefix="grow/bwd/unfixed/range-3/")

    # The two range-2 anchor runs face one another, so both directional grow
    # alternatives and every joining merge edge are present, including ends.
    for direction in ("fwd", "bwd"):
        ordinary = [gen for gen in by_path.values()
                    if gen["direction"] == direction and "unfixed" in gen["modes"]
                    and gen["range"] in {"range-2", None} and gen["path"].startswith("grow/")]
        assert len(ordinary) == 3
    fwd_merges = [gen for gen in by_path.values()
                  if gen["path"].startswith("merge/fwd/") and gen["range"] in {"range-2", None}
                  and "unfixed" in gen["modes"]]
    bwd_merges = [gen for gen in by_path.values()
                  if gen["path"].startswith("merge/bwd/") and gen["range"] in {"range-2", None}
                  and "unfixed" in gen["modes"]]
    assert {gen["final_pools"][0] for gen in fwd_merges} == {
        "frag7-merged-range2", "frag8-merged-range2", "frag9-merged", "frag10-merged",
    }
    assert {gen["final_pools"][0] for gen in bwd_merges} == {
        "frag6-merged-range2", "frag7-merged-range2", "frag8-merged-range2", "frag9-merged",
    }
    assert not any("range-3" in gen["path"] and gen["path"].startswith("merge/") for gen in by_path.values())


def test_fixed_scoring_orders_strictest_threshold_then_range_index(one_b_seven_f_snapshot):
    by_path = _by_path(_gens(one_b_seven_f_snapshot))
    gen = by_path["grow/fwd/fixed/gen-frag7"]
    filters = _actions_of(gen, "energy-filter")
    assert [action["parameters"]["protein"] for action in filters] == ["range-3", "range-2"]

    # frag6 has equal thresholds, so range-2 precedes range-3 by numeric index.
    gen = by_path["anchor/fwd/fixed/gen-frag3-6"]
    anchors = _actions_of(gen, "anchor")
    order = [(int(action["parameters"]["protein"].split("-")[1]),
              action["parameters"]["resid"]) for action in anchors]
    assert order == sorted(order)
    scored_frag6 = [action for action in _actions_of(gen, "energy-filter")
                    if action["parameters"]["fragment"] == 6]
    assert [action["parameters"]["protein"] for action in scored_frag6] == ["range-2", "range-3"]

    # Unfixed scoring remains in the chain's root frame, even where the
    # fragment itself interacts with both proteins.
    for range_id in ("range-2", "range-3"):
        gen = next(item for item in by_path.values()
                   if item["path"].startswith(f"grow/fwd/unfixed/{range_id}/")
                   and item["final_pools"][0].startswith("frag7"))
        assert [action["parameters"]["protein"] for action in _actions_of(gen, "energy-filter")] == [range_id]


def test_grow_and_merge_recipes_obey_r11_r12_and_score_only_when_needed(one_b_seven_f_snapshot):
    by_path = _by_path(_gens(one_b_seven_f_snapshot))
    grow = by_path["grow/fwd/fixed/gen-frag7"]
    kinds = [action["parameters"]["action"] for action in grow["actions"]]
    assert kinds == ["grow", "identity", "energy-filter", "energy-filter"]
    assert grow["actions"][0]["parameters"]["direction"] == "fwd"
    assert grow["actions"][1]["parameters"] == {
        "action": "identity", "input1": "frag7-fwd-grown",
    }

    # A one-stacking anchor restriction uses grow+identity, then no extra
    # energy filter; it is already a subset of the scored anchor pool.
    anchor = by_path["anchor/fwd/unfixed/range-2/gen-frag5-6"]
    assert any(action["name"] == "frag6-anchor-range2-grown" for action in anchor["actions"])
    assert any(action["name"] == "frag6-range2" and action["parameters"]["action"] == "identity"
               for action in anchor["actions"])

    # Grown-pool merges use restrict_input and have no redundant score step.
    merge = next(gen for gen in by_path.values()
                 if gen["path"].startswith("merge/fwd/") and gen["final_pools"] == ["frag9-merged"])
    grow_action = _actions_of(merge, "grow")[-1]
    assert "restrict_input" in grow_action["parameters"]
    assert not _actions_of(merge, "energy-filter")

    # Each R5 edge reads its original directional pools. Merge outputs are
    # not chained or redirected automatically.
    for candidate in by_path.values():
        if not candidate["path"].startswith("merge/"):
            continue
        assert all("-merged" not in external["pool"]
                   for external in candidate["external_inputs"])


def test_common_conformer_uses_two_mask_sides_and_preserves_range_outputs(one_b_seven_f_snapshot):
    by_path = _by_path(_gens(one_b_seven_f_snapshot))
    gen = by_path["common-conformer/gen-frag5-range-2-3"]
    mask = _actions_of(gen, "mask-common-conformer")
    applied = _actions_of(gen, "apply-mask")
    assert len(mask) == 1
    assert len(applied) == 2
    assert {action["parameters"]["mask"] for action in applied} == {1, 2}
    assert gen["final_pools"] == ["frag5-range2-common-conf", "frag5-range3-common-conf"]
    assert mask[0]["parameters"]["input1"] == "frag5-range2"
    assert mask[0]["parameters"]["input2"] == "frag5-range3"


def test_three_reachable_ranges_emit_pair_qualified_common_conformers(one_b_seven_f_snapshot):
    snapshot = deepcopy(one_b_seven_f_snapshot)
    range4 = {"chain":"D","start":1,"end":500}
    snapshot["range_specs"]["rows"].append(range4)
    snapshot["interaction_matrix"]["ranges"].append("range-4")
    snapshot["interaction_matrix"]["grid"] = [
        [*cells, 5 <= index <= 8]
        for index, cells in enumerate(snapshot["interaction_matrix"]["grid"], start=3)
    ]
    key = json.dumps(["U",6,4,10], separators=(",", ":"))
    snapshot["stacking_selection"].append(key)
    snapshot["stacking_selection"].sort()
    snapshot["stackings"].append({
        "key":key,"rna_resname":"U","rna_resid":6,"range":"range-4","protein_resid":10,
    })
    snapshot["stackings"].sort(key=lambda row:row["key"])
    frag5, frag6 = snapshot["kept_fragments"][2:4]
    snapshot["candidate_anchors"].append({
        "resid":6,"range":"range-4","stackings":[key],
        "locks":[{"fragment":frag5["key"],"nucleotide":2},
                 {"fragment":frag6["key"],"nucleotide":1}],
    })
    snapshot["candidate_anchors"].sort(key=lambda row:(row["resid"],int(row["range"].split("-")[1])))
    snapshot["class_table"].append({"key":key,"rrm":90,"tilt30":False,"disp45":False})
    snapshot["class_table"].sort(key=lambda row:row["key"])
    snapshot["energy_cells"].extend(
        {"fragment":row["key"],"range":"range-4","range_row":range4,
         "energy":-20.0,"threshold":"-8"}
        for row in snapshot["kept_fragments"] if 5 <= row["number"] <= 8
    )
    snapshot["energy_cells"].sort(key=lambda row:(int(row["fragment"].split("-")[1]),
                                                   int(row["range"].split("-")[1])))

    # Keep the already-normalized stage-5 hand-off consistent with the added
    # matrix column; the planner still takes its reachability from the matrix.
    snapshot["handoff"]["ranges"].append({"index":4,**range4})
    snapshot["handoff"]["runs"].append({
        "range":"range-4","first_fragment":5,"last_fragment":8,
        "anchors":[{"fragment":frag5["key"],"nucleotide":2},
                   {"fragment":frag6["key"],"nucleotide":1}],
    })
    for row in snapshot["kept_fragments"]:
        if 5 <= row["number"] <= 8:
            snapshot["handoff"]["interactions"][row["key"]].append("range-4")

    by_path = _by_path(_gens(snapshot))
    expected_pairs = {(2,3),(2,4),(3,4)}
    common_paths = {path for path in by_path if path.startswith("common-conformer/")}
    assert common_paths == {
        f"common-conformer/gen-frag{fragment}-range-{left}-{right}"
        for fragment in range(5,9) for left,right in expected_pairs
    }

    for fragment in range(5,9):
        expected_outputs = set()
        for left,right in expected_pairs:
            path = f"common-conformer/gen-frag{fragment}-range-{left}-{right}"
            gen = by_path[path]
            left_pool, right_pool = f"frag{fragment}-range{left}", f"frag{fragment}-range{right}"
            outputs = {f"frag{fragment}-range{left}-range{right}-common-conf",
                       f"frag{fragment}-range{right}-range{left}-common-conf"}
            expected_outputs.update(outputs)
            assert set(gen["final_pools"]) == outputs
            assert {row["pool"] for row in gen["external_inputs"]} == {left_pool,right_pool}
            actions = gen["actions"]
            mask = next(action for action in actions
                        if action["parameters"]["action"] == "mask-common-conformer")
            applications = [action for action in actions
                            if action["parameters"]["action"] == "apply-mask"]
            assert mask["parameters"]["input1"] == left_pool
            assert mask["parameters"]["input2"] == right_pool
            assert {action["parameters"]["mask"] for action in applications} == {1,2}
            assert {action["name"] for action in applications} == outputs
            by_side = {action["parameters"]["mask"]:action for action in applications}
            assert by_side[1]["name"] == f"frag{fragment}-range{left}-range{right}-common-conf"
            assert by_side[2]["name"] == f"frag{fragment}-range{right}-range{left}-common-conf"
            assert by_side[1]["parameters"]["input"] == left_pool
            assert by_side[2]["parameters"]["input"] == right_pool
        assert expected_outputs == {
            f"frag{fragment}-range{a}-range{b}-common-conf"
            for left,right in expected_pairs for a,b in ((left,right),(right,left))
        }


def test_reference_anchors_use_base_thresholds_and_sentinel_does_not_drop_descendants(one_b_seven_f_snapshot):
    snapshot = one_b_seven_f_snapshot
    by_path = _by_path(_gens(snapshot))
    reference_gen = by_path["anchor-refe/fwd/fixed/gen-frag3-6"]
    ref_actions = _actions_of(reference_gen, "anchor-refe")
    assert ref_actions
    assert all(action["parameters"]["reference"] == "reference.pdb" for action in ref_actions)
    assert all("protein" not in action["parameters"] for action in ref_actions)
    assert any(float(action["parameters"]["ovrmsd"]) == 0.75 for action in ref_actions)
    assert any(float(action["parameters"]["ovrmsd"]) == 0.5 for action in ref_actions)

    # A rounded reference over L_ov is emitted as the sentinel with an
    # explanatory YAML comment. The rest of the run recipe still exists.
    sentinel_snapshot = deepcopy(snapshot)
    base = next(row for row in sentinel_snapshot["threshold_bases"] if row["fragment"] == 3)
    base["nuc2_ovrmsd"] = "1.80"
    base["nuc2_ovrmsd_raw"] = 1.8
    sentinel = _by_path(_gens(sentinel_snapshot))["anchor-refe/fwd/fixed/gen-frag3-6"]
    sentinel_action = next(action for action in _actions_of(sentinel, "anchor-refe")
                           if action["parameters"]["fragment"] == 3
                           and action["parameters"]["nucleotide"] == "second")
    assert type(sentinel_action["parameters"]["ovrmsd"]) is int
    assert sentinel_action["parameters"]["ovrmsd"] == -1
    comment = " ".join(sentinel_action["comments"]).lower()
    assert "1.8" in comment and "1.80" in comment and "1.732" in comment
    assert "disabled" in comment and "positive" in comment
    assert sentinel["final_pools"] == ["frag3", "frag4", "frag5", "frag6"]
    assert any(action["parameters"]["action"] == "grow" for action in sentinel["actions"])


def test_multiple_stackings_fold_raw_anchors_and_change_r12_restriction(one_b_seven_f_snapshot):
    snapshot = deepcopy(one_b_seven_f_snapshot)
    new_keys = [
        json.dumps(["U", 6, 2, protein_resid], separators=(",", ":"))
        for protein_resid in (11, 12)
    ]
    snapshot["stacking_selection"].extend(new_keys)
    snapshot["stacking_selection"].sort()
    snapshot["stackings"].extend({
        "key":key,"rna_resname":"U","rna_resid":6,"range":"range-2",
        "protein_resid":protein_resid,
    } for key, protein_resid in zip(new_keys, (11, 12)))
    snapshot["stackings"].sort(key=lambda row: row["key"])
    anchor = next(row for row in snapshot["candidate_anchors"]
                  if row["resid"] == 6 and row["range"] == "range-2")
    anchor["stackings"].extend(new_keys)
    snapshot["class_table"].extend(
        {"key":key,"rrm":90,"tilt30":False,"disp45":False} for key in new_keys
    )
    snapshot["class_table"].sort(key=lambda row: row["key"])

    # The shared raw pools have one deterministic identity fold per fragment.
    gaussian = _by_path(_gens(snapshot))["anchor/fwd/unfixed/range-2/gen-frag5-6"]
    raw_actions = _actions_of(gaussian, "anchor")
    raw = [action["name"] for action in raw_actions]
    assert {f"frag5-anchor-range2-raw-{index}" for index in (1, 2, 3)} <= set(raw)
    assert {f"frag6-anchor-range2-raw-{index}" for index in (1, 2, 3)} <= set(raw)
    for fragment in (5,6):
        fragment_raw = [action for action in raw_actions if action["parameters"]["fragment"] == fragment]
        assert len(fragment_raw) == 3
        assert [action["parameters"]["resid"] for action in fragment_raw] == [10,11,12]
        raw_by_name = {action["name"]:action for action in fragment_raw}
        for index, protein_resid in enumerate((10,11,12),1):
            action_name = f"frag{fragment}-anchor-range2-raw-{index}"
            assert raw_by_name[action_name]["parameters"]["resid"] == protein_resid
    raw_order = [(int(action["parameters"]["protein"].split("-")[1]),
                  action["parameters"]["resid"]) for action in raw_actions]
    assert raw_order == sorted(raw_order)
    raw_folds = [action for action in _actions_of(gaussian, "identity")
                 if action["name"].endswith(("anchor-range2-ident-1", "anchor-range2-raw"))]
    expected_folds = {
        "frag5-anchor-range2-ident-1", "frag5-anchor-range2-raw",
        "frag6-anchor-range2-ident-1", "frag6-anchor-range2-raw",
    }
    assert {action["name"] for action in raw_folds} == expected_folds
    assert {action["name"] for action in _actions_of(gaussian, "identity")} == expected_folds
    for fragment in (5, 6):
        stem = f"frag{fragment}-anchor-range2-"
        fold_by_name = {action["name"]:action for action in raw_folds if action["name"].startswith(f"frag{fragment}-")}
        assert fold_by_name[stem+"ident-1"]["parameters"]["input1"] == stem+"raw-1"
        assert fold_by_name[stem+"ident-1"]["parameters"]["input2"] == stem+"raw-2"
        assert fold_by_name[stem+"raw"]["parameters"]["input1"] == stem+"ident-1"
        assert fold_by_name[stem+"raw"]["parameters"]["input2"] == stem+"raw-3"
        assert all("input2" in action["parameters"] for action in fold_by_name.values())
    assert not any(action["name"].startswith("frag6-anchor-range2-grown")
                   for action in gaussian["actions"])
    assert any("restrict_input" in action["parameters"]
               and action["parameters"].get("fragment") == 6
               for action in _actions_of(gaussian, "grow"))
    _assert_recipe_dag(gaussian)
    export_gens.validate_recipe(gaussian)

    # Reference actions are deduplicated by fragment/nucleotide candidate,
    # despite three Gaussian stackings for the same nucleotide. This leaves
    # one raw constraint per locked fragment, so R12 uses grow+identity here
    # even though the Gaussian recipe uses restrict_input.
    reference = _by_path(_gens(snapshot))["anchor-refe/fwd/unfixed/range-2/gen-frag5-6"]
    ref_actions = _actions_of(reference, "anchor-refe")
    ref_keys = [(action["parameters"]["fragment"], action["parameters"]["nucleotide"])
                for action in ref_actions]
    assert len(ref_keys) == len(set(ref_keys)) == 2
    ref_grows = _actions_of(reference, "grow")
    assert len(ref_grows) == 2
    assert all("restrict_input" not in action["parameters"] for action in ref_grows)
    assert {action["name"] for action in ref_grows} == {
        "frag5-anchor-range2-grown", "frag6-anchor-range2-grown",
    }
    ref_identity = {action["name"]:action for action in _actions_of(reference, "identity")}
    assert set(ref_identity) == {"frag5-range2", "frag6-range2"}
    assert set(ref_identity["frag6-range2"]["parameters"]) == {"action", "input1", "input2"}
    assert set(ref_identity["frag6-range2"]["parameters"][key] for key in ("input1", "input2")) == {
        "frag6-anchor-range2-scored", "frag6-anchor-range2-grown",
    }


def test_k2_inner_fragment_qualifies_each_multi_stacking_fold(one_b_seven_f_snapshot):
    """Two candidate fold chains on one K=2 inner fragment need distinct names."""
    snapshot = deepcopy(one_b_seven_f_snapshot)
    snapshot["candidate_anchors"] = [
        row for row in snapshot["candidate_anchors"]
        if row["range"] == "range-3" and row["resid"] in (4, 5)
    ]

    new_rows = []
    for resid, protein_resids in ((4, (94, 95)), (5, (136, 137))):
        anchor = next(row for row in snapshot["candidate_anchors"] if row["resid"] == resid)
        rna_resname = "G" if resid == 4 else "U"
        for protein_resid in protein_resids:
            key = json.dumps([rna_resname, resid, 3, protein_resid], separators=(",", ":"))
            anchor["stackings"].append(key)
            new_rows.append({
                "key": key,
                "rna_resname": rna_resname,
                "rna_resid": resid,
                "range": "range-3",
                "protein_resid": protein_resid,
            })
    for anchor in snapshot["candidate_anchors"]:
        anchor["stackings"].sort()

    retained_keys = {key for anchor in snapshot["candidate_anchors"] for key in anchor["stackings"]}
    snapshot["stacking_selection"] = sorted(retained_keys)
    snapshot["stackings"] = sorted(
        [row for row in snapshot["stackings"] if row["key"] in retained_keys] + new_rows,
        key=lambda row: row["key"],
    )
    snapshot["class_table"] = sorted(
        [row for row in snapshot["class_table"] if row["key"] in retained_keys]
        + [{"key": row["key"], "rrm": 90, "tilt30": False, "disp45": False}
           for row in new_rows],
        key=lambda row: row["key"],
    )

    gens = _gens(snapshot)
    gen = next(
        item for item in gens
        if item["path"].startswith("anchor/fwd/")
        and "unfixed" in item["modes"]
        and item["final_pools"] == ["frag3", "frag4", "frag5"]
    )
    actions = _action_index(gen)

    # Both adjacent nucleotides have three raw stackings on inner frag4.
    # Each binary fold chain gets its own candidate-qualified temporary.
    for ordinal in (4, 5):
        stem = f"frag4-anchor-from-nuc{ordinal}"
        assert {f"{stem}-{index}" for index in (1, 2, 3)} <= set(actions)
        fold = actions[f"{stem}-ident-1"]["parameters"]
        assert fold == {
            "action": "identity",
            "input1": f"{stem}-1",
            "input2": f"{stem}-2",
        }
        final_fold = actions[stem]["parameters"]
        assert final_fold == {
            "action": "identity",
            "input1": f"{stem}-ident-1",
            "input2": f"{stem}-3",
        }

    assert actions["frag4-anchor-ident"]["parameters"] == {
        "action": "identity",
        "input1": "frag4-anchor-from-nuc4",
        "input2": "frag4-anchor-from-nuc5",
    }
    assert len([action for action in _actions_of(gen, "anchor")
                if action["parameters"]["fragment"] == 4]) == 6
    _assert_recipe_dag(gen)
    export_gens.validate_recipe(gen)


def _add_stackings_at_nuc11(snapshot):
    anchor = next(row for row in snapshot["candidate_anchors"]
                  if row["resid"] == 11 and row["range"] == "range-2")
    added = []
    for protein_resid in (50, 51):
        key = json.dumps(["U", 11, 2, protein_resid], separators=(",", ":"))
        anchor["stackings"].append(key)
        added.append({
            "key": key,
            "rna_resname": "U",
            "rna_resid": 11,
            "range": "range-2",
            "protein_resid": protein_resid,
        })
    anchor["stackings"].sort()
    snapshot["stacking_selection"] = sorted(snapshot["stacking_selection"]
                                            + [row["key"] for row in added])
    snapshot["stackings"] = sorted(snapshot["stackings"] + added, key=lambda row: row["key"])
    snapshot["class_table"] = sorted(
        snapshot["class_table"]
        + [{"key": row["key"], "rrm": 90, "tilt30": False, "disp45": False}
           for row in added],
        key=lambda row: row["key"],
    )


def _grow_identity_pair(gen, *, restricted, restrict_pool):
    grow = _actions_of(gen, "grow")
    identity = _actions_of(gen, "identity")
    assert len(grow) == len(identity) == 1
    assert grow[0]["name"] == f"{restrict_pool}-merged-grown"
    if restricted:
        assert grow[0]["parameters"]["restrict_input"] == restrict_pool
    else:
        assert "restrict_input" not in grow[0]["parameters"]
    assert identity[0]["name"] == f"{restrict_pool}-merged"
    assert identity[0]["parameters"] == {
        "action": "identity",
        "input1": f"{restrict_pool}-merged-grown",
        "input2": restrict_pool,
    }
    return grow[0]


def test_single_anchor_merge_is_split_by_anchor_tree(one_b_seven_f_snapshot):
    by_path = _by_path(_gens(one_b_seven_f_snapshot))
    gaussian = by_path["merge/fwd/anchor/gen-frag10-merged"]
    reference = by_path["merge/fwd/anchor-refe/gen-frag10-merged"]
    assert "merge/fwd/gen-frag10-merged" not in by_path
    assert len(_actions_of(by_path["anchor/gen-frag10"], "anchor")) == 1

    # A single Gaussian anchor pool is too large for `restrict_input` (R12):
    # unrestricted grow + identity.  The reference pool is small, so that merge
    # keeps `restrict_input`, like a merge restricted by a grown pool.  Both
    # end with an identity with the restricting pool.
    gaussian_grow = _grow_identity_pair(gaussian, restricted=False, restrict_pool="frag10")
    reference_grow = _grow_identity_pair(reference, restricted=True, restrict_pool="frag10")
    assert gaussian_grow["parameters"]["input"] == reference_grow["parameters"]["input"] == "frag9"
    assert gaussian["final_pools"] == reference["final_pools"] == ["frag10-merged"]
    all_frag9_producers = sorted(
        gen["path"] for gen in by_path.values() if "frag9" in gen["final_pools"]
    )
    for gen, matching_tree in ((gaussian, "anchor"), (reference, "anchor-refe")):
        external_by_pool = {row["pool"]: row for row in gen["external_inputs"]}
        assert external_by_pool["frag10"]["producers"] == [f"{matching_tree}/gen-frag10"]
        # The grow-side pool remains an ordinary input with every permitted
        # producer; only the anchor restriction follows the qualified tree.
        assert external_by_pool["frag9"]["producers"] == all_frag9_producers
    _assert_recipe_dag(gaussian)
    _assert_recipe_dag(reference)
    export_gens.validate_recipe(gaussian)
    export_gens.validate_recipe(reference)

    # The bwd boundary merge at the other anchor is split the same way.
    for mode_level in ("fixed", "unfixed"):
        assert f"merge/bwd/{mode_level}/anchor/gen-frag6-merged" in by_path
        assert f"merge/bwd/{mode_level}/anchor-refe/gen-frag6-merged" in by_path

    # A merge ending on a grown pool has no anchor tree level.
    generic = [gen for path, gen in by_path.items()
               if path.startswith("merge/fwd/") and path.endswith("gen-frag9-merged")]
    assert generic
    assert all(not any(f"/{tree}/" in gen["path"] for tree in ("anchor", "anchor-refe"))
               for gen in generic)


def test_multi_stacking_anchor_merge_restricts_in_both_trees(one_b_seven_f_snapshot):
    snapshot = deepcopy(one_b_seven_f_snapshot)
    _add_stackings_at_nuc11(snapshot)
    by_path = _by_path(_gens(snapshot))
    assert len(_actions_of(by_path["anchor/gen-frag10"], "anchor")) == 3
    assert len(_actions_of(by_path["anchor-refe/gen-frag10"], "anchor-refe")) == 1

    # The combined Gaussian pool is cheap to restrict by, and the reference
    # pool is too: both trees give the same graph, so there is no split.
    assert not any(path.startswith("merge/fwd/anchor") for path in by_path)
    merge = by_path["merge/fwd/gen-frag10-merged"]
    _grow_identity_pair(merge, restricted=True, restrict_pool="frag10")
    external = {row["pool"]: row for row in merge["external_inputs"]}
    assert external["frag10"]["producers"] == ["anchor-refe/gen-frag10", "anchor/gen-frag10"]


def test_external_inputs_are_typed_ordered_contracts_and_recipes_are_dags(one_b_seven_f_snapshot):
    snapshot = one_b_seven_f_snapshot
    gens = _gens(snapshot)
    fragment_by_key = {row["key"]: row for row in snapshot["fragments"]}
    catalogue = _by_path(gens)
    for gen in gens:
        _assert_recipe_dag(gen)
        export_gens.validate_recipe(gen)
        for external in gen["external_inputs"]:
            assert external["fragment"] in fragment_by_key
            assert external["sequence"] == fragment_by_key[external["fragment"]]["sequence"]
            assert isinstance(external["frame"], str) and external["frame"]
            assert external["producers"], f"unproducible external pool {external['pool']} in {gen['path']}"
            for producer in external["producers"]:
                assert producer in catalogue
                assert external["pool"] in catalogue[producer]["final_pools"]


def test_recipe_validator_rejects_duplicate_and_unresolved_local_dependencies(one_b_seven_f_snapshot):
    gens = _gens(one_b_seven_f_snapshot)
    gen = next(item for item in gens if _actions_of(item, "identity"))

    duplicate = deepcopy(gen)
    duplicate["actions"][1]["name"] = duplicate["actions"][0]["name"]
    with pytest.raises(ValueError):
        export_gens.validate_recipe(duplicate)

    unresolved = deepcopy(gen)
    identity = next(action for action in unresolved["actions"]
                    if action["parameters"]["action"] == "identity")
    identity["parameters"]["input1"] = "frag999-unknown"
    with pytest.raises(ValueError):
        export_gens.validate_recipe(unresolved)


def test_action_recipes_use_only_supported_schema_keys_and_explicit_resolution(one_b_seven_f_snapshot):
    allowed = {
        "anchor": {"action", "fragment", "sequence", "exclude", "protein", "resid", "nucleotide",
                   "precomputed_filter_name", "precomputed_filter_threshold"},
        "anchor-refe": {"action", "fragment", "sequence", "exclude", "nucleotide", "ovrmsd", "reference"},
        "grow": {"action", "input", "restrict_input", "fragment", "sequence", "exclude", "direction",
                 "crmsd", "ovrmsd"},
        "energy-filter": {"action", "input", "fragment", "sequence", "exclude", "protein", "threshold",
                          "epsilon", "cascade", "min_cascade"},
        "identity": {"action", "input1", "input2"},
        "mask-common-conformer": {"action", "input1", "input2"},
        "apply-mask": {"action", "input", "mask_input", "mask"},
    }
    required = {
        "anchor": allowed["anchor"] - {"action"},
        "anchor-refe": allowed["anchor-refe"] - {"action", "reference"},
        "grow": allowed["grow"] - {"action", "restrict_input"},
        "energy-filter": allowed["energy-filter"] - {"action"},
        "identity": {"input1"},
        "mask-common-conformer": {"input1", "input2"},
        "apply-mask": {"input", "mask_input", "mask"},
    }
    for gen in _gens(one_b_seven_f_snapshot):
        for record in gen["actions"]:
            parameters = record["parameters"]
            action = parameters["action"]
            assert action in allowed
            assert set(parameters) <= allowed[action]
            assert required[action] <= set(parameters)
            if action in {"anchor", "anchor-refe", "grow", "energy-filter"}:
                assert type(parameters["fragment"]) is int and parameters["fragment"] > 0
                assert parameters["sequence"] == "auto"
                assert parameters["exclude"] == "auto"
            if action == "grow":
                assert parameters["direction"] in {"fwd", "bwd"}
            if action == "energy-filter":
                assert parameters["threshold"] == "auto"
                assert parameters["cascade"] == "auto"
                assert parameters["min_cascade"] == "auto"
                assert parameters["epsilon"] == 15
            if action == "apply-mask":
                assert parameters["mask"] in {1, 2}
            if action == "anchor":
                assert parameters["precomputed_filter_name"] == "gaussian-rrm"
                assert parameters["precomputed_filter_threshold"] == "rrm90"
            if action == "anchor-refe":
                assert parameters["reference"] == "reference.pdb"
            assert isinstance(record["comments"], list)


def test_single_range_recipes_are_common_to_fixed_and_unfixed_without_mode_path(one_b_seven_f_snapshot):
    snapshot = deepcopy(one_b_seven_f_snapshot)
    snapshot["range_specs"]["rows"] = snapshot["range_specs"]["rows"][:2]
    snapshot["interaction_matrix"]["ranges"] = ["range-2"]
    snapshot["interaction_matrix"]["grid"] = [[True] for _ in snapshot["kept_fragments"]]
    snapshot["energy_cells"] = [
        row for row in snapshot["energy_cells"] if row["range"] == "range-2"
    ]
    range2_row = snapshot["range_specs"]["rows"][1]
    snapshot["energy_cells"].extend(
        {"fragment":row["key"],"range":"range-2","range_row":range2_row,
         "energy":-20.0,"threshold":"-8"}
        for row in snapshot["kept_fragments"][:2]
    )
    snapshot["energy_cells"].sort(key=lambda row:(int(row["fragment"].split("-")[1]),
                                                   int(row["range"].split("-")[1])))
    anchor = deepcopy(snapshot["candidate_anchors"][0])
    anchor["range"] = "range-2"
    new_key = json.dumps(["G", 4, 2, 93], separators=(",", ":"))
    anchor["stackings"] = [new_key]
    snapshot["candidate_anchors"] = [anchor]
    snapshot["stackings"] = [{**snapshot["stackings"][0], "key":new_key, "range":"range-2"}]
    snapshot["stacking_selection"] = [new_key]
    snapshot["class_table"] = [{**snapshot["class_table"][0], "key":new_key}]

    gens = _gens(snapshot)
    assert gens
    assert all(gen["modes"] == ["fixed", "unfixed"] for gen in gens)
    assert all("/fixed/" not in gen["path"] and "/unfixed/" not in gen["path"] for gen in gens)
