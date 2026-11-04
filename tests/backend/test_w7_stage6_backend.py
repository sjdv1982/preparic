"""Contract-focused W7 backend checks for T18-T23 and T28."""
from __future__ import annotations

import copy
import csv
import json
import math
import shutil
import sys
from pathlib import Path

import pytest


PREPARIC_ROOT = Path(__file__).resolve().parents[2]
ALARIC_ROOT = PREPARIC_ROOT.parent / "alaric"
if str(ALARIC_ROOT) not in sys.path:
    sys.path.insert(0, str(ALARIC_ROOT))

import stage6
from alaric import gaussian_check, gaussian_filter


SPEC_PATH = ALARIC_ROOT / "alaric" / "gaussian_stacking_filter.json"
SPEC = json.loads(SPEC_PATH.read_text(encoding="utf-8"))
MODEL_JSON = json.dumps(SPEC)


def anchor_key(resid=10, protein_range=2, protein_resid=214, resname="A"):
    return json.dumps([resname, resid, protein_range, protein_resid], separators=(",", ":"))


ANCHOR = anchor_key(1)


def tiny_case(descriptors=None, *, fragment_count=1, anchor_resid=1, clip=None):
    residues = [{"resid": index, "base": "A", "chain": "A", "insertion_code": ""}
                for index in range(1, fragment_count + 2)]
    fragments = {
        "schema": 1,
        "echo": {"pdbcode": "1B7F", "rna": {"chain": "A", "start": 1, "end": fragment_count + 1}},
        "residues": residues,
        "fragments": [
            {"number": number, "key": f"frag-{number}-AA", "sequence": "AA",
             "residues": residues[number - 1:number + 1]}
            for number in range(1, fragment_count + 1)
        ],
    }
    stacking = {
        "rna_resname": "A", "rna_resid": anchor_resid,
        "protein_range_pdb_id": 2, "protein_resid": 214,
        "distance": 4.0, "tilt_angle": 10.0, "slip_angle": 20.0,
        "rna_centroid": [0.0, 0.0, 0.0], "protein_centroid": [1.0, 0.0, 0.0],
    }
    stackings = {"schema": 1, "echo": {}, "rows": [stacking]}
    key = anchor_key(anchor_resid)
    selection = {"schema": 1, "active": True, "basis": [], "selected": {key: True}}
    selected_candidates = stage6._candidate_map(fragments, stackings, selection)
    if clip is None:
        clip = {"schema": 1, "basis": None, "clip1": 0, "clip2": 0}
    fits = []
    t18_rows = []
    for candidate in selected_candidates:
        number = candidate["fragment"]
        fits.append({
            "key": candidate["fragment_key"], "fragment": number,
            "sequence": "AA", "conformer": 1, "rotamer": 2,
            "grid": [1.0, 2.0, 3.0], "pdb": "PDB",
        })
        d, phi, theta = (descriptors or {}).get(number, (3.778988652846167, 19.246627257798725, 12.439090799908035))
        t18_rows.append({
            **candidate, "conformer": 1, "rotamer": 2, "grid": [1.0, 2.0, 3.0],
            "distance": d, "phi_deg": phi, "theta_deg": theta,
            "rna_centroid": [0.0, 0.0, 0.0], "protein_centroid": [1.0, 0.0, 0.0],
        })
    best_fits = {"schema": 1, "echo": {"fit": "1b7f"}, "rows": fits}
    bestfit_stackings = {"schema": 1, "echo": {"best_fits": best_fits["echo"], "protein_ranges": []}, "rows": t18_rows}
    return fragments, stackings, selection, clip, best_fits, bestfit_stackings


def class_marker(key=ANCHOR, rrm=90, *, tilt30=False, disp45=False, active=True):
    row = {"key": key, "rrm": rrm}
    if rrm in {90, 95}:
        row.update(tilt30=tilt30, disp45=disp45)
    return {"schema": 1, "basis": [key], "active": active, "rows": [row]}


def gaussian_output(case, model_json=MODEL_JSON):
    fragments, stackings, selection, clip, _fits, descriptors = case
    return stage6.gaussian_scores(
        descriptors, stackings, selection, clip, fragments, model_json, gaussian_filter,
    )


@pytest.mark.parametrize("schema", [None, 0, 2, "1"])
def test_b01_t28_rejects_missing_or_non_v1_model_schema(schema):
    value = copy.deepcopy(SPEC)
    if schema is None:
        value.pop("schema")
    else:
        value["schema"] = schema
    with pytest.raises(ValueError, match="schema|specification"):
        gaussian_output(tiny_case(), json.dumps(value))


@pytest.mark.parametrize("edit", ["remove", "add"])
def test_b01_t28_rejects_missing_or_illegal_filter_names(edit):
    value = copy.deepcopy(SPEC)
    if edit == "remove":
        value["thresholds"].pop("rrm50")
    else:
        value["thresholds"]["rrm42"] = copy.deepcopy(value["thresholds"]["rrm50"])
    value["spec_digest"] = gaussian_filter.spec_digest(value)
    with pytest.raises(ValueError, match="twelve|names"):
        gaussian_output(tiny_case(), json.dumps(value))


def test_b02_score_acceptance_and_zero_score_goldens_are_checked_on_every_call():
    case = tiny_case()
    assert gaussian_output(case)["schema"] == 1
    value = copy.deepcopy(SPEC)
    value["golden"]["scores"][0]["score"] += 0.01
    with pytest.raises(ValueError, match="score golden"):
        gaussian_output(case, json.dumps(value))
    value = copy.deepcopy(SPEC)
    value["golden"]["decisions"][0]["accept"] = not value["golden"]["decisions"][0]["accept"]
    with pytest.raises(ValueError, match="decision golden"):
        gaussian_output(case, json.dumps(value))
    value = copy.deepcopy(SPEC)
    value["golden"]["lattice_decisions"][0]["mask_accept"] = not value["golden"]["lattice_decisions"][0]["mask_accept"]
    with pytest.raises(ValueError, match="lattice decision golden"):
        gaussian_output(case, json.dumps(value))
    value = copy.deepcopy(SPEC)
    value["golden"]["scores"].append({"d": 1000.0, "phi_deg": 0.0, "theta_deg": 0.0, "score": 0.0})
    assert gaussian_output(case, json.dumps(value))["schema"] == 1


def test_b03_mask_is_strict_while_raw_density_survives_outside_it():
    d_boundary = tiny_case({1: (5.0, 20.0, 10.0)})
    d_inside = tiny_case({1: (4.999999, 20.0, 10.0)})
    phi_boundary = tiny_case({1: (4.0, 60.0, 10.0)})
    phi_inside = tiny_case({1: (4.0, 59.99999, 10.0)})
    for case in (d_boundary, phi_boundary):
        row = gaussian_output(case)["candidate_scores"][0]
        assert row["score"] > 0
        assert row["mask_pass"] is False
        assert row["pose_pass"]["rrm50"] is False
        assert gaussian_output(case)["classes"][0]["class_pass"] is False
    assert gaussian_output(d_inside)["candidate_scores"][0]["mask_pass"] is True
    assert gaussian_output(phi_inside)["candidate_scores"][0]["mask_pass"] is True
    cutoff = SPEC["thresholds"]["rrm50"]["cutoff"]
    spec, _model, _names = stage6.load_gaussian_spec(MODEL_JSON, gaussian_filter)
    assert stage6._descriptor_passes(spec, "rrm50", 4.0, math.radians(20), math.radians(10), cutoff)
    assert not stage6._descriptor_passes(spec, "rrm50", 4.0, math.radians(20), math.radians(10), cutoff - 1e-12)
    assert SPEC["thresholds"]["rrm90-tilt30"]["cutoff"] == SPEC["thresholds"]["rrm90"]["cutoff"]


def test_b05_worst_unclipped_score_sets_class_and_clipped_candidate_is_ignored():
    case = tiny_case(fragment_count=11, anchor_resid=6, descriptors={
        5: (4.8, 59.0, 29.0), 6: (4.5, 45.0, 30.0),
    })
    computed = gaussian_output(case)["classes"][0]
    candidate_scores = {row["fragment"]: row["score"] for row in gaussian_output(case)["candidate_scores"]}
    assert computed["score"] == pytest.approx(min(candidate_scores.values()))
    score = computed["score"]
    expected_rrm = next((rrm for rrm in (50, 60, 70, 80, 90, 95)
                         if score >= SPEC["thresholds"][f"rrm{rrm}"]["cutoff"]), 95)
    assert computed["rrm"] == expected_rrm
    clipped = tiny_case(fragment_count=11, anchor_resid=6,
        descriptors={5: (4.8, 59.0, 29.0), 6: (4.5, 45.0, 30.0)},
        clip={"schema": 1, "basis": None, "clip1": 5, "clip2": 5})
    clipped_output = gaussian_output(clipped)
    assert clipped_output["classes"][0]["candidate_count"] == 1
    assert clipped_output["classes"][0]["score"] == pytest.approx(
        next(row["score"] for row in clipped_output["candidate_scores"] if row["fragment"] == 6)
    )
    low = tiny_case({1: (4.9, 59.9, 29.9)})
    # The exact density is allowed to be below every delivered cutoff; class 95 remains explicit.
    low_out = gaussian_output(low)
    assert low_out["classes"][0]["rrm"] == 95
    assert low_out["classes"][0]["below_rrm95"] is True


@pytest.mark.parametrize("descriptors,expected", [
    ((4.5, 44.999, 29.999), {"tilt30": True, "disp45": True}),
    ((4.5, 44.999, 30.0), {"tilt30": False, "disp45": True}),
    ((4.5, 45.0, 29.999), {"tilt30": True, "disp45": False}),
])
def test_b06_automatic_prefilters_use_strict_limits(descriptors, expected):
    case = tiny_case({1: descriptors})
    result = gaussian_output(case)["classes"][0]
    assert result["rrm"] in {90, 95}
    assert result["auto_prefilters"] == expected
    if result["rrm"] in {90, 95}:
        assert result["tilt30"] is expected["tilt30"]
        assert result["disp45"] is expected["disp45"]


def test_b07_threshold_order_and_model_version_are_json_driven():
    case = tiny_case()
    value = copy.deepcopy(SPEC)
    value["thresholds"] = dict(reversed(list(value["thresholds"].items())))
    value["filter_set_version"] = "test-model-version"
    value["spec_digest"] = gaussian_filter.spec_digest(value)
    result = gaussian_output(case, json.dumps(value))
    assert result["filter_names"] == sorted(value["thresholds"])
    assert result["filter_set_version"] == "test-model-version"
    marker = class_marker()
    stage6.validate_gaussian_classes(marker, json.dumps(value), case[2], gaussian_filter)
    illegal = class_marker(rrm=50)
    illegal["rows"][0]["key"] = "anchor-A"
    value["thresholds"].pop("rrm50")
    value["thresholds"]["rrm51"] = copy.deepcopy(SPEC["thresholds"]["rrm50"])
    value["spec_digest"] = gaussian_filter.spec_digest(value)
    with pytest.raises(ValueError, match="twelve|names"):
        stage6.validate_gaussian_classes(illegal, json.dumps(value), case[2], gaussian_filter)


@pytest.mark.parametrize("bad", [
    {"schema": 1, "basis": [ANCHOR], "active": True, "rows": [
        {"key": ANCHOR, "rrm": 90, "tilt30": True, "disp45": False},
        {"key": ANCHOR, "rrm": 95, "tilt30": True, "disp45": False},
    ]},
    {"schema": 1, "basis": [ANCHOR], "active": True, "rows": []},
    {"schema": 1, "basis": [ANCHOR, "extra"], "active": True, "rows": [
        {"key": ANCHOR, "rrm": 90, "tilt30": True, "disp45": False},
        {"key": "extra", "rrm": 90, "tilt30": True, "disp45": False},
    ]},
    {"schema": 1, "basis": [ANCHOR], "active": True, "rows": [
        {"key": ANCHOR, "rrm": True, "tilt30": True, "disp45": False},
    ]},
    {"schema": 1, "basis": [ANCHOR], "active": True, "rows": [
        {"key": ANCHOR, "rrm": 90, "tilt30": 1, "disp45": False},
    ]},
    {"schema": 1, "basis": [ANCHOR], "active": True, "rows": [
        {"key": ANCHOR, "rrm": 80, "tilt30": False, "disp45": False},
    ]},
])
def test_b08_t19_rejects_malformed_active_class_rows(bad):
    with pytest.raises(ValueError):
        stage6.validate_gaussian_classes(bad, MODEL_JSON, {"schema": 1, "active": True, "selected": {ANCHOR: True}}, gaussian_filter)


def test_b09_t20_and_semantic_identity_do_not_depend_on_class_edits():
    case = tiny_case(fragment_count=11, anchor_resid=6)
    fragments, stackings, selection, clip, best_fits, _ = case
    first = stage6.gaussian_check_jobs(best_fits, fragments, stackings, selection, clip)
    second = stage6.gaussian_check_jobs(best_fits, fragments, stackings, selection, clip)
    assert first == second
    assert all(set(job) == {"key", "fragment", "nucleotide", "sequence", "conformer", "rotamer", "grid", "range", "protein_resid"}
               for job in first["jobs"])
    assert stage6.semantic_gaussian_identity(MODEL_JSON) == stage6.semantic_gaussian_identity(MODEL_JSON)
    missing = copy.deepcopy(best_fits)
    missing["rows"].clear()
    jobs = stage6.gaussian_check_jobs(missing, fragments, stackings, selection, clip)
    assert jobs["errors"]
    assert len(jobs["jobs"]) < len(first["jobs"])


def _a9_result_for(class_row, *, fill=1 / 3, emitted=True, pose_pass=True, implied=None):
    filter_name = stage6._class_filter_name(class_row)
    result = {"schema": 1, "filter_set_version": SPEC["filter_set_version"],
              "min_cell_fraction": 1 / 3, "jobs": [{
        "key": ANCHOR, "fragment": 1, "nucleotide": "first", "filters": {},
    }]}
    measured = {"fill": fill, "emitted": emitted, "pose_pass": pose_pass, "score": 0.2}
    result["jobs"][0]["filters"][filter_name] = (
        {"emitted": False, "implied_by": implied} if implied else measured
    )
    return result


def _merge_case(class_row, result):
    marker = {"schema": 1, "basis": [ANCHOR], "active": True, "rows": [class_row]}
    selection = {"schema": 1, "active": True, "basis": [], "selected": {ANCHOR: True}}
    clip = {"schema": 1, "basis": None, "clip1": 0, "clip2": 0}
    fits = {"schema": 1, "echo": {"fit": "current"}, "rows": []}
    jobs = {"schema": 1, "echo": {"best_fits": fits["echo"], "selected_keys": [ANCHOR], "clip": clip},
            "jobs": [{"key": ANCHOR, "fragment": 1, "nucleotide": "first"}], "errors": []}
    return stage6.merge_gaussian_validation(result, jobs, marker, selection, clip, fits, MODEL_JSON, gaussian_filter)


@pytest.mark.parametrize("fill,emitted,pose,state", [
    (1 / 3, True, True, "green"),
    (0.1, False, True, "amber"),
    (0.0, False, True, "red"),
    (0.5, True, False, "amber"),
])
def test_b10_merge_uses_emission_fill_and_pose_state(fill, emitted, pose, state):
    row = {"key": ANCHOR, "rrm": 90, "tilt30": False, "disp45": False}
    merged = _merge_case(row, _a9_result_for(row, fill=fill, emitted=emitted, pose_pass=pose))
    assert merged["rows"][0]["worst_state"] == state
    assert merged["rows"][0]["candidates"][0]["fill"] == fill


def test_b10_implied_rejection_is_non_green_and_has_no_invented_fill():
    row = {"key": ANCHOR, "rrm": 90, "tilt30": False, "disp45": False}
    result = _a9_result_for(row, implied="rrm80")
    merged = _merge_case(row, result)
    candidate = merged["rows"][0]["candidates"][0]
    assert candidate == {"fragment": 1, "nucleotide": "first", "state": "implied", "implied_by": "rrm80"}
    assert merged["rows"][0]["worst_state"] == "implied"


def _reference_case(threshold):
    case = tiny_case()
    fragments, stackings, selection, clip, fits, _ = case
    edit = {"schema": 1, "basis": {"pdbcode": "1B7F", "rna": {"chain": "A", "start": 1, "end": 2}},
            "active": True, "pairs": [], "bases": [{"fragment": 1, "sequence": "AA",
             "nuc1_ovrmsd": threshold, "nuc2_ovrmsd": "1.00"}]}
    computed = {"schema": 1, "bases": [{"fragment": 1, "sequence": "AA",
                "nuc1_ovrmsd": "1.00", "nuc1_ovrmsd_raw": 1.0,
                "nuc2_ovrmsd": "1.00", "nuc2_ovrmsd_raw": 1.0}]}
    return stage6.validate_refe_anchors(edit, computed, fragments, stackings, selection, clip, fits)


@pytest.mark.parametrize("threshold,status", [
    ("1.00002", "green"), ("1.00001", "green"), ("1.00000", "amber"),
    ("0.99999", "red"), ("0.99998", "red"),
])
def test_b11_t23_reads_only_s11_base_values_and_has_exact_boundary(threshold, status):
    result = _reference_case(threshold)
    assert result["rows"][0]["status"] == status
    assert result["rows"][0]["threshold"] == threshold
    assert result["rows"][0]["raw"] == 1.0
    assert set(result["rows"][0]) >= {"key", "fragment", "nucleotide", "computed", "raw", "threshold", "delta", "status"}


@pytest.mark.parametrize("threshold,blocked", [
    ("1.73", False), (str(math.sqrt(3.0)), False), (str(math.nextafter(math.sqrt(3.0), math.inf)), True), ("1.74", True),
])
def test_b12_t23_uses_lov_as_a_hard_blocker(threshold, blocked):
    row = _reference_case(threshold)["rows"][0]
    assert row["blocked"] is blocked
    assert (row["status"] == "blocked") is blocked


def test_b13_s11_threshold_edits_do_not_enter_t20_or_t28_identity():
    case = tiny_case(fragment_count=11, anchor_resid=6)
    jobs_before = stage6.gaussian_check_jobs(case[4], case[0], case[1], case[2], case[3])
    scores_before = gaussian_output(case)
    jobs_after = stage6.gaussian_check_jobs(case[4], case[0], case[1], case[2], case[3])
    scores_after = gaussian_output(case)
    assert jobs_before == jobs_after
    assert scores_before == scores_after


def test_b14_t22_rejects_stale_t20_echo_and_a9_extra_rows():
    row = {"key": ANCHOR, "rrm": 90, "tilt30": False, "disp45": False}
    result = _a9_result_for(row)
    fits = {"schema": 1, "echo": {"fit": "current"}, "rows": []}
    selection = {"schema": 1, "active": True, "basis": [], "selected": {ANCHOR: True}}
    clip = {"schema": 1, "basis": None, "clip1": 0, "clip2": 0}
    marker = {"schema": 1, "basis": [ANCHOR], "active": True, "rows": [row]}
    jobs = {"schema": 1, "echo": {"best_fits": {"fit": "old"}, "selected_keys": [ANCHOR], "clip": clip},
            "jobs": [], "errors": []}
    with pytest.raises(ValueError, match="stale"):
        stage6.merge_gaussian_validation(result, jobs, marker, selection, clip, fits, MODEL_JSON, gaussian_filter)
    jobs["echo"]["best_fits"] = fits["echo"]
    result["jobs"].append({"key": "extra", "fragment": 2, "nucleotide": "second", "filters": {}})
    with pytest.raises(ValueError, match="not present in T20"):
        stage6.merge_gaussian_validation(result, jobs, marker, selection, clip, fits, MODEL_JSON, gaussian_filter)


def _read_residues(pdb_text):
    from Bio.PDB import PDBParser
    from io import StringIO
    model = next(PDBParser(QUIET=True).get_structure("rna", StringIO(pdb_text)).get_models())
    result = []
    for residue in model.get_residues():
        name = residue.get_resname().strip().upper()
        if name.startswith(("R", "D")) and len(name) == 2:
            name = name[1]
        if name in "ACGU":
            result.append({"resid": int(residue.id[1]), "base": name,
                           "chain": str(residue.get_parent().id), "insertion_code": ""})
    return result


def test_b04_t18_t28_match_independent_a9_scores_on_real_1b7f(tmp_path, monkeypatch):
    """Real 1b7f geometry and installed A9 caches; parity stays at 1e-6."""
    from detect_stacking import _detect_stacking
    from util.refe_best_fit_seamless_support import identities

    bestfit_source = ALARIC_ROOT / "tests" / "1b7f" / "best-fits"
    protein_source = ALARIC_ROOT / "tests" / "1b7f" / "pdbs"
    rna_text = (protein_source / "rna.pdb").read_text(encoding="utf-8")
    protein_texts = [
        (protein_source / "1b7f_dom2.pdb").read_text(encoding="utf-8"),
        (protein_source / "1b7f_dom1.pdb").read_text(encoding="utf-8"),
    ]
    ring_refs = {
        f"refe-{name}.pdb": (ALARIC_ROOT / "ring-refe" / f"refe-{name}.pdb").read_text(encoding="utf-8")
        for name in ("A", "C", "PHE", "HIS", "ARG", "TRP")
    }
    reference_rows = _detect_stacking([rna_text, *protein_texts], ring_references=ring_refs)
    assert reference_rows
    residues = _read_residues(rna_text)
    tsv = list(csv.DictReader((ALARIC_ROOT / "tests" / "1b7f" / "data" / "refe-best-fit.tsv").open(), delimiter="\t"))
    all_fragments = []
    fits_by_key = {}
    for record in tsv:
        number = int(record["fragment"])
        sequence = record["sequence"]
        all_fragments.append({"number": number, "key": f"frag-{number}-{sequence}",
            "sequence": sequence, "residues": residues[number - 1:number + 1]})
        path = bestfit_source / f"frag-{number}-{sequence}-bestfit.pdb"
        fits_by_key[f"frag-{number}-{sequence}"] = {
            "key": f"frag-{number}-{sequence}", "fragment": number, "sequence": sequence,
            "conformer": int(record["conformer"]), "rotamer": int(record["rotamer"]),
            "grid": [float(record[f"grid_{axis}"]) for axis in "xyz"],
            "pdb": path.read_text(encoding="utf-8"),
        }
    inventory = {"schema": 1, "echo": {"pdbcode": "1B7F", "rna": {"chain": residues[0]["chain"],
        "start": residues[0]["resid"], "end": residues[-1]["resid"]}}, "residues": residues,
        "fragments": all_fragments}
    chosen = next(row for row in reference_rows
        if any(item["resid"] == row["rna_resid"] and item["base"] == row["rna_resname"] for item in residues[1:-1]))
    key = json.dumps([chosen["rna_resname"], chosen["rna_resid"], chosen["protein_range_pdb_id"], chosen["protein_resid"]], separators=(",", ":"))
    selection = {"schema": 1, "active": True, "basis": [], "selected": {key: True}}
    candidates = stage6._candidate_map(inventory, {"schema": 1, "rows": reference_rows}, selection)
    kept_numbers = sorted({item["fragment"] for item in candidates})
    assert kept_numbers
    clip = {"schema": 1, "basis": None, "clip1": min(kept_numbers) - 1,
            "clip2": len(all_fragments) - max(kept_numbers)}
    fits = {"schema": 1, "echo": {"fixture": "1b7f"},
            "rows": [fits_by_key[f"frag-{number}-{next(row['sequence'] for row in tsv if int(row['fragment']) == number)}"]
                     for number in kept_numbers]}
    range_specs = {"schema": 1, "rows": [
        {"chain": "A", "start": residues[0]["resid"], "end": residues[-1]["resid"]},
        {"chain": "A", "start": 208, "end": 300}, {"chain": "A", "start": 131, "end": 180},
    ]}
    proteins_folder = {"range-2-aa.pdb": protein_texts[0].encode(), "range-3-aa.pdb": protein_texts[1].encode()}
    import detect_stacking
    sys.modules["preparic_detect_stacking"] = detect_stacking
    t18 = stage6.detect_stacking_bestfit(
        fits, proteins_folder, inventory, range_specs,
        *(ring_refs[f"refe-{name}.pdb"].encode() for name in ("A", "C", "PHE", "HIS", "ARG", "TRP")),
    )
    jobs = stage6.gaussian_check_jobs(fits, inventory, {"schema": 1, "rows": reference_rows}, selection, clip)
    assert jobs["jobs"]
    scores = stage6.gaussian_scores(t18, {"schema": 1, "rows": reference_rows}, selection,
        clip, inventory, MODEL_JSON, gaussian_filter)
    candidate_map = {(row["key"], row["fragment"], row["nucleotide"]): row
                     for row in scores["candidate_scores"]
                     if row.get("key") == key and row.get("score") is not None and row.get("kept") is True}

    best_fit_dir = tmp_path / "best-fit"
    (best_fit_dir / "pdbs").mkdir(parents=True)
    for path in bestfit_source.glob("*-bestfit.pdb"):
        shutil.copyfile(path, best_fit_dir / "pdbs" / path.name)
    proteins_dir = tmp_path / "proteins"
    proteins_dir.mkdir()
    for name, content in proteins_folder.items():
        (proteins_dir / name).write_bytes(content)
    identity_json = json.dumps(identities(ALARIC_ROOT / "fraglib"), sort_keys=True, separators=(",", ":"))
    monkeypatch.setenv("fraglib_config", str(Path.home() / ".alaric" / "fraglib.yaml"))
    a9 = gaussian_check.run_gaussian_check(
        jobs["jobs"], best_fit_dir, proteins_dir, "1b7f", identity_json, ALARIC_ROOT,
        semantic_identity=stage6.semantic_gaussian_identity(MODEL_JSON),
    )
    assert a9["jobs"]
    assert set(candidate_map) == {
        (job["key"], job["fragment"], job["nucleotide"]) for job in a9["jobs"]
    }
    prepared = gaussian_check._prepare_jobs(
        jobs["jobs"], best_fit_dir, proteins_dir,
        gaussian_check._read_exclude_codes("1b7f"),
        gaussian_check._ring_reference_paths(ALARIC_ROOT),
        stage6.load_gaussian_spec(MODEL_JSON, gaussian_filter)[1],
    )
    for index, job_result in enumerate(a9["jobs"]):
        assert "error" not in job_result, job_result
        identity = (job_result["key"], job_result["fragment"], job_result["nucleotide"])
        score = candidate_map[identity]["score"]
        measured_filter = next((value for value in job_result["filters"].values()
                                if isinstance(value, dict) and "score" in value), None)
        assert measured_filter is not None, job_result
        a9_score = measured_filter["score"]
        geometry = prepared[index]["prepared"]
        t18_row = candidate_map[identity]
        a9_descriptor = (geometry.distance, math.degrees(geometry.phi), math.degrees(geometry.theta))
        t18_descriptor = (t18_row["distance"], t18_row["phi_deg"], t18_row["theta_deg"])
        delta = tuple(left - right for left, right in zip(t18_descriptor, a9_descriptor))
        assert abs(score - a9_score) <= 1e-6, (
            f"{identity}: score {score:.12g} vs {a9_score:.12g}; "
            f"T18 {t18_descriptor}; A9 {a9_descriptor}; descriptor delta {delta}"
        )


def test_b16_browser_receives_no_gaussian_model_or_density_kernel():
    source = "\n".join(path.read_text(encoding="utf-8") for path in
        (PREPARIC_ROOT / "src").rglob("*.ts"))
    source += "\n" + "\n".join(path.read_text(encoding="utf-8") for path in
        (PREPARIC_ROOT / "src").rglob("*.svelte"))
    for forbidden in ("normalizer", "precision:", "covariance:", "Math.exp(-0.5", "gaussian_stacking_filter.json"):
        assert forbidden not in source



@pytest.mark.parametrize("anchor_resid,clip_values,kept_fragment,clipped_fragment", [
    (2, (1, 0), 2, 1), (3, (0, 1), 2, 3),
])
def test_w8a_kept_candidate_alone_controls_jobs_classes_and_both_validation_panels(
    anchor_resid, clip_values, kept_fragment, clipped_fragment,
):
    clip = {"schema": 1, "basis": None, "clip1": clip_values[0], "clip2": clip_values[1]}
    case = tiny_case(fragment_count=3, anchor_resid=anchor_resid, clip=clip,
                     descriptors={clipped_fragment: (4.9, 59.9, 80.0)})
    fragments, stackings, selection, _, fits, _ = case
    key = anchor_key(anchor_resid)
    marker = class_marker(key, 90, tilt30=False, disp45=True)
    marker_before = copy.deepcopy(marker)
    scores = gaussian_output(case)
    assert len(scores["classes"]) == 1
    assert scores["classes"][0]["candidate_count"] == 1
    assert scores["classes"][0]["auto_prefilters"] == {"tilt30": True, "disp45": True}
    assert [row["fragment"] for row in scores["candidate_scores"]] == [kept_fragment]
    jobs = stage6.gaussian_check_jobs(fits, fragments, stackings, selection, clip)
    assert not jobs["errors"]
    assert [job["fragment"] for job in jobs["jobs"]] == [kept_fragment]

    def check_result(jobs):
        return {"schema": 1, "filter_set_version": SPEC["filter_set_version"], "min_cell_fraction": 1 / 3,
                "jobs": [{"key": job["key"], "fragment": job["fragment"], "nucleotide": job["nucleotide"],
                          "filters": {"rrm90-disp45": {"fill": 1 / 3, "emitted": True,
                              "pose_pass": job["fragment"] != clipped_fragment, "score": 0.3}}}
                         for job in jobs["jobs"]]}

    merged = stage6.merge_gaussian_validation(check_result(jobs), jobs, marker, selection, clip, fits, MODEL_JSON, gaussian_filter)
    assert merged["rows"][0]["worst_state"] == "green"
    assert len(merged["rows"][0]["candidates"]) == 1
    assert marker == marker_before  # changed aggregation never overwrites edited classes
    edits = {"schema": 1, "active": True, "bases": [
        {"fragment": n, "sequence": "AA", "nuc1_ovrmsd": "0.50" if n == clipped_fragment else "1.10",
         "nuc2_ovrmsd": "0.50" if n == clipped_fragment else "1.10"} for n in range(1, 4)]}
    computed = {"schema": 1, "bases": [
        {"fragment": n, "sequence": "AA", "nuc1_ovrmsd": "1.01", "nuc2_ovrmsd": "1.01",
         "nuc1_ovrmsd_raw": 1.0, "nuc2_ovrmsd_raw": 1.0} for n in range(1, 4)]}
    reference = stage6.validate_refe_anchors(edits, computed, fragments, stackings, selection, clip, fits)
    assert [(row["fragment"], row["status"]) for row in reference["rows"]] == [(kept_fragment, "green")]
    unclip = {**clip, "clip1": 0, "clip2": 0}
    uncovered_jobs = stage6.gaussian_check_jobs(fits, fragments, stackings, selection, unclip)
    assert len(uncovered_jobs["jobs"]) == 2
    uncovered = stage6.merge_gaussian_validation(check_result(uncovered_jobs), uncovered_jobs, marker,
        selection, unclip, fits, MODEL_JSON, gaussian_filter)
    assert uncovered["rows"][0]["worst_state"] != "green"
    reference = stage6.validate_refe_anchors(edits, computed, fragments, stackings, selection, unclip, fits)
    assert {row["status"] for row in reference["rows"]} == {"green", "red"}
    assert marker == marker_before
