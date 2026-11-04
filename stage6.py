"""Backend helpers for Gaussian scoring and reference-anchor validation."""

from __future__ import annotations

import json
import math
import re
from pathlib import PurePosixPath
from typing import Any


_FRAGMENT_KEY = re.compile(r"^frag-([1-9][0-9]*)-([ACGU]{2})$")
_RRM_NAME = re.compile(r"^rrm(50|60|70|80|90|95)(?:-tilt30)?(?:-disp45)?$")
_REL_TOLERANCE = 1e-12
_MIN_CELL_FRACTION = 1.0 / 3.0


def _payload(value: Any) -> Any:
    """Unwrap Seamless byte buffers without changing ordinary values."""
    return getattr(value, "content", value)


def _text(value: Any, label: str) -> str:
    value = _payload(value)
    if isinstance(value, bytes):
        try:
            value = value.decode("utf-8")
        except UnicodeDecodeError as exc:
            raise ValueError(f"{label} is not UTF-8 text.") from exc
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} is empty or not text.")
    return value


def _json_object(value: Any, label: str) -> dict:
    if isinstance(value, dict):
        return value
    try:
        parsed = json.loads(_text(value, label))
    except (TypeError, json.JSONDecodeError) as exc:
        raise ValueError(f"{label} is not valid JSON.") from exc
    if type(parsed) is not dict:
        raise ValueError(f"{label} must be a JSON object.")
    return parsed


def _model(spec: dict, gaussian_filter) -> Any:
    try:
        model = spec["model"]
        return gaussian_filter.Model(
            mean=tuple(float(value) for value in model["mean"]),
            covariance=tuple(tuple(float(value) for value in row) for row in model["covariance"]),
            precision=tuple(tuple(float(value) for value in row) for row in model["precision"]),
            normalizer=float(model["normalizer"]),
        )
    except (KeyError, TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"Gaussian model parameters are invalid: {exc}") from exc


def _legal_filter_names(spec: dict, gaussian_filter) -> list[str]:
    thresholds = spec.get("thresholds")
    if type(thresholds) is not dict:
        raise ValueError("Gaussian specification thresholds must be an object.")
    try:
        expected = list(gaussian_filter.all_names())
    except Exception as exc:
        raise ValueError(f"Cannot read the delivered Gaussian filter names: {exc}") from exc
    names = list(thresholds)
    if len(names) != 12 or set(names) != set(expected):
        raise ValueError("Gaussian specification must contain exactly the twelve delivered filter names.")
    for name, entry in thresholds.items():
        if type(entry) is not dict or entry.get("accept") != ">=" \
            or not isinstance(entry.get("filters"), list):
            raise ValueError(f"Gaussian threshold entry {name!r} is malformed.")
        if any(value not in {"tilt30", "disp45"} for value in entry["filters"]):
            raise ValueError(f"Gaussian threshold entry {name!r} has an unknown prefilter.")
        try:
            cutoff = float(entry["cutoff"])
        except (KeyError, TypeError, ValueError, OverflowError) as exc:
            raise ValueError(f"Gaussian threshold {name!r} has no finite cutoff.") from exc
        if not math.isfinite(cutoff) or cutoff <= 0.0:
            raise ValueError(f"Gaussian threshold {name!r} has no finite positive cutoff.")
    return sorted(names)


def _descriptor_passes(spec: dict, name: str, d: float, phi: float, theta: float, score: float) -> bool:
    entry = spec["thresholds"][name]
    mask = spec["mask"]
    if not d < float(mask["distance_lt_A"]):
        return False
    if not phi < math.radians(float(mask["displacement_lt_deg"])):
        return False
    for flag in entry.get("filters", []):
        limit = float(spec["filters"][flag]["lt_deg"])
        value = theta if flag == "tilt30" else phi
        if not value < math.radians(limit):
            return False
    return score >= float(entry["cutoff"])


def _check_gaussian_goldens(spec: dict, gaussian_filter, model: Any) -> None:
    """Check score and mask/cutoff goldens on every T28 run."""
    golden = spec.get("golden")
    if type(golden) is not dict or type(golden.get("scores")) is not list \
        or type(golden.get("decisions")) is not list \
        or type(golden.get("lattice_decisions")) is not list:
        raise ValueError("Gaussian specification is missing its score and decision goldens.")
    tolerance = golden.get("tolerance_rel")
    if not isinstance(tolerance, (int, float)) or isinstance(tolerance, bool) \
        or float(tolerance) != _REL_TOLERANCE:
        raise ValueError("Gaussian score goldens must use relative tolerance 1e-12.")

    for index, item in enumerate(golden["scores"]):
        try:
            d = float(item["d"])
            phi = math.radians(float(item["phi_deg"]))
            theta = math.radians(float(item["theta_deg"]))
            expected = float(item["score"])
            actual = float(gaussian_filter.score(d, phi, theta, model=model))
        except (KeyError, TypeError, ValueError, OverflowError) as exc:
            raise ValueError(f"Gaussian score golden {index} is malformed: {exc}") from exc
        error = abs(actual - expected)
        allowed = _REL_TOLERANCE * abs(expected) if expected != 0.0 else 0.0
        if not math.isfinite(actual) or error > allowed:
            raise ValueError(
                f"Gaussian score golden {index} failed at relative tolerance 1e-12: "
                f"expected {expected:.17g}, got {actual:.17g}."
            )

    for index, item in enumerate(golden["decisions"]):
        try:
            name = item["name"]
            d = float(item["d"])
            phi = math.radians(float(item["phi_deg"]))
            theta = math.radians(float(item["theta_deg"]))
            expected = item["accept"]
            actual_score = float(gaussian_filter.score(d, phi, theta, model=model))
            actual = _descriptor_passes(spec, name, d, phi, theta, actual_score)
        except (KeyError, TypeError, ValueError, OverflowError) as exc:
            raise ValueError(f"Gaussian decision golden {index} is malformed: {exc}") from exc
        if type(expected) is not bool or actual is not expected:
            raise ValueError(f"Gaussian decision golden {index} failed for {name!r}.")

    mask = spec["mask"]
    for index, item in enumerate(golden["lattice_decisions"]):
        try:
            name = item["name"]
            x, y, z = (float(value) for value in item["xyz"])
            expected = item["mask_accept"]
            distance = math.sqrt(x * x + y * y + z * z)
            phi = math.acos(min(1.0, abs(z) / distance)) if distance else math.inf
            actual = distance < float(mask["distance_lt_A"]) \
                and phi < math.radians(float(mask["displacement_lt_deg"]))
            for flag in spec["thresholds"][name]["filters"]:
                if flag == "disp45":
                    actual = actual and phi < math.radians(float(spec["filters"][flag]["lt_deg"]))
                elif flag == "tilt30":
                    # The lattice golden's rotation is the identity orientation.
                    actual = actual and 0.0 < math.radians(float(spec["filters"][flag]["lt_deg"]))
        except (KeyError, TypeError, ValueError, OverflowError) as exc:
            raise ValueError(f"Gaussian lattice decision golden {index} is malformed: {exc}") from exc
        if type(expected) is not bool or actual is not expected:
            raise ValueError(f"Gaussian lattice decision golden {index} failed for {name!r}.")


def load_gaussian_spec(model_json: Any, gaussian_filter) -> tuple[dict, Any, list[str]]:
    """Validate the authoritative A7 document, delivered names and goldens."""
    spec = _json_object(model_json, "Gaussian model JSON")
    try:
        gaussian_filter.parse_spec(spec)
    except Exception as exc:
        raise ValueError(f"Gaussian specification is invalid: {exc}") from exc
    if type(spec.get("schema")) is not int or spec["schema"] != 1:
        raise ValueError("Gaussian specification schema must be 1.")
    if not isinstance(spec.get("filter_set_version"), str) or not spec["filter_set_version"].strip():
        raise ValueError("Gaussian specification filter_set_version is missing.")
    names = _legal_filter_names(spec, gaussian_filter)
    model = _model(spec, gaussian_filter)
    _check_gaussian_goldens(spec, gaussian_filter, model)
    return spec, model, names


def semantic_gaussian_identity(model_json: Any) -> str:
    """Canonical A9 key material: semantic spec plus decompressed array IDs."""
    spec = _json_object(model_json, "Gaussian model JSON")
    if type(spec.get("schema")) is not int or spec.get("schema") != 1:
        raise ValueError("Gaussian specification schema must be 1.")
    arrays: dict[str, dict[str, str]] = {}
    for name, entry in sorted(spec.get("thresholds", {}).items()):
        files = entry.get("precalc", {}).get("files", {})
        arrays[name] = {
            filename: str(files[filename]["npy_sha256"])
            for filename in ("filtered-concat", "filtered-sizes", "rotations", "translations")
        }
    identity = {"spec_digest": str(spec["spec_digest"]), "arrays": arrays}
    return json.dumps(identity, sort_keys=True, separators=(",", ":"))


def _range_rows(range_specs: Any) -> list[dict]:
    if type(range_specs) is not dict or range_specs.get("schema") != 1 \
        or not isinstance(range_specs.get("rows"), list):
        raise ValueError("Current protein range rows are not ready.")
    rows = range_specs["rows"]
    if not rows or any(type(row) is not dict for row in rows):
        raise ValueError("Current protein range rows are malformed.")
    return rows[1:]


def _bestfit_rows(best_fits: Any) -> list[dict]:
    if type(best_fits) is not dict or best_fits.get("schema") != 1 \
        or not isinstance(best_fits.get("rows"), list) or type(best_fits.get("echo")) is not dict:
        raise ValueError("Fresh best-fit rows are not ready.")
    rows = best_fits["rows"]
    result = []
    seen = set()
    for row in rows:
        if type(row) is not dict or not isinstance(row.get("key"), str):
            raise ValueError("A best-fit row is malformed.")
        match = _FRAGMENT_KEY.fullmatch(row["key"])
        if match is None or row["key"] in seen:
            raise ValueError(f"Invalid or duplicate best-fit key {row.get('key')!r}.")
        if row.get("sequence") != match.group(2) or row.get("fragment") != int(match.group(1)):
            raise ValueError(f"Best-fit identity does not match {row['key']}.")
        seen.add(row["key"])
        result.append(row)
    return result


def _fragment_rows(fragments: Any) -> list[dict]:
    if type(fragments) is not dict or fragments.get("schema") != 1 \
        or not isinstance(fragments.get("fragments"), list) or not isinstance(fragments.get("residues"), list):
        raise ValueError("Canonical fragment map is not ready.")
    rows = fragments["fragments"]
    if not rows:
        raise ValueError("Canonical fragment map is empty.")
    for row in rows:
        if type(row) is not dict or not isinstance(row.get("key"), str) \
            or type(row.get("number")) is not int or not isinstance(row.get("residues"), list):
            raise ValueError("Canonical fragment map has a malformed row.")
    return rows


def _stack_key(row: dict) -> str:
    try:
        parts = [row["rna_resname"], row["rna_resid"], row["protein_range_pdb_id"], row["protein_resid"]]
    except KeyError as exc:
        raise ValueError(f"A stacking row is missing {exc.args[0]!r}.") from exc
    if not isinstance(parts[0], str) or any(type(value) is not int for value in parts[1:]):
        raise ValueError("A stacking row has an invalid residue address.")
    return json.dumps(parts, separators=(",", ":"), ensure_ascii=False)


def _stacking_rows(stackings: Any) -> list[dict]:
    if type(stackings) is not dict or stackings.get("schema") != 1 \
        or not isinstance(stackings.get("rows"), list):
        raise ValueError("Fresh reference stacking rows are not ready.")
    return stackings["rows"]


def _selected_keys(selection: Any) -> list[str]:
    if type(selection) is not dict or selection.get("schema") != 1 \
        or selection.get("active") is not True or type(selection.get("selected")) is not dict:
        raise ValueError("An active keyed stacking selection is required.")
    if any(not isinstance(key, str) or type(value) is not bool for key, value in selection["selected"].items()):
        raise ValueError("Stacking selection rows must be boolean values.")
    return sorted(key for key, value in selection["selected"].items() if value)


def _candidate_map(fragments: Any, stackings: Any, selection: Any) -> list[dict]:
    fragment_rows = _fragment_rows(fragments)
    stacking_rows = _stacking_rows(stackings)
    selected = set(_selected_keys(selection))
    seen_stackings: set[str] = set()
    candidates: dict[tuple[str, int, str], dict] = {}
    for stacking in stacking_rows:
        if type(stacking) is not dict:
            raise ValueError("A reference stacking row is malformed.")
        key = _stack_key(stacking)
        if key in seen_stackings:
            raise ValueError(f"Duplicate reference stacking key {key}.")
        seen_stackings.add(key)
        if key not in selected:
            continue
        resid = stacking["rna_resid"]
        base = str(stacking["rna_resname"]).strip().upper()
        if len(base) == 2 and base[0] in {"R", "D"}:
            base = base[1]
        if base not in {"A", "C", "G", "U"}:
            continue
        for fragment in fragment_rows:
            for nucleotide, identity in enumerate(fragment["residues"], start=1):
                if type(identity) is not dict or identity.get("resid") != resid or identity.get("base") != base:
                    continue
                label = "first" if nucleotide == 1 else "second"
                address = (key, fragment["number"], label)
                candidates[address] = {
                    "key": key,
                    "fragment": fragment["number"],
                    "fragment_key": fragment["key"],
                    "nucleotide": label,
                    "sequence": fragment.get("sequence"),
                    "resid": resid,
                    "range": stacking["protein_range_pdb_id"],
                    "protein_resid": stacking["protein_resid"],
                    "rna_resname": stacking["rna_resname"],
                }
    missing = selected - seen_stackings
    if missing:
        raise ValueError(f"Selected stacking keys are missing from the reference table: {sorted(missing)}.")
    return [candidates[key] for key in sorted(candidates)]


def _clip_bounds(clip: Any, fragment_count: int) -> tuple[int, int]:
    if type(clip) is not dict:
        raise ValueError("The current fit clip is unavailable.")
    first, last = clip.get("clip1"), clip.get("clip2")
    if type(first) is not int or type(last) is not int or first < 0 or last < 0 \
        or first + last >= fragment_count:
        raise ValueError("The current fit clip is malformed.")
    return first, fragment_count - last


def _is_kept(fragment: int, bounds: tuple[int, int]) -> bool:
    return bounds[0] < fragment <= bounds[1]


def _compute_echo(best_fits: Any, range_specs: Any) -> dict:
    return {"best_fits": best_fits["echo"], "protein_ranges": _range_rows(range_specs)}


def _detector_residues(pdb_text: str) -> list[tuple[str, int]]:
    from io import StringIO
    from Bio.PDB import PDBParser

    structure = PDBParser(QUIET=True).get_structure("bestfit", StringIO(pdb_text))
    model = next(structure.get_models(), None)
    if model is None:
        return []
    allowed = {"A", "C", "G", "U", "RA", "RC", "RG", "RU", "DA", "DC", "DG", "DT"}
    rows = []
    for residue in model.get_residues():
        name = residue.get_resname().strip().upper()
        if name in allowed:
            rows.append((name, int(residue.id[1])))
    return rows


def detect_stacking_bestfit(
    best_fits: Any,
    proteins: Any,
    fragments: Any,
    range_specs: Any,
    ring_ref_a: Any,
    ring_ref_c: Any,
    ring_ref_phe: Any,
    ring_ref_his: Any,
    ring_ref_arg: Any,
    ring_ref_trp: Any,
) -> dict | None:
    """T18: run T4's Kabsch ring detector on each fitted fragment/range."""
    import preparic_detect_stacking

    if best_fits is None or proteins is None or fragments is None or range_specs is None:
        return None
    best_fit_rows = _bestfit_rows(best_fits)
    fragment_rows = _fragment_rows(fragments)
    fragment_by_key = {row["key"]: row for row in fragment_rows}
    ranges = _range_rows(range_specs)
    ring_references = {
        "refe-A.pdb": _text(ring_ref_a, "A ring reference"),
        "refe-C.pdb": _text(ring_ref_c, "C ring reference"),
        "refe-PHE.pdb": _text(ring_ref_phe, "PHE ring reference"),
        "refe-HIS.pdb": _text(ring_ref_his, "HIS ring reference"),
        "refe-ARG.pdb": _text(ring_ref_arg, "ARG ring reference"),
        "refe-TRP.pdb": _text(ring_ref_trp, "TRP ring reference"),
    }
    if type(proteins) is not dict:
        raise ValueError("Protein range PDBs are not a folder.")

    output = []
    seen = set()
    protein_cache: dict[int, str] = {}
    for fitted in best_fit_rows:
        fragment = fragment_by_key.get(fitted["key"])
        if fragment is None:
            continue
        try:
            pdb_residues = _detector_residues(fitted["pdb"])
        except Exception as exc:
            pdb_residues = []
            pdb_error = f"best-fit PDB could not be parsed: {exc}"
        else:
            pdb_error = "best-fit PDB has no complete RNA residue rows"
        actual_residues = fragment["residues"]
        residue_map: dict[int, int] = {}
        for index, identity in enumerate(actual_residues, start=1):
            if type(identity) is not dict or type(identity.get("resid")) is not int:
                continue
            if len(pdb_residues) >= index:
                residue_map[pdb_residues[index - 1][1]] = index
            residue_map.setdefault(identity["resid"], index)
        for range_number, range_row in enumerate(ranges, start=2):
            filename_options = (
                f"range-{range_number}-aa.pdb",
                f"range-{range_number}.pdb",
            )
            protein_text = protein_cache.get(range_number)
            if protein_text is None:
                source = next((proteins[name] for name in filename_options if name in proteins), None)
                if source is None:
                    raise ValueError(f"Protein PDB for range-{range_number} is missing.")
                protein_text = _text(source, f"range-{range_number} protein PDB")
                protein_cache[range_number] = protein_text
            if not pdb_residues:
                # Keep a named per-fragment diagnostic instead of dropping a candidate.
                output.append({
                    "key": None,
                    "fragment": fitted["fragment"],
                    "fragment_key": fitted["key"],
                    "nucleotide": None,
                    "sequence": fitted["sequence"],
                    "range": range_number,
                    "protein_resid": None,
                    "error": pdb_error,
                })
                continue
            try:
                found = preparic_detect_stacking._detect_stacking(
                    [fitted["pdb"], protein_text],
                    ring_references=ring_references,
                    include_all_pairs=True,
                )
            except Exception as exc:
                output.append({
                    "key": None,
                    "fragment": fitted["fragment"],
                    "fragment_key": fitted["key"],
                    "nucleotide": None,
                    "sequence": fitted["sequence"],
                    "range": range_number,
                    "protein_resid": None,
                    "error": f"ring geometry failed for {fitted['key']} / range-{range_number}: {exc}",
                })
                continue
            for row in found:
                nucleotide_index = residue_map.get(row["rna_resid"])
                if nucleotide_index is None or nucleotide_index > len(actual_residues):
                    continue
                identity = actual_residues[nucleotide_index - 1]
                if type(identity) is not dict:
                    continue
                key = json.dumps([
                    row["rna_resname"], identity["resid"], range_number, row["protein_resid"]
                ], separators=(",", ":"), ensure_ascii=False)
                address = (key, fitted["fragment"], nucleotide_index)
                if address in seen:
                    raise ValueError(f"T18 produced duplicate best-fit stacking candidate {address}.")
                seen.add(address)
                output.append({
                    "key": key,
                    "fragment": fitted["fragment"],
                    "fragment_key": fitted["key"],
                    "nucleotide": "first" if nucleotide_index == 1 else "second",
                    "sequence": fitted["sequence"],
                    "conformer": fitted["conformer"],
                    "rotamer": fitted["rotamer"],
                    "grid": fitted["grid"],
                    "range": range_number,
                    "protein_resid": row["protein_resid"],
                    "distance": float(row["distance"]),
                    "phi_deg": float(row["slip_angle"]),
                    "theta_deg": float(row["tilt_angle"]),
                    "rna_centroid": row["rna_centroid"],
                    "protein_centroid": row["protein_centroid"],
                })
    output.sort(key=lambda row: (
        row.get("key") or "", row.get("fragment", 0), row.get("nucleotide") or "",
        row.get("range", 0), row.get("protein_resid") or 0,
    ))
    return {"schema": 1, "echo": _compute_echo(best_fits, range_specs), "rows": output}


def gaussian_check_jobs(best_fits: Any, fragments: Any, stackings: Any, selection: Any, clip: Any) -> dict | None:
    """T20: make one A9 job for each selected, unclipped best-fit candidate."""
    if best_fits is None or fragments is None or stackings is None or selection is None or clip is None:
        return None
    fit_rows = _bestfit_rows(best_fits)
    fit_by_key = {row["key"]: row for row in fit_rows}
    fragment_rows = _fragment_rows(fragments)
    bounds = _clip_bounds(clip, len(fragment_rows))
    candidates = _candidate_map(fragments, stackings, selection)
    jobs = []
    errors = []
    for candidate in candidates:
        if not _is_kept(candidate["fragment"], bounds):
            continue
        fitted = fit_by_key.get(candidate["fragment_key"])
        if fitted is None:
            errors.append({**candidate, "error": f"Best-fit pose {candidate['fragment_key']} is missing."})
            continue
        try:
            job = {
                "key": candidate["key"],
                "fragment": candidate["fragment"],
                "nucleotide": candidate["nucleotide"],
                "sequence": fitted["sequence"],
                "conformer": int(fitted["conformer"]),
                "rotamer": int(fitted["rotamer"]),
                "grid": [float(value) for value in fitted["grid"]],
                "range": int(candidate["range"]),
                "protein_resid": int(candidate["protein_resid"]),
            }
        except (KeyError, TypeError, ValueError, OverflowError) as exc:
            errors.append({**candidate, "error": f"Best-fit pose {candidate['fragment_key']} is malformed: {exc}"})
            continue
        jobs.append(job)
    return {
        "schema": 1,
        "echo": {
            "best_fits": best_fits["echo"],
            "selected_keys": _selected_keys(selection),
            "clip": clip,
        },
        "jobs": jobs,
        "errors": errors,
    }


def _class_filter_name(row: dict) -> str:
    name = f"rrm{row['rrm']}"
    if row["rrm"] in {90, 95}:
        if row.get("tilt30") is True:
            name += "-tilt30"
        if row.get("disp45") is True:
            name += "-disp45"
    return name


def validate_gaussian_classes(marker: Any, model_json: Any, selection: Any, gaussian_filter) -> None:
    """T19: validate active S18 against names in the mounted A7 document."""
    if marker is None:
        return None
    if type(marker) is not dict or set(marker) != {"schema", "basis", "active", "rows"} \
        or type(marker.get("schema")) is not int or marker["schema"] != 1 \
        or type(marker.get("basis")) is not list or type(marker.get("active")) is not bool \
        or type(marker.get("rows")) is not list:
        raise ValueError("Gaussian classes must be a schema-1 marker or null.")
    if marker["active"] is not True:
        return None
    spec, _model_value, names = load_gaussian_spec(model_json, gaussian_filter)
    selected = _selected_keys(selection)
    if sorted(marker["basis"]) != selected:
        raise ValueError("Gaussian class basis does not match the selected stacking keys.")
    row_keys = []
    for row in marker["rows"]:
        if type(row) is not dict or not isinstance(row.get("key"), str) or type(row.get("rrm")) is not int:
            raise ValueError("A Gaussian class row is malformed.")
        if row["rrm"] not in {50, 60, 70, 80, 90, 95}:
            raise ValueError(f"Gaussian class {row.get('key')!r} has an invalid RRM value.")
        if row["rrm"] in {90, 95}:
            if set(row) != {"key", "rrm", "tilt30", "disp45"} \
                or type(row.get("tilt30")) is not bool or type(row.get("disp45")) is not bool:
                raise ValueError("RRM-90/95 rows require boolean tilt30 and disp45 flags.")
        elif set(row) != {"key", "rrm"}:
            raise ValueError("RRM-50/60/70/80 rows cannot contain prefilter flags.")
        if _class_filter_name(row) not in names or _class_filter_name(row) not in spec["thresholds"]:
            raise ValueError(f"Gaussian class row {row['key']!r} maps to an illegal filter name.")
        row_keys.append(row["key"])
    if len(row_keys) != len(set(row_keys)) or sorted(row_keys) != selected:
        raise ValueError("Gaussian classes must contain exactly one row per selected stacking key.")
    return None


def gaussian_scores(
    bestfit_stackings: Any,
    reference_stackings: Any,
    selection: Any,
    clip: Any,
    fragments: Any,
    model_json: Any,
    gaussian_filter,
) -> dict | None:
    """T28: validate A7, score reference and best-fit descriptors, derive classes."""
    if bestfit_stackings is None or reference_stackings is None or selection is None \
        or clip is None or fragments is None or model_json is None:
        return None
    spec, model, names = load_gaussian_spec(model_json, gaussian_filter)
    selected = _selected_keys(selection)
    fragment_rows = _fragment_rows(fragments)
    bounds = _clip_bounds(clip, len(fragment_rows))
    reference_rows = _stacking_rows(reference_stackings)
    bestfit_rows = bestfit_stackings.get("rows") if type(bestfit_stackings) is dict else None
    if type(bestfit_stackings) is not dict or bestfit_stackings.get("schema") != 1 \
        or type(bestfit_stackings.get("echo")) is not dict or not isinstance(bestfit_rows, list):
        raise ValueError("T18 best-fit descriptors are not ready.")
    expected = _candidate_map(fragments, reference_stackings, selection)
    by_identity = {}
    for row in bestfit_rows:
        if type(row) is not dict:
            raise ValueError("T18 returned a malformed descriptor row.")
        key = row.get("key")
        fragment = row.get("fragment")
        nucleotide = row.get("nucleotide")
        if not isinstance(key, str) or type(fragment) is not int or nucleotide not in {"first", "second"}:
            # T18 may include non-target/error rows; they cannot be scored or joined.
            continue
        identity = (key, fragment, nucleotide)
        if identity in by_identity:
            raise ValueError(f"T18 returned duplicate candidate descriptor {identity}.")
        by_identity[identity] = row

    expected_kept = [candidate for candidate in expected if _is_kept(candidate["fragment"], bounds)]
    expected_identities = {
        (candidate["key"], candidate["fragment"], candidate["nucleotide"])
        for candidate in expected_kept
    }
    missing = [candidate for candidate in expected_kept
               if (candidate["key"], candidate["fragment"], candidate["nucleotide"]) not in by_identity]
    if missing:
        first = missing[0]
        raise ValueError(
            f"T18 has no best-fit ring descriptor for selected anchor {first['key']} / "
            f"{first['fragment_key']} {first['nucleotide']}."
        )

    thresholds = spec["thresholds"]
    candidates_out = []
    score_by_identity = {}
    for row in bestfit_rows:
        if not isinstance(row, dict):
            raise ValueError("T18 returned a malformed descriptor row.")
        if type(row.get("fragment")) is int and not _is_kept(row["fragment"], bounds):
            continue
        if not isinstance(row.get("key"), str) or type(row.get("fragment")) is not int \
            or row.get("nucleotide") not in {"first", "second"}:
            if row.get("error"):
                candidates_out.append(dict(row))
            continue
        try:
            d = float(row["distance"])
            phi = math.radians(float(row["phi_deg"]))
            theta = math.radians(float(row["theta_deg"]))
            if not all(math.isfinite(value) for value in (d, phi, theta)) or d < 0:
                raise ValueError("descriptor values must be finite and distance nonnegative")
            score = float(gaussian_filter.score(d, phi, theta, model=model))
            if not math.isfinite(score):
                raise ValueError("Gaussian score is not finite")
        except (KeyError, TypeError, ValueError, OverflowError) as exc:
            identity = (row.get("key"), row.get("fragment"), row.get("nucleotide"))
            if identity in expected_identities:
                raise ValueError(f"Invalid selected best-fit descriptor {identity}: {exc}") from exc
            candidates_out.append({**row, "error": f"Invalid best-fit descriptor: {exc}"})
            continue
        mask = spec["mask"]
        mask_pass = d < float(mask["distance_lt_A"]) \
            and phi < math.radians(float(mask["displacement_lt_deg"]))
        pose_pass = {
            name: _descriptor_passes(spec, name, d, phi, theta, score)
            for name in names
        }
        result = {
            **row,
            "score": score,
            "mask_pass": mask_pass,
            "pose_pass": pose_pass,
            "kept": _is_kept(row["fragment"], bounds),
        }
        candidates_out.append(result)
        score_by_identity[(row["key"], row["fragment"], row["nucleotide"])] = result

    reference_scores = []
    reference_by_key = {}
    for row in reference_rows:
        key = _stack_key(row)
        try:
            d = float(row["distance"])
            phi = math.radians(float(row["slip_angle"]))
            theta = math.radians(float(row["tilt_angle"]))
            score = float(gaussian_filter.score(d, phi, theta, model=model))
        except (KeyError, TypeError, ValueError, OverflowError) as exc:
            raise ValueError(f"Reference stacking {key} has invalid Gaussian descriptors: {exc}") from exc
        reference = {"key": key, "distance": d, "phi_deg": math.degrees(phi),
                     "theta_deg": math.degrees(theta), "score": score}
        reference_scores.append(reference)
        reference_by_key[key] = reference

    computed_classes = []
    for key in selected:
        rows = [score_by_identity.get((item["key"], item["fragment"], item["nucleotide"]))
                for item in expected_kept if item["key"] == key]
        rows = [row for row in rows if row is not None and row.get("error") is None]
        if not rows:
            raise ValueError(f"Selected stacking {key} has no unclipped best-fit candidates to score.")
        worst_score = min(float(row["score"]) for row in rows)
        base_rrms = sorted(
            int(match.group(1)) for name in thresholds
            if (match := re.fullmatch(r"rrm(50|60|70|80|90|95)", name)) is not None
        )
        if base_rrms != [50, 60, 70, 80, 90, 95]:
            raise ValueError("Gaussian specification is missing one or more base RRM classes.")
        passing = []
        for rrm in base_rrms:
            base_name = f"rrm{rrm}"
            if base_name not in thresholds:
                continue
            if worst_score >= float(thresholds[base_name]["cutoff"]):
                passing.append(rrm)
        if passing:
            rrm = min(passing)
            below_last = False
        else:
            rrm = max(base_rrms)
            below_last = True
        mask_pass = all(row["mask_pass"] for row in rows)
        tilt30 = all(float(row["theta_deg"]) < float(spec["filters"]["tilt30"]["lt_deg"]) for row in rows)
        disp45 = all(float(row["phi_deg"]) < float(spec["filters"]["disp45"]["lt_deg"]) for row in rows)
        computed = {"key": key, "rrm": rrm}
        if rrm in {90, 95}:
            computed["tilt30"] = tilt30
            computed["disp45"] = disp45
        filter_name = _class_filter_name(computed)
        if filter_name not in thresholds:
            raise ValueError(f"Computed class {filter_name!r} is not in the Gaussian specification.")
        reference = reference_by_key.get(key)
        if reference is None:
            raise ValueError(f"Selected reference stacking {key} has no reference score.")
        class_pass = not below_last and mask_pass and all(
            _descriptor_passes(spec, filter_name, float(row["distance"]),
                               math.radians(float(row["phi_deg"])),
                               math.radians(float(row["theta_deg"])), float(row["score"]))
            for row in rows
        )
        computed_classes.append({
            **computed,
            "auto_prefilters": {"tilt30": tilt30, "disp45": disp45},
            "filter_name": filter_name,
            "score": worst_score,
            "reference_score": reference["score"],
            "below_rrm95": below_last,
            "mask_pass": mask_pass,
            "class_pass": class_pass,
            "candidate_count": len(rows),
        })

    return {
        "schema": 1,
        "echo": {
            "bestfit_stackings": bestfit_stackings["echo"],
            "selected_keys": selected,
            "clip": clip,
        },
        "filter_set_version": spec["filter_set_version"],
        "filter_names": names,
        "classes": computed_classes,
        "candidate_scores": candidates_out,
        "reference_scores": reference_scores,
    }


def _validate_a9(result: Any) -> dict:
    if type(result) is not dict or type(result.get("schema")) is not int or result.get("schema") != 1 \
        or not isinstance(result.get("filter_set_version"), str) or not isinstance(result.get("jobs"), list) \
        or type(result.get("min_cell_fraction")) not in {int, float} \
        or float(result["min_cell_fraction"]) != _MIN_CELL_FRACTION:
        raise ValueError("A9 Gaussian check result must use schema 1.")
    return result


def merge_gaussian_validation(
    check_result: Any,
    jobs_result: Any,
    class_marker: Any,
    selection: Any,
    clip: Any,
    best_fits: Any,
    model_json: Any,
    gaussian_filter,
) -> dict | None:
    """T22: merge class-independent A9 rows with current editable S18."""
    if check_result is None or jobs_result is None or class_marker is None \
        or selection is None or clip is None or best_fits is None:
        return None
    if type(best_fits) is not dict or not isinstance(best_fits.get("echo"), dict):
        raise ValueError("Fresh best-fit rows are not ready for Gaussian validation.")
    spec, _model_value, names = load_gaussian_spec(model_json, gaussian_filter)
    result = _validate_a9(check_result)
    if type(jobs_result) is not dict or jobs_result.get("schema") != 1 \
        or not isinstance(jobs_result.get("jobs"), list) or not isinstance(jobs_result.get("errors"), list):
        raise ValueError("T20 Gaussian job list is malformed.")
    expected_jobs_echo = {
        "best_fits": best_fits.get("echo") if type(best_fits) is dict else None,
        "selected_keys": _selected_keys(selection),
        "clip": clip,
    }
    if jobs_result.get("echo") != expected_jobs_echo:
        raise ValueError("T20 jobs are stale for the current best fits, selection or clip.")
    if result["filter_set_version"] != spec["filter_set_version"]:
        raise ValueError("A9 results use a different Gaussian filter-set version.")
    validate_gaussian_classes(class_marker, model_json, selection, gaussian_filter)
    class_rows = class_marker.get("rows") if type(class_marker) is dict else None
    if not isinstance(class_rows, list) or class_marker.get("active") is not True:
        raise ValueError("Active Gaussian class rows are unavailable.")
    jobs_by_identity = {}
    for job in jobs_result["jobs"]:
        if type(job) is not dict:
            continue
        identity = (job.get("key"), job.get("fragment"), job.get("nucleotide"))
        if identity in jobs_by_identity:
            raise ValueError(f"T20 returned duplicate job identity {identity}.")
        jobs_by_identity[identity] = job
    checked_by_identity = {}
    for row in result["jobs"]:
        if type(row) is not dict:
            continue
        identity = (row.get("key"), row.get("fragment"), row.get("nucleotide"))
        if identity in checked_by_identity:
            raise ValueError(f"A9 returned duplicate job identity {identity}.")
        checked_by_identity[identity] = row
    extra_results = set(checked_by_identity) - set(jobs_by_identity)
    if extra_results:
        raise ValueError(f"A9 returned jobs that were not present in T20: {sorted(extra_results)}.")

    priority = {"green": 0, "implied": 1, "amber": 2, "red": 3, "error": 4}
    merged_rows = []
    for class_row in class_rows:
        key = class_row["key"]
        filter_name = _class_filter_name(class_row)
        if filter_name not in names:
            raise ValueError(f"Editable class {filter_name!r} is not legal in the current model.")
        candidates = []
        for identity, job in jobs_by_identity.items():
            if identity[0] != key:
                continue
            checked = checked_by_identity.get(identity)
            if checked is None:
                candidates.append({
                    "fragment": identity[1], "nucleotide": identity[2], "state": "error",
                    "error": "A9 returned no row for this T20 job.",
                })
                continue
            if checked.get("error"):
                candidates.append({
                    "fragment": identity[1], "nucleotide": identity[2], "state": "error",
                    "error": checked["error"],
                })
                continue
            filter_result = checked.get("filters", {}).get(filter_name)
            if type(filter_result) is not dict:
                candidates.append({
                    "fragment": identity[1], "nucleotide": identity[2], "state": "error",
                    "error": f"A9 omitted filter {filter_name}.",
                })
                continue
            implied_by = filter_result.get("implied_by")
            if implied_by is not None:
                candidates.append({
                    "fragment": identity[1], "nucleotide": identity[2],
                    "state": "implied", "implied_by": implied_by,
                })
                continue
            try:
                fill = float(filter_result["fill"])
                emitted = filter_result["emitted"]
                pose_pass = filter_result["pose_pass"]
                if not math.isfinite(fill) or fill < 0 or type(emitted) is not bool or type(pose_pass) is not bool:
                    raise ValueError("A9 filter values are malformed.")
            except (KeyError, TypeError, ValueError, OverflowError) as exc:
                candidates.append({
                    "fragment": identity[1], "nucleotide": identity[2], "state": "error",
                    "error": f"A9 filter {filter_name} is malformed: {exc}",
                })
                continue
            minimum_fill = _MIN_CELL_FRACTION
            if emitted and pose_pass and fill >= minimum_fill:
                state = "green"
            elif fill == 0.0:
                state = "red"
            else:
                state = "amber"
            candidates.append({
                "fragment": identity[1], "nucleotide": identity[2],
                "state": state, "fill": fill, "emitted": emitted,
                "pose_pass": pose_pass, "score": filter_result.get("score"),
            })
        for error in jobs_result["errors"]:
            if isinstance(error, dict) and error.get("key") == key:
                candidates.append({
                    "fragment": error.get("fragment"), "nucleotide": error.get("nucleotide"),
                    "state": "error", "error": error.get("error", "The best-fit candidate is unavailable."),
                })
        if not candidates:
            candidates.append({"state": "error", "error": "No unclipped candidates were checked."})
        worst = max(candidates, key=lambda value: priority.get(value["state"], 4))["state"]
        merged_rows.append({"key": key, "filter_name": filter_name, "worst_state": worst, "candidates": candidates})

    return {
        "schema": 1,
        "echo": {
            "classes": class_rows,
            "clip": clip,
            "selected_keys": _selected_keys(selection),
            "best_fits": best_fits["echo"],
        },
        "filter_set_version": result["filter_set_version"],
        "filter_names": names,
        "rows": merged_rows,
    }


def read_gaussian_check_result(result_folder: Any) -> dict:
    """T21 adapter: parse the CLI's folder result into one plain cell."""
    if type(result_folder) is not dict:
        raise ValueError("A9 result folder is not available.")
    value = result_folder.get("gaussian-check.json")
    if value is None:
        value = result_folder.get("gaussian-check")
    return _validate_a9(_json_object(value, "A9 result"))

