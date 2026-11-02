"""Pure Python helpers for the preparation fragments and stage-4 hand-off."""

from __future__ import annotations

import json
import math
import re
from decimal import Decimal, InvalidOperation
from io import StringIO
from pathlib import PurePosixPath
from typing import Any


_FRAGMENT_FILE = re.compile(r"^frag-([1-9][0-9]*)-([ACGU]{2})-refe\.pdb$")
_FRAGMENT_KEY = re.compile(r"^frag-([1-9][0-9]*)-([ACGU]{2})$")
_PAIR_KEYS = {"down", "up", "sequence", "ovrmsd", "crmsd"}
_BASE_KEYS = {"fragment", "sequence", "nuc1_ovrmsd", "nuc2_ovrmsd"}


def _canonical_base(resname: str) -> str | None:
    name = resname.strip().upper()
    if len(name) == 2 and name[0] in {"R", "D"}:
        name = name[1]
    return name if name in {"A", "C", "G", "U"} else None


def _rna_residues(rna_pdb: str | bytes, rna_range: dict) -> list[dict]:
    from Bio.PDB import PDBParser

    if isinstance(rna_pdb, bytes):
        try:
            rna_pdb = rna_pdb.decode("utf-8")
        except UnicodeDecodeError as exc:
            raise ValueError("The prepared RNA PDB is not UTF-8 text.") from exc
    if not isinstance(rna_pdb, str) or not rna_pdb.strip():
        raise ValueError("The prepared RNA PDB is empty.")
    if type(rna_range) is not dict or set(rna_range) != {"chain", "start", "end"} \
        or not isinstance(rna_range["chain"], str) or not rna_range["chain"] \
        or type(rna_range["start"]) is not int or type(rna_range["end"]) is not int \
        or rna_range["start"] > rna_range["end"]:
        raise ValueError("The range-1 RNA address is invalid.")

    try:
        structure = PDBParser(QUIET=True).get_structure("rna", StringIO(rna_pdb))
    except Exception as exc:
        raise ValueError(f"The prepared RNA PDB is invalid: {exc}") from exc
    model = next(structure.get_models(), None)
    if model is None:
        raise ValueError("The prepared RNA PDB has no model.")

    # aareduce emits one range file and PDB chain ids are necessarily one
    # character. The requested assembly-copy name remains the S22 identity.
    chains = list(model.get_chains())
    if len(chains) != 1:
        raise ValueError("The prepared RNA PDB must contain exactly one chain.")
    residues = []
    start, end = rna_range["start"], rna_range["end"]
    for residue in chains[0].get_residues():
        hetflag, resid, insertion_code = residue.id
        if not start <= int(resid) <= end:
            continue
        base = _canonical_base(residue.get_resname())
        if base is None:
            raise ValueError(
                f"Prepared RNA residue {resid}{insertion_code.strip()} has unsupported base "
                f"{residue.get_resname()!r}."
            )
        residues.append({
            "chain": rna_range["chain"],
            "resid": int(resid),
            "insertion_code": insertion_code.strip(),
            "base": base,
        })
    residues.sort(key=lambda item: (item["resid"], item["insertion_code"]))
    if not residues or residues[0]["resid"] != start or residues[-1]["resid"] != end:
        observed = {item["resid"] for item in residues}
        missing = [value for value in (start, end) if value not in observed]
        missing_text = ", ".join(str(value) for value in missing) or "range endpoints"
        raise ValueError(f"Prepared RNA is missing residues {missing_text}.")
    return residues


def _fit_echo(rna_pdb: str | bytes, pdb_code: str, rna_range: dict) -> dict:
    """Build the shared fit echo from the exact T3 RNA identity list."""
    return {
        "pdb_code": str(pdb_code).strip().lower(),
        "rna": dict(rna_range),
        "resids": _rna_residues(rna_pdb, rna_range),
    }


def _parse_best_fits(best_fit: dict[str, bytes], fit_echo: dict) -> dict:
    """T9: parse T5's fit TSV and copy PDB texts into the plain S10 cell."""
    import csv
    from io import StringIO

    if type(best_fit) is not dict or not isinstance(best_fit.get("best-fits.tsv"), bytes):
        raise ValueError("The best-fit folder is missing best-fits.tsv.")
    try:
        tsv = best_fit["best-fits.tsv"].decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ValueError("The best-fit table is not UTF-8 text.") from exc
    expected_header = ["fragment", "sequence", "conformer", "rotamer", "grid_x", "grid_y", "grid_z", "rmsd"]
    reader = csv.DictReader(StringIO(tsv), delimiter="\t")
    if reader.fieldnames != expected_header:
        raise ValueError(f"The best-fit table columns must be {' '.join(expected_header)}.")
    rows = []
    seen = set()
    for source in reader:
        try:
            fragment = int(source["fragment"])
            sequence = source["sequence"]
            conformer = int(source["conformer"])
            rotamer = int(source["rotamer"])
            grid = [float(source[key]) for key in ("grid_x", "grid_y", "grid_z")]
            rmsd = float(source["rmsd"])
        except (TypeError, ValueError, KeyError) as exc:
            raise ValueError(f"Invalid row in best-fit table: {source!r}.") from exc
        if fragment < 1 or re.fullmatch(r"[ACGU]{2}", sequence or "") is None \
            or conformer < 1 or rotamer < 1 or not all(math.isfinite(value) for value in [*grid, rmsd]):
            raise ValueError(f"Invalid row in best-fit table: {source!r}.")
        key = f"frag-{fragment}-{sequence}"
        if key in seen:
            raise ValueError(f"Duplicate best-fit row {key}.")
        seen.add(key)
        pdb_path = f"pdbs/{key}-bestfit.pdb"
        pdb_bytes = best_fit.get(pdb_path)
        if not isinstance(pdb_bytes, bytes):
            raise ValueError(f"The best-fit folder is missing {pdb_path}.")
        try:
            pdb = pdb_bytes.decode("utf-8")
        except UnicodeDecodeError as exc:
            raise ValueError(f"The best-fit file {pdb_path} is not UTF-8 PDB text.") from exc
        rows.append({
            "key": key,
            "fragment": fragment,
            "sequence": sequence,
            "conformer": conformer,
            "rotamer": rotamer,
            "grid": grid,
            "rmsd": rmsd,
            "pdb": pdb,
        })
    rows.sort(key=lambda row: (row["fragment"], row["sequence"]))
    return {"schema": 1, "echo": fit_echo, "rows": rows}


def _list_fragments(
    refe_fragments: dict[str, bytes],
    rna_pdb: str | bytes,
    pdb_code: str,
    rna_range: dict,
) -> dict:
    """T24: publish the canonical fragment map before any best-fit work.

    The T6 files use template numbering (chain B, residues 1 and 2). This
    function uses their names and keys for ordinal windows, then maps those
    windows onto the real residue identities from T3's prepared RNA PDB.
    """
    if not isinstance(refe_fragments, dict):
        raise ValueError("The reference-fragment folder is invalid.")
    if not isinstance(pdb_code, str) or not pdb_code.strip():
        raise ValueError("A PDB code is required to list fragments.")
    residues = _rna_residues(rna_pdb, rna_range)
    by_number: dict[int, str] = {}
    for name in sorted(refe_fragments):
        if not isinstance(refe_fragments[name], bytes):
            raise ValueError(f"Reference fragment {name!r} is not bytes.")
        match = _FRAGMENT_FILE.fullmatch(PurePosixPath(name).name)
        if not match:
            raise ValueError(f"Unexpected reference-fragment file {name!r}.")
        number, sequence = int(match.group(1)), match.group(2)
        if number in by_number:
            raise ValueError(f"Duplicate reference fragment number {number}.")
        if number >= len(residues):
            raise ValueError(f"Reference fragment {name!r} exceeds the RNA residue count.")
        observed_sequence = residues[number - 1]["base"] + residues[number]["base"]
        if observed_sequence != sequence:
            raise ValueError(
                f"Reference fragment {name!r} has sequence {sequence}, but the RNA window "
                f"is {observed_sequence}."
            )
        by_number[number] = sequence

    expected = {
        index + 1: residues[index]["base"] + residues[index + 1]["base"]
        for index in range(len(residues) - 1)
        if residues[index + 1]["resid"] == residues[index]["resid"] + 1
    }
    if by_number != expected:
        missing = sorted(set(expected) - set(by_number))
        extra = sorted(set(by_number) - set(expected))
        details = []
        if missing:
            details.append("missing " + ", ".join(f"frag-{number}-{expected[number]}" for number in missing))
        if extra:
            details.append("unexpected " + ", ".join(f"frag-{number}-{by_number[number]}" for number in extra))
        raise ValueError("T6 reference fragments do not match the prepared RNA: " + "; ".join(details) + ".")

    fragments = []
    for number, sequence in sorted(by_number.items()):
        fragments.append({
            "number": number,
            "key": f"frag-{number}-{sequence}",
            "sequence": sequence,
            "residues": [dict(residues[number - 1]), dict(residues[number])],
        })
    adjacent_pairs = []
    fragments_by_number = {fragment["number"]: fragment for fragment in fragments}
    for number, fragment in fragments_by_number.items():
        next_fragment = fragments_by_number.get(number + 1)
        if next_fragment is None:
            continue
        if fragment["residues"][1] != next_fragment["residues"][0]:
            continue
        adjacent_pairs.append({
            "down": number,
            "up": number + 1,
            "sequence": fragment["sequence"] + next_fragment["sequence"][1],
        })

    return {
        "schema": 1,
        "echo": {"pdbcode": pdb_code.strip().upper(), "rna": dict(rna_range)},
        "residues": residues,
        "fragments": fragments,
        "adjacent_pairs": adjacent_pairs,
    }


def _candidate_fragments(fragments: dict, stackings: list[dict] | dict) -> list[dict]:
    """Map fresh T7 stackings onto T24's structural fragment identities.

    This stays a separate calculation because T24's identity depends only on
    the prepared RNA, range 1, and T6 output. Protein edits change T7 and this
    candidate map without making the structural S22 result stale.
    """
    if type(fragments) is not dict or type(fragments.get("schema")) is not int \
        or fragments.get("schema") != 1 or not isinstance(fragments.get("fragments"), list) \
        or not isinstance(fragments.get("residues"), list):
        raise ValueError("The canonical fragment list is not ready.")
    if isinstance(stackings, dict):
        if type(stackings.get("schema")) is not int or stackings.get("schema") != 1 \
            or not isinstance(stackings.get("rows"), list):
            raise ValueError("The stacking table must use schema 1.")
        stacking_rows = stackings["rows"]
    else:
        stacking_rows = stackings
    if not isinstance(stacking_rows, list):
        raise ValueError("The stacking table used for fragment candidates is invalid.")

    windows: dict[int, list[tuple[int, int, str]]] = {}
    for fragment in fragments["fragments"]:
        if type(fragment) is not dict or type(fragment.get("number")) is not int \
            or not isinstance(fragment.get("key"), str) or not isinstance(fragment.get("residues"), list):
            raise ValueError("The canonical fragment list contains an invalid row.")
        for nucleotide, identity in enumerate(fragment["residues"], start=1):
            if type(identity) is dict and type(identity.get("resid")) is int:
                windows.setdefault(identity["resid"], []).append(
                    (fragment["number"], nucleotide, fragment["key"])
                )

    residues = fragments["residues"]
    candidates: set[tuple[str, str, int]] = set()
    for row in stacking_rows:
        if type(row) is not dict:
            raise ValueError("A stacking row is invalid while building fragment candidates.")
        try:
            rna_resname = str(row["rna_resname"]).strip()
            rna_resid = row["rna_resid"]
            protein_range = row["protein_range_pdb_id"]
            protein_resid = row["protein_resid"]
        except KeyError as exc:
            raise ValueError(f"A stacking row is missing {exc.args[0]!r}.") from exc
        if not rna_resname or type(rna_resid) is not int \
            or type(protein_range) is not int or type(protein_resid) is not int:
            raise ValueError("Stacking residue addresses must be valid integers and names.")
        base = _canonical_base(rna_resname)
        identities = [identity for identity in residues
                      if type(identity) is dict and identity.get("resid") == rna_resid]
        if not base or not any(identity.get("base") == base for identity in identities):
            continue
        key = json.dumps([rna_resname, rna_resid, protein_range, protein_resid], separators=(",", ":"))
        for _number, nucleotide, fragment_key in windows.get(rna_resid, []):
            candidates.add((key, fragment_key, nucleotide))
    return [
        {"key": key, "fragment": fragment_key, "nucleotide": nucleotide}
        for key, fragment_key, nucleotide in sorted(candidates)
    ]


def _share_range_pdbs_aa(
    pdb_files: dict[str, bytes], range_specs: dict, exclude_pdb_code: str,
) -> dict:
    """T25: copy the aareduce PDB texts to a plain share without rewriting."""
    if not isinstance(pdb_files, dict):
        raise ValueError("The aareduce range-PDB folder is invalid.")
    if not isinstance(range_specs, dict) or type(range_specs.get("schema")) is not int \
        or range_specs.get("schema") != 1 \
        or not isinstance(range_specs.get("rows"), list) or not isinstance(range_specs.get("pdbcode"), str):
        raise ValueError("The current range specification is invalid.")
    if not isinstance(exclude_pdb_code, str) or not exclude_pdb_code.strip():
        raise ValueError("The exclusion PDB code is empty.")
    pdbs = []
    for index, _row in enumerate(range_specs["rows"], start=1):
        filename = f"range-{index}-aa.pdb"
        if filename not in pdb_files:
            raise ValueError(f"The aareduce folder is missing {filename}.")
        data = pdb_files[filename]
        if not isinstance(data, bytes):
            raise ValueError(f"The aareduce file {filename} is not bytes.")
        try:
            pdb = data.decode("utf-8")
        except UnicodeDecodeError as exc:
            raise ValueError(f"The aareduce file {filename} is not UTF-8 PDB text.") from exc
        pdbs.append({"range": f"range-{index}", "name": filename, "pdb": pdb})
    return {
        "schema": 1,
        "pdb_code": exclude_pdb_code.strip().lower(),
        "range_specs": range_specs,
        "pdbs": pdbs,
    }


def _reference_fragment_pdbs(refe_fragments: dict[str, bytes], fragments: dict) -> dict:
    """Publish exact T6 reference texts for stage-4 Mol* overlays."""
    if type(fragments) is not dict or type(fragments.get("schema")) is not int \
        or fragments.get("schema") != 1 \
        or not isinstance(fragments.get("fragments"), list):
        raise ValueError("The canonical fragment list is not ready.")
    if type(refe_fragments) is not dict:
        raise ValueError("The reference-fragment folder is invalid.")
    rows = []
    for fragment in fragments["fragments"]:
        if type(fragment) is not dict or not isinstance(fragment.get("key"), str):
            raise ValueError("The canonical fragment list contains an invalid row.")
        name = f"{fragment['key']}-refe.pdb"
        text = refe_fragments.get(name)
        if not isinstance(text, bytes):
            raise ValueError(f"The reference-fragment folder is missing {name}.")
        try:
            pdb = text.decode("utf-8")
        except UnicodeDecodeError as exc:
            raise ValueError(f"Reference fragment {name} is not UTF-8 PDB text.") from exc
        rows.append({"key": fragment["key"], "pdb": pdb})
    return {"schema": 1, "echo": fragments["echo"], "rows": rows}


def _valid_threshold_value(value: object) -> bool:
    if not isinstance(value, str) or re.fullmatch(r"(?:0|[1-9][0-9]*)(?:\.[0-9]{1,2})?", value) is None:
        return False
    try:
        number = Decimal(value)
    except InvalidOperation:
        return False
    return number.is_finite() and number > 0


def _valid_threshold_row(row: object, keys: set[str], value_keys: tuple[str, ...]) -> bool:
    return type(row) is dict and set(row) == keys and all(
        row[key] is None or _valid_threshold_value(row[key]) for key in value_keys
    )


def _validate_threshold_edits(edits: object, fragments: object) -> None:
    """T10: validate the edit marker against T24's canonical fragment map.

    Null and inactive decisions are valid sentinels. Blank individual values
    are allowed here; the 4→5 gate decides completeness for kept rows.
    """
    if edits is None:
        return None
    if type(edits) is not dict or type(edits.get("schema")) is not int \
        or edits.get("schema") != 1 or type(edits.get("active")) is not bool:
        raise ValueError("Threshold edits must be a schema-1 decision or null.")
    if edits["active"] is False:
        return None
    if type(fragments) is not dict or type(fragments.get("schema")) is not int \
        or fragments.get("schema") != 1:
        raise ValueError("The canonical fragment list is not ready.")
    if set(edits) != {"schema", "basis", "active", "pairs", "bases"}:
        raise ValueError("Threshold edits do not match schema 1.")
    if edits.get("basis") != fragments.get("echo"):
        raise ValueError("Threshold edits belong to a different RNA range.")
    expected_pairs = fragments.get("adjacent_pairs")
    expected_fragments = fragments.get("fragments")
    if not isinstance(expected_pairs, list) or not isinstance(expected_fragments, list):
        raise ValueError("The canonical fragment list is malformed.")
    pairs = edits.get("pairs")
    bases = edits.get("bases")
    if not isinstance(pairs, list) or not isinstance(bases, list):
        raise ValueError("Threshold pair and base rows must be lists.")
    expected_pair_rows = [(row.get("down"), row.get("up"), row.get("sequence")) for row in expected_pairs if isinstance(row, dict)]
    actual_pair_rows = []
    for row in pairs:
        if not _valid_threshold_row(row, _PAIR_KEYS, ("ovrmsd", "crmsd")):
            raise ValueError("A pair threshold row has invalid fields or values.")
        if type(row["down"]) is not int or type(row["up"]) is not int \
            or row["up"] != row["down"] + 1 or not isinstance(row["sequence"], str):
            raise ValueError("A pair threshold row has invalid residue addresses.")
        actual_pair_rows.append((row["down"], row["up"], row["sequence"]))
    if actual_pair_rows != expected_pair_rows:
        raise ValueError("Threshold pair rows do not match the current fragment list.")

    expected_base_rows = [
        (row.get("number"), row.get("sequence")) for row in expected_fragments if isinstance(row, dict)
    ]
    actual_base_rows = []
    for row in bases:
        if not _valid_threshold_row(row, _BASE_KEYS, ("nuc1_ovrmsd", "nuc2_ovrmsd")):
            raise ValueError("A base threshold row has invalid fields or values.")
        if type(row["fragment"]) is not int or row["fragment"] < 1 or not isinstance(row["sequence"], str):
            raise ValueError("A base threshold row has an invalid fragment address.")
        actual_base_rows.append((row["fragment"], row["sequence"]))
    if actual_base_rows != expected_base_rows:
        raise ValueError("Threshold base rows do not match the current fragment list.")
    return None


def _valid_clip(clip: object, fragment_count: int, selected_candidates: list[int] | tuple[int, ...] = ()) -> bool:
    if type(clip) is not dict or set(clip) != {"schema", "basis", "clip1", "clip2"} \
        or type(clip.get("schema")) is not int or clip.get("schema") != 1 \
        or type(clip.get("clip1")) is not int or type(clip.get("clip2")) is not int \
        or clip["clip1"] < 0 or clip["clip2"] < 0 \
        or type(fragment_count) is not int or fragment_count < 1 \
        or clip["clip1"] + clip["clip2"] >= fragment_count:
        return False
    if clip.get("basis") is not None and (
        type(clip["basis"]) is not dict or set(clip["basis"]) != {"pdbcode", "rna"}
        or not isinstance(clip["basis"].get("pdbcode"), str)
        or type(clip["basis"].get("rna")) is not dict
    ):
        return False
    if not isinstance(selected_candidates, (list, tuple)):
        return False
    candidates = []
    for number in selected_candidates:
        if type(number) is not int or number < 1 or number > fragment_count:
            return False
        candidates.append(number)
    return all(clip["clip1"] < number <= fragment_count - clip["clip2"] for number in candidates)


def _selected_candidate_numbers(fragments: dict, selection: object, stackings: object) -> list[int]:
    if selection is None or (type(selection) is dict and selection.get("active") is False):
        return []
    if type(selection) is not dict or type(selection.get("schema")) is not int \
        or selection.get("schema") != 1 or type(selection.get("active")) is not bool \
        or set(selection) != {"schema", "basis", "active", "selected"} \
        or type(selection.get("selected")) is not dict \
        or any(type(flag) is not bool for flag in selection["selected"].values()):
        raise ValueError("The stacking selection must be a schema-1 decision.")
    if type(fragments) is not dict or type(fragments.get("schema")) is not int \
        or fragments.get("schema") != 1:
        raise ValueError("The canonical fragment list is not ready.")
    if selection.get("active") is not True:
        return []
    selected = {key for key, value in selection["selected"].items() if value is True}
    stacking_rows = stackings.get("rows") if type(stackings) is dict else stackings
    if type(stackings) is dict and (type(stackings.get("schema")) is not int
        or stackings.get("schema") != 1 or not isinstance(stacking_rows, list)):
        raise ValueError("The fresh stacking table is malformed.")
    if not isinstance(stacking_rows, list):
        raise ValueError("The fresh stacking table is not ready.")
    available = set()
    for row in stacking_rows:
        if type(row) is not dict:
            raise ValueError("The fresh stacking table contains an invalid row.")
        try:
            available.add(json.dumps([
                str(row["rna_resname"]).strip(), row["rna_resid"],
                row["protein_range_pdb_id"], row["protein_resid"],
            ], separators=(",", ":")))
        except KeyError as exc:
            raise ValueError(f"A stacking row is missing {exc.args[0]!r}.") from exc
    if selected - available:
        raise ValueError("The stacking selection contains keys outside the fresh stacking table.")
    numbers = set()
    candidates = _candidate_fragments(fragments, stackings)
    for candidate in candidates:
        if type(candidate) is not dict or candidate.get("key") not in selected:
            continue
        match = _FRAGMENT_KEY.fullmatch(str(candidate.get("fragment", "")))
        if match is not None:
            numbers.add(int(match.group(1)))
    return sorted(numbers)


def _clamp_clip(clip: dict, fragment_count: int, selected_candidates: list[int] | tuple[int, ...]) -> dict:
    if type(clip) is not dict or type(fragment_count) is not int or fragment_count < 1:
        raise ValueError("Cannot clamp a malformed clip or fragment list.")
    if type(clip.get("clip1")) is not int or type(clip.get("clip2")) is not int \
        or clip["clip1"] < 0 or clip["clip2"] < 0:
        raise ValueError("Clip counts must be nonnegative integers.")
    candidates = list(selected_candidates)
    if any(type(number) is not int or not 1 <= number <= fragment_count for number in candidates):
        raise ValueError("Selected anchor candidates contain an invalid fragment number.")
    result = dict(clip)
    if candidates:
        result["clip1"] = min(clip["clip1"], min(candidates) - 1)
        result["clip2"] = min(clip["clip2"], fragment_count - max(candidates))
    if result["clip1"] + result["clip2"] >= fragment_count:
        result["clip1"] = min(result["clip1"], fragment_count - 1)
        result["clip2"] = min(result["clip2"], fragment_count - result["clip1"] - 1)
    return result


def _clip_fragments(
    best_fits: dict,
    refe_fragments: dict[str, bytes],
    thresholds: dict,
    clip: dict,
    fragments: dict,
    selected_candidates: list[int] | tuple[int, ...] = (),
) -> dict[str, bytes]:
    """T11: make the kept-only B11 folder without changing source PDB bytes."""
    if thresholds is None or (type(thresholds) is dict and thresholds.get("active") is False):
        return {}
    _validate_threshold_edits(thresholds, fragments)
    if type(fragments) is not dict or type(fragments.get("schema")) is not int \
        or fragments.get("schema") != 1:
        raise ValueError("The canonical fragment list must use schema 1.")
    fragment_rows = fragments.get("fragments")
    if not isinstance(fragment_rows, list) or not fragment_rows:
        raise ValueError("The canonical fragment list is empty or invalid.")
    fragment_count = len(fragment_rows)
    if not _valid_clip(clip, fragment_count, selected_candidates):
        raise ValueError("The clip is invalid or would remove a selected stacking candidate.")
    if clip.get("basis") != fragments.get("echo"):
        raise ValueError("The clip belongs to a different RNA range.")
    if type(best_fits) is not dict or type(best_fits.get("schema")) is not int \
        or best_fits.get("schema") != 1 or not isinstance(best_fits.get("rows"), list):
        raise ValueError("The best-fit result is not ready.")
    if type(refe_fragments) is not dict:
        raise ValueError("The reference-fragment folder is invalid.")
    kept_rows = fragment_rows[clip["clip1"]:fragment_count - clip["clip2"]]
    kept_keys = {row["key"] for row in kept_rows}
    fit_by_key = {}
    for row in best_fits["rows"]:
        if type(row) is not dict or not isinstance(row.get("key"), str):
            raise ValueError("The best-fit table contains an invalid row.")
        if row["key"] in fit_by_key:
            raise ValueError(f"Duplicate best-fit row {row['key']!r}.")
        fit_by_key[row["key"]] = row

    output: dict[str, bytes] = {}
    for row in kept_rows:
        key = row["key"]
        match = _FRAGMENT_KEY.fullmatch(key)
        if not match:
            raise ValueError(f"Invalid canonical fragment key {key!r}.")
        number, sequence = int(match.group(1)), match.group(2)
        reference_name = f"frag-{number}-{sequence}-refe.pdb"
        reference = refe_fragments.get(reference_name)
        if not isinstance(reference, bytes):
            raise ValueError(f"The reference-fragment folder is missing {reference_name}.")
        fit = fit_by_key.get(key)
        if fit is None or not isinstance(fit.get("pdb"), str):
            raise ValueError(f"The best-fit result is missing PDB text for {key}.")
        output[f"reference/{reference_name}"] = reference
        output[f"bestfit/frag-{number}-{sequence}-bestfit.pdb"] = fit["pdb"].encode("utf-8")

    pairs = thresholds["pairs"]
    # Build membership from exact fragment numbers to avoid sequence-key ambiguity.
    kept_numbers = {int(_FRAGMENT_KEY.fullmatch(key).group(1)) for key in kept_keys}
    kept_pairs = [row for row in pairs if row["down"] in kept_numbers and row["up"] in kept_numbers]
    kept_bases = [row for row in thresholds["bases"] if row["fragment"] in kept_numbers]
    kept_residue_ids = {json.dumps(identity, sort_keys=True) for row in kept_rows for identity in row["residues"]}
    kept_residues = [identity for identity in fragments["residues"] if json.dumps(identity, sort_keys=True) in kept_residue_ids]
    output["thresholds.json"] = (json.dumps({"schema": 1, "pairs": kept_pairs, "bases": kept_bases}, sort_keys=True) + "\n").encode("utf-8")
    output["clip.json"] = (json.dumps({
        "schema": 1,
        "clip1": clip["clip1"],
        "clip2": clip["clip2"],
        "kept_residues": kept_residues,
    }, sort_keys=True) + "\n").encode("utf-8")
    return output
