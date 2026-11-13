"""Pure E4 export-tree renderer tests using committed normalized data.

The fixture is deliberately independent of the live Preparic mounts. It adds
the exact prepared range-1 bytes and kept-fragment PDB payloads omitted by the
planner-only fixture, then sends the normalized snapshot and recipes directly
to ``render_tree``.
"""

from __future__ import annotations

from copy import deepcopy
import base64
import hashlib
import json
from pathlib import Path

import pytest
import yaml

import export_gens
import export_snapshot
import export_tree


_DATA = Path(__file__).parent / "data" / "export-gens"


def _restore_fixture_values(value):
    if isinstance(value, dict):
        if set(value) == {"$bytes_base64"}:
            return base64.b64decode(value["$bytes_base64"], validate=True)
        return {key: _restore_fixture_values(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_restore_fixture_values(item) for item in value]
    return value


def _complete_tree_snapshot():
    snapshot = _restore_fixture_values(json.loads((_DATA / "1b7f.snapshot.json").read_text()))

    # E3 needs only the protein B1 files. E4 renders every prepared range, so
    # supply the raw RNA range as well. B4 is exactly B3's all-atom range 1.
    snapshot["range_pdb_files"]["range-1.pdb"] = (
        b"REMARK exact prepared raw RNA range\r\nEND\r\n"
    )
    # anchor-refe-bases.pdb reads base atoms from the reference, so B4/B3's
    # synthetic range 1 carries complete nucleobases for every RNA position.
    snapshot["range_pdb_files_aa"]["range-1-aa.pdb"] = _reference_pdb(snapshot["fragments"])
    snapshot["reference_pdb"] = snapshot["range_pdb_files_aa"]["range-1-aa.pdb"]
    snapshot["reference_sha256"] = hashlib.sha256(snapshot["reference_pdb"]).hexdigest()

    # The E3 fixture records only planner inputs. Add exact, synthetic B11
    # payloads for every kept fragment, along with the source-only metadata
    # members that the renderer must not copy into the archive.
    clipped = {
        "clip.json": json.dumps(snapshot["clip"], sort_keys=True).encode("utf-8") + b"\n",
        "thresholds.json": b'{"schema":1}\n',
    }
    for fragment in snapshot["kept_fragments"]:
        number, sequence = fragment["number"], fragment["sequence"]
        clipped[f"reference/frag-{number}-{sequence}-refe.pdb"] = (
            f"REMARK exact kept reference fragment {number}\r\nEND\r\n".encode("ascii")
        )
        clipped[f"bestfit/frag-{number}-{sequence}-bestfit.pdb"] = (
            f"REMARK exact kept bestfit fragment {number}\r\nEND\r\n".encode("ascii")
        )
    snapshot["clipped_fragments"] = clipped
    snapshot["best_fit_statistics"] = _best_fit_statistics(snapshot)
    return snapshot


def _rna_residues(fragments):
    residues = [row["residues"][0] for row in fragments]
    residues.append(fragments[-1]["residues"][1])
    return residues


def _reference_pdb(fragments):
    lines = ["REMARK exact prepared RNA with complete nucleobases"]
    serial = 0
    for residue in _rna_residues(fragments):
        for name in ("P", "C1'", *export_snapshot.ANCHOR_REFE_BASE_ATOMS[residue["base"]], "O3'"):
            serial += 1
            atom_name = f" {name:<3}"
            lines.append(
                f"ATOM  {serial:5d} {atom_name} {'R' + residue['base']:>3} A{residue['resid']:4d}    "
                f"{serial * 0.5:8.3f}{serial * -0.25:8.3f}{serial * 0.125:8.3f}  1.00  0.00"
                f"           {name[0]}"
            )
    lines.append("END")
    return ("\r\n".join(lines) + "\r\n").encode("ascii")


def _best_fit_statistics(snapshot):
    """Synthetic normalised statistics; normalisation itself is tested in E2."""
    bases = {row["fragment"]: row for row in snapshot["threshold_bases"]}
    energies = {(row["fragment"], row["range"]): row["energy"] for row in snapshot["energy_cells"]}
    protein_ranges = [f"range-{index}" for index in range(2, len(snapshot["range_specs"]["rows"]) + 1)]
    fragments = []
    for fragment in snapshot["kept_fragments"]:
        base = bases[fragment["number"]]
        fragment_energies = {}
        for range_id in protein_ranges:
            energy = energies.get((fragment["key"], range_id))
            fragment_energies[range_id] = (
                {"energy": energy, "error": None} if energy is not None
                else {"energy": None, "error": "fixture scoring error"}
            )
        fragments.append({
            "fragment": fragment["key"], "rmsd": fragment["number"] / 10,
            **{name: base[name] for name in ("nuc1_ovrmsd_raw", "nuc1_ovrmsd",
                                             "nuc2_ovrmsd_raw", "nuc2_ovrmsd")},
            "energies": fragment_energies,
        })
    pairs = [
        {"down": pair["down"], "up": pair["up"], "sequence": pair["sequence"],
         "crmsd_raw": pair["crmsd_raw"], "crmsd_computed": pair["crmsd"], "crmsd": pair["crmsd"],
         "ovrmsd_raw": pair["ovrmsd_raw"], "ovrmsd_computed": pair["ovrmsd"], "ovrmsd": pair["ovrmsd"]}
        for pair in snapshot["kept_pairs"]
    ]
    return {"schema": 1, "fragments": fragments, "pairs": pairs}


@pytest.fixture
def tree_inputs():
    snapshot = _complete_tree_snapshot()
    recipes = export_gens.plan_gens(snapshot)
    digest = export_snapshot.snapshot_digest(snapshot)
    return snapshot, recipes, digest


def _render(snapshot, recipes, digest):
    return export_tree.render_tree(snapshot, recipes, digest)


def _json_member(tree, path):
    value = tree[path]
    assert isinstance(value, bytes)
    assert value.endswith(b"\n")
    return json.loads(value.decode("utf-8"))


def _expected_gen_inventory(recipes):
    return [
        {
            "path": recipe["path"],
            "modes": recipe["modes"],
            "range": recipe["range"],
            "direction": recipe["direction"],
            "final_pools": recipe["final_pools"],
            "external_inputs": recipe["external_inputs"],
            "actions": [action["name"] for action in recipe["actions"]],
        }
        for recipe in recipes
    ]


def _expected_data_paths(snapshot):
    paths = {
        "DATA/reference.pdb",
        "DATA/pdbcode.txt",
        "DATA/constraints.json",
        "DATA/energy-thresholds.json",
        "DATA/score-filter.yaml",
        "DATA/export-manifest.json",
        "DATA/best-fit-statistics.json",
        "DATA/anchor-refe-bases.pdb",
    }
    paths.update(f"DATA/{name}" for name in snapshot["range_pdb_files"])
    paths.update(f"DATA/{name}" for name in snapshot["range_pdb_files_aa"])
    for fragment in snapshot["kept_fragments"]:
        number, sequence = fragment["number"], fragment["sequence"]
        paths.add(f"DATA/clipped-fragments/reference/frag-{number}-{sequence}-refe.pdb")
        paths.add(f"DATA/clipped-fragments/bestfit/frag-{number}-{sequence}-bestfit.pdb")
    return paths


def test_render_tree_has_exact_members_and_preserves_prepared_bytes(tree_inputs):
    snapshot, recipes, digest = tree_inputs
    tree = _render(snapshot, recipes, digest)

    expected = _expected_data_paths(snapshot)
    expected.update(
        f"{recipe['path']}/{action['name']}/alaric.yaml"
        for recipe in recipes for action in recipe["actions"]
    )
    assert set(tree) == expected
    assert all(isinstance(path, str) and isinstance(value, bytes)
               for path, value in tree.items())
    assert all(not path.startswith("alaric_project/") for path in tree)

    assert tree["DATA/reference.pdb"] == snapshot["reference_pdb"]
    assert tree["DATA/pdbcode.txt"] == (snapshot["exclude_pdb_code"] + "\n").encode("utf-8")
    for name, value in snapshot["range_pdb_files"].items():
        assert tree[f"DATA/{name}"] == value
    for name, value in snapshot["range_pdb_files_aa"].items():
        assert tree[f"DATA/{name}"] == value

    kept_numbers = {row["number"] for row in snapshot["kept_fragments"]}
    for fragment in snapshot["kept_fragments"]:
        number, sequence = fragment["number"], fragment["sequence"]
        ref_source = f"reference/frag-{number}-{sequence}-refe.pdb"
        fit_source = f"bestfit/frag-{number}-{sequence}-bestfit.pdb"
        assert tree[f"DATA/clipped-fragments/{ref_source}"] == snapshot["clipped_fragments"][ref_source]
        assert tree[f"DATA/clipped-fragments/{fit_source}"] == snapshot["clipped_fragments"][fit_source]
    assert all(f"/frag-{number}-" not in path
               for path in tree if path.startswith("DATA/clipped-fragments/")
               for number in set(range(1, 12)) - kept_numbers)
    assert "DATA/clipped-fragments/clip.json" not in tree
    assert "DATA/clipped-fragments/thresholds.json" not in tree


def test_best_fit_statistics_member_is_the_normalised_kept_table(tree_inputs):
    snapshot, recipes, digest = tree_inputs
    tree = _render(snapshot, recipes, digest)
    statistics = _json_member(tree, "DATA/best-fit-statistics.json")
    assert statistics == snapshot["best_fit_statistics"]
    assert [row["fragment"] for row in statistics["fragments"]] \
        == [row["key"] for row in snapshot["kept_fragments"]]
    assert [(row["down"], row["up"]) for row in statistics["pairs"]] \
        == [(row["down"], row["up"]) for row in snapshot["kept_pairs"]]

    stale = deepcopy(snapshot)
    stale["best_fit_statistics"]["fragments"].pop()
    with pytest.raises(ValueError, match="kept fragments"):
        _render(stale, recipes, export_snapshot.snapshot_digest(stale))


def test_anchor_refe_bases_pdb_holds_verbatim_reference_bases_of_every_candidate(tree_inputs):
    snapshot, recipes, digest = tree_inputs
    tree = _render(snapshot, recipes, digest)
    text = tree["DATA/anchor-refe-bases.pdb"].decode("utf-8")
    assert text.endswith("END\n") and "\r" not in text

    candidates = {}
    for recipe in recipes:
        for action in recipe["actions"]:
            parameters = action["parameters"]
            if parameters["action"] == "anchor-refe":
                candidates[(parameters["fragment"], parameters["nucleotide"])] = parameters["ovrmsd"]
    assert candidates
    positions = sorted({number + (nucleotide == "second") for number, nucleotide in candidates})

    residues = _rna_residues(snapshot["fragments"])
    reference_lines = snapshot["reference_pdb"].decode("ascii").split("\r\n")
    expected_atoms = []
    for position in positions:
        residue = residues[position - 1]
        for name in export_snapshot.ANCHOR_REFE_BASE_ATOMS[residue["base"]]:
            expected_atoms.extend(
                line for line in reference_lines
                if line.startswith("ATOM") and int(line[22:26]) == residue["resid"]
                and line[12:16].strip() == name
            )
    assert [line for line in text.splitlines() if line.startswith("ATOM")] == expected_atoms

    remarks = [line for line in text.splitlines() if line.startswith("REMARK   2 ")]
    assert [int(line.split()[4]) for line in remarks] == positions
    keys = {row["number"]: row["key"] for row in snapshot["fragments"]}
    for (number, nucleotide), threshold in candidates.items():
        position = number + (nucleotide == "second")
        remark = remarks[positions.index(position)]
        assert f"{keys[number]} {nucleotide} ovrmsd {threshold}" in remark


def test_anchor_refe_bases_pdb_rejects_a_reference_missing_a_base_atom(tree_inputs):
    snapshot, recipes, _digest = tree_inputs
    changed = deepcopy(snapshot)
    reference = changed["reference_pdb"].decode("ascii").split("\r\n")
    position = next(
        action["parameters"]["fragment"]
        for recipe in recipes for action in recipe["actions"]
        if action["parameters"]["action"] == "anchor-refe"
        and action["parameters"]["nucleotide"] == "first"
    )
    resid = _rna_residues(changed["fragments"])[position - 1]["resid"]
    reference = [line for line in reference
                 if not (line.startswith("ATOM") and int(line[22:26]) == resid
                         and line[12:16].strip() == "C5")]
    changed["reference_pdb"] = "\r\n".join(reference).encode("ascii")
    changed["range_pdb_files_aa"]["range-1-aa.pdb"] = changed["reference_pdb"]
    changed["reference_sha256"] = hashlib.sha256(changed["reference_pdb"]).hexdigest()
    with pytest.raises(ValueError, match="lacks base atom C5"):
        _render(changed, recipes, export_snapshot.snapshot_digest(changed))


def test_constraints_energy_and_score_filter_members_follow_schema(tree_inputs):
    snapshot, recipes, digest = tree_inputs
    tree = _render(snapshot, recipes, digest)

    constraints = _json_member(tree, "DATA/constraints.json")
    expected_constraints = {"schema": 1}
    for fragment in snapshot["kept_fragments"]:
        expected_constraints[f"frag{fragment['number']}"] = {"sequence": fragment["sequence"]}
    expected_constraints["pairs"] = [
        {
            "down": pair["down"],
            "up": pair["up"],
            "sequence": pair["sequence"],
            "crmsd": pair["crmsd"],
            "ovrmsd": pair["ovrmsd"],
        }
        for pair in snapshot["kept_pairs"]
    ]
    assert constraints == expected_constraints
    assert "frag1" not in constraints and "frag2" not in constraints
    assert set(constraints) == {"schema", "pairs"} | {
        f"frag{fragment['number']}" for fragment in snapshot["kept_fragments"]
    }

    energy = _json_member(tree, "DATA/energy-thresholds.json")
    expected_energy = {"schema": 1}
    for cell in snapshot["energy_cells"]:
        expected_energy.setdefault(cell["fragment"], {})[cell["range"]] = cell["threshold"]
    assert energy == expected_energy
    assert all(isinstance(value, str) for fragment in energy.values() if isinstance(fragment, dict)
               for value in fragment.values())

    score_filter = yaml.safe_load(tree["DATA/score-filter.yaml"].decode("utf-8"))
    assert score_filter == {
        "min_cascade": 1000,
        "cascade": [
            {"plateau_distance": 5, "neighbour_distance": 7, "max_error": 1.0},
            {"plateau_distance": 10, "neighbour_distance": 12, "max_error": 0.05},
        ],
    }
    assert tree["DATA/score-filter.yaml"].endswith(b"\n")


def test_manifest_is_compact_ordered_and_keeps_source_qualified_producers(tree_inputs):
    snapshot, recipes, digest = tree_inputs
    tree = _render(snapshot, recipes, digest)
    manifest = _json_member(tree, "DATA/export-manifest.json")

    assert manifest["schema"] == 1
    assert manifest["snapshot_sha256"] == digest
    assert manifest["reference_sha256"] == snapshot["reference_sha256"]
    assert manifest["range_specs"] == snapshot["range_specs"]
    assert manifest["handoff"] == snapshot["handoff"]
    assert manifest["class_table"] == snapshot["class_table"]
    assert manifest["gens"] == _expected_gen_inventory(recipes)
    assert manifest["disabled_reference_anchors"] == []

    # The manifest carries the compact recipe inventory, not action mappings.
    assert all(set(gen) == {
        "path", "modes", "range", "direction", "final_pools", "external_inputs", "actions"
    } for gen in manifest["gens"])
    assert all(isinstance(action, str) for gen in manifest["gens"] for action in gen["actions"])
    assert all("parameters" not in gen for gen in manifest["gens"])

    # Source-qualified R12 merge inputs remain tied to their matching anchor
    # tree in the manifest; the adjacent grown-side input keeps all producers.
    for source in ("anchor", "anchor-refe"):
        gen = next(item for item in recipes
                   if item["path"] == f"merge/fwd/{source}/gen-frag10-merged")
        record = next(item for item in manifest["gens"] if item["path"] == gen["path"])
        target = next(row for row in record["external_inputs"] if row["pool"] == "frag10")
        assert target["producers"]
        assert target["producers"] == [f"{source}/gen-frag10"]
        grown = next(row for row in record["external_inputs"] if row["pool"] == "frag9")
        assert set(grown["producers"]) == {
            path
            for item in recipes
            if "frag9" in item["final_pools"]
            for path in [item["path"]]
        }


def test_renderer_output_is_deterministic_and_text_has_final_newline(tree_inputs):
    snapshot, recipes, digest = tree_inputs
    first = _render(snapshot, recipes, digest)
    second = _render(snapshot, recipes, digest)
    assert first == second

    for path, value in first.items():
        if path.endswith((".json", ".yaml", ".txt", ".pdb")) and not path.endswith(".pdb"):
            value.decode("utf-8")
            assert value.endswith(b"\n"), path


def test_sentinel_comments_are_real_yaml_comments_and_manifested(tree_inputs):
    snapshot, _, _ = tree_inputs
    snapshot = deepcopy(snapshot)
    base = next(row for row in snapshot["threshold_bases"] if row["fragment"] == 3)
    base["nuc2_ovrmsd"] = "1.80"
    base["nuc2_ovrmsd_raw"] = 1.8
    recipes = export_gens.plan_gens(snapshot)
    digest = export_snapshot.snapshot_digest(snapshot)
    tree = _render(snapshot, recipes, digest)
    manifest = _json_member(tree, "DATA/export-manifest.json")

    sentinels = []
    for recipe in recipes:
        for action in recipe["actions"]:
            params = action["parameters"]
            if params["action"] == "anchor-refe" and params["fragment"] == 3 \
                    and params["nucleotide"] == "second":
                sentinels.append((recipe, action))
    assert sentinels
    listed_paths = set()
    for recipe, action in sentinels:
        path = f"{recipe['path']}/{action['name']}/alaric.yaml"
        listed_paths.add(f"{recipe['path']}/{action['name']}/alaric.yaml")
        emitted = tree[path].decode("utf-8")
        comment_lines = [line for line in emitted.splitlines() if line.startswith("#")]
        first_mapping_line = next(i for i, line in enumerate(emitted.splitlines())
                                  if line and not line.startswith("#"))
        assert comment_lines
        assert all(i < first_mapping_line and line.startswith("#")
                   for i, line in enumerate(emitted.splitlines()) if line.startswith("#"))
        assert "comment:" not in emitted
        assert "disabled" in " ".join(comment_lines).lower()
        parsed = yaml.safe_load(emitted)
        assert parsed["ovrmsd"] == -1

    disabled = manifest["disabled_reference_anchors"]
    record = next(row for row in disabled
                  if row["candidate"] == {"fragment": "frag-3-AG", "nucleotide": "second"})
    assert record["raw"] == 1.8
    assert record["rounded"] == "1.80"
    assert record["limit"] == pytest.approx(3 ** 0.5)
    assert set(record["action_paths"]) == listed_paths


def test_render_tree_rejects_a_digest_that_does_not_match_snapshot(tree_inputs):
    snapshot, recipes, digest = tree_inputs
    wrong = "0" * 64 if digest != "0" * 64 else "1" * 64
    with pytest.raises(ValueError):
        _render(snapshot, recipes, wrong)


@pytest.mark.parametrize("bad_path", [
    "../outside/gen-frag3", "/absolute/gen-frag3", "grow//fwd/gen-frag3", "grow\\fwd\\gen-frag3",
])
def test_render_tree_rejects_unsafe_generated_paths(tree_inputs, bad_path):
    snapshot, recipes, digest = tree_inputs
    changed = deepcopy(recipes)
    changed[0]["path"] = bad_path
    with pytest.raises(ValueError):
        _render(snapshot, changed, digest)


@pytest.mark.parametrize(("source_path", "bad_path"), [
    # A source-tree name is only a top-level anchor tree, or the optional
    # merge source immediately before its terminal gen folder.
    ("grow/fwd/fixed/gen-frag7", "grow/fwd/anchor/gen-frag7"),
    ("anchor/fwd/fixed/gen-frag3-6", "anchor/fwd/fixed/anchor/gen-frag3-6"),
    ("merge/fwd/anchor/gen-frag10-merged", "merge/fwd/anchor/fixed/gen-frag10-merged"),
    # Grow and merge require direction before their mode and range levels.
    ("grow/fwd/fixed/gen-frag7", "grow/fixed/fwd/gen-frag7"),
    ("grow/fwd/fixed/gen-frag7", "grow/fixed/gen-frag7"),
    ("merge/fwd/anchor/gen-frag10-merged", "merge/fixed/fwd/anchor/gen-frag10-merged"),
    ("merge/fwd/anchor/gen-frag10-merged", "merge/fixed/anchor/gen-frag10-merged"),
    # A range is valid only after `unfixed`; it cannot be attached to fixed.
    ("anchor/fwd/unfixed/range-2/gen-frag5-6",
     "anchor/fwd/range-2/unfixed/gen-frag5-6"),
    ("anchor/fwd/unfixed/range-2/gen-frag5-6", "anchor/fwd/range-2/gen-frag5-6"),
    ("anchor/fwd/fixed/gen-frag3-6", "anchor/fwd/fixed/range-2/gen-frag3-6"),
    ("merge/fwd/anchor/gen-frag10-merged", "merge/fwd/fixed/range-2/anchor/gen-frag10-merged"),
    # Common-conformer has no direction, mode, or range classifier levels.
    ("common-conformer/gen-frag5-range-2-3", "common-conformer/fwd/gen-frag5-range-2-3"),
])
def test_render_tree_rejects_misclassified_generator_paths(tree_inputs, source_path, bad_path):
    snapshot, recipes, digest = tree_inputs
    changed = deepcopy(recipes)
    recipe = next(item for item in changed if item["path"] == source_path)
    assert bad_path not in {item["path"] for item in recipes}
    recipe["path"] = bad_path

    with pytest.raises(ValueError):
        _render(snapshot, changed, digest)


def test_render_tree_rejects_multiple_merge_source_qualifiers(tree_inputs):
    snapshot, recipes, digest = tree_inputs

    # The single-stacking nuc11 anchor yields valid source-qualified merges.
    changed = deepcopy(recipes)
    recipe = next(item for item in changed
                  if item["path"] == "merge/fwd/anchor/gen-frag10-merged")
    bad_path = "merge/fwd/anchor/anchor-refe/gen-frag10-merged"
    assert bad_path not in {item["path"] for item in recipes}
    recipe["path"] = bad_path

    with pytest.raises(ValueError):
        _render(snapshot, changed, digest)


def test_render_tree_rejects_duplicate_generated_members(tree_inputs):
    snapshot, recipes, digest = tree_inputs
    changed = deepcopy(recipes)
    changed.append(deepcopy(changed[0]))
    with pytest.raises(ValueError):
        _render(snapshot, changed, digest)


def test_render_tree_rejects_prefixing_generator_paths(tree_inputs):
    snapshot, recipes, digest = tree_inputs
    changed = deepcopy(recipes[:2])
    changed[0]["path"] = "anchor/gen-frag3"
    changed[1]["path"] = "anchor/gen-frag3/child/gen-frag99"
    with pytest.raises(ValueError, match="prefix"):
        _render(snapshot, changed, digest)


def test_render_tree_rejects_generated_member_file_directory_collision(tree_inputs, monkeypatch):
    snapshot, recipes, digest = tree_inputs
    changed = deepcopy(recipes[:2])
    first_action_name = changed[0]["actions"][0]["name"]
    changed[0]["path"] = "anchor/gen-frag3"
    changed[1]["path"] = f"anchor/gen-frag3/{first_action_name}/alaric.yaml/gen-frag99"

    # Each recipe is independently well-formed and its action names are safe.
    # The catalogue path hierarchy normally catches this earlier, so bypass
    # only that uniqueness pass to exercise the final generated-key check.
    for recipe in changed:
        export_gens.validate_recipe(recipe)
    monkeypatch.setattr(export_tree, "_validate_recipes", lambda recipes: list(recipes))
    with pytest.raises(ValueError, match="directory prefix"):
        _render(snapshot, changed, digest)


def test_emitted_action_yamls_are_accepted_by_alaric_schema(tree_inputs):
    schema = pytest.importorskip("alaric.middle.schema")
    snapshot, recipes, digest = tree_inputs
    tree = _render(snapshot, recipes, digest)
    from pathlib import PurePosixPath

    for path, value in tree.items():
        if not path.endswith("/alaric.yaml"):
            continue
        action_name = PurePosixPath(path).parent.name
        mapping = yaml.safe_load(value.decode("utf-8"))
        assert isinstance(mapping, dict), path
        schema.normalize_action(action_name, Path(path), mapping)


def test_emitted_constraints_are_accepted_by_alaric_resolver(tree_inputs, tmp_path):
    resolve = pytest.importorskip("alaric.middle.resolve")
    Project = pytest.importorskip("alaric.middle.project").Project
    snapshot, recipes, digest = tree_inputs
    tree = _render(snapshot, recipes, digest)

    data_dir = tmp_path / "DATA"
    data_dir.mkdir()
    constraints_path = data_dir / "constraints.json"
    constraints_path.write_bytes(tree["DATA/constraints.json"])
    resolved = resolve._constraints(Project(tmp_path))
    emitted = _json_member(tree, "DATA/constraints.json")

    expected = {
        "schema": emitted["schema"],
        "pairs": [
            {
                "down": pair["down"],
                "up": pair["up"],
                "sequence": pair["sequence"].upper(),
                "crmsd": float(pair["crmsd"]),
                "ovrmsd": float(pair["ovrmsd"]),
            }
            for pair in emitted["pairs"]
        ],
    }
    expected.update({
        key: {"sequence": value["sequence"].upper()}
        for key, value in emitted.items()
        if key.startswith("frag")
    })
    assert resolved == expected

