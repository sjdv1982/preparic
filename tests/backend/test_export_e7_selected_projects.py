"""E7: extract real renderer output, choose gens and use Alaric consumers.

Run in the Alaric env with both checkouts on PYTHONPATH. These are graph and
deployment tests; executing scientific kernels is a separate acceptance gate.
"""

import json
import re

import pytest

import export_gens
import export_snapshot
import export_tree
from pathlib import Path

from test_export_tree import _complete_tree_snapshot


graph_module = pytest.importorskip("alaric.middle.graph")
Project = pytest.importorskip("alaric.middle.project").Project
deploy = pytest.importorskip("alaric.middle.deploy").deploy
resolve = pytest.importorskip("alaric.middle.resolve")


def _write_tree(root, tree):
    for member, value in tree.items():
        target = root / member
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(value)


def _select(root, tree, recipes, paths):
    selected = [next(recipe for recipe in recipes if recipe["path"] == path) for path in paths]
    writes = {}
    for recipe in selected:
        for action in recipe["actions"]:
            name = action["name"]
            target = root / name / "alaric.yaml"
            if name in writes or target.exists():
                raise FileExistsError(f"selected gen alternatives collide at {name}")
            writes[name] = tree[f"{recipe['path']}/{name}/alaric.yaml"]
    for name, value in writes.items():
        (root / name).mkdir()
        (root / name / "alaric.yaml").write_bytes(value)
    return selected


def _catalogue():
    snapshot = _complete_tree_snapshot()
    # The normalized planner fixture uses contact IDs 121 below the original
    # Alaric domain PDBs. Supply committed 1B7F coordinates at those fixture
    # addresses, replacing the planner's placeholder protein bytes, so the
    # real Gaussian resolver can look up TYR/PHE and bind the ring identities.
    alaric_root = Path(resolve.__file__).resolve().parents[2]
    for number, source in ((2, "1b7f_dom1.pdb"), (3, "1b7f_dom2.pdb")):
        lines = []
        for line in (alaric_root / "tests/1b7f/pdbs" / source).read_text().splitlines():
            if line.startswith(("ATOM  ", "HETATM")):
                line = line[:22] + f"{int(line[22:26]) - 121:4d}" + line[26:]
                lines.append(line)
        snapshot["range_pdb_files"][f"range-{number}.pdb"] = (
            "\n".join(lines) + "\nEND\n"
        ).encode("ascii")
    # The planner oracle uses integer-valued decimal text for readability.
    # Real S24/export threshold values use exactly two fractional digits.
    # Pad those fixture integers without changing their numerical values.
    for cell in snapshot["energy_cells"]:
        whole, separator, fraction = cell["threshold"].partition(".")
        assert len(fraction) <= 2
        cell["threshold"] = whole + "." + fraction.ljust(2, "0")
    recipes = export_gens.plan_gens(snapshot)
    tree = export_tree.render_tree(snapshot, recipes, export_snapshot.snapshot_digest(snapshot))
    return snapshot, recipes, tree


FIXED_ANCHORS = ["anchor/fwd/fixed/gen-frag3-6", "anchor-refe/gen-frag10"]
CASES = [
    pytest.param(FIXED_ANCHORS + [
        "grow/fwd/fixed/gen-frag7", "grow/fwd/fixed/gen-frag8",
        "grow/fwd/fixed/gen-frag9", "merge/fwd/gen-frag10",
    ], "fixed", "frag10-merged", id="fixed-anchor-target-r12"),
    pytest.param(FIXED_ANCHORS + [
        "grow/fwd/fixed/gen-frag7", "grow/bwd/fixed/gen-frag8",
        "grow/bwd/gen-frag9", "merge/fwd/fixed/gen-frag8",
    ], "fixed", "frag8-merged", id="fixed-grown-target-restrict"),
    pytest.param([
        "anchor-refe/fwd/unfixed/range-3/gen-frag3-5",
        "anchor-refe/fwd/unfixed/range-2/gen-frag5-6", "anchor-refe/gen-frag10",
        "grow/fwd/unfixed/range-2/gen-frag7", "grow/fwd/unfixed/range-2/gen-frag8",
        "grow/fwd/unfixed/range-2/gen-frag9", "grow/fwd/unfixed/range-3/gen-frag6",
        "grow/fwd/unfixed/range-3/gen-frag7", "grow/fwd/unfixed/range-3/gen-frag8",
        "merge/fwd/gen-frag10",
        *[f"common-conformer/gen-frag{number}-range-2-3" for number in range(5, 9)],
    ], "unfixed", "frag10-merged", id="unfixed-two-frames-common-conformer"),
]


@pytest.mark.parametrize("paths,mode,target", CASES)
def test_e7_selected_catalogues_resolve_and_deploy(tmp_path, monkeypatch, paths, mode, target):
    monkeypatch.delenv("ALARIC_ANCHOR_MIN_CELL_FRACTION", raising=False)
    snapshot, recipes, tree = _catalogue()
    root = tmp_path / "alaric_project"
    _write_tree(root, tree)
    project = Project(root)
    # Extraction alone is a catalogue, with no executable action at its root.
    assert project.action_dirs() == []
    assert graph_module.ActionGraph(project).build() == {}
    selected = _select(root, tree, recipes, paths)
    resolved = graph_module.ActionGraph(project).build()
    assert target in resolved
    assert len(resolved) == sum(len(recipe["actions"]) for recipe in selected)
    for recipe in selected:
        for action in recipe["actions"]:
            name = action["name"]
            current = resolved[name]
            assert (root / name / "alaric.yaml").read_bytes() == tree[
                f"{recipe['path']}/{name}/alaric.yaml"]
            if current.action in {"anchor", "anchor-refe", "grow", "energy-filter"}:
                assert isinstance(current.params["fragment"], int)
                assert current.params["sequence"] != "auto"
                assert current.params["exclude"] == ["1b7f"]
            if current.action == "anchor":
                assert current.params["precomputed_filter_name"] == "gaussian-rrm"
                assert current.params["precomputed_filter_identity"]["spec_digest"]
                assert current.params["anchor_reference_identity"]["nucleotide_sha256"]
                assert current.params["anchor_reference_identity"]["protein_sha256"]
            if current.action == "grow":
                assert current.params["direction"] in {"forward", "backward"}
                assert isinstance(current.params["crmsd"], float)
                assert isinstance(current.params["ovrmsd"], float)
            if current.action == "energy-filter":
                assert isinstance(current.params["threshold"], float)
                assert current.params["epsilon"] == 15
                assert current.params["min_cascade"] == 1000
                if mode == "unfixed" and recipe["range"] is not None:
                    assert current.params["protein"] == recipe["range"]
            deploy("local", root / name)
            script = (root / name / "run.sh").read_text()
            # Deployment's shell prologue deliberately preserves this literal
            # environment expression; action-template fields must be rendered.
            assert re.search(r"\{\{\s*[a-z_][a-z_0-9]*\s*\}\}", script) is None
            assert "ALARIC_ANCHOR_MIN_CELL_FRACTION" not in script
    merged = resolved[target]
    if target == "frag10-merged":
        assert merged.action == "identity"
        assert merged.params["input1"].action == "grow"
        assert merged.params["input2"].name == "frag10"
    else:
        assert merged.action == "grow" and merged.params["restrict_input"].name == "frag8"
    if mode == "fixed":
        assert {"anchor", "anchor-refe"} <= {action.action for action in resolved.values()}
        # A scored pool at a two-range fragment is filtered in tandem.
        filters = [action for action in resolved.values() if action.action == "energy-filter"]
        assert any(action.params["input"].action == "energy-filter"
                   and action.params["input"].params["protein"] != action.params["protein"]
                   for action in filters)
    else:
        applications = [action for action in resolved.values() if action.action == "apply-mask"]
        assert len(applications) == 8
        for action in applications:
            pair = action.params["mask_input"]
            assert action.params["input"].name == pair.params[f"input{action.params['mask']}"].name
            assert resolve.final_fragment(action) in range(5, 9)
    manifest = json.loads(tree["DATA/export-manifest.json"])
    assert manifest["snapshot_sha256"] == export_snapshot.snapshot_digest(snapshot)


def test_e7_opposite_direction_alternatives_cannot_silently_overwrite(tmp_path):
    _, recipes, tree = _catalogue()
    _write_tree(tmp_path, tree)
    _select(tmp_path, tree, recipes, ["anchor-refe/fwd/fixed/gen-frag3-6"])
    before = {path: path.read_bytes() for path in tmp_path.glob("*/alaric.yaml")}
    with pytest.raises(FileExistsError, match="collide"):
        _select(tmp_path, tree, recipes, ["anchor-refe/bwd/fixed/gen-frag3-6"])
    assert {path: path.read_bytes() for path in tmp_path.glob("*/alaric.yaml")} == before
