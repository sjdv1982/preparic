# Stage 3 — stacking interaction implementation plan

## Goal and fixed contracts

Stage 3 consumes the **pre-aareduce** `range_pdbs` output, detects protein–RNA
aromatic stacking interactions in Seamless, exposes the result and its selection
state as mounted/shared cells, and lets the user inspect and select interactions in
both a table and the Mol* viewer.

The implementation must preserve these contracts:

- `range_pdbs[0]` is RNA; every later entry is a protein range PDB.
- Only residues with every atom in the configured ring atom set participate.
- A candidate is emitted when centroid distance is `<= 6.0 Å` and the acute
  plane/plane tilt is `<= 45.0°`. Cutoffs are inclusive.
- Slip is the acute angle between the centroid displacement vector and the
  **protein ring plane** (`0°` in-plane, `90°` normal to the plane).
- Results are deterministic and sorted by `(protein_range_pdb_id,
  protein_resid, rna_resid)`.
- `protein_range_pdb_id` is the one-based ID of the original `range_pdbs`
  element/file. Consequently the first possible protein ID is `2`.
- A stacking row has the JSON shape below. Centroids are retained as rendering
  metadata so the viewer draws the exact geometry measured by the detector:

  ```json
  {
    "rna_resname": "G",
    "rna_resid": 4,
    "protein_range_pdb_id": 3,
    "protein_resid": 93,
    "distance": 4.123,
    "tilt_angle": 12.345,
    "slip_angle": 71.234,
    "rna_centroid": [1.0, 2.0, 3.0],
    "protein_centroid": [4.0, 5.0, 6.0]
  }
  ```

- Row identity for preserving a selection is the participant tuple
  `(rna_resname, rna_resid, protein_range_pdb_id, protein_resid)`. Measurement
  changes update the row but do not make the same interaction “new”.
- Newly seen rows default to selected exactly when `distance <= 5.0` and
  `tilt_angle <= 30.0`.
- `stackings` is a read-only shared `plain` cell mounted at
  `deployment-mounts/stackings.json`. `stacking_selection` is a writable shared
  `plain` cell mounted read/write at
  `deployment-mounts/stacking-selection.json` and contains only `list[bool]`.
- Any non-initial range edit while the app is in stage 2 clears the app-held
  stacking rows and selection. The old backend result must not immediately
  repopulate the cleared UI; stage 3 waits for the next backend stacking revision.

For all Python/Seamless commands, run in the `preparic` conda environment.
For source exploration, use the activated Serena project and prefer symbol/
reference queries. Before changing Seamless wrapping or graph wiring, consult
`~/seamless1/seamless/docs/agent/contracts/direct-delayed-and-transformation.md`,
`modules-and-closures.md`, and `identity-and-caching.md`. Re-check the contract
after wiring the detector into the graph.

## Phase 1 — pure detector and isolated backend tests

**Active Seamless concern:** DAG/plain step behavior. The detector must be a
referentially transparent ordinary Python step before it is placed in the graph.

1. Add `detect_stacking.py` with a top-level `_detect_stacking(range_pdbs)`
   transformer function and small private helpers for canonical residue names,
   complete-ring extraction, SVD plane fitting, and stable angle calculation.
2. Reuse and adapt ring definitions from
   `~/alaric-benchmark/detect-stacking.py`. Do not retain its old normal/lateral
   separation filters: stage 3 uses only the requested distance and tilt filters.
3. Validate the boundary value as `list[str]`; return `[]` for fewer than two
   inputs or empty PDB text, and raise clear `ValueError`s for malformed boundary
   types/PDBs. Read only the first model, as the range extractor currently emits
   one selected structure.
4. Clamp dot products before `acos`/`asin`, use `abs(dot)` for unoriented plane
   normals, skip zero-length centroid displacement safely, and emit ordinary
   Python numbers/lists (never numpy values or class instances) across the
   Seamless boundary.
5. Add `tests/backend/test_detect_stacking.py` in the same isolated style as
   `test_range_pdbs.py`. Cover a positive synthetic pair with independently
   checked distance/tilt/slip; inclusive and just-outside distance/tilt cutoffs;
   absent ring atoms; non-ring residues; RNA aliases; multiple protein files and
   one-based IDs; deterministic ordering; empty/invalid input; and the repository
   `1b7f` fixture (RNA range 1 against protein ranges 2/3).

**Exit check:**

```bash
conda run -n preparic pytest -q tests/backend/test_detect_stacking.py
```

The unwrapped function passes all isolated tests and returns only JSON/plain-cell
compatible data.

## Phase 2 — Seamless graph, mount/share contract, and integration tests

**Active Seamless concern:** minimal wrapping followed immediately by Seamless
contract validation.

### Phase 2A — isolated Seamless transformation

1. Import `_detect_stacking` into `backend.py` rather than inlining it.
2. First build `_detect_stacking` as an isolated Seamless transformation in the
   style already used for `_ranges_to_pdb` in `tests/backend/test_seamless.py`.
   Prove successful, empty, validation-error, and deterministic execution without
   constructing `build_context()` or any application workflow. Keep
   `test_seamless.py` limited to isolated transformation behavior.

### Phase 2B — application workflow wiring

3. Add `context.detect_stacking = _detect_stacking`, wire its sole input pin to
   `context.range_pdbs` (not `range_pdb_files_aa`), and expose its result through
   the plain cell `context.stackings`.
4. Mount `stackings` as a cell-authoritative write mount at the agreed JSON path
   and share it read-only.
5. Add an `_initial_stacking_selection` loader which accepts only a JSON list of
   booleans and otherwise falls back safely to `[]`. Create
   `context.stacking_selection`, initialize it, mount it read/write with file
   authority, and share it writable.
6. Extend state-graph error reporting later consumed by the browser so detector
   failures are distinguishable from range extraction failures.
7. Put all workflow assertions in a separate integration module such as
   `tests/backend/test_backend_workflow.py`; do not add `build_context()` workflow
   tests to `test_seamless.py`. Prove the graph dependency is pre-aareduce,
   `stackings` computes the same value as the isolated function, the empty
   pipeline is stable, both mounts have the intended celltypes/modes, and
   selection initialization rejects non-boolean or non-list persisted JSON.
8. Re-check explicit inputs, JSON-compatible outputs, imported-module capture,
   deterministic result identity, and persistent cell-backed caching against the
   Seamless contract documents. No remote execution, scratch/fingertip, or HPC
   behavior belongs in this UI graph phase.

**Exit check:**

```bash
conda run -n preparic pytest -q tests/backend/test_detect_stacking.py tests/backend/test_seamless.py tests/backend/test_backend_workflow.py
```

The small wrapped detector and the extraction-to-detection graph both pass, and
mounted output bytes are valid canonical JSON.

## Phase 3 — stage-3 domain state and read-only table

**Active Seamless concern:** scaling the verified result into the application data
model without selection synchronization yet.

1. Extend the frontend types with `Stage = 1 | 2 | 3`, `Stacking`, and the small
   viewer input types needed later. Add runtime parsing/normalization helpers in a
   dedicated `src/lib/stackings.ts`; do not trust shared JSON merely because the
   TypeScript type is correct.
2. Extend app state with stacking rows and a parallel boolean selection array.
   Provide explicit helpers to enter stage 3, return to stage 2, clear stacking
   state, compute stable row keys, merge a new table while preserving selections
   for existing keys, and auto-select genuinely new rows using the `5 Å / 30°`
   rule.
3. Add `Stage3Stackings.svelte` and `StackingTable.svelte`. Show all requested
   statistics in a read-only, accessible table; format numbers for display without
   changing stored values. In this phase the selection column may be displayed
   from local state, but backend synchronization is deferred to phase 4.
4. Make stage 2 Complete enter stage 3, add the third progress indicator, and add
   stage-3 Back navigation. Stage 1 and 2 rendering/picking behavior must remain
   unchanged.
5. Detect range-row changes while stage 2 is active. Treat the first assignment on
   entering stage 2 as initialization; every subsequent content change (local or
   shared) clears rows/selection and records the current stacking revision as
   stale so the prior backend table cannot flash back into the UI.
6. Add focused Vitest coverage for schema rejection, stable identities, merge and
   auto-selection, reordering, statistic-only updates, stage transitions, and the
   initialization-versus-edit reset rule. Keep complex state logic in plain `.ts`
   helpers so it can be tested without mounting Svelte components.

**Exit check:**

```bash
npm test -- --run
npm run check
npm run build
```

Stage 3 renders an accessible read-only stacking table from normalized data, new
rows receive deterministic local defaults, and range edits invalidate it.

## Phase 4 — bidirectional Seamless selection synchronization

**Active Seamless concern:** explicit browser/backend dataflow and feedback-loop
control for two writable/readable shared values.

1. Extend the browser share interfaces and share-list attachment logic for
   read-only `stackings` and writable `stacking_selection`. Parse plain-cell JSON
   defensively, expose independent table/selection revisions and validation
   errors, and extend detector errors from the state graph.
2. Use content signatures and pending-write signatures to distinguish echoes from
   genuine backend updates. Never compare arrays by object identity.
3. On a backend stacking-table revision, merge by stable row key, preserving old
   selections and auto-assigning only new rows. Honor the stage-2 stale-revision
   barrier from phase 3.
4. On a backend selection revision:
   - if it is `list[bool]` with exactly the current table length, apply it;
   - if its only defect is length mismatch, retain it as the pending backend
     selection and try it against the next **non-empty backend-modified** stacking
     table;
   - for any other defect, reject it with a visible/non-fatal sync error and do not
     overwrite valid app state.
5. When a new table is merged and no valid backend selection supersedes it, push
   the post-auto-assignment list immediately. When a checkbox changes, push the
   complete updated boolean list immediately. Do not send writes when state was
   merely applying a valid backend value, and do not create an echo loop.
6. Add a Select checkbox per table row, with an accessible label identifying both
   residues. Add Vitest tests using a fake Seamless connector for exact-length
   remote application, queued length mismatch, next-empty-table behavior,
   invalid-type rejection, auto-assignment push, user push, write rejection, and
   echo suppression.

**Exit check:**

```bash
npm test -- --run
npm run check
npm run build
```

The table and `stacking_selection` converge in both directions under normal,
out-of-order, invalid, and delayed-table updates without a feedback loop.

## Phase 5 — Mol* stacking rendering and end-to-end browser acceptance

**Active Seamless concern:** application realism and final contract regression
after all dataflow consumers are present.

1. Extend `RenderSpec` and the stage-3 render effect so only selected stacking
   interactions are passed to Mol*. Stage 3 keeps the full structure cartoon for
   context, but replaces stage-2 sidechain-stick behavior with exactly the RNA and
   protein residues participating in selected rows.
2. Resolve a protein range PDB ID through the corresponding stage-2 range row and
   chain metadata. Build residue expressions for both participants, de-duplicate
   overlapping participants, and render one `ball-and-stick` representation with
   intra-residue bonds only. No unselected residue may receive sticks.
3. Add a Mol* tagged shape/representation for centroid connectors. Render each
   connector from the detector-provided centroid pair as a dashed line (short
   cylinder segments with gaps is acceptable), use a clear contrasting color,
   and tag all created state nodes with the existing render cleanup tag so updates
   and stage navigation cannot leak representations.
4. Serialize viewer operations through the existing render token/operation queue,
   preserve stage-1/stage-2 rendering, and ensure empty selection removes all
   stacking sticks/connectors.
5. Add unit tests for participant-to-expression mapping, de-duplication, and dash
   segment generation. Run all backend/frontend tests, type checks, and the
   production build.
6. Start `backend.py` and Vite in the `preparic` environment, then use the Chrome
   DevTools MCP server for acceptance testing. Exercise a known `1b7f` flow into
   stage 3; verify the table values and default checkbox; toggle a row and inspect
   DOM plus Mol* state/visual changes; mutate `stacking_selection` from the backend
   and verify the checkbox/viewer update; exercise the queued length-mismatch case;
   go back to stage 2, edit a range, and verify immediate table/selection clearing.
   Inspect console errors and failed network requests. Capture a screenshot only
   as supporting evidence after semantic assertions pass.

**Exit checks:**

```bash
conda run -n preparic pytest -q
npm test -- --run
npm run check
npm run build
```

Chrome acceptance passes with no unexpected console errors, selected interactions
alone have sticks and dashed centroid connectors, backend edits flow to the UI,
and range edits invalidate stage 3.

## Final integration review

After phase 5, inspect the combined diff rather than trusting phase-local tests.
Confirm that every new mount/share name matches between Python and TypeScript,
that all generated/runtime mount contents remain untracked, and that no unrelated
working-tree files were modified. Re-run the complete exit-check matrix once more
after any integration fix.

### Seamless capability target

- Implemented: explicit pure detector DAG step, Python transformation wrapping,
  deterministic content-addressed result identity, persistent `plain` cells,
  extraction-to-detection incremental recomputation, mounted results, and browser
  sharing.
- Intentionally omitted: remote jobserver/daskserver execution, HPC execution,
  scratch/fingertip, and extra parallelism. The detector is a small interactive UI
  graph step, and this stage provides no configured local cluster or observed
  pressure that would justify the operational complexity. These remain follow-up
  options if real structures make detection latency material.
