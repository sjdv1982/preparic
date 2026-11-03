# Focused W3 Chrome DevTools acceptance

Prepared, **not executed**. This resumes `tests/W3_BROWSER_ACCEPTANCE.md`;
it does not convert that deferred runbook into a historical pass. Primary oracle:
`/home/agent/alaric/plans/final/oracle-contract.md` P4.1–4.4, P5, P10, P11;
secondary: CB T0–T11/T24/T25 and §8.6.5, ST §4/§7. W3 exit is implementation-plan §6.

Use the Chrome DevTools MCP skill. No Playwright, browser CLI or Node execution
of observation files. Browser interaction uses fresh `take_snapshot` UIDs and
MCP `click`, `fill`, `hover`, `handle_dialog`. Read-only `evaluate_script` may
inspect live Vite state. Never assign reactive state to manufacture a result.

## Launch and supporting tests

Use the disposable-copy launch instructions in `../w1/README.md`, replacing
the temporary directory prefix with `preparic-w3-acceptance`. Verify its proxy
and share listener before browser actions; do not use production mounts.
Launch command from that disposable directory:

```sh
SEAMLESS_SHARE_HOST=127.0.0.1 SEAMLESS_SHARE_PORT=5947 conda run --live-stream -n preparic ./preparic-start --log "$task_dir/service.log"
```

The following are exact supporting commands, separate from browser validation:

```sh
cd /home/agent/preparic
conda run --live-stream -n preparic python -m pytest tests/backend/test_w3_preparation.py tests/backend/test_w3_stage4.py tests/backend/test_w3_composite_assembly.py
conda run --live-stream -n preparic npm exec -- vitest run tests/frontend/w3-stage4.test.ts
conda run --live-stream -n preparic npm run check
```

No `seamless1` environment. Do not report these commands as executed until a
separate execution agent runs them.

At each checkpoint save `take_snapshot`, screenshot, relevant GET/PUT responses,
the output below, and mounted decision bytes. Substitute the actual page ID:

```js
await tools.mcp__chrome_devtools__evaluate_script({
  pageId: A,
  sourcePath: '/home/agent/preparic/tests/browser/w3/observe.js',
  waitForStableDom: true
});
```

`pass` checks only local invariants; the numbered scenarios require independent
temporal/network/backend evidence. S23 hashes attest received texts, **not** the
atoms actually rendered. Never substitute a title or canvas screenshot for atom
identity evidence. Capture results outside tracked test files.

## Required cases

1. **W3-B1 early stage 4 with held fit (P4.1/P4.3/P10).** Use mounted valid
   1B7F or another scientifically supported project; record exact range rows,
   preparation readiness and library identity. At stage 2 verify no active S8
   and no fit library execution. Continue to stage 3, select at least one anchor,
   and establish a real pending T5 window. If needed inspect launcher descendants
   (`ps -eo pid,ppid,args --forest`) and SIGSTOP only the identified disposable
   fit child, recording PID, command, parent and time; release with SIGCONT in a
   cleanup/finally step. Do not pause the share server or the whole backend.
   If the fitter cannot be isolated, use an actual slow cold fit and record the
   window; an already finished fit does not establish this case.

   With fresh S22/reference fragments/S23 and T5 still pending, click Continue
   to best fits. Assert stage 4, keyed empty skeleton, visible computed column,
   editable inputs, Copy disabled and fitting/pending explanation. Type `2.34`
   into `Fragment 1 nucleotide 1 ovRMSD`, blur using a real UI action and wait
   for its threshold PUT acknowledgment. Confirm mounted value and observer
   retain `2.34` while fit remains held. Stage 4→5 stays closed with scoring
   unavailable. Release T5; observe fresh S9, Copy enabling and stage staying 4
   with the exact typed value unchanged. Retain state-graph/GET timing proving
   fresh results arrived after the edit, rather than merely a disabled button.

2. **W3-B2 Copy, validation, revert (P4.3).** While fresh computed values are
   visibly alongside edits, click Copy; dismiss the overwrite dialog and assert
   the typed value and decision bytes remain. Repeat and accept; all available
   cells equal S9, and diff marks/count agree. Type a value different from its
   computed cell, assert mark/count, then use that cell's Revert and assert only
   that cell returns to the computed string. Type `10.01`, assert above-10 warning.
   Try `1.234`, `0`, `-1`, `NaN`, `Infinity`: each must be refused/explained with
   no invalid or rounded value committed. Restore a valid value between cases;
   capture PUT bodies and decision snapshots, not only input DOM strings.

3. **W3-B3 independent protein edit and cold-fit restoration (P4.3/P5/P10).**
   Record the edited table, T5 identities and fresh S9, Back to stage 2 through
   stage 3, change only a protein range while retaining RNA range 1, and Continue
   normally. Assert T5/S9 identity and freshness remain unchanged: protein-only
   edits do not refit RNA or invalidate S9. Per-range energetics rescoring belongs
   to W4/W5 and must not be claimed established by this W3 browser case.

   Separately, at stage 4 save the acknowledged typed threshold and all mounted
   decisions. Cold-restart only the disposable backend with identical mounts,
   PDB code and RNA basis, and immediately reload every browser tab (O6).
   Establish pending T5 as in B1 if necessary. Assert restored stage 4 and exact
   typed threshold, editable inputs, and Copy closed until this connection reads
   fresh S9. Release T5 and capture the new S9 reread; Copy enables, stage stays
   4, and the user's threshold stays unchanged. Preserve before/after S9 echoes,
   decision bytes, graph execution and GET/PUT timing. Matching echoes alone do
   not establish current-connection freshness. A completed cache hit without an
   observed stale interval does not prove freshness blocking. A changed RNA
   basis must instead create a new empty skeleton when Continue activates S11.

4. **W3-B4 prepared bytes, visible atoms and overlays (P11/P4.2).** With no row
   hovered, match each observer S23 SHA-256 against the corresponding disposable
   `deployment-mounts/pdbs-aa/range-N-aa.pdb` (`sha256sum`); retain the share GET
   response and files. Obtain the actual live Mol* state through available
   DevTools inspection/debugger access, without creating a second viewer. For
   each visible prepared structure enumerate atom name, residue identity,
   coordinates and range label; compare the multiset with the matching backend
   PDB ATOM/HETATM fixed columns (coordinates within PDB precision). Inspect
   visibility state to exclude hidden/staged old structures and require no
   downloaded CIF atoms among visible representations. Two files with the same
   internal chain ID must still map to their distinct `range-N` labels. If the
   fixture lacks this condition, record that subcase unestablished.

   Hover a base row via MCP `hover`; assert observer `overlayFragments` equals
   `[fragment]`, reference/best-fit structures visible with original coordinates,
   and unchanged base structures. Then hover a pair row: assert the observer
   array equals `[down, up]` in that order and capture the live visible Mol*
   labels/atoms proving reference and best-fit overlays for **both** fragments
   appear together. No overlay from the previously hovered row may remain.
   Move away and confirm the array is `[]` and all fragment overlays disappear.
   State arrays alone do not prove that overlays rendered. Inspect table headings/cells to
   exclude raw threshold values. S23 hashes and screenshots alone cannot pass
   viewer atom or overlay checks. If live Mol* introspection is unavailable,
   report the byte check separately and mark atom equality unestablished.

5. **W3-B5 stale S23 hides former structure (P11/P5).** Establish stage 4, then
   restart only the disposable backend with identical mounts and immediately
   reload the browser (O6). Hold a real preparation child if needed as in B1.
   Before fresh S23 reread, stage remains 4, threshold inputs remain editable,
   viewer shows “preparing all-atom structure…”, and neither old prepared nor
   downloaded structure is visible. After release/reread, repeat B4. Also cover
   range edit → reentry when a stale S23 delivery can be held independently;
   if 3→4 correctly blocks on missing S23, record that gate instead of claiming
   it exercises stale S23 while already at stage 4.

6. **W3-B6 clip/anchors/completeness (P4.4/P10/D2/D4).** With unanchored outer
   fragments, set a legal integer clip; verify clipped rows grey but editable.
   Leave only clipped-row thresholds blank, fill kept rows, and verify gate
   depends on fresh energetics/T10, not Copy-source freshness. Unclip a blank
   kept row and verify gate closes. Refuse a clip cutting any candidate of a
   selected anchor. Retain legal clip, Back to stage 3, select an anchor whose
   candidate conflicts, Continue: network acknowledgments must show clip clamp
   accepted before threshold marker activation. Backend task identities/counters
   before/after clip must prove no A1–A6 or gridless rerun; browser alone cannot
   establish that. Capture graph evidence or mark this subcase unestablished.

Stages 5–6 editability/read-only clip and viewer checks remain downstream
requirements. Placeholder panels or unavailable Continue 4→5/5→6 routes must be
reported as downstream/unavailable, never an acceptance pass for those stages. Keep
the assembly build/disconnect cases in the original deferred runbook applicable.
