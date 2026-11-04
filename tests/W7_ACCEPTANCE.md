# W7 A1 acceptance set

Authored independently of W7 implementation; **not executed**. Baseline: Preparic
`b1b7aea`, Alaric `63972dc`. Authority: implementation-plan §9; OC P6.1–P6.6;
ST 4.2/4.4/6; CB stage 6, S18–S19/S25. Run Python checks only in the
`preparic` Conda environment (`conda run -n preparic python -m pytest ...`).
Browser acceptance must use **Chrome DevTools MCP**, not Playwright.

W8b amendment (2026-11-04): reference thresholds come directly from computed
S9 bases. S11 contains only editable pairs. T23/S20 and the reference panel are
removed; 4→5 additionally requires fresh S9. Over-limit computed reference
thresholds become disabled controls in the export bundle, with associated
comments, and contribute no Proceed-anyway condition.

## Backend and graph cases

Implement these as focused Preparic backend acceptance tests once W7's callable
interfaces exist; do not silently skip absent interfaces. Use copied JSON and tiny
candidate descriptors for pure cases, plus one real 1b7f integration fixture.

| ID | Setup / action | Required assertion |
|---|---|---|
| B01 | T28 runs with delivered A7 JSON; mutate schema to absent, 0, 2, string `1`; remove or add an illegal threshold name | Correct schema and all twelve model names required; malformed input raises; no fresh S25 and 5→6 remains closed. Names/cutoffs originate in JSON. |
| B02 | Independently corrupt one density golden and one acceptance/boundary golden; run T28 twice with second-run corruption | Every run checks relative `1e-12` score tolerance and acceptance goldens; second run fails rather than trusting a prior successful load. Include zero-valued vectors without dividing by zero. |
| B03 | Descriptor at distance 5 and phi 60 degrees, and immediately inside each; score at cutoff and below | Strict JSON mask rejects equality; `score >= cutoff` accepts equality inside mask. Raw density remains available outside mask but cannot yield passing auto class. Cutoffs are delivered five-decimal values, including variant parent cutoff. |
| B04 | Run T18 on 1b7f best-fit PDBs and matching protein ranges, score via T28 and independent A9 job results | Same selected candidate identities and score absolute error <= `1e-6`. Pin the six ring-reference inputs; Kabsch normals agree with T4. Diagnose three-decimal PDB rounding if failing; never relax tolerance. |
| B05 | One anchor has candidate scores spanning multiple class cutoffs, a lower clipped candidate, and a high reference score | Lowest **unclipped** score picks smallest passing N from 50,60,70,80,90,95; reference is display-only. Test equality at each cutoff. Below last cutoff retains 95 with below-RRM-95 diagnostic. Empty/error candidate set is explicit, never silently drops anchor. |
| B06 | Class 90/95 with candidate theta 29.999/30 and phi 44.999/45, including clipped failures | Auto flags require every unclipped candidate to satisfy strict limits. Class <=80 has no flags. Mask failure cannot establish acceptance even at high raw density. |
| B07 | Permute JSON threshold insertion order and vary a legal cutoff/model version | T19/T28 use JSON keys; S25 carries all twelve names and filter_set_version. No order-dependent filter mapping or separate browser legal-name list. Class-table row selecting a missing name is rejected. |
| B08 | S18 rows: duplicate, missing/extra selected key, invalid rrm, bool-as-number, absent/nonboolean flags at 90/95, flags at <=80 | Frontend and T19 reject malformed active tables; exact selected-key set, one row each. Null/inactive marker is sentinel. Backend never writes S18. |
| B09 | Compute T20/T21 identity then edit each class/flag, revert, Copy and deactivate/reactivate S18 | T20 job list and T21 transformation identity stay identical; no new A9 calls. One job per existing kept best-fit candidate, full fragment/nucleotide/conformer/rotamer/grid/range/protein identity; missing pose yields diagnostic, not invented job or discarded anchor. |
| B10 | Feed T22 fixed all-filter oracle results; select another legal class | Merge changes immediately with no T21 rerun. Worst candidate governs row: green emitted+passing with fill >=1/3; amber 0<fill<1/3 or positive fill+pose fail; red zero; explicit implied result remains non-green and never invents 0% fill. Per-job error survives merge. Include exact fill 1/3. |
| B11 | Inspect stage-6 backend graph and browser shares | No T23, S20, refe-validation mount/share or reference-anchor panel remains. |
| B12 | Export computed candidate thresholds 1.73, 1.74 and raw 1.731 rounded 1.74 | 1.73 is copied; above sqrt(3), action has numeric -1 plus raw/rounded/coverage comment. Independent export conditions remain effective. |
| B13 | Change S11 pairs while holding fit, clip, descriptors and classes fixed | Fitting/T18/T28/T20/T21 identities stay unchanged. S24 changes trigger only T27. Neither overwrites S18. |
| B14 | Delay each computed share, then deliver old content with successful node status | Freshness checks compare content echoes, including optimistic S18 rows for S19, S12, selected keys, S10 echo; S25 includes S17 echo/selection/clip. Old success remains grey/pending; failed T18/T28 closes entry and displays banner; failed T21 displays installation/error banner. |
| B15 | Back 6→5 then Continue, with same selected basis, changed fit/model; repeat with changed selected basis | Same basis restores edited rows, updates computed/diff marks; changed basis seeds fresh computed rows. Continue requires fresh S17/S25, valid readable matrix, complete true S24 cells and nonfailed T27. No automatic write repairs on failed reads. |
| B16 | Inspect backend graph and browser dependency/build inputs | No T17 or browser Gaussian model share; B15 scoring module is attached to T28; source/bundle contain no Gaussian model mean/precision/normalizer/cutoff/golden evaluation or density kernel. Text labels and filter-name formatting are allowed. Verify import graph as well as token search; Molstar geometry alone is not Gaussian scoring. |

## Frontend pure cases

Use existing Vitest patterns in tests/frontend and fake timers for debounce.

- F01: invalid S18 shape cannot establish stage 6; stale computed shares never
  change derived stage from valid edit markers alone.
- F02: <=80 drops flags; returning to 90/95 seeds current auto flags, not stale
  edited flags. Revert changes only chosen row; diff follows current computed row.
- F03: Copy cancellation queues no write; confirmation queues one replacement.
- F04: local edit updates table, diff and both relevant pending indicators
  immediately; at 499 ms no PUT, at 500 ms one latest-value S18 write. Multiple
  edits coalesce. Backend response must not flash the previous green result.
- F05: remote S18 update or Back during debounce discards deferred stale action;
  409 adopts server value. Equal-value notifications must follow existing client
  retry/authority semantics; no assumed cancellation of an already queued PUT.
- F06: incoming computed/model changes retain edited rows; version changes display
  `model changed`. Invalid/stale names yield actionable state, not a fallback class.

## A4 browser scenarios (Chrome DevTools MCP)

Use two tabs against an isolated mounted project copy. Record snapshots, console
errors, network writes, S18/S11 mounted values, and transformer identities before
and after. Keep real A9 results for numeric acceptance; synthetic share delays
are permissible for race scenarios and must be labelled as injected.

1. **C01 entry:** navigate 1b7f through stages 3–5 with all required S24 cells
   filled. Observe disabled Continue while S17/S25 stale, then enter stage 6.
   Snapshot candidate/reference scores, computed/editable classes, twelve names,
   the Gaussian panel; count one row per selected stacking. Fail a copied golden and
   verify T28 banner/closed gate; restore fixture afterwards.
2. **C02 editor:** change a 90/95 flag, set <=80 then back, revert one row, cancel
   Copy then confirm Copy. Inspect network for immediate optimistic display and
   exactly latest S18 PUT after >=500 ms quiet. Capture unchanged T20/T21
   identities and changed T22 merge output.
3. **C03 validation:** choose filters exercising emitted-pass, low fill, zero,
   positive-fill pose failure and implied result; compare to A9 outputs and
   confirm worst-candidate row. Delay S19 while editing and verify grey
   `validating…`; inject job/cache installation error and verify visible banner.
4. **C04 computed reference controls:** verify no stage-6 reference panel or base
   inputs exist. Stage 4 still shows fragment residues, sequence, conformer and
   fit RMSD read-only. Check pairs-only S11 and fresh-S9 admission at 4→5.
   Sentinel/export comment coverage belongs to the W8b bundle tests.
5. **C05 persistence/drift:** edit classes, Back to 5, Continue; reload and open
   second tab. Applicable edits survive. Change attached model version or refit
   fixture: computed values/diff/model-changed update without replacing S18.
6. **C06 races:** while S18 debounce pending, remote-edit same cell in tab B;
   repeat with Back before timeout and with conflict response. A must adopt remote
   state and never resurrect stale edit; while S11/S18 writes are optimistic,
   delayed old S19 must remain grey. Restore isolated mounts after scenario.
7. **C07 no browser scoring:** inspect loaded modules/network response bodies
   with DevTools and corroborate B16 source/import/build audit. Browser receives
   S25 scores/classes/names, not JSON Gaussian parameters or a scoring module.
   No runtime browser-vs-A9 numerical comparison is introduced.

A4 evidence must identify which cases ran, fixture/source revisions and every
injection; missing cases remain unverified rather than marked passed.

## Initial executable set and A2 handoff

Added without running: `tests/backend/test_w7_class_decision.py` uses the existing
`backend.load_marker_payload`; `tests/frontend/w7-class-state.test.ts` uses existing
`markerIsWellFormed`, `echoIsFresh`, and `deferredWriteIsCurrent`. These exercise
strict S18 schemas and retained decisions, per-flag validation, exact keys,
optimistic content mismatch and invalidation of deferred writes. They do **not**
claim to test real debounce timing or T28, which require the W7 implementation.

From Preparic, A2/A4 must run:

```sh
conda run -n preparic python -m pytest tests/backend/test_w7_class_decision.py
conda run -n preparic npm test -- tests/frontend/w7-class-state.test.ts
conda run -n preparic npm run check
```

A2 must additionally provide executable coverage for B01–B15 against its actual
new interfaces and fake-timer coverage of F02–F05 against its actual editor/write
code. In particular the `500 ms` boundary, real T20/T21 identity stability,
T28 golden failure on a second run, real 1b7f A9 parity and W8b sentinel behavior
cannot be discharged by the initial state-helper vectors. No new callable
signature is imposed by this test handoff.
