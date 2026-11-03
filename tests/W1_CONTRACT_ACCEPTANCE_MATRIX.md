# W1 contract acceptance matrix — A1 preparation

Status: A1 prepared the matrix; A2 ran the focused backend/frontend commands and a corrected post-fix D10 Chrome check; A4 repeated the feasible pre-fix browser checks with an isolated Vite proxy. Earlier custom-port browser runs whose Vite proxy still targeted 5813 are invalidated. An additional independent post-fix browser rerun was unavailable because the shared Chrome DevTools profile remained locked. See the execution ledger below; several contract cases remain unverified.
Authority: Alaric `plans/final/oracle-contract.md` P1–P11; `state-transitions.md` §§3–8; `contract-boundaries.md` §§8.1–8.6. The implementation plan's W1 section determines scope, not the oracle.

## Existing executable checks

Run from `/home/agent/preparic`, in the `preparic` conda environment. A2 rerun outcomes are recorded in the execution ledger below:

```sh
conda run -n preparic python -m pytest tests/backend/test_w1_project.py tests/backend/test_w1_spawn_runtime.py
conda run -n preparic npm exec -- vitest run tests/frontend/w1-project-state.test.ts tests/frontend/seamless-stacking-sync.test.ts
```

The backend cases cover mounted canonical/undecided values, preserved malformed-file diagnostics, T0→T1 envelope compatibility, file-authoritative decisions/CIF, and the eight-spawn-worker runtime without TCP services. Frontend cases cover stage prefix/basis restoration, threshold blanks, cross-project lineage, addressed matrix extent, read validation, full echoes, deferred-intent predicates and error classification. Connector mocks test stacking synchronization. These checks cannot establish browser wiring, real acknowledgments, mounted-file propagation, viewer identity, restart restoration, or follower-tab behavior.

## Browser execution protocol

Use Chrome DevTools **MCP**, following `/home/agent/.codex/skills/chrome-devtools/SKILL.md`: `list_pages`/`new_page`, navigate, `wait_for`, `take_snapshot`, then interact using fresh snapshot UIDs and the correct `pageId`. Inspect real requests with network tools and DOM/runtime observations with `evaluate_script`. Do not substitute Playwright, a CLI browser, component mocks, or injected state for user navigation. Fault injection may hold network traffic or modify disposable mount inputs; record precisely what was injected and restore it afterwards.

Use an isolated checkout/copy with disposable `deployment-mounts/` (the backend resolves that directory relative to its source). `vite.config.ts` hardcodes both backend proxy routes to port 5813. When using a unique `SEAMLESS_SHARE_PORT`, update only the disposable copy's two Vite proxy targets to the same port; record the exact diff and verify the browser response headers show the isolated backend. Launch there with `SEAMLESS_SHARE_HOST=127.0.0.1 SEAMLESS_SHARE_PORT=<unique-port> conda run --no-capture-output -n preparic ./preparic-start --log <unique-log-path>`; use the URL actually printed by Vite. Do not stop or mutate somebody else's active service. Stop/restart only this isolated service. Snapshot mounts before/after each scenario. Never use `seamless1` as the test environment.

Complete the remaining A, B and C steps in `tests/W1_BROWSER_ACCEPTANCE.md`. Corrected A4 restored stage 3 with matching panel/footer; C4 now describes a P9-compliant live edit followed by disconnect, not typing while offline. Add the following missing contract cases. A checkbox below means evidence has been captured, not merely that an automated unit check exists.

## Additional browser cases absent from the existing runbook

- [ ] D1 — **No armed Continue** (P4; ST 4.1–4.2). Hold S7 delivery at stage 2 after a completed range edit. Click/attempt Continue while disabled, then release S7. Assert no S8 PUT and no stage change until a new explicit click. Repeat stage 3 with S22/S23 withheld individually. Zero S7 rows and zero selected keys must each close their gate. Record requests, footer and panel.
- [ ] D2 — **Actual local freshness** (P5; CB 8.4). Withhold computation after a committed protein-range edit, keeping fragment keys unchanged. Assert old dependent result cannot enable Continue; retain the old state-graph `complete` status to demonstrate that status is not freshness. Revert exactly and verify a previously delivered matching result can become fresh without requiring another `update`. Verify an optimistic pending edit closes gates before its HTTP acknowledgment. Later-stage energy UI checks remain unavailable until that UI is delivered.
- [ ] D3 — **External-file changes and orphan cleanup** (P3/P4; ST 3.5, 7.2a). From stage 4 edit the disposable mounted range file to a valid changed RNA/protein row while tabs are idle. Both tabs fall to the valid prefix without browser repair PUTs. Restore exact bytes and confirm restoration. Repeat the change, then explicitly Continue through stage 3: inspect request ordering; later active markers must deactivate and acknowledge before next-marker activation. Withhold/refuse a cleanup acknowledgment: action must abort, never skip an entry gate.
- [ ] D4 — **Real read barriers and invalid values** (P7–P9; ST 8). Delay one required decision read at a time (S1/S2/S4/S8/S11/S12/S14/S24/S18); initial page stays read-only until all are read, including 204/null. Delay the first state graph and assert no failure banner. Introduce malformed/schema-unsupported mounted decisions on the disposable service: browser refuses invalid reads, backend preserves startup diagnostics before normalizing, and no observer silently rewrites decisions in response to missing computed results. Explicit null differs from an empty keyed selection. Capture backend diagnostic and browser explanation.
- [ ] D5 — **Scoped failures** (P6; CB 8.4). Inject a real T0/T1 failure, then T4, T24/T25, and T5 or T10 failure separately. Assert the banner first appears at stages 2, 3, 3, and 4 respectively, persists in later reached stages, and is absent from earlier stages. Confirm a conversion-blocked node with no exception shows `blocked`; pending relevant input writes suppress old failures. This case tests the W1 central error mapper and suppression barriers; scientific validation behavior belongs to its delivery package.
- [ ] D6 — **Stage-1 pair ownership** (W1 arbiter; ST 7.1). Failed RCSB download leaves both code and mounted CIF unchanged. Delay successful download, advance stage in the other idle tab, then release download: completion drops and reports it. Independently induce CIF conflict while the code commit still wins: retry matching CIF only while code/stage ownership holds; 1→2 remains closed until stored S21 matches accepted S1 and the viewer matches. This narrowly exercises the required code/CIF exception; it does not assert arbitrary racing-transition correctness.
- [ ] D7 — **Commit granularity** (P4). Type several characters in a range row: no ranges PUT until row completion. A checkbox action emits one keyed selection write. Stage-4 clip edits coalesce by debounce and abort on a remote revision/stage/connection change. Class debounce is deferred until stage 6 exists.
- [ ] D8 — **Fit scheduling identity** (P10; ST 5.3, 7.2a). With S8 null/inactive, stage-2 committed RNA edits launch no T5 library work. Continue 2→3 starts it. Toggle stage-3 keys without changing RNA: no repeated numerical fit identity/launch. A subsequent external range edit while S8 remains stale-active lowers stage to 2 and restarts the fit; do not require an active stage-3 prefix. Record transformation identity/execution logs rather than inferring fit behavior from the spinner. Check worker queuing/supersession under the configured eight-worker pool without introducing persistent services.
- [ ] D9 — **Persisted viewer bytes** (P2/P11; CB 8.6.5). For stages 1–3 explicitly inspect persisted CIF share reads and visible entry/assembly identity; an independent RCSB fetch is insufficient restoration evidence. At a restored stage 4, hold S23 after reconnect/restart: show `preparing all-atom structure…`, no old range atoms or CIF fallback. Detailed all-atom fidelity and overlay acceptance belongs to W3 and is excluded here.
- [ ] D10 — **URL stage-1 boundaries** (P9; ST §7). On empty stage-1 mounts, `?pdb=1B7F` must be applied and load the structure; on nonempty stage-1 mounts a conflicting `?pdb=2JDQ` must not overwrite; beyond stage 1 neither `?pdb=` nor `?continue=1` changes decisions/stage. Verify mounted files, panel, footer and hidden stage-1 input, not just URL/title.

## Execution evidence — W1 A2, A4, and A2 post-fix

### A2 executable checks

- `conda run -n preparic python -m pytest tests/backend/test_w1_project.py tests/backend/test_w1_spawn_runtime.py` — **22 passed** in 5.42s.
- `conda run -n preparic npm exec -- vitest run tests/frontend/w1-project-state.test.ts tests/frontend/seamless-stacking-sync.test.ts` — **63 passed** in 600ms (53 + 10).
- After the startup URL fix, the same two A1 commands were rerun in `preparic`: `conda run -n preparic python -m pytest tests/backend/test_w1_project.py tests/backend/test_w1_spawn_runtime.py` — **22 passed** in 8.46s; `conda run -n preparic npm exec -- vitest run tests/frontend/w1-project-state.test.ts tests/frontend/seamless-stacking-sync.test.ts` — **63 passed** in 1.00s.

### Browser target correction and A4 environment

- The A2 browser run used Vite `http://localhost:5174/` with `SEAMLESS_SHARE_PORT=5814`, but the disposable `vite.config.ts` still proxied `/preparic` and `/seamless-client.js` to 5813. A later isolated-copy attempt had the same issue with its custom port. Those runs could reach the default service, so their browser evidence and statuses are excluded below. The A2 Python/Vitest results above remain valid.
- One attempted A4 launch logged in `/tmp/preparic-w1-resume-a4.log` is also excluded: Vite auto-selected 5174 because 5173 was occupied, and the backend failed to bind 5927 (`EADDRINUSE`). Corrected evidence came from later starts after confirming both listeners were free. Those runs used `/tmp/preparic-w1-resume-a4`, Vite `http://localhost:5173/`, backend share port 5927, and Chrome DevTools pages 2 and 3 in the same browser context. The command was `SEAMLESS_SHARE_HOST=127.0.0.1 SEAMLESS_SHARE_PORT=5927 conda run --no-capture-output -n preparic ./preparic-start --log <unique-log-path>`. Valid logs: `/tmp/preparic-w1-resume-a4-cold-stage2.log`, `...-cold-stage3.log`, and `...-empty-url.log`.
- The only disposable-copy source diff was in `vite.config.ts`: both targets changed from `http://127.0.0.1:5813` to `http://127.0.0.1:5927`. Browser share responses showed `server: Python/3.15 aiohttp/3.14.5`, `X-Seamless-Marker: 1`; the unique backend listened on 127.0.0.1:5927. No service from `/home/agent/preparic` was started or written.
- For the pre-fix A4 empty-mount `?pdb=1B7F` check, the copy's PDB/CIF files were set to zero bytes before service start. Startup normalized both to one LF (SHA-256 `01ba4719c80b6fe911b091a7c05124b64eeece964e09c058ef8f9805daca546b`). Other decisions matched the source baseline: `ranges.yaml` `2dbd7bd6…742b2b68`, `stacking-selection.json` `60f0e766…2731fd2c`, `stackings.json` `c585be4c…7a4a11dc`, `biounit.json` `4355a46b…460dd865`, and `best-fit-clip.json` `3dbc498d…399562d`. Range PDB and all-atom payloads also matched the source snapshot. The test files were retained under `deployment-mounts-empty-url-test/`; after shutdown, the standard disposable `deployment-mounts/` was restored from the source baseline. The observed pre-fix failure is superseded by the post-fix A2 result below.
- The source baseline remained unchanged: PDB code `1B7F` SHA-256 `83852e3b30e7a887942565365a437cdf656b427bea820bc266b40405c109bc3f`, CIF SHA-256 `5abe0fa6bbf9cea80a7ee3ad1310e42b26e30aabff37b594a73f5a85830b8a07`. At completion the disposable copy's PDB/CIF/range/selection/stacking hashes matched that baseline. The isolated backend/Vite processes were stopped and ports 5173/5927 were free.
- The earlier `_derive_all` `OverflowError`/503 and `META__FILE__fraglib.yaml` messages were not reproduced by valid corrected-target logs and came from runs with an unverified/wrong target; do not attribute them to this isolated browser run. The corrected empty-CIF fixture did produce the expected `structure_id` failure explaining that the stored CIF was empty.

### A, B, and C browser cases

| Case | Status | Evidence and limit |
|---|---|---|
| A1 | PASS | Corrected A4 cold restart/reload restored stages 1, 2, and 3 with the same derived panel/footer and viewer identity `1B7F`; every open tab was reloaded after each restart. |
| A2 | UNAVAILABLE | The earlier overlay observation used the mismatched proxy target; a correct-target per-share readiness-delay check was not repeated. |
| A3 | UNAVAILABLE | Persisted CIF/viewer identity `1B7F` was observed at stages 1–3; stage 4+ all-atom restoration was unavailable because later W2 UI/data was not delivered. |
| A4 | UNAVAILABLE | Fresh downstream execution/result echoes after restart were not captured. |
| A5 | UNAVAILABLE | A stage-3 URL boundary was observed during the excluded pre-isolation run; the corrected-target stage-3 URL case was not repeated. |
| B1 | UNAVAILABLE | Corrected A4 confirmed idle follower propagation for 1→2 and 2→3. The 3→4 transition was unavailable because W2 UI/data was not delivered; full reachable-transition acceptance remains incomplete. |
| B2 | UNAVAILABLE | The earlier external `ranges.yaml` edit used the mismatched proxy target. A corrected-target external-file change/revert was not repeated. |
| B3 | UNAVAILABLE | A remote code/assembly change during a pending viewer load was not completed. |
| B4 | UNAVAILABLE | Stage 4 was not reached on the corrected target because W2 UI/data was unavailable; marker deactivation/reactivation was not tested. |
| C1 | UNAVAILABLE | Earlier disconnect evidence used the mismatched proxy target; the corrected-target disconnect overlay case was not repeated. |
| C2 | UNAVAILABLE | Earlier reconnect evidence used the mismatched proxy target; concurrent upstream change during same-process reconnect was not captured on the corrected target. |
| C3 | UNAVAILABLE | Deferred range edit before timer expiry and handed-off/ambiguous HTTP corrective-write paths were not exercised. |
| C4 | UNAVAILABLE | The old offline-typing instruction was invalid under P9 and was replaced with a live edit followed by disconnect before deferred download/commit. This corrected PDB-specific deferred-write case has not been run. A prior synthetic input through a disabled UI is invalid evidence. |
| C5 | UNAVAILABLE | Browser equal-value write/acknowledgment behavior was not exercised. |

### D browser cases

| Case | Status | Evidence and limit |
|---|---|---|
| D1 | UNAVAILABLE | Corrected A4 toggled all five initially selected keys off on page 2; after all acknowledgments, all seven keys were false, the stage-3 Continue control stayed disabled on both tabs, and the mounted marker hash was `2a22924b9db47d3b57af0142204ad8c15b54ff3566995e12e0366b551e434526`. Zero rows and withheld S7/S22/S23 computations were not tested. |
| D2 | UNAVAILABLE | Computation withholding, old-result freshness, and optimistic edit gate checks were not exercised. |
| D3 | UNAVAILABLE | The earlier external range change/revert used the mismatched proxy target and is excluded. Corrected-target prefix changes, marker cleanup request ordering, and withheld-ack abort behavior were not tested. |
| D4 | UNAVAILABLE | Per-share delayed reads, malformed decision values, and null-vs-empty read checks were not fault-injected. |
| D5 | UNAVAILABLE | Scoped T0/T1, T4, T24/T25, and T5/T10 error-stage behavior was not tested. The earlier W2 `META__FILE__fraglib.yaml` observation came from an unverified target and requires a separate isolated rerun. |
| D6 | UNAVAILABLE | The successful manual PDB entry was observed before isolation was confirmed and is excluded. Failed download, delayed-download/stage-advance drop, and CIF conflict paths were not exercised on the corrected target. |
| D7 | UNAVAILABLE | Range-row commit granularity and stage-4 clip coalescing were not exercised. |
| D8 | UNAVAILABLE | Fit identity and execution-log scheduling checks were not exercised. |
| D9 | UNAVAILABLE | Persisted CIF/viewer identity passed for stages 1–3. Stage-4 S23-held preparing state was not tested. |
| D10 | PASS | A2 post-fix, correctly isolated target: an empty-mount `/?pdb=1B7F` populated S1/CIF and displayed `Ready · 1B7F`; retained valid markers then restored stage 3. `/?pdb=2JDQ&continue=1` at stage 3 stayed at stage 3 with code/viewer `1B7F`. After Back to stage 1, the same conflicting URL kept stage 1 and code/input/viewer `1B7F`. Details follow. |

### A2 D10 revalidation after the startup URL fix

- Fresh copy `/tmp/preparic-w1-startup-url-fix`; Vite `http://127.0.0.1:5185/`, unique backend share port `5929`. Only the disposable copy's `vite.config.ts` was changed: proxy targets use 5929 and Vite binds 5185. Browser pages 2 and 3 used the same isolated namespace. The original checkout/service ports were not started or changed.
- Before startup, disposable `pdbcode.txt` and `cif.cif` were both zero bytes (SHA-256 `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`). After opening `/?pdb=1B7F`, the app loaded and committed the structure; S1 returned 200 with marker 2 and final mounted `pdbcode.txt` contained `1B7F` (SHA-256 `83852e3b30e7a887942565365a437cdf656b427bea820bc266b40405c109bc3f`). CIF returned 200 with marker 2, length 379599, and SHA-256 `5abe0fa6bbf9cea80a7ee3ad1310e42b26e30aabff37b594a73f5a85830b8a07`; both final hashes match the original baseline. The viewer showed `Ready · 1B7F`, and the retained matching ranges/selection markers re-derived stage 3. No explicit RCSB URL appeared in the DevTools request list; the mounted pair and rendered structure did.
- With the project at stage 3, page 3 opened `/?pdb=2JDQ&continue=1`; it stayed on the Stackings panel/footer, with `projectStage=3`, `appState.stage=3`, shared/local PDB code `1B7F`, viewer `Ready · 1B7F`, and no 2JDQ request/write.
- After Back 3→2→1, page 3 reloaded `/?pdb=2JDQ&continue=1`. It stayed at stage 1 with the `1B7F` input and viewer ready; PDB/CIF share reads returned the existing `1B7F` pair. No conflicting code was applied and Continue was not auto-fired.
- After shutdown, the disposable `deployment-mounts/` was restored from `/home/agent/preparic/deployment-mounts/`; `diff -qr` was clean. PDB/CIF hashes returned to `83852e3b30e7a887942565365a437cdf656b427bea820bc266b40405c109bc3f` and `5abe0fa6bbf9cea80a7ee3ad1310e42b26e30aabff37b594a73f5a85830b8a07`. Vite/backend services were stopped; ports 5185/5929 were free.

### Independent post-fix A4 attempt

- A fresh copy at `/tmp/preparic-w1-a4-postfix` was prepared with empty stage-1 PDB/CIF mounts and its two proxy targets set to the unique share port 5931. The service started in `preparic` with Vite at `http://localhost:5173/` and backend `127.0.0.1:5931`.
- Chrome DevTools MCP `new_page` and, after cleanup, `list_pages` both failed because the shared Chrome profile was already running in another agent session. No page navigation, snapshot, network inspection, or browser write was performed in this attempt. The service was stopped; ports 5173/5931 are free; the disposable copy's standard mounts were restored from its baseline snapshot (PDB `83852e3b…c3f`, CIF `5abe0fa6…8a07`, ranges `2dbd7bd6…2b2b68`, selection `60f0e766…31fd2c`). The independent recheck is therefore unavailable; D10 remains PASS on the implementation agent's corrected A2 Chrome evidence above.

### W1 fixes and focused browser evidence

- `src/lib/seamless.svelte.ts`: `setSharedPdbCode` now rejects writes unless the share is live and the project is at stage 1.
- `src/components/Stage1Structure.svelte`: PDB input disables only when the viewer/share is unavailable or `appState.loading` is true. The pre-fix A4 corrected-proxy run showed the input enabled but did not trigger a load; the post-fix A2 run populated the requested structure.
- `src/App.svelte`: the sidebar remains interactive at stage 1 while no structure is loaded. The pre-fix A4 corrected-proxy run showed `sidebar.inert=false`; no stage-2/3 empty-CIF inertness claim is carried forward from the mismatched-target run.
- `src/App.svelte` and `src/components/Stage1Structure.svelte`: stage-1 URL intent now remains local until Stage1 queues the normal structure load; the pending PDB intent is marked after `queuePdbLoad` so `cancelSharedCifWrite()` cannot erase it. The existing code/CIF arbiter still commits the pair after load success. A2's corrected D10 rerun passed.
- Post-fix A2 Chrome evidence: on empty stage-1 mounts, `/?pdb=1B7F` populated the shared PDB/CIF pair and showed `Ready · 1B7F`; the retained valid decision markers then restored stage 3. S1 and CIF share responses were 200 with marker 2, and the CIF response length was 379599. The DevTools request list did not explicitly record an external RCSB URL request, so this evidence establishes the observed mounted pair and viewer result without claiming that external request was captured. A2 also confirmed the conflicting URL stayed at stage 3 and, after Back to stage 1, did not replace the existing code/input/viewer.
- A clean-close reconnect and the earlier `META__FILE__fraglib.yaml` report were observed only in pre-isolation runs; neither is counted as corrected-target browser evidence.

## Evidence ledger and completion rule

For every A–D scenario record: case ID; backend/mount directory; initial and final decision files; two page IDs/URLs; snapshots before/after; pertinent share requests and acknowledgments (including equal-value success, 409 and ambiguous response where applicable); result echoes/state graph; logs; actual versus expected outcome. Mark each `PASS`, `FAIL`, or `UNAVAILABLE` with reason. Existing pure tests are supporting evidence, never replacements for browser evidence. A prepared checkbox/runbook is not a passed test.

Current exact gap: D10 passed on the implementation agent's corrected post-fix Chrome run; this agent's independent browser rerun could not proceed because the shared Chrome profile was locked. Cases marked `UNAVAILABLE` remain unverified, including the full D6 failed-download/stale-stage sequence and the corrected P9-compliant C4 deferred-write sequence. The two A1 executable commands passed after the fix. W4–W9-only UI remains `UNAVAILABLE`. Arbitrary simultaneous tab transition races and tabs kept alive across backend restart are excluded by O7/O6.
