# Focused W1 browser acceptance set

Prepared for the current W1 cycle; **not executed by the test author**. This set
resumes the interrupted browser validation. Existing historical PASS entries in
`tests/W1_CONTRACT_ACCEPTANCE_MATRIX.md` are not results for this cycle.

Authority: `alaric/plans/final/oracle-contract.md` P1–P11, with secondary detail
from `contract-boundaries.md` §8 and `state-transitions.md` §§3–8. W1 exit criteria
are in `implementation-plan.md` §4. No arbitrary racing transition guarantee is
tested (O7); backend restarts require immediate browser reload (O6).

## Prerequisites and launch

Use Chrome DevTools MCP and its skill. Node dependencies must already be installed
in `/home/agent/preparic/node_modules`; `preparic` must contain the configured
Seamless installation and scientific inputs used by the current checkout.
Chrome MCP `list_pages({})` was the only browser capability probe performed by
the author; it succeeded with one `about:blank` page. No application was opened.

Create a disposable source/mount copy. Do not reuse a previous run's mounts.
These commands are preparation commands for the execution agent:

```sh
task_dir=$(mktemp -d /tmp/preparic-w1-acceptance.XXXXXX)
rsync -a --exclude=.git --exclude=node_modules /home/agent/preparic/ "$task_dir/"
ln -s /home/agent/preparic/node_modules "$task_dir/node_modules"
cp -a "$task_dir/deployment-mounts" "$task_dir/mounts-before"
sed -i 's/127.0.0.1:5813/127.0.0.1:5947/g' "$task_dir/vite.config.ts"
```

Confirm 5947 is unused (`ss -ltnp 'sport = :5947'`). If occupied choose a free
port and replace it in **both** proxies and the environment below. Start from
the disposable copy, recording the actual Vite URL, launcher session and log:

```sh
cd "$task_dir"
SEAMLESS_SHARE_HOST=127.0.0.1 SEAMLESS_SHARE_PORT=5947 \
  conda run --live-stream -n preparic ./preparic-start --log "$task_dir/service.log"
```

Use Ctrl-C on this launcher for restarts/shutdown. Do not use a global stop
command against other services. Before tests, inspect a real share response via
MCP `get_network_request`; verify the disposable backend listener and proxy diff,
response marker and mounted bytes. A wrong proxy invalidates the entire run.

Supporting executable tests, run separately from the browser cases:

```sh
cd /home/agent/preparic
conda run --live-stream -n preparic python -m pytest tests/backend/test_w1_project.py tests/backend/test_w1_spawn_runtime.py
conda run --live-stream -n preparic npm exec -- vitest run tests/frontend/w1-project-state.test.ts tests/frontend/seamless-stacking-sync.test.ts tests/frontend/seamless-debounced-clip.test.ts tests/frontend/w1-reconnect-regression.test.ts
```

## MCP commands and observation assertions

Discover tool schemas, then invoke `new_page` twice with the recorded Vite URL
in the same default browser context. Save the returned page IDs as A/B. For each
navigation/action: `wait_for`, `take_snapshot`, then use fresh snapshot UIDs for
`click`/`fill`. A UI action must use MCP interaction, never an application setter.

At every checkpoint execute this exact tool call, substituting A or B's page ID:

```js
await tools.mcp__chrome_devtools__evaluate_script({
  pageId: A,
  sourcePath: '/home/agent/preparic/tests/browser/w1/observe.js',
  waitForStableDom: true
});
```

These files are function expressions accepted by `evaluate_script`'s default
`function` format. They dynamically import the **existing live Vite modules** and
read state/DOM; they do not mutate state, dispatch events, fetch share values,
or repair the implementation. A result with `pass: false` fails that checkpoint.
Capture the full returned snapshot in evidence. A passing observation alone does
not pass a scenario: verify network requests, CIF bytes and mounts independently.
Do not run these scripts with Node or treat them as Vitest modules.

Use `list_network_requests({pageId: A, includePreservedRequests: true})` and
`get_network_request({pageId: A, reqid: ...})` to retain PUT bodies, response
codes/markers and reread responses. Record mount hashes with:

```sh
find "$task_dir/deployment-mounts" -type f -exec sha256sum {} +
```

## Required focused cases, in execution order

1. **W1-B1 idle follower (P1/P3/P4, O7).** Restore a valid `1B7F` project at
   stage 1 using Back actions. Record both observers and mount hashes. Continue
   in A, leaving B idle. Assert both show stage 2 and matching ranges marker;
   Continue A to stage 3 and assert both show the stackings panel/footer. Repeat
   3→4 if scientific inputs make its gate ready. Back in A must deactivate the
   current marker while retaining its basis/data and move B back without reload.
   Compare decisions in both observers and mounts. Record unavailable stages
   separately; a missing stage does not pass that transition.

2. **W1-C4 interrupted stale PDB intent (P4/P9, ST 4.1).** Return both tabs to
   stage 1 with accepted `1B7F`, ready viewer and no pending writes. Save code/CIF
   bytes and marker. In B apply MCP `emulate({pageId: B,
   networkConditions: 'Slow 3G'})`; type `2JDQ` into the enabled PDB input using
   its current snapshot UID. Before its download/commit completes apply
   `emulate({pageId: B, networkConditions: 'Offline'})`. Confirm B actually reaches
   `status: disconnected`, `ready: false`, overlay present, and no accepted
   `2JDQ` code/CIF PUT. If the offline emulation does not close the websocket or
   the write already committed, this attempt is **not established**: reset by
   normal UI actions and repeat with stronger throttling (CPU throttle plus Slow
   3G) or a recorded transport fault. Never type into a disabled offline control
   or manufacture disconnect by assigning reactive state.

   In live A press Continue to stage 2, await its acknowledgment, then restore B
   with `emulate({pageId: B})`. Wait for real reconnect/rereads and viewer restore.
   Run the common observer and this specialized assertion:

   ```js
   await tools.mcp__chrome_devtools__evaluate_script({
     pageId: B,
     sourcePath: '/home/agent/preparic/tests/browser/w1/assert-c4.js'
   });
   ```

   Assert `pass: true`, unchanged mounted `1B7F` code/CIF, and no deferred `2JDQ`
   commit after reread. A hidden input, when mounted, must say `1B7F`; if absent,
   return B to stage 1 and verify the newly visible input says `1B7F`. Reload B
   and verify it remains `1B7F`, with no repair PUT. Keep snapshots before fault,
   disconnected, after reconnect, after Back and after reload. The specialized
   assertion expects stage 2 and must run **before** Back.

3. **W1-C2 distrust until reread (P5/P9).** With B offline in the same backend
   process, change a completed upstream row in A. After reconnect, observe B's
   read-only barrier and cleared freshness before reread, then matching decisions
   and stage after reread. Network evidence must contain fresh required decision
   reads and computed rereads before those results become usable. A previously
   held matching echo or graph `complete` state alone cannot establish freshness.
   The common observer asserts overlay/inertness and absence of trusted results
   whenever `ready` is false. Transport faults must be established as in C4.

4. **W1-D6 failed download (P1/P4/P7).** At stage 1 snapshot code/CIF. Type
   `0ZZZ` while live, wait for a confirmed failed RCSB download (or inject and
   record an actual download failure). Assert an explanation and unchanged code
   and CIF bytes; no PUT for the failed structure. The reserved-looking code is
   not assumed permanently absent: verify the request actually failed. Reload
   and assert accepted structure/input restored without a browser repair write.
   CIF-conflict ownership/retry is still covered by the full D6 runbook and must
   not be reported passed by this failed-download subcase.

5. **W1-A cold restart (P2/P3/P5, O6).** At every reached stage, save all decision
   files and CIF, stop/restart only the disposable launcher with identical mounts,
   immediately reload **both** tabs and take early/settled observations. Before
   results return the same footer/panel stage must restore from decisions. Mounted
   decisions/CIF retain values; no browser reset/repair writes. Inspect logs and
   result identities to establish cold fit/downstream execution at stage ≥3.
   Missing results close freshness-dependent gates. Stage 4's blank skeleton
   remains stage 4 and editable; stale S23 shows all-atom preparing wording.
   Viewer atoms must be checked independently: title/structureId is insufficient
   evidence that persisted bytes or fresh S23 were rendered.

6. **W1-D1 no armed Continue (P4/P5).** At stage 2 while a real S7 result is
   pending after a completed row edit, capture disabled Continue and attempt the
   user action; let computation complete. Verify no S8 PUT or stage advance until
   a second explicit Continue. At stage 3 uncheck every stacking, await all
   acknowledgments and verify Continue disabled in both tabs. Repeat with S22 or
   S23 pending when attainable; record unavailable fault windows explicitly.

7. **W1-D10 URL boundary (P9).** At populated stage 1 and at stage 3 open
   `/?pdb=2JDQ&continue=1` in a new idle tab. Observe accepted code/input/viewer
   `1B7F` and unchanged stage/decisions without unwanted PUTs. Run empty-stage-1
   `/?pdb=1B7F` using a **separate disposable mount copy**, asserting successful
   stored code/CIF identity. Do not erase the main run's mounts to create it.

The full `tests/W1_BROWSER_ACCEPTANCE.md` A–C and matrix D1–D10 remain applicable,
especially per-share read barriers, malformed reads, cleanup acknowledgment
ordering, optimistic/revert freshness, equal-value HTTP acknowledgment, ambiguous
in-flight corrective writes, T10/T27 failure vetoes and fit identity. This focused
set prioritizes the interrupted feasible browser paths; it does not replace those
obligations or mark unavailable later-stage UI complete.

## Focused regressions for current W1 fixes

**W1-C2b retained result with omitted/failed reread (P5/P9).** Retain old S22
`fragments` while its current-connection reread fails (503 for all retries) or
it is omitted from the reconnect share list. Let all decision reads and S21
succeed. After the decision barrier opens, the retained matching S22 echo must
remain stale and 3→4 must remain closed. A real newly delivered S22 can restore
freshness. If per-share transport faults are unavailable, mark this browser
subcase unavailable and use the deterministic connector regression
`tests/frontend/w1-reconnect-regression.test.ts` as supporting evidence. It
exercises both faults, retained payload, decision-barrier release, zero repair
writes, and subsequent real result delivery. It is included in the exact Vitest
command above and has not been executed by the author.

**W1-D4 malformed/baseless decision explanations (P7/P8, I7).** Snapshot the
exact original file before each injection, using disposable mounts only. While
the service runs at stage 3 replace `best-fit-threshold-edits.json` with
`{"schema":1,"active":true,"pairs":[],"bases":[]}` (missing basis). Await
file-authoritative propagation. Both tabs must refuse that marker, remain/fall
to stage 3 and show an explanation for the next-stage threshold decision. The
observer captures `decisionErrors` and visible `alerts`; the explanation must
persist through unrelated share deliveries. A transient error immediately cleared
by recording the refused value fails. There must be no browser repair PUT.
Restore valid original bytes and check a valid read clears the explanation.

At stage 1 repeat with malformed distant `gaussian-classes.json`: the diagnostic
may be stored, but no Gaussian error banner may leak into stage 1. Restore its
bytes. Then replace `ranges.yaml` with `[]` (legacy/baseless value). Both tabs
must fall to stage 1 and show a visible current/next-stage ranges explanation
there, even though the range editor is absent. A message hidden exclusively in
the stage-2 editor fails. Restore original bytes and verify a valid read clears
the error. Use actual mount filenames from the inventory if they differ and
record the mapping. Separately restart with malformed stored content and capture
the preserved backend startup diagnostic before normalization; a startup reset
to null does not demonstrate browser handling of a malformed received value.

## Evidence and completion

Write one row per case/subcase: PASS/FAIL/UNAVAILABLE, initial/final mount hashes,
page IDs/URLs, MCP snapshots/assertion results, relevant share request IDs and
acknowledgments, stored CIF/S23 byte evidence, backend log path and exact fault.
Capture both tabs after propagation. Preserve intermediate reconnect evidence;
checking only a settled page cannot prove the read barrier. Stop the disposable
launcher when done. Leave application files untouched and preserve the initial
mount copy for comparison. No prepared runbook or successful capability probe
counts as executed acceptance.
