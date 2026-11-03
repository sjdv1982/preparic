# W4 pending energetics retest — 2026-11-03

## Scope

Focused Chrome DevTools validation of the W4 defect where pending energetics could be reported as an anchor scoring error. This checks the stage-5 pending UI after a backend restart; it does not validate numerical energetics or A5 scoring correctness.

## Setup

- Ran the disposable Preparic copy at `/tmp/preparic-w4-acceptance.x68Rcz` with the `preparic` conda environment (Vite 5283, share/backend 5983).
- Entered stage 5 through the UI with fresh energetics, then restarted the disposable service and immediately reloaded the stage-5 page.
- The initial energy and threshold data requests returned 404 while the restarted graph was pending (`/preparic/energetics` requests 590 and 606; `/preparic/energy_thresholds` requests 591 and 607).
- The disposable deployment mounts were compared with the saved pre-run snapshot and matched. The service is stopped and its browser tab is closed.

## Result

The immediate post-restart snapshot remained on stage 5 and showed:

- Energetics values as `—` and computed limits as `— / —`.
- “Waiting for current gridless energetics…” status.
- Recompute enabled and Continue to Gaussian filters disabled.
- No “Selected anchor … has a scoring error” message.

The graph subsequently completed; stage 5 stayed active, values appeared, and no scoring error was shown.

I clicked Recompute near the end of the pending interval, but completion raced with the click response. The specific Recompute wait-message branch was therefore not independently observed. The automatic false scoring-error condition while energetics were pending was observed and did not recur.

No test suites or regression tests were run. Source fix details are in [W4-source-fix-2026-11-03.md](W4-source-fix-2026-11-03.md).
