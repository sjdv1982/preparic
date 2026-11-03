# W5 stage-6 reactive-loop retest — 2026-11-03

## Scope

Focused browser validation of the W5 failure where Continue to stage 6 accepted the Gaussian marker but the page then entered an endless Svelte reactive update loop. This covers the stage-5→6 transition and initial stage-6 render; it does not validate Gaussian scoring or scientific results.

## Setup and steps

- Ran the disposable Preparic copy at `/tmp/preparic-w5-acceptance.d6GWgV` in the `preparic` conda environment (Vite 5319, share/backend 6019), using Chrome DevTools MCP.
- Resumed from the saved stage-5 fixture. After backend restart, the page initially showed the energetics graph as blocked and disabled Continue while current outputs recomputed. After the graph settled, energy thresholds were current and Continue became enabled.
- Clicked **Continue to Gaussian filters** once.

## Result

- The marker write `PUT /preparic/gaussian_classes?marker=1` (request 142) and its readback (request 143) both returned HTTP 200.
- The page rendered stage 6 with the “Review stacking classes” heading, editable Gaussian class table, anchor validation panels, and stacking scores. The stage-5 page was no longer shown.
- The browser console contained no `effect_update_depth_exceeded` message. It did report eight generic 404 resource-load errors during initial resynchronization; the page recovered, the graph settled, and the stage transition succeeded.
- The component had compiled in the `preparic` environment with zero warnings before this browser check.

The disposable service is stopped, its tab is closed, and ports 5319/6019 are free. The disposable deployment mounts were restored to the pre-run snapshot; the before/after hash files match. Main project deployment mounts were not used.

No regression suite was run. The original broader W5 findings remain in [W5.md](W5.md). Mount snapshots and hashes for this retest are under [evidence/W5-loop-fix-20261103](evidence/W5-loop-fix-20261103/).
