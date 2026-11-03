# W5 S24 editor visibility retest — 2026-11-03

## Scope

Focused Chrome DevTools validation of the W5 defect where the editable energy-threshold table collapsed to about 0 px and could not be used. This checks visible layout, typing and S24 persistence; it does not revalidate Copy, per-cell revert, per-range reset, or scientific energy results.

## Setup and steps

- Ran the disposable Preparic copy at `/tmp/preparic-w5-acceptance.d6GWgV` in the `preparic` conda environment (Vite 5319, share/backend 6019), using Chrome DevTools MCP.
- Entered stage 5 through the UI from the prepared 1B7F stage-3 fixture. Copied fresh stage-4 thresholds and waited for the stage-5 energetics/thresholds to settle.
- Measured `.stage5-energy-wrap` and its S24 input with a read-only DOM probe, then filled the `frag-5-UU / range-2` threshold with `-12.34` through the UI and pressed Tab to commit.

## Result

- The energy wrapper measured **328×881 px** (the stage panel was 694 px high), compared with the prior 328×2 px observation. The S24 input measured **78×24 px**, was enabled and writable, and was within the viewport after the browser scrolled it into view.
- The typed value appeared in the field. S24 PUT request 167 and readback request 168 both returned HTTP 200; the readback contained `frag-5-UU: "-12.34"` for range 2.
- Continue remained disabled because other true cells still had empty thresholds. The completeness message no longer listed the edited range-2 cell.

The disposable service is stopped, its tab is closed, and ports 5319/6019 are free. The disposable mounts were restored to their pre-run snapshot; before/after hash files match. Main project deployment mounts were not used.

No regression suite was run. The original W5 results are in [W5.md](W5.md); evidence and mount hashes for this retest are under [evidence/W5-s24-visibility-20261103](evidence/W5-s24-visibility-20261103/).
