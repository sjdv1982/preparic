# W4 pending energetics source fix

Date: 2026-11-03 (Europe/Paris)

## Defect

The W4 browser record showed an absent S13 result after backend restart, while the Stage 5 matrix diagnostic called an anchor lock a scoring error. The same message could be produced by pressing Recompute while energies were pending. An absent or stale computed share is not row-level evidence of a ligand scoring error.

## Contract basis

- `plans/final/state-transitions.md` §2.2 distinguishes an absent S13 (404) and a stale/pending result from a fresh result; the Stage 5 absent/computing behavior is to wait for S13.
- `plans/final/implementation-plan.md` §7 specifies that Stage 5 matrix read validation is against fresh energetics, with per-ligand errors kept row-local.
- `plans/final/oracle-contract.md` P5.1–P5.6 says per-ligand errors are displayed on their rows, and the anchor deselection diagnostic applies when a fresh scoring result establishes an error on an anchor lock.

## Change

In `src/components/Stage5Energetics.svelte`, Stage 5 defers its row-error/matrix validation while energetics are not fresh. Recompute now reports an anchor scoring error only after fresh energetics are present; while S13 is absent, stale or pending, it uses the existing wait message. Matrix controls remain unavailable because `matrixUsable` still requires fresh energetics, and the Stage 5 gate continues to explain that current energetics are pending.

## Verification

No tests or regression suite were run, per instruction. Chrome DevTools retest is assigned to the coordinator on a disposable app: reload Stage 5 after backend restart, confirm absent/pending S13 shows the waiting gate without an anchor scoring error, then confirm fresh S13 restores matrix validation. This report records the source change only; browser retest results are not yet available.
