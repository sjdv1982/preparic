# Preparic

Browser app for preparing an Alaric run. One shared project goes through six
stages: PDB entry and biological assembly, residue ranges, stacking selection,
best-fit thresholds and clipping, energetics and the interaction matrix, and the
Gaussian stacking filter. The export stage is not implemented yet. The frontend
is Svelte 5 + Mol*; the backend is a Seamless `Context` (`backend.py`) that runs
preparic's own steps and Alaric's sandbox-tested steps.

What each stage does, its gates and the files it keeps are described in
[STAGES.md](STAGES.md).

## Requirements

- The `preparic` conda environment (`environment.yml`) for the app, the
  backend and the tests. The backend also needs the `seamless` and
  `seamless_workflow` packages, which `environment.yml` does not install.
- An Alaric checkout with its `alaric` conda environment (see Alaric's
  `INSTALLATION.md`). The fitter, reference cutter, best-fit tables,
  gridless energetics and the Gaussian check run in that environment.
- `~/.alaric/fraglib.yaml`, pointing at the fragment library.
- The twelve Gaussian stacking-filter caches, installed once per Alaric
  checkout under `precomputed-filters/gaussian/` (Alaric
  `documents/gaussian-filter.md` §6). A missing or corrupt cache is reported as
  an installation error in stage 6 that names the generator command; preparic
  never builds the caches itself.
- A browser with WebGL enabled.

## Start the app

```sh
mamba env create -f environment.yml
mamba activate preparic
npm install
./preparic-install-hook   # once, in the activated env: puts this directory on PATH on activation
preparic-start            # backend + Vite as subprocesses, output to the terminal
```

Launcher options:

```sh
preparic-start --alaric-dir PATH  # Alaric's alaric/ package directory
preparic-start --force    # also force Vite to rebuild optimized dependencies
preparic-start --log      # output of both servers goes to preparic.log (--log FILE for another file)
preparic-start --log --shell  # only Vite logs; the backend runs in `ipython -i` (else `python -i`) on the terminal, with `ctx` available. Exit the shell to stop both
preparic-stop             # stops any running preparic-start, and with it both servers
preparic-restart [--log] [--force]  # preparic-stop, then preparic-start with the same options
```

All Alaric paths derive from one setting: `--alaric-dir`, else `ALARIC_DIR`,
else `~/alaric/alaric`. It names the `alaric/` package directory, and the
repository root is its parent. The backend never imports the installed
`alaric` package; it attaches Alaric modules or stages Alaric CLIs with their
declared dependency lists inside Seamless transformations.

Use `--force` when Vite's optimized dependency cache is stale or inconsistent.
`npm run dev` alone serves only the frontend, without a backend.

For a quick stage-2 test, open the app with `?pdb=1B7F&continue=1`; it loads
the first biological unit and advances after the viewer is ready. URL
parameters apply only in stage 1 with an empty PDB code.

## Execution model

- **Eight local spawn workers.** The backend starts exactly eight Seamless
  spawn workers on this machine. Every transformation runs in one of them.
  There are no remote workers.
- **No persistent result caching.** The backend loads no Seamless
  configuration and needs no hashserver or database; it starts and runs with
  neither installed. Results are reused only within one live backend.
- **Durable decisions.** Every decision the user makes is a file-authoritative
  mount under `deployment-mounts/` (the table in [STAGES.md](STAGES.md)),
  including the downloaded CIF. A missing or malformed file loads as
  undecided. For a malformed file, the backend prints the diagnostic to its log
  before it resets the file. A stored decision is never migrated silently.
- **Computed files are output only.** The other files in `deployment-mounts/`
  (`stackings.json`, `best-fit-thresholds.json`, `pdbs-aa/`,
  `clipped-fragments/`, `interactions.json`, …) are written by the backend for
  inspection and the later export. They are never read back, so an old file
  never counts as a fresh result.
- **One stage for all tabs.** The stage is derived from the decision files,
  never from computation state. When one tab presses Continue while another tab
  is not editing, the other tab follows to the same stage. Racing edits in
  two tabs are not guarded against.

## Restarting the backend

After `preparic-restart` (or any backend restart), **reload every open
browser tab**. A tab kept open across a backend restart has undefined
behavior.

The restarted backend restores the decisions and the CIF from
`deployment-mounts/`, and the project reopens at the same stage. Everything
else is recomputed cold: from stage 3 on that includes the best fit, the
threshold tables, gridless energetics, Gaussian scores and the Gaussian check.
Loading the fragment library can take minutes. Meanwhile the stage stays
where it was, results show as pending, and Continue stays closed where it
needs a fresh result. Stage 4's threshold table stays editable.

Long-lived folder results (the reference fragments and the best fits) are
also written to temporary write-only backup folders for the backend's
lifetime. They are removed when the backend stops and never adopted on start.

## Tests

Run tests in the `preparic` conda environment (`pytest` is included in
`environment.yml`; update an existing environment if necessary):

```sh
conda activate preparic
npm run check                   # svelte-check / TypeScript
npm test                        # tests/frontend: Vitest
npm run test:backend            # tests/backend: pytest, including native Bash tests
python -m pytest -m 'not integration'  # direct Python and Bash checks
python -m pytest -m integration        # standalone Seamless transformations
bash tests/backend/test_range_pdbs_aareduce.sh  # native Bash suite by itself
npm run build
```

The frontend tests cover the stage derivation, marker validation, freshness,
deferred writes, and each stage's pure logic (`w1-*`, `w3-stage4`,
`w45-stage5`, `w6-anchor-constants`, `w7-*`, and the clip and selection
tests). The backend tests cover the decision loader, the eight-worker spawn
runtime with no services (`test_w1_spawn_runtime.py`), preparation
(`test_w3_*`), stage 4, stage 5 (`test_w45_stage5_acceptance.py`) and
stage 6 (`test_w7_*`).

The integration tests build standalone Seamless `Transformer`/`Transformation`
objects with the code, celltypes and conda environments used by `backend.py`.
They test individual steps and the extraction-to-reduction chain without
workflow Contexts, mounts, share servers or remote services. Persistent
caching is disabled so previous deployments cannot supply cached results.

Fixtures are checked into `tests/backend/data/1b7f`; tests do not read live
deployment files or download data. One strict expected failure records an
upstream Seamless error when materializing an empty-list transformation result
(`Checksum.resolution` creates a trivial Buffer without its checksum). The
direct empty-input checks and the empty Bash transformation still run normally.

Browser acceptance scenarios (two tabs, restart, reconnect) are described in
`tests/W1_BROWSER_ACCEPTANCE.md`, `tests/W3_BROWSER_ACCEPTANCE.md`,
`tests/W7_ACCEPTANCE.md` and `tests/browser/*/README.md`. They are run by hand
with Chrome DevTools.

## Background documents

The contract that defines the stages is in the Alaric repository under
`plans/final/`: `oracle-contract.md` (requirements and their checks),
`state-transitions.md` (the stage model) and `contract-boundaries.md` (the
transformer and shared-cell inventory). `outline-phase-1.md`,
`implementation-phase-1.md`, `stage-3-implementation-plan.md` and
`preparic-on-seamless.md` in this directory are earlier design notes, kept as
history.
