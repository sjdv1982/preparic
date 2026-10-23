# Preparic · Phase 1

Browser app for choosing a PDB biological assembly and defining residue ranges with Mol*.

## Start the app

```sh
mamba env create -f environment.yml
mamba activate preparic
npm install
npm run dev
```

`npm run dev` serves the app only. To run it together with the Seamless backend (`backend.py`), use the scripts in this directory:

```sh
./preparic-install-hook   # once, in the activated env: puts this directory on PATH on activation
preparic-start            # backend + Vite as subprocesses, output to the terminal
preparic-start --force    # also force Vite to rebuild optimized dependencies
preparic-start --log      # same, but output of both goes to preparic.log (--log FILE for another file)
preparic-start --log --shell  # only Vite logs; the backend runs in `ipython -i` (else `python -i`) on the terminal, with `ctx` available. Exit the shell to stop both
preparic-stop             # stops any running preparic-start, and with it both servers
preparic-restart [--log] [--force]  # preparic-stop, then preparic-start with the same options
```

Use `--force` when Vite's optimized dependency cache is stale or inconsistent.

The backend needs the `seamless` and `seamless_workflow` packages, which `environment.yml` does not install.

The viewer needs a browser with WebGL enabled. Use `npm run check`, `npm test`, and `npm run build` for local verification.

For a quick Stage 2 test, open the app with `?pdb=1B7F&continue=1`; it loads the first biological unit and advances after the viewer is ready.

## Tests

Run tests in the `preparic` conda environment (`pytest` is included in
`environment.yml`; update an existing environment if necessary):

```sh
conda activate preparic
npm test                       # tests/frontend: Vitest
npm run test:backend            # tests/backend: pytest, including native Bash tests
python -m pytest -m 'not integration'  # direct Python and Bash checks
python -m pytest -m integration        # standalone Seamless transformations
bash tests/backend/test_range_pdbs_aareduce.sh  # native Bash suite by itself
```

The direct Python tests are unit tests. The Bash suite includes isolated dispatch
tests with a fake reducer, plus a script test using the real aareduce/PDB2PQR tools.
The integration tests build standalone Seamless `Transformer`/`Transformation`
objects with the code, celltypes and `preparic` conda environment used by
`backend.py`. They test individual steps and the extraction-to-reduction chain
without workflow Contexts, mounts, share servers or remote services. Persistent
caching is disabled so previous deployments cannot supply cached results.

Fixtures are checked into `tests/backend/data/1b7f`; tests do not read live
deployment files or download data. One strict expected failure records an
upstream Seamless error when materializing an empty-list transformation result
(`Checksum.resolution` creates a trivial Buffer without its checksum). The
direct empty-input checks and the empty Bash transformation still run normally.
