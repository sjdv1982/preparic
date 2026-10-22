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
preparic-start --log      # same, but output of both goes to preparic.log (--log FILE for another file)
preparic-start --log --shell  # only Vite logs; the backend runs in `ipython -i` (else `python -i`) on the terminal, with `ctx` available. Exit the shell to stop both
preparic-stop             # stops any running preparic-start, and with it both servers
preparic-restart [--log]  # preparic-stop, then preparic-start with the same options
```

The backend needs the `seamless` and `seamless_workflow` packages, which `environment.yml` does not install.

The viewer needs a browser with WebGL enabled. Use `npm run check`, `npm test`, and `npm run build` for local verification.

For a quick Stage 2 test, open the app with `?pdb=1B7F&continue=1`; it loads the first biological unit and advances after the viewer is ready.
