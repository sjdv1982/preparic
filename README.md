# Preparic · Phase 1

Browser app for choosing a PDB biological assembly and defining residue ranges with Mol*.

## Start the app

```sh
mamba env create -f environment.yml
mamba activate preparic
npm install
npm run dev
```

The viewer needs a browser with WebGL enabled. Use `npm run check`, `npm test`, and `npm run build` for local verification.

For a quick Stage 2 test, open the app with `?pdb=1B7F&continue=1`; it loads the first biological unit and advances after the viewer is ready.
