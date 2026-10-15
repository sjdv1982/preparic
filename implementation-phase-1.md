# Phase 1 plan: molecular-selection web app (Svelte + Mol*)

## Context
`outline-phase-1.md` describes the first phase of a web app for defining selections on a molecular structure. There is no server yet, and later phases will add more stages. Phase 1 has a camera-and-picking-only Mol* viewer and two stages:
1. **Stage 1:** pick a PDB code and a biological unit (biounit).
2. **Stage 2:** build a table of residue ranges `(chain, start, end)` by typing values or clicking in the viewer.

`/home/agent/preparic` contains only the outline and `molstar-bio-unit.md`, which holds verified Mol* findings and is reused below. This is a greenfield build.

**Decisions made with you:**
- Stack: Svelte 5 (runes) + TypeScript + Vite.
- Numbering: **label** numbering (`label_asym_id`, `label_seq_id`).
- When stage 2 is finished: nothing happens yet. The OK button calls a no-op `onComplete` hook that the next phase will fill in.

**Defaults I chose (easy to change):**
- The range table has no external grid library. Its requirements are too specific to Mol* picking.
- Only **polymer** chains appear in the chain choices. Ligands and waters have no `label_seq_id` range, so they cannot be part of a range.
- Molecules are downloaded from RCSB `https://models.rcsb.org/{id}.bcif`, which has CORS and is what Mol* itself uses.

## Project layout
```
package.json          svelte@5, vite, @sveltejs/vite-plugin-svelte, typescript, svelte-check, vitest, molstar@5.13.x (pinned)
vite.config.ts, tsconfig.json, index.html
src/main.ts
src/App.svelte                     layout: side panel + viewer; renders the current stage
src/lib/state.svelte.ts            app state ($state runes) + stage transitions
src/lib/types.ts                   ChainInfo, RangeRow, ActiveCell, OtherChainsMode, AssemblyInfo
src/lib/chains.ts                  pure: units → named ChainInfo[] (A, A-1, A-2 …)
src/lib/ranges.ts                  pure: validation, start/end swap, row ops (add/remove/move)
src/lib/molstar/viewer.ts          MolViewer class wrapping PluginContext (the only Mol*-aware module)
src/lib/molstar/expressions.ts     MolScript builders (chain copy, residue range, union, complement)
src/components/Viewer.svelte       hosts canvas, owns MolViewer instance
src/components/Stage1Structure.svelte
src/components/Stage2Ranges.svelte
src/components/RangeTable.svelte
tests/chains.test.ts, tests/ranges.test.ts
```

## Mol* viewer (`src/lib/molstar/viewer.ts`)
- **No Mol* UI.** Create `new PluginContext(spec)` from `molstar/lib/mol-plugin/context`, then call `await plugin.init()` and `plugin.initViewerAsync(canvas, container)`, falling back to `initViewer` depending on the installed version. Without the React UI there are no buttons at all; only the camera and picking work.
- **Behaviours.** Start from `DefaultPluginSpec()` and remove the behaviours that act on a click, `FocusLoci` (zoom on click) and `SelectLoci` (default selection), so a click only does what we define. Keep `HighlightLoci` for hover feedback. Check the behaviour names against the installed version.
- **Click handling.** Subscribe to `plugin.behaviors.interaction.click`. Read the clicked loci (`StructureElement.Loci`) and resolve it to `{labelAsymId, instanceId, labelSeqId}` with `StructureProperties` and `unit.conformation.operator.instanceId`. Pass that to the app through an `onPick` callback.
- **Pick granularity.** Call `plugin.managers.interactivity.setProps({ granularity: 'chain' | 'residue' })` depending on the active cell.
- **API used by the app:**
  - **`loadPdb(id)`:** `builders.data.download({url, isBinary:true})` → `parseTrajectory(data,'mmcif')` → `createModel(trajectory)`, with model index 0 for NMR entries. Returns `AssemblyInfo[]` built from `ModelSymmetry.Provider.get(model).assemblies`.
    - For each assembly, the participating polymer chains are collected from `operatorGroups` (asymIds × operator count). The polymer filter uses `model.entities` and the atomic hierarchy.
    - If an entry has no assemblies, the list holds a single "Asymmetric unit" entry mapped to `{name:'model'}`.
    - Uses a request token so that if an older load finishes after a newer one has started, its result is thrown away.
  - **`showAssembly(id)`:** deletes the previous structure subtree and calls `createStructure(model, {name:'assembly', params:{id}})`. It validates the id first, because Mol* silently falls back to the model for unknown ids.
  - **`listChains()`:** walks `structure.units` and collects polymer units as `(label_asym_id, instanceId, isIdentity, label_seq_id min/max/set)`, then passes them to `chains.ts` for naming.
  - **`render(spec)`:** `spec = { focus: ChainKey | null, alwaysVisible: RangeSpec[], others: 'visible'|'transparent'|'invisible' }`.
    - With `focus` null, everything is shown as a polymer cartoon.
    - Otherwise it builds two components with `tryCreateComponentFromExpression`:
      - **focus set:** the focused chain plus the always-visible ranges, drawn as a normal cartoon.
      - **rest:** the complement of the focus set, drawn as a cartoon with `alpha≈0.2` (transparent), or with no representation at all (invisible).
    - Calls are serialized so the latest one wins. Old components are removed through a tagged group.
- **Expressions (`expressions.ts`).** Follow the verified pattern in `molstar-bio-unit.md:97-106`: `atomGroups` with a `chain-test` on `label_asym_id` + `instanceId`, an optional `residue-test` `inRange(label_seq_id, s, e)`, and `exceptBy` for the complement.

## Chain naming (`src/lib/chains.ts`, pure, unit-tested)
- Group units by `label_asym_id`.
- A chain that occurs only once is named `A`.
- When a chain occurs several times, the identity-operator copy is `A` and the other copies are `A-1`, `A-2`, … in unit order.
  - If no copy uses the identity operator, the copies are named `A-1`, `A-2`, ….
- `instanceId` is the stable key. Do not use `operator.name` (see `molstar-bio-unit.md:96`).
- `ChainInfo` holds `{ name, labelAsymId, instanceId, seqIds:Set<number>, seqMin, seqMax }`.

## App state & stages (`state.svelte.ts`)
- **State:** `stage: 1|2`, `pdbId`, `assemblies`, `assemblyId`, `chains`, `rows: RangeRow[]`, `activeCell`, `otherChainsMode` (default `'transparent'`).
- **Stage 1 → 2:** OK is enabled once a structure is shown. It calls `listChains()` and starts with one empty row.
- **Stage 2 → 1 (Back):** clears `rows`, `activeCell` and `chains`, and calls `render({focus:null})`. The PDB code and biounit belong to stage 1, so they are kept.
- **Extensibility:** stages are a simple array/union, so later phases can add stages 3 and up. Each stage has a "reset my data" function that Back calls.

## Stage 1 (`Stage1Structure.svelte`)
- **PDB code field:** a text input checked against `^[0-9][A-Za-z0-9]{3}$`. When the code is valid, `loadPdb` runs after a ~300 ms debounce.
  - Shows a loading state while downloading, and errors such as "entry not found" or a network failure.
- **Biounit menu:** appears after a successful load. Entries are labelled with index and participating chains, e.g. `1 — A, B ×2`, truncated for very large assemblies, with the full list in a tooltip. The first biounit is selected by default and shown. Changing the selection calls `showAssembly`.
- **Edits after loading:** changing the code clears the menu and the viewer.

## Stage 2 (`Stage2Ranges.svelte` + `RangeTable.svelte`)
- **Table:** columns are a row checkbox ("always visible"), Chain, Start and End, plus row controls (move up, move down, remove). An "Add row" button inserts below the active row and copies the active row's chain.
- **Excel-style cells:**
  - A single click activates a cell (`activeCell`). Typing or a double-click edits it; Enter or Tab commits; arrow keys and Tab move between cells; Delete clears.
  - The Chain cell offers a list of the valid chain names. Start and End are integers.
  - Invalid values are marked red: an unknown chain, a residue not in that chain's `seqIds`, or an empty cell when OK is pressed.
- **Swap rule (`ranges.ts`):**
  - When `activeCell` moves to another row, or to null, the row being left is normalised: if `start > end`, the two values swap.
  - A click outside the table counts as moving to null. A click in the viewer does not, because that is how values get picked.
  - OK normalises all rows too.
- **What the active cell does in the viewer:**

  | Active cell | Granularity | `render` call | A click in the viewer |
  |---|---|---|---|
  | Chain cell, or none | `chain` | `render({focus:null})` | sets the row's chain to the clicked chain name |
  | Start/End cell with a valid chain | `residue` | `render({focus: row chain, alwaysVisible: checked rows, others: otherChainsMode})` | sets Start/End to the clicked `label_seq_id`, but only if it is in the focused chain; otherwise ignored |
  | Start/End cell, chain invalid | `chain` | `render({focus:null})` | ignored |

- **Other-chains menu:** offers visible, transparent and invisible, and drives `others`.
- **Always-visible checkbox:** a ticked row with a valid range is added to `alwaysVisible`. It is always drawn normally, even when other chains are transparent or hidden.
- **OK:** validates and normalises, then calls `onComplete(state)`, which is a no-op in phase 1. **Back:** as described in the stage section.

## Verification
1. **Unit tests:** `npx vitest run`.
   - `chains.test.ts`: a single copy; identity plus copies giving `A`, `A-1`, `A-2`; copies with no identity giving `A-1`, `A-2`; several chains mixed.
   - `ranges.test.ts`: the swap happens when leaving the row and not on a viewer pick; validation; add, remove and move rows.
2. **Static checks:** `npx svelte-check` and `npx tsc --noEmit`.
3. **Manual run** with `npm run dev`, driving the browser with the run / claude-in-chrome skill:
   - **2JDQ:** the menu shows 2 biounits, the first is shown by default, and switching works.
   - **2BTV:** biounit 3 (chain A ×5) shows chains `A, A-1…A-4`. Picking a residue on a copy fills Start/End from that copy only. Transparent and invisible modes work, and the always-visible checkbox keeps a range drawn.
   - **1CRN:** a single chain.
   - **An NMR entry, e.g. 1D3Z:** model 0 is used.
   - **An invalid or unknown code:** an error is shown.
   - **Back from stage 2:** the ranges are erased, and the PDB code and biounit are kept.
   - **Swap:** Start=50 and End=20 swap to 20/50 when moving to another row.
