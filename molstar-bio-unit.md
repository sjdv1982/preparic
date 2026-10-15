# Mol* and biological-unit (assembly) selection

## Verdict

**Yes.** Mol* natively supports biological units ("assemblies"). It reads `_pdbx_struct_assembly`, `_pdbx_struct_assembly_gen` and `_pdbx_struct_oper_list` from the mmCIF, lists the assemblies, and builds any one of them on demand.

Checked against `molstar@5.13.1` (type declarations and compiled source), then run headlessly on real mmCIF files:

- **2JDQ**: 2 assemblies, each building a 3-unit structure.
- **2BTV**: 6 assemblies. Assembly 1 expands 15 asymmetric-unit chains with 60 operators into 900 units (about 2.9M atoms).

## How it works

An assembly is **not a separate model**. It is a `Structure` built from the deposited `Model`: for each generator, the chains in `asym_id_list` are copied once per operator.

- Chains keep their `label_asym_id`. Copies are told apart by the unit's symmetry operator (`unit.conformation.operator`, with `.name`, `.instanceId` and `.assembly = {id, operId, operList}`).
- `Assembly` objects expose only `id`, `details` and `operatorGroups`.
- The assembly id is the mmCIF `_pdbx_struct_assembly.id` string (`'1'`, `'2'`, …). Lookup is case-insensitive.

## API, by layer

**1. Enumerate the assemblies of a loaded model**

```ts
import { ModelSymmetry } from 'molstar/lib/mol-model-formats/structure/property/symmetry';
const asms = ModelSymmetry.Provider.get(model)?.assemblies ?? [];
// [{ id: '1', details: 'author_and_software_defined_assembly' }, ...]
```

`Assembly` does not carry `oligomeric_details`, `oligomeric_count` or `method_details`. Read them from the raw category (verified):

```ts
MmcifFormat.is(model.sourceData) && model.sourceData.data.db.pdbx_struct_assembly; // Table: .oligomeric_details.value(i) ...
```

**2. Build one (plugin state, the normal route)**

```ts
const trajectory = await plugin.builders.structure.parseTrajectory(data, 'mmcif');
const model      = await plugin.builders.structure.createModel(trajectory);
const structure  = await plugin.builders.structure.createStructure(model, {
  name: 'assembly', params: { id: '2' }          // or { name: 'model' } | { name: 'auto' }
});
```

- `createStructure` wraps the `StateTransforms.Model.StructureFromModel` transform, whose param is `type: RootStructureDefinition.Params`.
- To switch assemblies on a live structure, update that transform's params; `canAutoUpdate` allows this for `assembly`.
- The same param is available through presets:
  ```ts
  plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default', {
    structure: { name: 'assembly', params: { id: '2' } } });
  ```
- Without a `structure` param, the `default` preset ("Default (Assembly)") shows the **first assembly in file order**, or the model if there are none.

**3. Low level (no plugin, no rendering)**

```ts
const base = Structure.ofModel(model);
const s = await StructureSymmetry.buildAssembly(base, '2').run();   // throws if id unknown
```

**4. Declarative (MolViewSpec)**

The `structure` node takes `type: 'assembly'` and `assembly_id`. In the builder, use `.assemblyStructure({ assembly_id: '2' })`. `null` means the first assembly.

**Not offered**

- `Viewer.loadPdb` / `loadStructureFromUrl` accept only `LoadStructureOptions { representationParams }`, so they cannot select an assembly.
- The stock viewer's URL parameters (`pdb=`, `structure-url=`, …) have no assembly parameter.
- Stock Mol* therefore does not give you a "PDB code → assembly" UI; you build it from layers 1–2.

## Suggested design for your app

1. **PDB code field.** Download the mmCIF (RCSB or PDBe) through the plugin and run `parseTrajectory` and `createModel`.
2. **Assembly dropdown.** Populate it from layer 1, labelling each entry `id: details (oligomeric_details)`. Add a "Deposited asymmetric unit" entry mapped to `{ name: 'model' }`.
3. **On change.** Rebuild or update the `StructureFromModel` node.

Mol* does not have to be loaded first just to list assemblies. Lightweight JSON APIs are an alternative:

- RCSB: `https://data.rcsb.org/rest/v1/core/entry/{ID}` → `rcsb_entry_container_identifiers.assembly_ids` (verified: `['1','2']` for 2JDQ).
- PDBe: `https://www.ebi.ac.uk/pdbe/api/pdb/entry/assembly/{id}` (verified; returns entities, composition and `details` per assembly).

Browser CORS behaviour of these two was not tested.

## Chains occurring several times (no `A-1`/`A-2`)

In memory every copy keeps the plain id `A`. Each copy is a separate *unit*, told apart only by its operator. 2BTV assembly 3 (chain A, 5 copies):

| unit id | `label_asym_id` | `operator.name` | `instance_id` | export `suffix` |
|---|---|---|---|---|
| 0 | A | ASM_1 | ASM-1 | *(empty)* |
| 15 | A | ASM_2 | ASM-2 | `_2` |
| 30 | A | ASM_3 | ASM-3 | `_3` |

- **`A_2` naming exists only on export.** `to_mmCIF` writes `label_asym_id`/`auth_asym_id` as `A`, `A_2`, `A_3`, … (identity operator unsuffixed), drops the assembly categories, and adds `_molstar_atom_site_operator_mapping` to record the mapping (verified).
- **Use `instance_id` as the stable key.** `operator.name` is a running index, but with product operators (2BTV assembly 6) the 29th copy of `A` has name `ASM_29` and instance id `ASM-X0-54`.
- **Selecting one copy, or a residue range in it** (verified on 2BTV assembly 3; chain A is 845 residues, `label_seq_id` 57–901):
  ```ts
  const P = MS.struct.atomProperty;
  MS.struct.generator.atomGroups({
    'chain-test': MS.core.logic.and([
      MS.core.rel.eq([P.macromolecular.label_asym_id(), 'A']),
      MS.core.rel.eq([P.core.instanceId(), 'ASM-2'])]),
    'residue-test': MS.core.rel.inRange([P.macromolecular.label_seq_id(), 100, 110]), // omit for whole chain
  });
  ```
  Results: all copies → 5 units / 34120 atoms; `ASM-2` → 1 unit / 6824; `ASM-2` + residues 100–110 → 1 unit / 101 atoms (11 residues). Without the instance test, the range hits all 5 copies (505 atoms). `auth_seq_id` gave the same atoms here; add `pdbx_PDB_ins_code` if the range crosses insertion codes.
- **Applying it in the plugin** (signatures checked, not run): `Script.getStructureSelection(expr, structure)` → `StructureSelection.toLociWithSourceUnits` → `plugin.managers.structure.selection.fromLoci('set', loci)`. For a persistent component use `plugin.builders.structure.tryCreateComponentFromExpression(structureRef, expr, key)`.
- **For your app,** identify a copy as `(label_asym_id, instance_id)`. Display it as `A_2` if you want the export convention.

## Other gotchas for selections

- **Store selections against the model's atoms, not assembly instances,** if you want them to survive an assembly switch. Chain/residue-based MolScript queries are re-evaluated on the new structure; instance-level constraints are not portable between assemblies.
- **Unknown id falls back silently** to the plain model (Mol* only logs a warning), so validate the id against the enumerated list.
- **Single model only.** `buildAssembly` throws for multi-model structures (NMR). Select a `modelIndex` in `ModelFromTrajectory` first.
- **Assemblies can be very large or software-predicted.** Show `details` and the count, and expect heavy rendering. Consider `mergeBySymmetry: true` (a `StructureFromModel` param) and an LOD representation preset (`auto-lod`, `coarse-surface`).
- **No `'deposited'` id** exists in this version. Use `{ name: 'model' }` for the asymmetric unit.
- **Pre-expanded files** such as `https://files.rcsb.org/download/2JDQ-assembly1.cif` exist (HTTP 200 verified), but they rename chains and lose the mapping to the deposited model. Prefer on-the-fly assembly building.

## Verification and sources

- Source files read: `mol-plugin-state/helpers/root-structure`, `mol-model-formats/structure/property/assembly`, `mol-model/structure/structure/symmetry`, `mol-plugin-state/builder/structure*`, `extensions/plugin/loaders`, `extensions/mvs/tree/mvs/*`.
- The runtime check used Node scripts against RCSB mmCIFs.
