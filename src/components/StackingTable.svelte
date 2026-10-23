<script lang="ts">
  import { appState } from '../lib/state.svelte';
  import { setSharedStackingSelection } from '../lib/seamless.svelte';
  import { formatStackingNumber } from '../lib/stackings';
</script>

{#if appState.stackings.length === 0}
  <p class="stacking-empty" role="status">No aromatic stacking interactions are available yet.</p>
{:else}
  <div class="stacking-table-wrap">
    <table class="stacking-table" aria-label="Detected aromatic stacking interactions">
      <caption class="sr-only">Aromatic stacking interactions detected between RNA residues and protein range PDBs.</caption>
      <thead>
        <tr>
          <th scope="col">Selection</th>
          <th scope="col">RNA residue</th>
          <th scope="col">Protein range PDB</th>
          <th scope="col">Protein residue</th>
          <th scope="col">Distance</th>
          <th scope="col">Plane tilt</th>
          <th scope="col">Slip angle</th>
        </tr>
      </thead>
      <tbody>
        {#each appState.stackings as stacking, index (`${stacking.protein_range_pdb_id}:${stacking.protein_resid}:${stacking.rna_resname}:${stacking.rna_resid}`)}
          <tr>
            <td class="stacking-selection-cell">
              <input
                type="checkbox"
                checked={appState.stackingSelection[index] ?? false}
                aria-label={`Select RNA ${stacking.rna_resname} ${stacking.rna_resid} with protein residue ${stacking.protein_resid} in range PDB ${stacking.protein_range_pdb_id}`}
                onchange={event => {
                  const next = [...appState.stackingSelection];
                  next[index] = event.currentTarget.checked;
                  appState.stackingSelection = next;
                  setSharedStackingSelection(next);
                }}
              />
            </td>
            <th scope="row">{stacking.rna_resname} {stacking.rna_resid}</th>
            <td>{stacking.protein_range_pdb_id}</td>
            <td>{stacking.protein_resid}</td>
            <td>{formatStackingNumber(stacking.distance)} Å</td>
            <td>{formatStackingNumber(stacking.tilt_angle)}°</td>
            <td>{formatStackingNumber(stacking.slip_angle)}°</td>
          </tr>
        {/each}
      </tbody>
    </table>
  </div>
{/if}
