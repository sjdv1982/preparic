<script lang="ts">
  import { seamlessState } from '../lib/seamless.svelte';
  import { appState } from '../lib/state.svelte';
  import StackingTable from './StackingTable.svelte';

  let { onBack }: { onBack: () => void } = $props();
</script>

<section class="stage-content stage-three-content">
  <div class="eyebrow">Stage 03 <span>Stackings</span></div>
  <h1>Review aromatic stackings</h1>
  <p class="stage-description">Inspect the RNA–protein aromatic interactions detected for the current residue ranges.</p>

  {#if seamlessState.stackingTransformError}
    <p class="stacking-sync-error" role="alert"><strong>Stacking detection failed.</strong> {seamlessState.stackingTransformError}</p>
  {/if}
  {#if seamlessState.stackingsError}
    <p class="stacking-sync-error" role="alert"><strong>Could not read stackings.</strong> {seamlessState.stackingsError}</p>
  {/if}
  {#if seamlessState.stackingSelectionError}
    <p class="stacking-sync-error" role="alert"><strong>Selection sync issue.</strong> {seamlessState.stackingSelectionError}</p>
  {/if}

  <div class="stacking-summary" aria-live="polite">
    <strong>{appState.stackings.length}</strong>
    <span>{appState.stackings.length === 1 ? 'interaction detected' : 'interactions detected'}</span>
  </div>

  <StackingTable />

  <div class="stage-footer stage-footer-split">
    <button class="button button-secondary" onclick={onBack}><span aria-hidden="true">←</span> Back to ranges</button>
  </div>
</section>
