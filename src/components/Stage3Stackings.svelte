<script lang="ts">
  import { seamlessState } from '../lib/seamless.svelte';
  import { appState } from '../lib/state.svelte';
  import { canEnterStage4 } from '../lib/stage4';
  import StackingTable from './StackingTable.svelte';

  let { onBack, onContinue }: { onBack: () => void; onContinue: () => void } = $props();
  let canContinue = $derived(canEnterStage4({
    selection: seamlessState.projectSnapshot.stacking_selection,
    stackings: seamlessState.computed.stackings,
    stackingsFresh: seamlessState.computedFreshness.stackings === true,
    fragmentsFresh: seamlessState.computedFreshness.fragments === true,
    rangePdbsFresh: seamlessState.computedFreshness.range_pdbs_aa === true,
  }));
</script>

<section class="stage-content stage-three-content">
  <div class="eyebrow">Stage 03 <span>Stackings</span></div>
  <h1>Review aromatic stackings</h1>
  <p class="stage-description">Inspect the RNA–protein aromatic interactions detected for the current residue ranges.</p>

  {#if seamlessState.stackingsError}
    <p class="stacking-sync-error" role="alert"><strong>Could not read stackings.</strong> {seamlessState.stackingsError}</p>
  {/if}
  {#if seamlessState.stackingSelectionError}
    <p class="stacking-sync-error" role="alert"><strong>Selection sync issue.</strong> {seamlessState.stackingSelectionError}</p>
  {/if}

  {#if !seamlessState.computedFreshness.stackings
    && (seamlessState.status === 'live' || !seamlessState.hasStackingsSnapshot)}
    <p class="stacking-empty" role="status">Waiting for stacking results that match the current structure and ranges…</p>
  {:else}
    <div class="stacking-summary" aria-live="polite">
      <strong>{appState.stackings.length}</strong>
      <span>{appState.stackings.length === 1 ? 'interaction detected' : 'interactions detected'}</span>
    </div>
    <StackingTable />
  {/if}

  <div class="stage-footer stage-footer-split">
    <button class="button button-secondary" onclick={onBack}><span aria-hidden="true">←</span> Back to ranges</button>
    <button class="button button-primary" onclick={onContinue} disabled={!canContinue}>Continue to best fits <span aria-hidden="true">→</span></button>
  </div>
</section>
