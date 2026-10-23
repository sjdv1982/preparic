<script lang="ts">
  import { appState } from '../lib/state.svelte';
  import { normalizeRow, validateRows } from '../lib/ranges';
  import RangeTable from './RangeTable.svelte';

  let { onBack, onComplete }: { onBack: () => void; onComplete: () => void } = $props();
  let showAllErrors = $state(false);
  let issues = $derived(validateRows(appState.rows, appState.chains, showAllErrors));

  function finish() {
    appState.rows = appState.rows.map(normalizeRow);
    showAllErrors = true;
    if (validateRows(appState.rows, appState.chains).length > 0) return;
    onComplete();
  }
</script>

<section class="stage-content stage-two-content">
  <div class="eyebrow">Stage 02 <span>Ranges</span></div>
  <h1>Define residue ranges</h1>
  <p class="stage-description">Add ranges by typing residue numbers or picking directly from the structure.</p>

  <div class="viewer-instructions"><span class="instruction-icon">⌖</span><span>Select a cell, then click a chain or residue in the viewer to fill it.</span></div>

  <RangeTable issues={issues} />

  <div class="range-settings">
    <label class="field-label" for="other-chains">Other chains</label>
    <select id="other-chains" class="select-input" bind:value={appState.otherChainsMode}>
      <option value="visible">Visible</option>
      <option value="transparent">Transparent</option>
      <option value="invisible">Invisible</option>
    </select>
    <div class="field-note">Controls chains outside the focused range.</div>
    <button type="button" class="button button-secondary sidechain-toggle" aria-pressed={appState.showSidechainSticks} onclick={() => appState.showSidechainSticks = !appState.showSidechainSticks}>
      <span>Sidechain sticks</span>
      <span class="sidechain-toggle-state">{appState.showSidechainSticks ? 'On' : 'Off'}</span>
    </button>
    <div class="field-note">Show or hide sticks for every defined range.</div>
  </div>

  <div class="stage-footer stage-footer-split">
    <button class="button button-secondary" onclick={onBack}><span aria-hidden="true">←</span> Back</button>
    <button class="button button-primary" onclick={finish}>Complete <span aria-hidden="true">✓</span></button>
  </div>
</section>
