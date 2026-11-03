<script lang="ts">
  import { seamlessState, setProjectDecision, setProjectDecisionDebounced } from '../lib/seamless.svelte';
  import { appState } from '../lib/state.svelte';
  import {
    copyComputedThresholds,
    canContinueStage4,
    selectedCandidateFragments,
    thresholdsCompleteForClip,
    validateThresholdEdits,
    validateThresholdValue,
    type ClipMarker,
    type FragmentInventory,
    type ThresholdMarker,
  } from '../lib/stage4';

  let { onBack, onContinue, embedded = false }: {
    onBack: () => Promise<boolean>; onContinue: () => Promise<boolean>; embedded?: boolean;
  } = $props();
  let edits = $state<ThresholdMarker | null>(null);
  let message = $state('');
  let clip1Draft = $state('0');
  let clip2Draft = $state('0');

  const record = (value: unknown): value is Record<string, any> =>
    value !== null && typeof value === 'object' && !Array.isArray(value);

  $effect(() => {
    const marker = seamlessState.projectSnapshot.threshold_edits;
    if (record(marker) && marker.active === true) {
      edits = structuredClone(marker) as ThresholdMarker;
    }
  });

  $effect(() => {
    const clip = seamlessState.projectSnapshot.best_fit_clip;
    if (record(clip)) {
      clip1Draft = String(clip.clip1 ?? 0);
      clip2Draft = String(clip.clip2 ?? 0);
    }
  });

  let inventory = $derived(seamlessState.computed.fragments as FragmentInventory | undefined);
  let computedThresholds = $derived(seamlessState.computed.best_fit_thresholds as Record<string, any> | undefined);
  let bestFits = $derived(seamlessState.computed.best_fits as Record<string, any> | undefined);
  let clip = $derived(seamlessState.projectSnapshot.best_fit_clip as ClipMarker);
  let computedFresh = $derived(seamlessState.computedFreshness.best_fit_thresholds === true);
  let stackingRows = $derived(seamlessState.computed.stackings);
  let candidatesFresh = $derived(seamlessState.computedFreshness.stackings === true
    && seamlessState.computedFreshness.fragments === true);
  let energeticsFresh = $derived(seamlessState.computedFreshness.energetics === true);
  let selectedKeys = $derived.by(() => {
    const marker = seamlessState.projectSnapshot.stacking_selection;
    if (!record(marker) || !record(marker.selected)) return [] as string[];
    return Object.entries(marker.selected).filter(([, selected]) => selected === true).map(([key]) => key);
  });
  let candidates = $derived(candidatesFresh && inventory
    ? selectedCandidateFragments(inventory, selectedKeys, stackingRows) : []);
  let validatorFailed = $derived(seamlessState.stateGraphErrors.some(issue => issue.name === 'validate_threshold_edits' && issue.state === 'failed'));
  let fragmentCount = $derived(inventory?.fragments?.length ?? 0);
  let complete = $derived(!!edits && thresholdsCompleteForClip(edits, fragmentCount, clip));
  let canContinue = $derived(canContinueStage4({
    edits,
    fragmentCount,
    clip,
    selectedCandidates: candidates,
    energeticsFresh,
    thresholdValidatorFailed: validatorFailed,
  }) && !!inventory && validateThresholdEdits(edits, inventory)
    && seamlessState.ready && seamlessState.status === 'live' && seamlessState.projectStage === 4
    && seamlessState.computedFreshness.fragments === true && seamlessState.computedFreshness.stackings === true
    && !seamlessState.pendingInputs.includes('threshold_edits')
    && !seamlessState.pendingInputs.includes('best_fit_clip'));
  let diffCount = $derived.by(() => {
    if (!edits || !computedFresh || !computedThresholds) return 0;
    let count = 0;
    for (const row of edits.pairs) {
      const source = computedThresholds.pairs?.find((item: any) => item.down === row.down && item.up === row.up);
      if (source && (row.ovrmsd !== source.ovrmsd || row.crmsd !== source.crmsd)) {
        if (row.ovrmsd !== source.ovrmsd) count += 1;
        if (row.crmsd !== source.crmsd) count += 1;
      }
    }
    for (const row of edits.bases) {
      const source = computedThresholds.bases?.find((item: any) => item.fragment === row.fragment);
      if (source) {
        if (row.nuc1_ovrmsd !== source.nuc1_ovrmsd) count += 1;
        if (row.nuc2_ovrmsd !== source.nuc2_ovrmsd) count += 1;
      }
    }
    return count;
  });

  function isKept(fragment: number): boolean {
    return fragment > (clip?.clip1 ?? 0) && fragment <= fragmentCount - (clip?.clip2 ?? 0);
  }

  function computedValue(row: Record<string, any>, key: string): string | null {
    if (!computedFresh || !computedThresholds) return null;
    if ('fragment' in row) {
      const source = computedThresholds.bases?.find((item: any) => item.fragment === row.fragment);
      return typeof source?.[key] === 'string' ? source[key] : null;
    }
    const source = computedThresholds.pairs?.find((item: any) => item.down === row.down && item.up === row.up);
    return typeof source?.[key] === 'string' ? source[key] : null;
  }

  function save(next: ThresholdMarker): boolean {
    if (!inventory || !validateThresholdEdits(next, inventory)) {
      message = 'The threshold edit does not match the current fragment list.';
      return false;
    }
    if (!setProjectDecision('threshold_edits', next)) {
      message = 'The shared threshold table rejected this edit.';
      return false;
    }
    edits = structuredClone(next);
    message = '';
    return true;
  }

  function editValue(group: 'pairs' | 'bases', index: number, key: string, value: string): void {
    if (!edits) return;
    const normalized = value.trim();
    const nextValue = normalized === '' ? null : normalized;
    if (nextValue !== null && !validateThresholdValue(nextValue)) {
      message = 'Enter a finite positive value with at most two decimal places.';
      return;
    }
    const next = structuredClone(edits);
    Object.assign(next[group][index], { [key]: nextValue });
    save(next);
  }

  function revertCell(group: 'pairs' | 'bases', index: number, key: string): void {
    if (!edits || !computedFresh) return;
    const source = computedValue(edits[group][index] as Record<string, any>, key);
    if (source === null) return;
    const next = structuredClone(edits);
    Object.assign(next[group][index], { [key]: source });
    save(next);
  }

  function copyAll(): void {
    if (!edits || !computedFresh || !computedThresholds) return;
    const hasUserValues = [...edits.pairs.flatMap(row => [row.ovrmsd, row.crmsd]),
      ...edits.bases.flatMap(row => [row.nuc1_ovrmsd, row.nuc2_ovrmsd])].some(value => value !== null);
    if (hasUserValues && !window.confirm('Replace the current threshold edits with the computed values?')) return;
    const copied = copyComputedThresholds(edits, computedThresholds, true);
    save(copied);
  }

  function saveClip(): void {
    if (seamlessState.pendingInputs.includes('best_fit_clip')) {
      message = 'Wait for the current clip change to finish before applying another.';
      return;
    }
    if (!candidatesFresh) {
      message = 'Wait for fresh stacking and fragment results before changing the clip.';
      return;
    }
    const first = Number(clip1Draft);
    const last = Number(clip2Draft);
    if (!Number.isSafeInteger(first) || !Number.isSafeInteger(last) || first < 0 || last < 0
      || first + last >= fragmentCount) {
      message = 'Clip counts must be nonnegative integers and leave at least one fragment.';
      return;
    }
    if (candidates.some(number => number <= first || number > fragmentCount - last)) {
      message = 'The clip must keep every candidate fragment of each selected stacking.';
      return;
    }
    const basis = inventory?.echo;
    if (!basis) return;
    const next: ClipMarker = { schema: 1, basis, clip1: first, clip2: last };
    message = 'Applying the clip change…';
    void setProjectDecisionDebounced('best_fit_clip', next, 500).then(accepted => {
      message = accepted ? '' : 'The clip change was dropped because the shared project changed. Review the accepted clip and try again.';
    });
  }

  async function leaveStage(): Promise<void> {
    if (!await onBack()) message = 'The threshold table could not be deactivated before leaving stage 4.';
  }

  async function continueToEnergetics(): Promise<void> {
    message = '';
    if (!await onContinue()) message = 'The shared project could not enter stage 5. Wait for fresh best-fit and range results, then try again.';
  }

  function isWarning(value: string | null): boolean {
    return value !== null && Number(value) > 10;
  }
</script>

<section class="stage-content stage-four-content" class:stage4-embedded={embedded}>
  {#if embedded}
    <h2 class="stage4-embedded-title">Fit RMSD thresholds</h2>
    <p class="stage-description">Edit the RMSD limits used for the kept fragments.</p>
  {:else}
    <div class="eyebrow">Stage 04 <span>Best fits</span></div>
    <h1>Set fit thresholds</h1>
    <p class="stage-description">Review each dinucleotide's best fit and set the RMSD limits used for the kept fragments.</p>
  {/if}

  {#if seamlessState.stateGraphErrors.length}
    {#each seamlessState.stateGraphErrors as issue (`${issue.name}:${issue.state}`)}
      {#if ['best_fit', 'refe_fragments', 'best_fit_tables', 'best_fit_thresholds', 'parse_best_fits', 'validate_threshold_edits', 'clip_fragments'].includes(issue.name)}
        <p class="stage4-error" role="alert"><strong>{issue.name} {issue.state}.</strong> {issue.message}</p>
      {/if}
    {/each}
  {/if}

  {#if seamlessState.computedFreshness.fragments !== true || !inventory}
    <p class="stage4-pending" role="status">Preparing reference fragments…</p>
  {:else if !edits}
    <p class="stage4-error" role="alert">The threshold decision is unavailable. Return to stackings and continue again.</p>
  {:else}
    <div class="stage4-toolbar">
      {#if !embedded}<div class="stage4-clip-controls">
        <label>Clip 5′ <input aria-label="Clip fragments from 5 prime end" inputmode="numeric" value={clip1Draft} onchange={event => clip1Draft = event.currentTarget.value} disabled={seamlessState.pendingInputs.includes('best_fit_clip')} /></label>
        <label>Clip 3′ <input aria-label="Clip fragments from 3 prime end" inputmode="numeric" value={clip2Draft} onchange={event => clip2Draft = event.currentTarget.value} disabled={seamlessState.pendingInputs.includes('best_fit_clip')} /></label>
        <button class="button button-secondary" onclick={saveClip} disabled={seamlessState.pendingInputs.includes('best_fit_clip')}>Apply clip</button>
      </div>{/if}
      <div class="stage4-copy-controls">
        <span>{computedFresh ? `${diffCount} cells differ from computed` : 'Computed thresholds pending'}</span>
        <button class="button button-secondary" disabled={!computedFresh} onclick={copyAll}>Copy computed values</button>
      </div>
    </div>

    {#if !computedFresh}
      <p class="stage4-pending" role="status">Computing reference thresholds… Copy is unavailable until fresh values arrive.</p>
    {/if}

    <div class="stage4-table-wrap">
      <table class="stage4-table" aria-label="Editable fragment and pair RMSD thresholds">
        <thead><tr><th>Fragment or pair</th><th>Sequence</th><th>Residues / pose</th><th>Fit RMSD</th><th>Computed limits</th><th>First limit</th><th>Second limit</th></tr></thead>
        <tbody>
          {#each edits.bases as row, index (`base-${row.fragment}`)}
            {@const comp1 = computedValue(row, 'nuc1_ovrmsd')}
            {@const comp2 = computedValue(row, 'nuc2_ovrmsd')}
            {@const fragmentInfo = inventory.fragments.find(item => item.number === row.fragment)}
            {@const fitInfo = seamlessState.computedFreshness.best_fits === true ? bestFits?.rows?.find((item: any) => item.key === fragmentInfo?.key) : null}
            {@const residueText = fragmentInfo?.residues?.map(item => `${item.chain}:${item.resid}${item.insertion_code}`).join(' / ') ?? '—'}
            <tr class:clipped={!isKept(row.fragment)} onmouseenter={() => appState.stage4OverlayFragments = [row.fragment]} onmouseleave={() => appState.stage4OverlayFragments = []}>
              <th scope="row">frag-{row.fragment}</th><td>{row.sequence}</td><td>{residueText}{fitInfo ? ` · conformer ${fitInfo.conformer}` : ''}</td>
              <td>{fitInfo && Array.isArray(fitInfo.grid) ? Number(fitInfo.rmsd ?? NaN).toFixed(2) + ' Å' : '—'}</td>
              <td>{comp1 ?? '—'} / {comp2 ?? '—'}</td>
              <td><input aria-label={`Fragment ${row.fragment} nucleotide 1 ovRMSD`} value={row.nuc1_ovrmsd ?? ''} onchange={event => editValue('bases', index, 'nuc1_ovrmsd', event.currentTarget.value)} />
                {#if row.nuc1_ovrmsd !== comp1 && comp1 !== null}<span class="stage4-diff" title="Differs from computed">●</span>{/if}
                {#if isWarning(row.nuc1_ovrmsd)}<small class="stage4-warning">Above 10 Å</small>{/if}
                <button class="stage4-revert" aria-label={`Revert fragment ${row.fragment} nucleotide 1 to computed`} disabled={!computedFresh || row.nuc1_ovrmsd === comp1} onclick={() => revertCell('bases', index, 'nuc1_ovrmsd')}>↶</button></td>
              <td><input aria-label={`Fragment ${row.fragment} nucleotide 2 ovRMSD`} value={row.nuc2_ovrmsd ?? ''} onchange={event => editValue('bases', index, 'nuc2_ovrmsd', event.currentTarget.value)} />
                {#if row.nuc2_ovrmsd !== comp2 && comp2 !== null}<span class="stage4-diff" title="Differs from computed">●</span>{/if}
                {#if isWarning(row.nuc2_ovrmsd)}<small class="stage4-warning">Above 10 Å</small>{/if}
                <button class="stage4-revert" aria-label={`Revert fragment ${row.fragment} nucleotide 2 to computed`} disabled={!computedFresh || row.nuc2_ovrmsd === comp2} onclick={() => revertCell('bases', index, 'nuc2_ovrmsd')}>↶</button></td>
            </tr>
          {/each}
          {#each edits.pairs as row, index (`pair-${row.down}-${row.up}`)}
            {@const comp1 = computedValue(row, 'ovrmsd')}
            {@const comp2 = computedValue(row, 'crmsd')}
            <tr class:clipped={!isKept(row.down) || !isKept(row.up)} onmouseenter={() => appState.stage4OverlayFragments = [row.down, row.up]} onmouseleave={() => appState.stage4OverlayFragments = []}>
              <th scope="row">{row.down} → {row.up}</th><td>{row.sequence}</td><td>Dinucleotide pair</td><td>—</td>
              <td>{comp1 ?? '—'} / {comp2 ?? '—'}</td>
              <td><input aria-label={`Pair ${row.down}-${row.up} ovRMSD`} value={row.ovrmsd ?? ''} onchange={event => editValue('pairs', index, 'ovrmsd', event.currentTarget.value)} />
                {#if row.ovrmsd !== comp1 && comp1 !== null}<span class="stage4-diff" title="Differs from computed">●</span>{/if}
                {#if isWarning(row.ovrmsd)}<small class="stage4-warning">Above 10 Å</small>{/if}
                <button class="stage4-revert" aria-label={`Revert pair ${row.down}-${row.up} ovRMSD to computed`} disabled={!computedFresh || row.ovrmsd === comp1} onclick={() => revertCell('pairs', index, 'ovrmsd')}>↶</button></td>
              <td><input aria-label={`Pair ${row.down}-${row.up} cRMSD`} value={row.crmsd ?? ''} onchange={event => editValue('pairs', index, 'crmsd', event.currentTarget.value)} />
                {#if row.crmsd !== comp2 && comp2 !== null}<span class="stage4-diff" title="Differs from computed">●</span>{/if}
                {#if isWarning(row.crmsd)}<small class="stage4-warning">Above 10 Å</small>{/if}
                <button class="stage4-revert" aria-label={`Revert pair ${row.down}-${row.up} cRMSD to computed`} disabled={!computedFresh || row.crmsd === comp2} onclick={() => revertCell('pairs', index, 'crmsd')}>↶</button></td>
            </tr>
          {/each}
        </tbody>
      </table>
    </div>
    {#if message}<p class="stage4-error" role="alert">{message}</p>{/if}
    {#if validatorFailed}<p class="stage4-error" role="alert">Threshold validation failed. Correct the table before continuing.</p>{/if}
    <p class="stage4-gate" aria-live="polite">{complete ? 'Kept fragment thresholds are complete.' : 'Fill every threshold for kept fragments and pairs before continuing.'}</p>
  {/if}

  {#if !embedded}<div class="stage-footer stage-footer-split">
    <button class="button button-secondary" onclick={leaveStage}><span aria-hidden="true">←</span> Back to stackings</button>
    <button class="button button-primary" onclick={continueToEnergetics}
      disabled={!canContinue}>
      Continue to energetics <span aria-hidden="true">→</span>
    </button>
  </div>{/if}
</section>
