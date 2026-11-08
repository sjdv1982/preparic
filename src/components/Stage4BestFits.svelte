<script lang="ts">
  import { seamlessState, setProjectDecision, setProjectDecisionDebounced, setStage4ActiveNucleotide } from '../lib/seamless.svelte';
  import {
    activeDinucleotidePairs,
    cloneThresholdMarker,
    canContinueStage4,
    fragmentIdentityKey,
    selectedCandidateGroups,
    thresholdsCompleteForClip,
    validateThresholdEdits,
    validateThresholdValue,
    type ActiveDinucleotidePair,
    type ClipMarker,
    type FragmentInventory,
    type ThresholdMarker,
  } from '../lib/stage4';

  let { onBack, onContinue, embedded = false, readOnly = false }: {
    onBack: () => Promise<boolean>; onContinue: () => Promise<boolean>; embedded?: boolean; readOnly?: boolean;
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
      edits = cloneThresholdMarker(marker as ThresholdMarker);
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
  let bestFitsFresh = $derived(seamlessState.computedFreshness.best_fits === true);
  let stackingRows = $derived(seamlessState.computed.stackings);
  let candidatesFresh = $derived(seamlessState.computedFreshness.stackings === true
    && seamlessState.computedFreshness.fragments === true);
  let energeticsFresh = $derived(seamlessState.computedFreshness.energetics === true);
  let selectedKeys = $derived.by(() => {
    const marker = seamlessState.projectSnapshot.stacking_selection;
    if (!record(marker) || !record(marker.selected)) return [] as string[];
    return Object.entries(marker.selected).filter(([, selected]) => selected === true).map(([key]) => key);
  });
  let fragmentCount = $derived(inventory?.fragments?.length ?? 0);
  let candidateGroups = $derived(candidatesFresh && inventory
    ? selectedCandidateGroups(inventory, selectedKeys, stackingRows) : []);
  let activePairs = $derived(activeDinucleotidePairs(inventory, clip));
  // In stage 4 the active nucleotide is the shared S26 value, for collaborative
  // viewing. Embedded in stages 5–6, where S26 has no meaning, the choice is per tab.
  let localActiveKey = $state<string | null>(null);
  let activeKey = $derived.by(() => {
    if (embedded) return localActiveKey;
    const shared = seamlessState.stage4ActiveNucleotide;
    return shared ? fragmentIdentityKey(shared) : null;
  });
  // Null, or a nucleotide without both fragments, shows the first one without writing.
  let activePairIndex = $derived(Math.max(0, activePairs.findIndex(item => item.key === activeKey)));
  let activePair = $derived(activePairs[activePairIndex] ?? null);
  let activeThreshold = $derived.by(() => {
    if (!edits || !activePair) return null;
    const index = edits.pairs.findIndex(row => row.down === activePair.down && row.up === activePair.up);
    return index >= 0 ? { index, row: edits.pairs[index] } : null;
  });
  let activeComputedPair = $derived.by(() => computedFresh && computedThresholds && activePair
    ? computedThresholds.pairs?.find((item: any) => item.down === activePair.down && item.up === activePair.up) ?? null
    : null);
  let canCopyActive = $derived(typeof activeComputedPair?.ovrmsd === 'string'
    && typeof activeComputedPair?.crmsd === 'string'
    && validateThresholdValue(activeComputedPair.ovrmsd)
    && validateThresholdValue(activeComputedPair.crmsd));
  let activeFits = $derived.by(() => {
    const rows = Array.isArray(bestFits?.rows) ? bestFits.rows : [];
    return {
      pre: bestFitsFresh && activePair ? rows.find((item: any) => item.key === activePair.pre.key) : null,
      post: bestFitsFresh && activePair ? rows.find((item: any) => item.key === activePair.post.key) : null,
    };
  });
  let clipLimits = $derived.by(() => {
    const fallback = Math.max(0, fragmentCount - 1);
    if (!candidateGroups.length) return { clip1: fallback, clip2: fallback };
    if (candidateGroups.some(group => group.length === 0)) return { clip1: 0, clip2: 0 };
    return {
      clip1: Math.min(fallback, ...candidateGroups.map(group => Math.max(...group) - 1)),
      clip2: Math.min(fallback, ...candidateGroups.map(group => fragmentCount - Math.min(...group))),
    };
  });
  let validatorFailed = $derived(seamlessState.stateGraphErrors.some(issue => issue.name === 'validate_threshold_edits' && issue.state === 'failed'));
  let complete = $derived(!!edits && thresholdsCompleteForClip(edits, fragmentCount, clip));
  let canContinue = $derived(canContinueStage4({
    edits,
    fragmentCount,
    clip,
    selectedCandidates: candidateGroups,
    thresholdsFresh: computedFresh,
    energeticsFresh,
    thresholdValidatorFailed: validatorFailed,
  }) && !!inventory && validateThresholdEdits(edits, inventory)
    && seamlessState.ready && seamlessState.status === 'live' && seamlessState.projectStage === 4
    && seamlessState.computedFreshness.fragments === true && seamlessState.computedFreshness.stackings === true
    && !seamlessState.pendingInputs.includes('threshold_edits')
    && !seamlessState.pendingInputs.includes('best_fit_clip'));

  function computedValue(row: Record<string, any>, key: string): string | null {
    if (!computedFresh || !computedThresholds) return null;
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
    edits = cloneThresholdMarker(next);
    message = '';
    return true;
  }

  function editValue(group: 'pairs', index: number, key: string, value: string): void {
    if (!edits) return;
    const normalized = value.trim();
    const nextValue = normalized === '' ? null : normalized;
    if (nextValue !== null && !validateThresholdValue(nextValue)) {
      message = 'Enter a finite positive value with at most two decimal places.';
      return;
    }
    const next = cloneThresholdMarker(edits);
    Object.assign(next[group][index], { [key]: nextValue });
    save(next);
  }

  function explainInvalidDraft(value: string): void {
    const normalized = value.trim();
    message = normalized !== '' && !validateThresholdValue(normalized)
      ? 'Enter a finite positive value with at most two decimal places.'
      : '';
  }

  function copyActive(): void {
    if (!edits || !activeThreshold || !canCopyActive || !activeComputedPair) return;
    const source = activeComputedPair;
    if (typeof source?.ovrmsd !== 'string' || typeof source?.crmsd !== 'string'
      || !validateThresholdValue(source.ovrmsd) || !validateThresholdValue(source.crmsd)) return;
    const next = cloneThresholdMarker(edits);
    Object.assign(next.pairs[activeThreshold.index], { ovrmsd: source.ovrmsd, crmsd: source.crmsd });
    save(next);
  }

  function chooseActivePair(index: string | number): void {
    const selectedIndex = Number(index);
    const selected = activePairs[selectedIndex];
    if (!selected) return;
    if (embedded) localActiveKey = selected.key;
    else setStage4ActiveNucleotide(selected.residue);
  }

  function formatRmsd(value: unknown): string {
    return typeof value === 'number' && Number.isFinite(value) ? `${value.toFixed(2)} Å` : '—';
  }

  function residueLabel(residue: ActiveDinucleotidePair['residue']): string {
    return `${residue.chain}:${residue.resid}${residue.insertion_code}`;
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
    if (candidateGroups.some(group => group.length === 0
      || !group.some(number => number > first && number <= fragmentCount - last))) {
      message = 'The clip must keep at least one candidate fragment of every selected stacking.';
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
    <h2 class="stage4-embedded-title">Dinucleotide fit comparison</h2>
    <p class="stage-description">Compare both fits around the active RNA nucleotide.</p>
  {:else}
    <div class="eyebrow">Stage 04 <span>Best fits</span></div>
    <h1>Compare local best fits</h1>
    <p class="stage-description">Choose an RNA nucleotide to compare its pre and post dinucleotide fits and set the overlap RMSD thresholds.</p>
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
        <label>Clip 5′ <input aria-label="Clip fragments from 5 prime end" inputmode="numeric" max={clipLimits.clip1} value={clip1Draft} onchange={event => clip1Draft = event.currentTarget.value} disabled={seamlessState.pendingInputs.includes('best_fit_clip')} /></label>
        <label>Clip 3′ <input aria-label="Clip fragments from 3 prime end" inputmode="numeric" max={clipLimits.clip2} value={clip2Draft} onchange={event => clip2Draft = event.currentTarget.value} disabled={seamlessState.pendingInputs.includes('best_fit_clip')} /></label>
        <button class="button button-secondary" onclick={saveClip} disabled={seamlessState.pendingInputs.includes('best_fit_clip')}>Apply clip</button>
      </div>{/if}
      <span class="stage4-threshold-status">{computedFresh ? 'Computed overlap thresholds ready' : 'Computed thresholds pending'}</span>
    </div>

    {#if activePairs.length === 0}
      <p class="stage4-pending" role="status">This RNA range has no internal nucleotide with both a pre and post dinucleotide fragment.</p>
    {:else if activePair}
      <div class="stage4-active-picker">
        <div class="stage4-active-picker-heading">
          <div><span class="stage4-kicker">Active nucleotide</span><strong>{activePair.residue.base}</strong><span>{residueLabel(activePair.residue)}</span>{#if activePair.clipped}<small class="stage4-clipped-note">Outside clip</small>{/if}</div>
          <label for="stage4-active-nucleotide">Choose nucleotide</label>
          <select id="stage4-active-nucleotide" class="select-input" aria-label="Active nucleotide" value={activePairIndex} onchange={event => chooseActivePair(event.currentTarget.value)}>
            {#each activePairs as item, index (`${item.key}-${index}`)}
              <option value={index}>{item.clipped ? `(${item.residue.chain}:${item.residue.resid}${item.residue.insertion_code} · ${item.residue.base})` : `${item.residue.chain}:${item.residue.resid}${item.residue.insertion_code} · ${item.residue.base}`}</option>
            {/each}
          </select>
        </div>
        <input class="stage4-active-slider" type="range" aria-label="Active nucleotide slider" min="0" max={activePairs.length - 1} step="1" value={activePairIndex} oninput={event => chooseActivePair(event.currentTarget.value)} />
        <div class="stage4-slider-ends"><span>{residueLabel(activePairs[0].residue)}</span><span>{residueLabel(activePairs[activePairs.length - 1].residue)}</span></div>
      </div>

      <div class="stage4-fit-grid" aria-label="Pre and post best-fit fragments">
        {#each [
          { position: 'Pre', fragment: activePair.pre, fit: activeFits.pre, color: '#e879f9' },
          { position: 'Post', fragment: activePair.post, fit: activeFits.post, color: '#2dd4bf' },
        ] as fitCard (fitCard.position)}
          <article class="stage4-fit-card" style={`--fit-color: ${fitCard.color}`}>
            <header><span class="stage4-fit-swatch" aria-hidden="true"></span><span>{fitCard.position} fragment</span><strong>{fitCard.fragment.key}</strong></header>
            <div class="stage4-fit-card-body">
              <div class="stage4-fit-sequence">{fitCard.fragment.sequence}</div>
              <div class="stage4-fit-rmsd"><span>Overall fit RMSD</span><strong>{formatRmsd(fitCard.fit?.rmsd)}</strong></div>
              <div class="stage4-fit-residues">{residueLabel(fitCard.fragment.residues[0])} + {residueLabel(fitCard.fragment.residues[1])}{fitCard.fit ? ` · conformer ${fitCard.fit.conformer}` : ''}</div>
            </div>
          </article>
        {/each}
      </div>

      <section class="stage4-overlap-card" aria-label="Shared nucleotide overlap thresholds">
        <header class="stage4-overlap-header">
          <div><span class="stage4-kicker">Shared overlap · {activePair.residue.base} at {residueLabel(activePair.residue)}</span><h2>frag-{activePair.down} ↔ frag-{activePair.up}</h2></div>
          {#if !readOnly && !activePair.clipped}<button class="button button-secondary" disabled={!canCopyActive || !activeThreshold} onclick={copyActive}>Copy cRMSD / ovRMSD to thresholds</button>{/if}
        </header>
        <div class="stage4-overlap-values">
          <div><span>Computed cRMSD</span><strong>{typeof activeComputedPair?.crmsd === 'string' ? activeComputedPair.crmsd : '—'}{typeof activeComputedPair?.crmsd === 'string' ? ' Å' : ''}</strong></div>
          <div><span>Computed ovRMSD</span><strong>{typeof activeComputedPair?.ovrmsd === 'string' ? activeComputedPair.ovrmsd : '—'}{typeof activeComputedPair?.ovrmsd === 'string' ? ' Å' : ''}</strong></div>
        </div>
        {#if activePair.clipped}
          <p class="stage4-pending" role="status">This nucleotide is outside the clip. The computed overlap values remain visible; threshold inputs and copying are unavailable.</p>
        {:else if activeThreshold}
          <div class="stage4-threshold-fields">
            <label>cRMSD threshold <span>Å</span>{#if readOnly}<output aria-label="Overlap cRMSD threshold">{activeThreshold.row.crmsd ?? '—'}</output>{:else}<input aria-label="Overlap cRMSD threshold" inputmode="decimal" value={activeThreshold.row.crmsd ?? ''} oninput={event => explainInvalidDraft(event.currentTarget.value)} onchange={event => editValue('pairs', activeThreshold!.index, 'crmsd', event.currentTarget.value)} />{/if}
              {#if activeThreshold.row.crmsd !== computedValue(activeThreshold.row as Record<string, any>, 'crmsd') && computedValue(activeThreshold.row as Record<string, any>, 'crmsd') !== null}<span class="stage4-diff" title="Differs from computed">●</span>{/if}
              {#if isWarning(activeThreshold.row.crmsd)}<small class="stage4-warning">Above 10 Å</small>{/if}</label>
            <label>ovRMSD threshold <span>Å</span>{#if readOnly}<output aria-label="Overlap ovRMSD threshold">{activeThreshold.row.ovrmsd ?? '—'}</output>{:else}<input aria-label="Overlap ovRMSD threshold" inputmode="decimal" value={activeThreshold.row.ovrmsd ?? ''} oninput={event => explainInvalidDraft(event.currentTarget.value)} onchange={event => editValue('pairs', activeThreshold!.index, 'ovrmsd', event.currentTarget.value)} />{/if}
              {#if activeThreshold.row.ovrmsd !== computedValue(activeThreshold.row as Record<string, any>, 'ovrmsd') && computedValue(activeThreshold.row as Record<string, any>, 'ovrmsd') !== null}<span class="stage4-diff" title="Differs from computed">●</span>{/if}
              {#if isWarning(activeThreshold.row.ovrmsd)}<small class="stage4-warning">Above 10 Å</small>{/if}</label>
          </div>
        {:else}
          <p class="stage4-error" role="alert">The active fragment overlap is missing from the threshold decision.</p>
        {/if}
      </section>

      {#if !bestFitsFresh}<p class="stage4-pending" role="status">Best-fit structures are still being prepared.</p>{/if}
      {#if !computedFresh}<p class="stage4-pending" role="status">Computing overlap RMSD values… They can be copied after they are ready.</p>{/if}
    {/if}

    {#if message}<p class="stage4-error" role="alert">{message}</p>{/if}
    {#if validatorFailed}<p class="stage4-error" role="alert">Threshold validation failed. Correct the table before continuing.</p>{/if}
    <p class="stage4-gate" aria-live="polite">{complete ? 'Kept pair thresholds are complete.' : 'Fill both RMSD thresholds for every kept pair before continuing.'}</p>
  {/if}

  {#if !embedded}<div class="stage-footer stage-footer-split">
    <button class="button button-secondary" onclick={leaveStage}><span aria-hidden="true">←</span> Back to stackings</button>
    <button class="button button-primary" onclick={continueToEnergetics}
      disabled={!canContinue}>
      Continue to energetics <span aria-hidden="true">→</span>
    </button>
  </div>{/if}
</section>
