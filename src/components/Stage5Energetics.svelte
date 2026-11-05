<script lang="ts">
  import { seamlessState, setEnergyThresholdEdits, setProjectDecisionAwaited } from '../lib/seamless.svelte';
  import {
    candidateAnchors,
    cloneStage5Decision,
    makeInteractionMatrix,
    recomputeInteractionMatrix,
    stage5CanAdvance,
    stage5RangeIds,
    stage5Rows,
    setStage5UntilNextAnchor,
    stage5RangePalette,
    toggleStage5MatrixCell,
    untilNextAnchorState,
    validateStage5InteractionMatrix,
    validEnergyEdit,
    type Stage5S14,
    type Stage5S24,
    type Stage5Range,
  } from '../lib/stage5';

  let {
    onBack,
    onContinue,
    disabledRangeIds = [],
    hoveredFragmentNumbers = [],
    onRangeVisibilityChange,
    onHoverFragment,
  }: {
    onBack: () => Promise<boolean>;
    onContinue: () => Promise<boolean>;
    disabledRangeIds?: string[];
    hoveredFragmentNumbers?: number[];
    onRangeVisibilityChange?: (ids: string[]) => void;
    onHoverFragment?: (number: number | null) => void;
  } = $props();
  let matrixDraft = $state<Stage5S14 | null>(null);
  let energyEdits = $state<Stage5S24 | null>(null);
  let message = $state('');
  let matrixBusy = $state(false);
  let locallyHoveredFragment = $state<number | null>(null);

  const isRecord = (value: unknown): value is Record<string, any> =>
    value !== null && typeof value === 'object' && !Array.isArray(value);
  const sameAddress = (left: unknown, right: unknown): boolean => isRecord(left) && isRecord(right)
    && left.chain === right.chain && left.start === right.start && left.end === right.end;

  $effect(() => {
    const stored = seamlessState.projectSnapshot.interaction_matrix;
    if (isRecord(stored) && stored.active === true) matrixDraft = cloneStage5Decision(stored as Stage5S14);
  });
  $effect(() => {
    const stored = seamlessState.energyThresholdEdits;
    if (isRecord(stored)) energyEdits = cloneStage5Decision(stored as Stage5S24);
  });

  let fragments = $derived(seamlessState.computed.fragments);
  let rangeMarker = $derived(seamlessState.projectSnapshot.ranges);
  let rangeRows = $derived(isRecord(rangeMarker) && Array.isArray(rangeMarker.rows) ? rangeMarker.rows : []);
  let matrixRanges = $derived(stage5RangeIds(rangeRows) ?? []);
  let displayRanges = $derived(rangeRows.slice(1).map((row: any, index: number) => ({
    row,
    index,
    rangeId: `range-${index + 2}`,
    color: stage5RangePalette[index % stage5RangePalette.length],
  })).filter((item: any) => !disabledRangeIds.includes(item.rangeId)));
  let matrixStructureMatches = $derived(!!matrixDraft
    && JSON.stringify(matrixDraft.ranges) === JSON.stringify(matrixRanges));
  let clip = $derived(seamlessState.projectSnapshot.best_fit_clip);
  let selection = $derived(seamlessState.projectSnapshot.stacking_selection);
  let stackings = $derived(seamlessState.computed.stackings);
  let energetics = $derived(seamlessState.computed.energetics);
  let computedThresholds = $derived(seamlessState.computed.energy_thresholds);
  let rows = $derived(stage5Rows(fragments, clip));
  let anchors = $derived.by(() => rows
    ? candidateAnchors(fragments, selection, stackings, rows) : null);
  let selectedAnchors = $derived(anchors ?? []);
  let energeticsFresh = $derived(seamlessState.computedFreshness.energetics === true);
  let thresholdsFresh = $derived(seamlessState.computedFreshness.energy_thresholds === true);
  let canEdit = $derived(seamlessState.ready && seamlessState.status === 'live' && seamlessState.projectStage >= 5);
  let validatorFailed = $derived(seamlessState.stateGraphErrors.some(issue =>
    issue.name === 'validate_energy_threshold_edits' && issue.state === 'failed'));
  let gaussianBackendFailed = $derived(seamlessState.stateGraphErrors.some(issue =>
    ['detect_stacking_bestfit', 'gaussian_scores_transformer'].includes(issue.name) && issue.state === 'failed'));
  let matrixError = $derived.by(() => {
    if (!matrixDraft || !rows || !anchors) return 'The current clip or selected anchors are not ready.';
    if (!energeticsFresh) return null;
    for (const anchor of anchors) {
      for (const lock of anchor.locks) {
        const pair = rowEnergy(anchor.range, lock.fragment, 'reference');
        const bestfit = rowEnergy(anchor.range, lock.fragment, 'bestfit');
        if (![pair, bestfit].every(value => isRecord(value) && value.error == null
          && typeof value.energy === 'number' && Number.isFinite(value.energy))) {
          return `Selected anchor ${anchor.resid} in ${anchor.range} has a scoring error at ${lock.fragment}. Deselect it in stage 3 to continue.`;
        }
      }
    }
    return validateStage5InteractionMatrix(matrixDraft, rows, anchors, energetics, matrixRanges);
  });
  let matrixUsable = $derived(matrixError === null && energeticsFresh);
  let diffCount = $derived.by(() => {
    if (!energyEdits || !thresholdsFresh || !isRecord(computedThresholds) || !isRecord(computedThresholds.values)) return 0;
    let count = 0;
    for (const entry of energyEdits.ranges) {
      for (const [fragment, value] of Object.entries(entry.values)) {
        const rangeIndex = rangeRows.slice(1).findIndex((row: any) => sameAddress(row, entry.row)) + 2;
        const computed = computedThresholds.values[fragment]?.[`range-${rangeIndex}`];
        if (typeof computed === 'string' && value !== computed) count += 1;
      }
    }
    return count;
  });
  let gate = $derived(stage5CanAdvance(seamlessState.projectSnapshot, {
    ready: seamlessState.ready,
    stage: seamlessState.projectStage,
    energeticsFresh,
    bestFitsFresh: seamlessState.computedFreshness.best_fits === true,
    geometryFresh: seamlessState.computedFreshness.bestfit_stackings === true,
    gaussianFresh: seamlessState.computedFreshness.gaussian_scores === true,
    gaussianBackendFailed,
    validatorFailed,
    fragments,
    rangeRows,
    stackings,
    anchors: anchors ?? undefined,
    energetics,
    energyEdits,
  }));

  function rangeIdFor(row: Stage5Range): string {
    const index = rangeRows.slice(1).findIndex((item: any) => sameAddress(item, row)) + 2;
    return `range-${index}`;
  }

  function thresholdValue(row: Stage5Range, fragment: string): string | null {
    const entry = energyEdits?.ranges.find(item => sameAddress(item.row, row));
    const value = entry?.values[fragment];
    return typeof value === 'string' ? value : null;
  }

  function computedValue(row: Stage5Range, fragment: string): string | null {
    if (!thresholdsFresh || !isRecord(computedThresholds) || !isRecord(computedThresholds.values)) return null;
    const value = computedThresholds.values[fragment]?.[rangeIdFor(row)];
    return typeof value === 'string' ? value : null;
  }

  function saveEdits(next: Stage5S24): void {
    if (!setEnergyThresholdEdits(next)) {
      message = 'The shared energy-threshold table rejected this edit.';
      return;
    }
    energyEdits = cloneStage5Decision(next);
    message = '';
  }

  function editThreshold(row: Stage5Range, fragment: string, raw: string): void {
    if (!energyEdits) return;
    const value = raw.trim() === '' ? null : raw.trim();
    if (value !== null && !validEnergyEdit(value)) {
      message = 'Enter a finite signed energy threshold with at most two decimal places.';
      return;
    }
    const next = cloneStage5Decision(energyEdits);
    const entry = next.ranges.find(item => sameAddress(item.row, row));
    if (!entry) return;
    entry.values[fragment] = value;
    saveEdits(next);
  }

  function revertThreshold(row: Stage5Range, fragment: string): void {
    const value = computedValue(row, fragment);
    if (!energyEdits || value === null) return;
    const next = cloneStage5Decision(energyEdits);
    const entry = next.ranges.find(item => sameAddress(item.row, row));
    if (!entry) return;
    entry.values[fragment] = value;
    saveEdits(next);
  }

  function copyThresholds(): void {
    if (!energyEdits || !thresholdsFresh || !isRecord(computedThresholds) || !isRecord(computedThresholds.values)) return;
    const hasValues = energyEdits.ranges.some(entry => Object.values(entry.values).some(value => value !== null));
    if (hasValues && !window.confirm('Replace the current energy-threshold edits with the computed values?')) return;
    const next = cloneStage5Decision(energyEdits);
    for (const entry of next.ranges) {
      const rangeId = rangeIdFor(entry.row);
      for (const fragment of Object.keys(entry.values)) {
        const computed = computedThresholds.values[fragment]?.[rangeId];
        entry.values[fragment] = typeof computed === 'string' ? computed : null;
      }
    }
    saveEdits(next);
  }

  function rowEnergy(range: string, key: string, kind: 'reference' | 'bestfit'): any {
    if (!isRecord(energetics) || !Array.isArray(energetics.ranges)) return null;
    const data = energetics.ranges.find((entry: any) => entry?.range === range);
    return Array.isArray(data?.rows) ? data.rows.find((entry: any) => entry?.key === key)?.[kind] ?? null : null;
  }

  function energyLabel(value: unknown): string {
    return typeof value === 'number' && Number.isFinite(value) ? value.toFixed(3) : '—';
  }

  function energyColor(value: unknown): string {
    if (typeof value !== 'number' || !Number.isFinite(value)) return 'transparent';
    const capped = Math.max(-20, Math.min(20, value));
    const position = (20 - capped) * 2.5;
    return `linear-gradient(90deg, #324761 ${position}%, #142337 ${position}%)`;
  }

  function setRangeDisabled(rangeId: string, disabled: boolean): void {
    const next = disabled
      ? [...new Set([...disabledRangeIds, rangeId])]
      : disabledRangeIds.filter(item => item !== rangeId);
    onRangeVisibilityChange?.(next);
  }

  function hoverFragment(number: number | null): void {
    locallyHoveredFragment = number;
    onHoverFragment?.(number);
  }

  function isKept(key: string): boolean { return !!rows?.includes(key); }

  function isMatrixCellActive(range: string, fragment: string): boolean {
    if (!matrixDraft || !matrixStructureMatches) return false;
    const row = matrixDraft.rows.indexOf(fragment);
    if (row < 0) return false;
    const column = matrixDraft.ranges.indexOf(range);
    return column >= 0 && matrixDraft.grid[row][column];
  }

  async function persistMatrix(next: Stage5S14): Promise<void> {
    matrixBusy = true;
    if (!await setProjectDecisionAwaited('interaction_matrix', next)) {
      message = 'The shared interaction matrix rejected this change.';
      matrixBusy = false;
      return;
    }
    matrixDraft = cloneStage5Decision(next);
    message = '';
    matrixBusy = false;
  }

  async function recompute(): Promise<void> {
    if (!rows || !anchors || !energeticsFresh || !energetics) {
      const failed = energeticsFresh && isRecord(energetics) ? anchors?.find(anchor => anchor.locks.some(lock => {
        const reference = rowEnergy(anchor.range, lock.fragment, 'reference');
        const bestfit = rowEnergy(anchor.range, lock.fragment, 'bestfit');
        return ![reference, bestfit].every(value => isRecord(value) && value.error == null
          && typeof value.energy === 'number' && Number.isFinite(value.energy));
      })) : undefined;
      message = failed
        ? `Selected anchor ${failed.resid} in ${failed.range} has a scoring error. Deselect it in stage 3 to continue.`
        : 'Wait for fresh energetics and the current selected anchors before recomputing.';
      return;
    }
    const skeleton = makeInteractionMatrix(fragments, clip, rangeRows, selection, stackings);
    const computed = skeleton && recomputeInteractionMatrix(skeleton, anchors, energetics);
    if (!computed) { message = 'The current matrix or scoring results are malformed.'; return; }
    await persistMatrix(computed);
  }

  async function toggleCell(row: number, column: number): Promise<void> {
    if (!matrixDraft || !anchors || !energetics || !energeticsFresh) return;
    const next = toggleStage5MatrixCell(matrixDraft, row, column, anchors, energetics);
    if (!next) {
      message = 'That change would remove an anchor lock, create an unanchored run, or include a scoring-error row.';
      return;
    }
    await persistMatrix(next);
  }

  async function toggleUntilNextAnchor(column: number, enabled: boolean): Promise<void> {
    if (!matrixDraft || !anchors || !energetics || !energeticsFresh) return;
    const next = setStage5UntilNextAnchor(matrixDraft, column, enabled, anchors, energetics);
    if (!next) {
      const state = untilNextAnchorState(matrixDraft, column, anchors);
      message = state.forced
        ? 'This join is required because the selected anchors touch or overlap.'
        : 'That join is unavailable for the next anchor on this protein range.';
      return;
    }
    await persistMatrix(next);
  }

  async function leaveStage(): Promise<void> {
    if (!await onBack()) message = 'The interaction matrix could not be deactivated before leaving stage 5.';
  }

  async function continueToGaussian(): Promise<void> {
    if (!gate.allowed) return;
    if (!await onContinue()) message = 'The shared project could not enter the Gaussian stage.';
  }
</script>

<section class="stage-content stage-five-content">
  <div class="eyebrow">Stage 05 <span>Energetics and interactions</span></div>
  <h1>ATTRACT scores and interactions</h1>
  <p class="stage-description">Negative ATTRACT scores are better. Clipped fragments stay visible; their interactions remain unavailable.</p>

  <div class="stage5-toolbar">
    <button class="button button-secondary" disabled={!thresholdsFresh} onclick={copyThresholds}>Copy computed thresholds</button>
    <span>{thresholdsFresh ? `${diffCount} edited cells differ from computed values` : 'Computed thresholds are updating.'}</span>
    <button class="button button-secondary" onclick={recompute} disabled={matrixBusy}>Recompute matrix</button>
  </div>
  {#if !energeticsFresh}
    <p class="stage5-status" role="status">Waiting for current gridless energetics. Best-fit changes score every range; a single range edit rescores that range.</p>
  {/if}

  <section class="stage5-range-visibility" aria-label="Protein range visibility">
    <strong>Protein ranges</strong>
    {#each rangeRows.slice(1) as proteinRange, rangeIndex (`${proteinRange.chain}:${proteinRange.start}-${proteinRange.end}`)}
      {@const rangeId = `range-${rangeIndex + 2}`}
      {@const color = stage5RangePalette[rangeIndex % stage5RangePalette.length]}
      <label class="stage5-range-toggle"><input type="checkbox" checked={!disabledRangeIds.includes(rangeId)}
        aria-label={`Show ${rangeId}, chain ${proteinRange.chain} ${proteinRange.start} to ${proteinRange.end}`}
        onchange={event => setRangeDisabled(rangeId, !event.currentTarget.checked)} />
        <span class="stage5-range-swatch" style={`background:${color.css}`}></span>
        Range {rangeIndex + 2} · {proteinRange.chain} {proteinRange.start}–{proteinRange.end}
      </label>
    {/each}
  </section>
  <p class="stage5-matrix-note">Interaction columns sit with each protein range. Recompute grows out from anchor locks while either absolute ATTRACT score reaches 3; clips and scoring errors stop growth.</p>
  {#if matrixError}<p class="stage5-error" role="alert">{matrixError}. Recompute when current energetics are ready.</p>{/if}
  {#if !anchors}<p class="stage5-error" role="alert">Selected anchors could not be mapped to the kept fragment rows.</p>{/if}
  <div class="stage5-energy-wrap">
    <table class="stage5-energy-table" aria-label="ATTRACT scores, thresholds, and interaction matrix by protein range">
      <thead>
        <tr><th scope="col" rowspan="2">Fragment</th>{#each displayRanges as item (item.rangeId)}
          {@const anchorsInRange = selectedAnchors.map((anchor, anchorIndex) => ({ anchor, anchorIndex }))
            .filter(entry => entry.anchor.range === item.rangeId)}
          <th scope="colgroup" colspan="4" class="stage5-range-heading" style={`background:${item.color.css};color:#0b111b`}>
            <strong>Range {item.index + 2} · chain {item.row.chain} {item.row.start}–{item.row.end}</strong>
            {#if anchorsInRange.length === 0}<small>No selected anchor locks</small>{/if}
            {#each anchorsInRange as entry, anchorPosition (`${entry.anchor.resid}:${entry.anchorIndex}`)}
              {@const anchor = entry.anchor}
              {@const join = matrixDraft && anchors ? untilNextAnchorState(matrixDraft, entry.anchorIndex, anchors) : null}
              <div class="stage5-anchor-locks">
                <small>Anchor {anchor.resid}: {anchor.locks.map(lock => `${lock.fragment} (nt ${lock.nucleotide})`).join(', ')}</small>
                {#if anchorPosition + 1 < anchorsInRange.length}
                  {@const nextAnchor = anchorsInRange[anchorPosition + 1].anchor}
                  <label class="stage5-until-next"><input type="checkbox"
                    aria-label={`Join anchor ${anchor.resid} to anchor ${nextAnchor.resid} in ${item.rangeId}`}
                    checked={join?.joined || join?.forced || false}
                    disabled={!matrixDraft || !matrixStructureMatches || !matrixUsable || !canEdit || matrixBusy || !!join?.forced}
                    onchange={event => void toggleUntilNextAnchor(entry.anchorIndex, event.currentTarget.checked)} />Until next anchor</label>
                {/if}
              </div>
            {/each}
          </th>
        {/each}</tr>
        <tr>{#each displayRanges as item (`headers-${item.rangeId}`)}
          <th>Reference ATTRACT</th><th>Best-fit ATTRACT</th><th>Threshold</th><th>Interaction</th>
        {/each}</tr>
      </thead>
      <tbody>
        {#each (isRecord(fragments) && Array.isArray(fragments.fragments) ? fragments.fragments : []) as fragment (fragment.key)}
          <tr class:clipped={!isKept(fragment.key)} class:viewer-highlighted={hoveredFragmentNumbers.includes(fragment.number)}
            class:table-highlighted={locallyHoveredFragment === fragment.number}
            onpointerenter={() => hoverFragment(fragment.number)} onpointerleave={() => hoverFragment(null)}>
            <th scope="row">{fragment.key}</th>
            {#each displayRanges as item (item.rangeId)}
              {@const proteinRange = item.row}
              {@const rangeId = item.rangeId}
              {@const matrixColumn = matrixDraft?.ranges.indexOf(rangeId) ?? -1}
              {@const matrixRow = rows?.indexOf(fragment.key) ?? -1}
              {@const reference = rowEnergy(rangeId, fragment.key, 'reference')}
              {@const bestfit = rowEnergy(rangeId, fragment.key, 'bestfit')}
              {@const threshold = thresholdValue(proteinRange, fragment.key)}
              {@const computed = computedValue(proteinRange, fragment.key)}
              <td class="stage5-energy" style={`background:${energyColor(reference?.energy)}`} title={reference?.error ?? 'Capped at ±20; negative is better'}>{reference?.error ? `Error: ${reference.error}` : energyLabel(reference?.energy)}</td>
              <td class="stage5-energy" style={`background:${energyColor(bestfit?.energy)}`} title={bestfit?.error ?? 'Capped at ±20; negative is better'}>{bestfit?.error ? `Error: ${bestfit.error}` : energyLabel(bestfit?.energy)}</td>
              <td class="stage5-threshold-cell" class:inactive={!isMatrixCellActive(rangeId, fragment.key)}>
                <input aria-label={`${fragment.key} ${rangeId} energy threshold`} value={threshold ?? ''} disabled={!canEdit} onchange={event => editThreshold(proteinRange, fragment.key, event.currentTarget.value)} />
                {#if computed !== null}<small class="stage5-computed-threshold">Computed: {computed}</small>{/if}
                {#if computed !== null && threshold !== computed}<span class="stage5-diff" title="Differs from computed">●</span>{/if}
                <button class="stage4-revert" aria-label={`Revert ${fragment.key} ${rangeId} to computed`} disabled={!thresholdsFresh || computed === null || threshold === computed} onclick={() => revertThreshold(proteinRange, fragment.key)}>↶</button>
              </td>
              <td class="stage5-matrix-cell">
                <label><input type="checkbox" aria-label={`${fragment.key} / ${rangeId} interaction`}
                  checked={matrixRow >= 0 && matrixColumn >= 0 && matrixDraft?.grid?.[matrixRow]?.[matrixColumn] === true}
                  disabled={!matrixDraft || !matrixStructureMatches || !matrixUsable || !canEdit || matrixBusy || matrixRow < 0 || matrixColumn < 0}
                  onchange={() => void toggleCell(matrixRow, matrixColumn)} />
                  <span>Interact</span></label>
              </td>
            {/each}
          </tr>
        {/each}
      </tbody>
    </table>
  </div>

  {#if message}<p class="stage5-error" role="alert">{message}</p>{/if}
  {#if validatorFailed}<p class="stage5-error" role="alert">Energy-threshold validation failed. Correct the addressed S24 table before continuing.</p>{/if}
  {#if gaussianBackendFailed}<p class="stage5-error" role="alert">Best-fit stacking geometry or Gaussian scoring failed. Review the stage-6 error before continuing.</p>{/if}
  <p class="stage5-gate" aria-live="polite">{#if gate.allowed}Stage 5 requirements are fresh.{:else}{gate.reasons.join(' ')}{/if}</p>

  <div class="stage-footer stage-footer-split">
    <button class="button button-secondary" onclick={leaveStage}><span aria-hidden="true">←</span> Back to best fits</button>
    <button class="button button-primary" onclick={continueToGaussian} disabled={!gate.allowed || matrixBusy}>Continue to Gaussian filters <span aria-hidden="true">→</span></button>
  </div>
</section>
