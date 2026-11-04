<script lang="ts">
  import { untrack } from 'svelte';
  import Stage4BestFits from './Stage4BestFits.svelte';
  import {
    seamlessState,
    deactivateProjectDecisionAndWait,
    setEnergyThresholdEdits,
    setProjectDecisionDebounced,
  } from '../lib/seamless.svelte';
  import { cloneStage5Decision, stage5Rows, validEnergyEdit, type Stage5S24, type Stage5Range } from '../lib/stage5';
  import {
    editableGaussianRows,
    automaticGaussianPrefilters,
    copyComputedGaussianClasses,
    gaussianCandidateRows,
    gaussianClassDiffCount,
    gaussianFilterName,
    gaussianReferenceRows,
    gaussianRrmOptions,
    reconcileGaussianClassDraft,
    revertGaussianClassRow,
    selectedGaussianKeys,
    setGaussianFlag,
    switchGaussianClass,
    type GaussianClassMarker,
    type GaussianClassRow,
    type GaussianRrm,
  } from '../lib/stage6';

  let { onBack }: { onBack: () => Promise<boolean> } = $props();
  let classDraft = $state<GaussianClassMarker | null>(null);
  let energyEdits = $state<Stage5S24 | null>(null);
  let message = $state('');
  let classWriteRevision = 0;
  let initialFilterVersion = $state<string | null>(null);
  let modelChanged = $state(false);

  const isRecord = (value: unknown): value is Record<string, any> =>
    value !== null && typeof value === 'object' && !Array.isArray(value);
  const sameAddress = (left: unknown, right: unknown): boolean => isRecord(left) && isRecord(right)
    && left.chain === right.chain && left.start === right.start && left.end === right.end;

  $effect(() => {
    const marker = seamlessState.projectSnapshot.gaussian_classes;
    const writePending = seamlessState.pendingInputs.includes('gaussian_classes');
    // Reconcile shared updates without subscribing to the draft we replace.
    classDraft = reconcileGaussianClassDraft(marker, untrack(() => classDraft), writePending);
  });

  $effect(() => {
    const stored = seamlessState.energyThresholdEdits;
    if (isRecord(stored)) energyEdits = cloneStage5Decision(stored as Stage5S24);
  });

  let scores = $derived(seamlessState.computed.gaussian_scores);
  let validation = $derived(seamlessState.computed.gaussian_validation);
  let selected = $derived(selectedGaussianKeys(seamlessState.projectSnapshot.stacking_selection));
  let computedFresh = $derived(seamlessState.computedFreshness.gaussian_scores === true);
  let geometryFresh = $derived(seamlessState.computedFreshness.bestfit_stackings === true);
  let validationFresh = $derived(seamlessState.computedFreshness.gaussian_validation === true);
  let filterNames = $derived(isRecord(scores) && Array.isArray(scores.filter_names)
    ? scores.filter_names.filter((name: unknown): name is string => typeof name === 'string') : []);
  let rrmOptions = $derived(gaussianRrmOptions(filterNames));
  let computedClasses = $derived(isRecord(scores) && Array.isArray(scores.classes) ? scores.classes : []);
  let classDiffCount = $derived(gaussianClassDiffCount(classDraft, scores));
  let candidateScores = $derived(gaussianCandidateRows(scores, selected));
  let referenceScores = $derived(gaussianReferenceRows(scores, selected));
  let rangeMarker = $derived(seamlessState.projectSnapshot.ranges);
  let rangeRows = $derived(isRecord(rangeMarker) && Array.isArray(rangeMarker.rows) ? rangeMarker.rows : []);
  let fragments = $derived(seamlessState.computed.fragments);
  let clip = $derived(seamlessState.projectSnapshot.best_fit_clip);
  let keptRows = $derived(stage5Rows(fragments, clip));
  let computedThresholds = $derived(seamlessState.computed.energy_thresholds);
  let thresholdsFresh = $derived(seamlessState.computedFreshness.energy_thresholds === true);
  let canEdit = $derived(seamlessState.ready && seamlessState.status === 'live' && seamlessState.projectStage === 6);
  let filterVersion = $derived(isRecord(scores) && typeof scores.filter_set_version === 'string'
    ? scores.filter_set_version : null);
  let classRows = $derived(classDraft?.rows ?? []);
  let clipSummary = $derived(isRecord(clip) ? `${clip.clip1} 5′ / ${clip.clip2} 3′ fragments` : 'pending');

  $effect(() => {
    if (!computedFresh || !filterVersion) return;
    if (initialFilterVersion === null) initialFilterVersion = filterVersion;
    else if (initialFilterVersion !== filterVersion) modelChanged = true;
  });

  function classComputed(key: string): GaussianClassRow | null {
    const value = computedClasses.find((row: any) => row?.key === key);
    return editableGaussianRows(value === undefined ? [] : [value])?.[0] ?? null;
  }

  function updateClass(key: string, update: (row: GaussianClassRow, computed: unknown) => GaussianClassRow): void {
    if (!classDraft || !canEdit || !computedFresh) return;
    const current = classDraft.rows.find(row => row.key === key);
    const computed = computedClasses.find((row: any) => row?.key === key);
    if (!current || !computed || !classComputed(key)) return;
    const nextRow = update(current, computed);
    const next: GaussianClassMarker = {
      ...classDraft,
      rows: classDraft.rows.map(row => row.key === key ? nextRow : row),
    };
    classDraft = next;
    message = '';
    const revision = ++classWriteRevision;
    void setProjectDecisionDebounced('gaussian_classes', next, 500).then(accepted => {
      if (revision !== classWriteRevision) return;
      message = accepted ? '' : 'This class edit was dropped because the shared project changed. Review the accepted value and try again.';
    });
  }

  function changeRrm(key: string, raw: string): void {
    const rrm = Number(raw) as GaussianRrm;
    if (!rrmOptions.includes(rrm)) return;
    updateClass(key, (current, computed) => {
      const automatic = automaticGaussianPrefilters(computed);
      return automatic ? switchGaussianClass(current, rrm, automatic) : current;
    });
  }

  function changeFlag(key: string, flag: 'tilt30' | 'disp45', value: boolean): void {
    updateClass(key, (current) => setGaussianFlag(current, flag, value));
  }

  function revertClass(key: string): void {
    if (!classDraft || !canEdit || !computedFresh) return;
    const next = revertGaussianClassRow(classDraft, key, scores);
    if (!next) return;
    classDraft = next;
    message = '';
    const revision = ++classWriteRevision;
    void setProjectDecisionDebounced('gaussian_classes', next, 500).then(accepted => {
      if (revision !== classWriteRevision) return;
      message = accepted ? '' : 'The class revert was dropped because the shared project changed.';
    });
  }

  function copyComputedClasses(): void {
    if (!classDraft || !computedFresh) return;
    const confirmed = classDiffCount === 0
      || window.confirm('Replace the current Gaussian class edits with the computed classes?');
    if (!confirmed) return;
    const next = copyComputedGaussianClasses(classDraft, scores, true);
    if (!next || next.rows.length !== selected.length || next.rows.some((row, index) => row.key !== selected[index])) {
      message = 'Fresh computed classes do not match the selected stacking keys.';
      return;
    }
    classDraft = next;
    message = '';
    const revision = ++classWriteRevision;
    void setProjectDecisionDebounced('gaussian_classes', next, 500).then(accepted => {
      if (revision !== classWriteRevision) return;
      message = accepted ? '' : 'Copy was dropped because the shared project changed.';
    });
  }

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

  function saveEnergyEdits(next: Stage5S24): void {
    if (!setEnergyThresholdEdits(next)) {
      message = 'The shared energy-threshold table rejected this edit.';
      return;
    }
    energyEdits = cloneStage5Decision(next);
    message = '';
  }

  function editEnergy(row: Stage5Range, fragment: string, raw: string): void {
    if (!energyEdits || !canEdit) return;
    const value = raw.trim() === '' ? null : raw.trim();
    if (value !== null && !validEnergyEdit(value)) {
      message = 'Enter a finite signed energy threshold with at most two decimal places.';
      return;
    }
    const next = cloneStage5Decision(energyEdits);
    const entry = next.ranges.find(item => sameAddress(item.row, row));
    if (!entry) return;
    entry.values[fragment] = value;
    saveEnergyEdits(next);
  }

  function revertEnergy(row: Stage5Range, fragment: string): void {
    const value = computedValue(row, fragment);
    if (!energyEdits || value === null) return;
    const next = cloneStage5Decision(energyEdits);
    const entry = next.ranges.find(item => sameAddress(item.row, row));
    if (!entry) return;
    entry.values[fragment] = value;
    saveEnergyEdits(next);
  }

  function copyEnergyThresholds(): void {
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
    saveEnergyEdits(next);
  }

  function rowEnergy(range: string, key: string, kind: 'reference' | 'bestfit'): any {
    if (!isRecord(seamlessState.computed.energetics) || !Array.isArray(seamlessState.computed.energetics.ranges)) return null;
    const data = seamlessState.computed.energetics.ranges.find((entry: any) => entry?.range === range);
    return Array.isArray(data?.rows) ? data.rows.find((entry: any) => entry?.key === key)?.[kind] ?? null : null;
  }

  function energyLabel(value: unknown): string {
    return typeof value === 'number' && Number.isFinite(value) ? value.toFixed(3) : '—';
  }

  function energyColor(value: unknown): string {
    if (typeof value !== 'number' || !Number.isFinite(value)) return 'transparent';
    const capped = Math.max(-20, Math.min(20, value));
    const position = (capped + 20) * 2.5;
    return `linear-gradient(90deg, #324761 ${position}%, #142337 ${position}%)`;
  }

  function leaveStage(): void {
    void onBack().then(accepted => {
      if (!accepted) message = 'The class table could not be set aside before leaving stage 6.';
    });
  }

  function errorFor(...names: string[]): string {
    return seamlessState.stateGraphErrors
      .filter(issue => names.includes(issue.name))
      .map(issue => `${issue.name} ${issue.state}: ${issue.message}`)
      .join(' ');
  }

  let gaussianErrors = $derived(errorFor('detect_stacking_bestfit', 'gaussian_scores_transformer',
    'validate_gaussian_classes', 'gaussian_check_transformer', 'gaussian_check', 'merge_gaussian_validation'));
  let bestfitStackings = $derived(seamlessState.computed.bestfit_stackings);
  let gaussianValidationRows = $derived(validationFresh && isRecord(validation) && Array.isArray(validation.rows)
    ? validation.rows : []);
</script>

<section class="stage-content stage-six-content">
  <div class="eyebrow">Stage 06 <span>Gaussian and reference validation</span></div>
  <h1>Review stacking classes</h1>
  <p class="stage-description">Scores and classes come from the current backend model. Editing a class only changes validation results, not the best-fit jobs.</p>
  <p class="stage6-clip">Fit clip · {clipSummary} (read only)</p>

  {#if gaussianErrors}
    <p class="stage6-error" role="alert">{gaussianErrors}</p>
  {/if}
  {#if !geometryFresh || !computedFresh}
    <p class="stage6-pending" role="status">Computing best-fit stacking geometry and Gaussian scores…</p>
  {/if}
  {#if modelChanged}
    <p class="stage6-model-changed" role="status">Model changed · current filter set {filterVersion}</p>
  {/if}

  <section class="stage6-section" aria-label="Computed and editable Gaussian classes">
    <div class="stage6-heading">
      <div><h2>Gaussian class table</h2><p>{filterNames.length} legal filters · model {filterVersion ?? 'pending'}</p></div>
      <div class="stage6-toolbar">
        <span>{classDiffCount} rows differ from computed</span>
        <button class="button button-secondary" disabled={!computedFresh || !canEdit} onclick={copyComputedClasses}>Copy computed → editable</button>
      </div>
    </div>
    {#if !classDraft}
      <p class="stage6-error" role="alert">The editable class table is unavailable. Return to stage 5 and continue again.</p>
    {:else}
      <div class="stage6-table-wrap">
        <table class="stage6-table" aria-label="Editable Gaussian classes by selected stacking">
          <thead><tr><th>Selected stacking</th><th>Computed class</th><th>Editable filter</th><th>Automatic / editable prefilters</th><th>Score · reference score</th><th>Actions</th></tr></thead>
          <tbody>
            {#each classRows as row (row.key)}
              {@const computed = computedClasses.find((item: any) => item?.key === row.key)}
              {@const differs = gaussianClassDiffCount({ ...(classDraft as GaussianClassMarker), rows: [row], basis: [row.key] }, {
                ...(isRecord(scores) ? scores : {}), classes: computed ? [computed] : [],
              }) > 0}
              {@const hasFilter = filterNames.includes(gaussianFilterName(row))}
              <tr class:stage6-different={differs}>
                <th scope="row"><code>{row.key}</code></th>
                <td>
                  {#if computed}
                    <strong>{computed.filter_name}</strong>
                    {#if computed.below_rrm95}<span class="stage6-below">below RRM-95</span>{/if}
                    {#if computed.class_pass && computed.mask_pass}<span class="stage6-pass">passes mask and score</span>
                    {:else}<span class="stage6-fail">does not pass applied mask / class</span>{/if}
                  {:else}—{/if}
                </td>
                <td>
                  <label class="stage6-visually-hidden" for={`rrm-${row.key}`}>Editable RRM class for {row.key}</label>
                  <select id={`rrm-${row.key}`} value={row.rrm} disabled={!canEdit || !computedFresh}
                    onchange={event => changeRrm(row.key, event.currentTarget.value)}>
                    {#if !rrmOptions.includes(row.rrm)}<option value={row.rrm}>RRM-{row.rrm} · unavailable</option>{/if}
                    {#each rrmOptions as rrm}<option value={rrm}>RRM-{rrm}</option>{/each}
                  </select>
                  {#if !hasFilter}<small class="stage6-fail">Saved filter {gaussianFilterName(row)} is not in the current model.</small>{/if}
                </td>
                <td class="stage6-flags">
                  {#if row.rrm === 90 || row.rrm === 95}
                    <label><input type="checkbox" checked={row.tilt30} disabled={!canEdit || !computedFresh || !filterNames.includes(`rrm${row.rrm}-tilt30`)}
                      onchange={event => changeFlag(row.key, 'tilt30', event.currentTarget.checked)} />tilt &lt; 30°</label>
                    <label><input type="checkbox" checked={row.disp45} disabled={!canEdit || !computedFresh || !filterNames.includes(`rrm${row.rrm}-disp45`)}
                      onchange={event => changeFlag(row.key, 'disp45', event.currentTarget.checked)} />displacement &lt; 45°</label>
                  {:else}<span>—</span>{/if}
                </td>
                <td>
                  {#if computed}<span>{Number(computed.score).toFixed(6)}</span><small>ref {Number(computed.reference_score).toFixed(6)}</small>{:else}—{/if}
                </td>
                <td><button class="stage4-revert" aria-label={`Revert ${row.key} to computed class`} disabled={!computedFresh || !differs || !canEdit} onclick={() => revertClass(row.key)}>↶ Revert row</button></td>
              </tr>
            {/each}
          </tbody>
        </table>
      </div>
    {/if}
  </section>

  <section class="stage6-section" aria-label="Gaussian validation results">
    <div class="stage6-heading"><div><h2>Gaussian anchor validation</h2><p>One A9 job list is shared by every editable filter choice.</p></div></div>
    {#if !validationFresh}
      <p class="stage6-pending" role="status">validating…</p>
    {:else}
      <div class="stage6-validation-list">
        {#each gaussianValidationRows as row (row.key)}
          <article class="stage6-validation-card">
            <h3><code>{row.key}</code> · {row.filter_name} · {row.worst_state}</h3>
            <ul>
              {#each row.candidates as candidate, index (`${candidate.fragment}:${candidate.nucleotide}:${index}`)}
                <li class={`stage6-state-${candidate.state}`}>
                  <strong>frag-{candidate.fragment} · {candidate.nucleotide}</strong>
                  {#if candidate.state === 'implied'}<span>implied by {candidate.implied_by}</span>
                  {:else if candidate.state === 'error'}<span>{candidate.error}</span>
                  {:else}<span>{candidate.state} · fill {Number(candidate.fill).toFixed(4)}{candidate.pose_pass ? '' : ' · pose fails'}</span>{/if}
                </li>
              {/each}
            </ul>
          </article>
        {/each}
      </div>
    {/if}
  </section>

  <section class="stage6-section" aria-label="Best-fit and reference stacking scores">
    <div class="stage6-heading"><div><h2>Stacking scores</h2><p>Raw Gaussian densities remain visible when geometry lies outside the acceptance mask.</p></div></div>
    <div class="stage6-score-grid">
      <div><h3>Best-fit candidates</h3><div class="stage6-table-wrap"><table class="stage6-table"><thead><tr><th>Anchor / fragment</th><th>distance</th><th>φ</th><th>θ</th><th>score</th><th>mask</th></tr></thead><tbody>
        {#each candidateScores as row, index (`${row.key}:${row.fragment}:${row.nucleotide}:${index}`)}<tr><td><code>{row.key}</code><small>{row.fragment_key} · {row.nucleotide}</small></td><td>{Number(row.distance).toFixed(3)} Å</td><td>{Number(row.phi_deg).toFixed(2)}°</td><td>{Number(row.theta_deg).toFixed(2)}°</td><td>{Number(row.score).toFixed(6)}</td><td>{row.mask_pass ? 'pass' : 'outside'}</td></tr>{/each}
      </tbody></table></div></div>
      <div><h3>Reference anchors</h3><div class="stage6-table-wrap"><table class="stage6-table"><thead><tr><th>Anchor</th><th>distance</th><th>φ</th><th>θ</th><th>score</th></tr></thead><tbody>
        {#each referenceScores as row (row.key)}<tr><td><code>{row.key}</code></td><td>{Number(row.distance).toFixed(3)} Å</td><td>{Number(row.phi_deg).toFixed(2)}°</td><td>{Number(row.theta_deg).toFixed(2)}°</td><td>{Number(row.score).toFixed(6)}</td></tr>{/each}
      </tbody></table></div></div>
    </div>
  </section>

  <Stage4BestFits embedded={true} onBack={onBack} onContinue={async () => false} />

  <section class="stage6-section" aria-label="Editable energy thresholds">
    <div class="stage6-heading"><div><h2>Energy thresholds</h2><p>The same S24 cells remain editable in stages 5 and 6.</p></div>
      <div class="stage6-toolbar"><button class="button button-secondary" disabled={!thresholdsFresh || !canEdit} onclick={copyEnergyThresholds}>Copy computed thresholds</button></div>
    </div>
    {#if !energyEdits}<p class="stage6-error" role="alert">The shared energy-threshold table is unavailable.</p>
    {:else}
      <div class="stage5-energy-wrap"><table class="stage5-energy-table" aria-label="Editable energy thresholds by protein range">
        <thead><tr><th>Fragment</th>{#each rangeRows.slice(1) as proteinRange, rangeIndex (`${proteinRange.chain}:${proteinRange.start}-${proteinRange.end}`)}<th colspan="3">Range {rangeIndex + 2} · chain {proteinRange.chain} {proteinRange.start}–{proteinRange.end}</th>{/each}</tr>
          <tr><th>Reference / best-fit energy, threshold</th>{#each rangeRows.slice(1) as _range, rangeIndex (`r-${rangeIndex}`)}<th>Reference</th><th>Best fit</th><th>Threshold</th>{/each}</tr></thead>
        <tbody>{#each (isRecord(fragments) && Array.isArray(fragments.fragments) ? fragments.fragments : []) as fragment (fragment.key)}
          <tr class:clipped={!keptRows?.includes(fragment.key)}><th scope="row">{fragment.key}</th>
            {#each rangeRows.slice(1) as proteinRange, rangeIndex (`${proteinRange.chain}:${proteinRange.start}-${proteinRange.end}`)}
              {@const rangeId = `range-${rangeIndex + 2}`}{@const reference = rowEnergy(rangeId, fragment.key, 'reference')}{@const bestfit = rowEnergy(rangeId, fragment.key, 'bestfit')}{@const threshold = thresholdValue(proteinRange, fragment.key)}{@const computed = computedValue(proteinRange, fragment.key)}
              <td class="stage5-energy" style={`background:${energyColor(reference?.energy)}`} title={reference?.error ?? 'Capped at ±20'}>{reference?.error ? `Error: ${reference.error}` : energyLabel(reference?.energy)}</td>
              <td class="stage5-energy" style={`background:${energyColor(bestfit?.energy)}`} title={bestfit?.error ?? 'Capped at ±20'}>{bestfit?.error ? `Error: ${bestfit.error}` : energyLabel(bestfit?.energy)}</td>
              <td class="stage5-threshold-cell"><input aria-label={`${fragment.key} ${rangeId} energy threshold`} value={threshold ?? ''} disabled={!canEdit} onchange={event => editEnergy(proteinRange, fragment.key, event.currentTarget.value)} />
                {#if computed !== null}<small class="stage5-computed-threshold">Computed: {computed}</small>{/if}
                {#if computed !== null && threshold !== computed}<span class="stage5-diff" title="Differs from computed">●</span>{/if}
                <button class="stage4-revert" aria-label={`Revert ${fragment.key} ${rangeId} to computed`} disabled={!thresholdsFresh || computed === null || threshold === computed} onclick={() => revertEnergy(proteinRange, fragment.key)}>↶</button></td>
            {/each}
          </tr>
        {/each}</tbody>
      </table></div>
    {/if}
  </section>

  {#if message}<p class="stage6-error" role="alert">{message}</p>{/if}
  <div class="stage-footer stage-footer-split">
    <button class="button button-secondary" onclick={leaveStage}><span aria-hidden="true">←</span> Back to interactions</button>
  </div>
</section>
