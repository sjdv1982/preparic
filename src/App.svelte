<script lang="ts">
  import Stage1Structure from './components/Stage1Structure.svelte';
  import Stage2Ranges from './components/Stage2Ranges.svelte';
  import Stage3Stackings from './components/Stage3Stackings.svelte';
  import Stage4BestFits from './components/Stage4BestFits.svelte';
  import Stage5Energetics from './components/Stage5Energetics.svelte';
  import Viewer from './components/Viewer.svelte';
  import { clearStructureDependentState, enterStage2, enterStage3, appState, invalidateStackingState, returnToStage1, returnToStage2 } from './lib/state.svelte';
  import { cancelSharedCifWrite, deactivateProjectDecisionAndWait, prepareProjectTransition, seamlessState, setEnergyThresholdEditsAwaited, setProjectDecisionAwaited, setSharedBiounitIndex, setSharedCif, setSharedPdbCode, setSharedRangeRows, setSharedStackingSelectionAndWait, waitForShareWrite } from './lib/seamless.svelte';
  import { findChain } from './lib/chains';
  import { cloneRangeRows, updateRow, validateRows } from './lib/ranges';
  import { mergeStackings, observeStage2RangeRows, resolveSelectedStackings } from './lib/stackings';
  import { canContinueStage4, canEnterStage4, clipForSelectedAnchors, makeThresholdSkeleton, selectedCandidateFragments, validateThresholdEdits } from './lib/stage4';
  import { adaptInteractionMatrix, candidateAnchors, makeInteractionMatrix, recomputeInteractionMatrix, reconcileEnergyThresholdEdits, stage5CanAdvance, stage5Rows, type Stage5S14 } from './lib/stage5';
  import { markerIsWellFormed } from './lib/project-state';
  import type { ChainKey, PickInfo, RangeRow, RangeSpec } from './lib/types';
  import type { MolViewer } from './lib/molstar/viewer';

  const startupQuery = new URLSearchParams(window.location.search);
  const startupPdbId = startupQuery.get('pdb');
  const startupAutoContinue = startupQuery.get('continue') === '1';
  let startupPdbIdForRestore = $state<string | null>(null);

  const decisionErrorStage: Record<string, number> = {
    ranges: 2,
    stacking_selection: 3,
    threshold_edits: 4,
    best_fit_clip: 4,
    interaction_matrix: 5,
    energy_threshold_edits: 5,
    gaussian_classes: 6,
  };

  function decisionErrorIsVisible(key: string): boolean {
    if (key === 'energy_threshold_edits') return appState.stage === 5 || appState.stage === 6;
    const stage = decisionErrorStage[key];
    return stage !== undefined && stage <= appState.stage + 1;
  }

  let seamlessStateInitialized = false;
  let startupStageRestoreHandled = false;
  let startupPdbApplied = false;
  let appliedRangeRevision = 0;
  let rangeObservation = { active: false, signature: null as string | null };

  $effect(() => {
    if (!seamlessState.ready) return;

    const requestedStartupPdbId = startupPdbId?.trim().toUpperCase() ?? '';
    if (!startupPdbApplied && seamlessState.status === 'live' && seamlessState.projectStage === 1 && !seamlessState.pdbCode
      && !seamlessState.pdbCodePendingWrite && /^[0-9][A-Za-z0-9]{3}$/.test(requestedStartupPdbId)) {
      startupPdbApplied = true;
      startupPdbIdForRestore = requestedStartupPdbId;
      appState.pdbId = requestedStartupPdbId;
    }
    const startupPdbAwaitingCommit = startupPdbApplied && seamlessState.projectStage === 1
      && !seamlessState.pdbCode && appState.pdbId === startupPdbIdForRestore;

    if (!seamlessStateInitialized) {
      seamlessStateInitialized = true;
      if (!seamlessState.pdbCodePendingWrite && !startupPdbAwaitingCommit) {
        appState.pdbId = seamlessState.pdbCode;
      }
      if (!seamlessState.biounitIndexPendingWrite) {
        appState.biounitIndex = seamlessState.biounitIndex;
      }
      return;
    }

    if (seamlessState.pdbCode !== appState.pdbId && !seamlessState.pdbCodePendingWrite && !startupPdbAwaitingCommit) {
      cancelSharedCifWrite();
      appState.pdbId = seamlessState.pdbCode;
      appState.assemblies = [];
      appState.assemblyId = '';
      appState.structureReady = false;
      appState.loading = false;
      appState.error = '';
      clearStructureDependentState();
      return;
    }

    if (seamlessState.biounitIndex !== appState.biounitIndex) {
      appState.biounitIndex = seamlessState.biounitIndex;
      appState.assemblies = [];
      appState.assemblyId = '';
      appState.structureReady = false;
      clearStructureDependentState();
    }
  });

  $effect(() => {
    if (!seamlessState.ready || !seamlessStateInitialized) return;
    setSharedBiounitIndex(appState.biounitIndex);
  });

  $effect(() => {
    if (appState.stage < 2 || !seamlessState.rangesReady) return;
    if (seamlessState.rangeRevision > appliedRangeRevision) {
      appliedRangeRevision = seamlessState.rangeRevision;
      appState.rows = cloneRangeRows(seamlessState.rangeRows);
      return;
    }
  });

  $effect(() => {
    const result = observeStage2RangeRows(rangeObservation, appState.stage, appState.rows);
    rangeObservation = result.observation;
    if (result.changed) invalidateStackingState(seamlessState.stackingRevision);
  });

  let viewer = $state<MolViewer | null>(null);
  let viewerError = $state('');
  let preparedViewerKey = '';
  let preparedViewerActive = false;

  function handleViewerReady(created: MolViewer) {
    viewer = created;
    created.setOnCifDownloaded(setSharedCif);
    viewerError = '';
  }

  function rangeForRow(row: RangeRow | undefined): RangeSpec | null {
    if (!row) return null;
    const chain = appState.chains.find(item => item.name === row.chain);
    const rawStart = row.start.trim();
    const rawEnd = row.end.trim();
    const start = Number(rawStart);
    const end = Number(rawEnd);
    if (!chain || !rawStart || !rawEnd || !Number.isInteger(start) || !Number.isInteger(end)) return null;
    if (!chain.seqIds.has(start) || !chain.seqIds.has(end)) return null;
    return {
      labelAsymId: chain.labelAsymId,
      instanceId: chain.instanceId,
      start: Math.min(start, end),
      end: Math.max(start, end),
    };
  }

  function handlePick(pick: PickInfo) {
    if (!seamlessState.ready || appState.stage !== 2) return;
    const cell = appState.activeCell;
    if (!cell) return;
    const row = appState.rows.find(item => item.id === cell.rowId);
    if (!row) return;
    const selectedRange = rangeForRow(row);
    if (selectedRange) {
      if (selectedRange.labelAsymId !== pick.labelAsymId || selectedRange.instanceId !== pick.instanceId) return;
      if (pick.labelSeqId < selectedRange.start) {
        appState.rows = updateRow(appState.rows, row.id, { start: String(pick.labelSeqId), end: String(selectedRange.end) });
      } else if (pick.labelSeqId > selectedRange.end) {
        appState.rows = updateRow(appState.rows, row.id, { start: String(selectedRange.start), end: String(pick.labelSeqId) });
      } else if (cell.column !== 'chain') {
        appState.rows = updateRow(appState.rows, row.id, { [cell.column]: String(pick.labelSeqId) });
      }
      setSharedRangeRows(appState.rows);
      return;
    }
    if (cell.column === 'chain') {
      const chain = findChain(appState.chains, pick.labelAsymId, pick.instanceId);
      if (chain) {
        appState.rows = updateRow(appState.rows, row.id, { chain: chain.name });
        setSharedRangeRows(appState.rows);
      }
      return;
    }
    const focusChain = appState.chains.find(chain => chain.name === row.chain);
    if (!focusChain || focusChain.labelAsymId !== pick.labelAsymId || focusChain.instanceId !== pick.instanceId) return;
    appState.rows = updateRow(appState.rows, row.id, { [cell.column]: String(pick.labelSeqId) });
    setSharedRangeRows(appState.rows);
  }

  async function proceedToRanges() {
    if (!viewer || !seamlessState.ready || !seamlessState.structureIdReady
      || seamlessState.structureId !== seamlessState.pdbCode) {
      appState.error = 'The shared structure is still being checked. Try Continue again when it is ready.';
      return;
    }
    const chains = viewer.listChains();
    if (chains.length === 0) {
      appState.error = 'This unit has no polymer chains with residue numbering.';
      appState.stage = 1;
      return;
    }
    appState.chains = chains;
    if (!appState.rows.length && seamlessState.rangeRows.length > 0) {
      appState.rows = cloneRangeRows(seamlessState.rangeRows);
    }
    if (!appState.rows.length) {
      const firstChain = chains[0];
      appState.rows = [{ id: `range-start-${Date.now()}`, alwaysVisible: false, chain: firstChain.name,
        start: String(firstChain.seqMin), end: String(firstChain.seqMax) }];
    }
    if (seamlessState.projectStage === 1) {
      if (!await prepareProjectTransition(2)) {
        appState.error = 'The project changed while preparing this transition. Review it and try again.';
        return;
      }
      if (!seamlessState.ready || !seamlessState.structureIdReady
        || seamlessState.structureId !== seamlessState.pdbCode
        || seamlessState.pdbCode !== appState.pdbId.trim().toUpperCase()) {
        appState.error = 'The shared structure changed while preparing this transition. Check it and try again.';
        return;
      }
      if (!setSharedRangeRows(appState.rows) || !await waitForShareWrite('ranges')) {
        appState.error = 'Save a complete range decision before continuing.';
        return;
      }
    }
    appliedRangeRevision = seamlessState.rangeRevision;
    enterStage2(seamlessState.rangesReady ? seamlessState.rangeRows : undefined);
  }

  $effect(() => {
    if (startupStageRestoreHandled || !seamlessState.ready || !seamlessStateInitialized
      || !appState.structureReady || appState.loading || !seamlessState.structureIdReady) return;
    if (seamlessState.projectStage === 1 && startupPdbIdForRestore === seamlessState.pdbCode && startupAutoContinue
      && seamlessState.structureId === seamlessState.pdbCode) proceedToRanges();
    startupStageRestoreHandled = true;
  });

  async function goBack() {
    if (!await deactivateProjectDecisionAndWait('ranges')) {
      seamlessState.rangeError = 'The range decision could not be deactivated. Review the shared project and try again.';
      return;
    }
    returnToStage1();
    void viewer?.render({ focus: null });
  }

  async function complete(): Promise<void> {
    if (!seamlessState.computedFreshness.stackings || appState.stackings.length === 0) {
      seamlessState.stackingSelectionError = 'Wait for fresh stacking results before continuing.';
      return;
    }
    const rows = seamlessState.projectSnapshot.ranges;
    if (!rows || typeof rows !== 'object' || !('active' in rows) || rows.active !== true
      || seamlessState.pendingInputs.includes('ranges') || validateRows(appState.rows, appState.chains).length > 0) {
      seamlessState.stackingSelectionError = 'Save and validate the current ranges before continuing.';
      return;
    }
    const currentRows = (rows as { rows?: Array<{ chain: string; start: number; end: number }> }).rows;
    const localRows = appState.rows.map(row => ({ chain: row.chain, start: Number(row.start), end: Number(row.end) }));
    if (!Array.isArray(currentRows) || JSON.stringify(currentRows) !== JSON.stringify(localRows)
      || currentRows.some(row => {
        const chain = appState.chains.find(item => item.name === row.chain);
        if (!chain) return true;
        let observed = 0;
        for (const resid of chain.seqIds) if (resid >= row.start && resid <= row.end) observed += 1;
        return observed !== row.end - row.start + 1;
      })) {
      seamlessState.stackingSelectionError = 'Commit valid ranges with no unobserved residues before continuing.';
      return;
    }
    appState.stackingSelection = mergeStackings(appState.stackings, appState.stackingSelection, appState.stackings).stackingSelection;
    if (!await prepareProjectTransition(3)) {
      seamlessState.stackingSelectionError = 'The project changed while preparing this transition. Review it and try again.';
      return;
    }
    const confirmedRanges = seamlessState.projectSnapshot.ranges;
    const confirmedRows = confirmedRanges && typeof confirmedRanges === 'object' && 'rows' in confirmedRanges
      ? confirmedRanges.rows : null;
    const confirmedLocalRows = appState.rows.map(row => ({ chain: row.chain, start: Number(row.start), end: Number(row.end) }));
    if (!seamlessState.ready || !seamlessState.computedFreshness.stackings || !Array.isArray(confirmedRows)
      || !confirmedRanges || typeof confirmedRanges !== 'object' || !('active' in confirmedRanges)
      || confirmedRanges.active !== true || seamlessState.pendingInputs.includes('ranges')
      || JSON.stringify(confirmedRows) !== JSON.stringify(confirmedLocalRows)
      || validateRows(appState.rows, appState.chains).length > 0
      || confirmedRows.some(row => {
        const chain = appState.chains.find(item => item.name === row.chain);
        if (!chain) return true;
        let observed = 0;
        for (const resid of chain.seqIds) if (resid >= row.start && resid <= row.end) observed += 1;
        return observed !== row.end - row.start + 1;
      })) {
      seamlessState.stackingSelectionError = 'The range or stacking results changed while preparing this transition. Review them and try again.';
      return;
    }
    appState.stackingSelection = mergeStackings(appState.stackings, appState.stackingSelection, appState.stackings).stackingSelection;
    if (!await setSharedStackingSelectionAndWait(appState.stackingSelection)) {
      seamlessState.stackingSelectionError = 'The stacking selection was not accepted by the shared project.';
      return;
    }
    enterStage3();
  }

  async function backToRanges(): Promise<void> {
    if (!await deactivateProjectDecisionAndWait('stacking_selection')) {
      seamlessState.stackingSelectionError = 'The stacking selection could not be deactivated. Review the shared project and try again.';
      return;
    }
    returnToStage2();
  }

  async function backToStackings(): Promise<boolean> {
    if (!await deactivateProjectDecisionAndWait('threshold_edits')) {
      seamlessState.syncNotice = 'The threshold decision could not be deactivated. Review the shared project and try again.';
      return false;
    }
    return true;
  }

  async function continueToStage4(): Promise<void> {
    const stackings = seamlessState.computed.stackings;
    const selection = seamlessState.projectSnapshot.stacking_selection;
    const inventory = seamlessState.computed.fragments;
    if (!canEnterStage4({
      selection,
      stackings,
      stackingsFresh: seamlessState.computedFreshness.stackings === true,
      fragmentsFresh: seamlessState.computedFreshness.fragments === true,
      rangePdbsFresh: seamlessState.computedFreshness.range_pdbs_aa === true,
    })) {
      seamlessState.stackingSelectionError = 'Stage 4 needs the current selection, reference fragments, and all-atom PDBs to be fresh.';
      return;
    }
    if (!inventory || typeof inventory !== 'object' || !('echo' in inventory)
      || !('fragments' in inventory)) {
      seamlessState.stackingSelectionError = 'The canonical fragment list is not ready.';
      return;
    }
    const inventoryFragments = inventory.fragments;
    if (!Array.isArray(inventoryFragments)) {
      seamlessState.stackingSelectionError = 'The canonical fragment list is not ready.';
      return;
    }
    if (!await prepareProjectTransition(4)) {
      seamlessState.stackingSelectionError = 'Later-stage decisions could not be deactivated before opening stage 4.';
      return;
    }
    const selected = selection && typeof selection === 'object' && 'selected' in selection
      && selection.selected && typeof selection.selected === 'object'
      ? selection.selected as Record<string, boolean> : {};
    const selectedKeys = Object.keys(selected).filter(key => selected[key] === true);
    const candidateNumbers = selectedCandidateFragments(inventory, selectedKeys, stackings);
    const basis = inventory.echo;
    const currentClip = seamlessState.projectSnapshot.best_fit_clip as any;
    const clipBasis = currentClip?.basis;
    const basisMatches = clipBasis && clipBasis.pdbcode === (basis as any).pdbcode
      && clipBasis.rna?.chain === (basis as any).rna.chain
      && clipBasis.rna?.start === (basis as any).rna.start
      && clipBasis.rna?.end === (basis as any).rna.end;
    let nextClip = basisMatches
      ? clipForSelectedAnchors(currentClip, inventory, selectedKeys, stackings)
      : { schema: 1 as const, basis, clip1: 0, clip2: 0 };
    if (!Number.isSafeInteger(nextClip.clip1) || !Number.isSafeInteger(nextClip.clip2)
      || nextClip.clip1 < 0 || nextClip.clip2 < 0
      || nextClip.clip1 + nextClip.clip2 >= inventoryFragments.length
      || candidateNumbers.some(number => number <= nextClip.clip1 || number > inventoryFragments.length - nextClip.clip2)) {
      seamlessState.stackingSelectionError = 'The existing clip is invalid for this fragment list and selected anchors.';
      return;
    }
    if (!basisMatches || nextClip.clip1 !== currentClip.clip1 || nextClip.clip2 !== currentClip.clip2) {
      if (!await setProjectDecisionAwaited('best_fit_clip', nextClip)) {
        seamlessState.stackingSelectionError = 'The clip could not be reset or clamped before opening stage 4.';
        return;
      }
    }
    const existing = seamlessState.projectSnapshot.threshold_edits;
    const retained = validateThresholdEdits(existing, inventory)
      ? { ...(existing as any), active: true }
      : makeThresholdSkeleton(inventory, seamlessState.pdbCode, true);
    if (!retained || !await setProjectDecisionAwaited('threshold_edits', retained)) {
      seamlessState.stackingSelectionError = 'The stage-4 threshold table could not be saved.';
      return;
    }
    seamlessState.stackingSelectionError = '';
  }

  async function continueToStage5(): Promise<boolean> {
    const snapshot = seamlessState.projectSnapshot;
    const inventory = seamlessState.computed.fragments;
    const clip = snapshot.best_fit_clip;
    const selection = snapshot.stacking_selection;
    const stackings = seamlessState.computed.stackings;
    const rangeMarker = snapshot.ranges;
    const rangeRows = rangeMarker && typeof rangeMarker === 'object' && 'rows' in rangeMarker
      && Array.isArray(rangeMarker.rows) ? rangeMarker.rows : null;
    const edits = snapshot.threshold_edits;
    const fragmentCount = inventory && typeof inventory === 'object' && 'fragments' in inventory
      && Array.isArray(inventory.fragments) ? inventory.fragments.length : 0;
    const selectedKeys = selection && typeof selection === 'object' && 'selected' in selection
      && selection.selected && typeof selection.selected === 'object'
      ? Object.keys(selection.selected as Record<string, boolean>)
        .filter(key => (selection.selected as Record<string, boolean>)[key] === true) : [];
    const selectedCandidates = selectedCandidateFragments(inventory, selectedKeys, stackings);
    const thresholdValidatorFailed = seamlessState.stateGraphErrors.some(issue =>
      issue.name === 'validate_threshold_edits' && issue.state === 'failed');
    if (!seamlessState.ready || seamlessState.status !== 'live' || seamlessState.projectStage !== 4
      || seamlessState.pendingInputs.includes('threshold_edits') || seamlessState.pendingInputs.includes('best_fit_clip')
      || seamlessState.computedFreshness.fragments !== true || seamlessState.computedFreshness.energetics !== true
      || !rangeRows || !inventory || !selection || !clip || !validateThresholdEdits(edits, inventory)
      || thresholdValidatorFailed) return false;
    const currentEcho = (inventory as any).echo;
    const clipBasis = (clip as any).basis;
    if (!currentEcho || !clipBasis || clipBasis.pdbcode !== currentEcho.pdbcode
      || !clipBasis.rna || !currentEcho.rna || clipBasis.rna.chain !== currentEcho.rna.chain
      || clipBasis.rna.start !== currentEcho.rna.start || clipBasis.rna.end !== currentEcho.rna.end) return false;
    if (!canContinueStage4({
      edits,
      fragmentCount,
      clip,
      selectedCandidates,
      energeticsFresh: seamlessState.computedFreshness.energetics === true,
      thresholdValidatorFailed,
    })) return false;
    const matrixSkeleton = makeInteractionMatrix(inventory, clip, rangeRows, selection, stackings);
    const anchors = matrixSkeleton ? candidateAnchors(inventory, selection, stackings, matrixSkeleton.rows) : null;
    if (!matrixSkeleton || !anchors) return false;

    const reconciled = reconcileEnergyThresholdEdits(
      seamlessState.energyThresholdEdits,
      inventory,
      rangeRows,
    );
    if (!reconciled) return false;
    if (JSON.stringify(seamlessState.energyThresholdEdits) !== JSON.stringify(reconciled)
      && !await setEnergyThresholdEditsAwaited(reconciled)) return false;

    if (!await prepareProjectTransition(5)) return false;
    const existing = snapshot.interaction_matrix;
    let nextMatrix: Stage5S14 | null = null;
    if (markerIsWellFormed('interaction_matrix', existing)
      && existing && typeof existing === 'object' && 'basis' in existing && 'ranges' in existing
      && 'rows' in existing && 'grid' in existing
      && JSON.stringify(existing.basis) === JSON.stringify(matrixSkeleton.basis)
      && JSON.stringify(existing.ranges) === JSON.stringify(matrixSkeleton.ranges)) {
      nextMatrix = adaptInteractionMatrix(existing, matrixSkeleton.rows, anchors, seamlessState.computed.energetics);
    }
    if (!nextMatrix) {
      nextMatrix = recomputeInteractionMatrix(matrixSkeleton, anchors, seamlessState.computed.energetics);
    }
    // A lock-row scoring error is a stage-5 matrix diagnostic, not a veto on
    // entering stage 5. The user must be able to see the affected anchor and
    // return to stage 3 to deselect it.
    if (!nextMatrix) nextMatrix = matrixSkeleton;
    nextMatrix = { ...nextMatrix, active: true };
    return await setProjectDecisionAwaited('interaction_matrix', nextMatrix);
  }

  async function backToStage4(): Promise<boolean> {
    return await deactivateProjectDecisionAndWait('interaction_matrix');
  }

  async function continueToStage6(): Promise<boolean> {
    const snapshot = seamlessState.projectSnapshot;
    const inventory = seamlessState.computed.fragments;
    const rangeMarker = snapshot.ranges;
    const rangeRows = rangeMarker && typeof rangeMarker === 'object' && 'rows' in rangeMarker
      && Array.isArray(rangeMarker.rows) ? rangeMarker.rows : null;
    const rows = stage5Rows(inventory, snapshot.best_fit_clip);
    const anchors = rows ? candidateAnchors(
      inventory, snapshot.stacking_selection, seamlessState.computed.stackings, rows,
    ) : null;
    const classes = snapshot.gaussian_classes;
    const selected = snapshot.stacking_selection && typeof snapshot.stacking_selection === 'object'
      && 'selected' in snapshot.stacking_selection && snapshot.stacking_selection.selected
      && typeof snapshot.stacking_selection.selected === 'object'
      ? Object.keys(snapshot.stacking_selection.selected as Record<string, boolean>)
        .filter(key => (snapshot.stacking_selection as any).selected[key] === true).sort() : [];
    const classRowsMatch = markerIsWellFormed('gaussian_classes', classes)
      && classes && typeof classes === 'object' && 'basis' in classes
      && JSON.stringify([...(classes.basis as string[])].sort()) === JSON.stringify(selected);
    const validatorFailed = seamlessState.stateGraphErrors.some(issue =>
      issue.name === 'validate_energy_threshold_edits' && issue.state === 'failed');
    const gate = stage5CanAdvance(snapshot, {
      ready: seamlessState.ready && seamlessState.status === 'live',
      stage: seamlessState.projectStage,
      energeticsFresh: seamlessState.computedFreshness.energetics === true,
      bestFitsFresh: seamlessState.computedFreshness.best_fits === true,
      geometryFresh: seamlessState.computedFreshness.bestfit_stackings === true,
      gaussianFresh: seamlessState.computedFreshness.gaussian_scores === true,
      gaussianClassesValid: !!classRowsMatch,
      validatorFailed,
      fragments: inventory,
      rangeRows,
      stackings: seamlessState.computed.stackings,
      anchors: anchors ?? undefined,
      energetics: seamlessState.computed.energetics,
      energyEdits: seamlessState.energyThresholdEdits,
    });
    if (!gate.allowed || !classRowsMatch || !classes || typeof classes !== 'object') return false;
    return await setProjectDecisionAwaited('gaussian_classes', { ...classes, active: true });
  }

  $effect(() => {
    if (!viewer) return;
    if (appState.stage >= 4 && !seamlessState.ready) {
      preparedViewerActive = true;
      return;
    }
    if (!seamlessState.ready) return;
    if (appState.stage >= 4) {
      preparedViewerActive = true;
      const share = seamlessState.rangePdbsAa as any;
      if (seamlessState.computedFreshness.range_pdbs_aa !== true || !share
        || share.schema !== 1 || !Array.isArray(share.pdbs)) {
        if (preparedViewerKey !== 'pending') {
          preparedViewerKey = 'pending';
          void viewer.clearPreparedPdbs();
        }
        return;
      }
      const inventory = seamlessState.computed.fragments as any;
      const referencePdbs = seamlessState.computed.reference_fragment_pdbs as any;
      const fits = seamlessState.computed.best_fits as any;
      const overlayNumbers = appState.stage === 4 ? appState.stage4OverlayFragments : [];
      const fragments = Array.isArray(inventory?.fragments)
        ? overlayNumbers.map(number => inventory.fragments.find((item: any) => item.number === number))
          .filter((item: any) => item && typeof item.key === 'string')
        : [];
      const overlays: { label: string; pdb: string }[] = [];
      for (const fragment of fragments) {
        if (seamlessState.computedFreshness.reference_fragment_pdbs === true) {
          const reference = referencePdbs?.rows?.find((item: any) => item.key === fragment.key)?.pdb;
          if (typeof reference === 'string') overlays.push({ label: `Reference ${fragment.key}`, pdb: reference });
        }
        if (seamlessState.computedFreshness.best_fits === true) {
          const bestfit = fits?.rows?.find((item: any) => item.key === fragment.key)?.pdb;
          if (typeof bestfit === 'string') overlays.push({ label: `Best fit ${fragment.key}`, pdb: bestfit });
        }
      }
      const loadKey = JSON.stringify([
        seamlessState.pdbCode.toLowerCase(), share.range_specs,
        share.pdbs.map((item: any) => [item.range, item.name, item.pdb]),
        overlayNumbers, overlays,
      ]);
      if (loadKey === preparedViewerKey) return;
      preparedViewerKey = loadKey;
      viewer.setGranularity('residue');
      const isCurrentTarget = () => seamlessState.ready && appState.stage >= 4 && preparedViewerKey === loadKey;
      void viewer.showPreparedPdbs(share.pdbs, overlays, isCurrentTarget).catch(error => {
        viewerError = error instanceof Error ? error.message : String(error);
      });
      return;
    }
    if (preparedViewerActive) {
      preparedViewerActive = false;
      preparedViewerKey = '';
      void viewer.restoreLastVisibleView().then(() => {
        if (appState.stage === 3 && seamlessState.computedFreshness.stackings) {
          const selectedStackings = resolveSelectedStackings(
            appState.stackings,
            appState.stackingSelection,
            appState.rows,
            appState.chains,
          );
          void viewer?.render({ focus: null, stackings: selectedStackings });
        }
      }).catch(error => { viewerError = error instanceof Error ? error.message : String(error); });
      return;
    }
    if (appState.stage === 3 && seamlessState.computedFreshness.stackings) {
      viewer.setGranularity('residue');
      const selectedStackings = resolveSelectedStackings(
        appState.stackings,
        appState.stackingSelection,
        appState.rows,
        appState.chains,
      );
      void viewer.render({ focus: null, stackings: selectedStackings });
      return;
    }
    if (appState.stage !== 2) {
      viewer.setGranularity('chain');
      void viewer.render({ focus: null });
      return;
    }

    const cell = appState.activeCell;
    const activeRow = cell ? appState.rows.find(item => item.id === cell.rowId) : undefined;
    const activeRowHasChain = activeRow && appState.chains.some(item => item.name === activeRow.chain);
    const row = activeRowHasChain
      ? activeRow
      : appState.rows.find(item => appState.chains.some(chain => chain.name === item.chain));
    const chain = row ? appState.chains.find(item => item.name === row.chain) : undefined;
    if (!row || !chain) {
      viewer.setGranularity('chain');
      void viewer.render({ focus: null });
      return;
    }
    const selectedRange = rangeForRow(row);

    const alwaysVisible: RangeSpec[] = [];
    for (const pinned of appState.rows) {
      if (!pinned.alwaysVisible) continue;
      const pinnedRange = rangeForRow(pinned);
      if (pinnedRange) alwaysVisible.push(pinnedRange);
    }
    const sidechainRanges = appState.showSidechainSticks
      ? appState.rows.flatMap(item => {
          const range = rangeForRow(item);
          return range ? [range] : [];
        })
      : [];
    const focus: ChainKey = { labelAsymId: chain.labelAsymId, instanceId: chain.instanceId };
    const isRangeCell = !!cell && cell.rowId === row.id && (cell.column === 'start' || cell.column === 'end');
    viewer.setGranularity(selectedRange !== null || isRangeCell ? 'residue' : 'chain');
    void viewer.render({ focus, range: selectedRange ?? undefined, alwaysVisible, sidechainRanges, others: appState.otherChainsMode });
  });

</script>

<svelte:head>
  <title>Molecular Selection · Preparic</title>
  <meta name="description" content="Define chain and residue-range selections on a molecular structure." />
</svelte:head>

<main class="app-shell">
  <header class="topbar">
    <a class="brand" href="/" aria-label="Preparic home">
      <span class="brand-mark"><span></span><span></span><span></span><span></span></span>
      <span>preparic</span>
    </a>
    <div class="topbar-center"><span class="topbar-dot"></span> Molecular selection workspace</div>
    <div class="topbar-right"><span class="phase-pill">PHASE 01</span></div>
  </header>

  <div class="workspace">
    <aside class="sidebar" inert={!seamlessState.ready || (appState.stage > 1 && appState.stage <= 3 && !appState.structureReady)}>
      <div class="stage-progress" aria-label="Workflow progress">
        <div class:current={appState.stage === 1} class:complete={appState.stage > 1} class="progress-step">
          <span class="progress-number">{appState.stage > 1 ? '✓' : '01'}</span>
          <span><strong>Structure</strong><small>Choose a PDB unit</small></span>
        </div>
        <div class="progress-rail"></div>
        <div class:current={appState.stage === 2} class:complete={appState.stage > 2} class="progress-step">
          <span class="progress-number">{appState.stage > 2 ? '✓' : '02'}</span>
          <span><strong>Ranges</strong><small>Select residues</small></span>
        </div>
        <div class="progress-rail"></div>
        <div class:current={appState.stage === 3} class:complete={appState.stage > 3} class="progress-step">
          <span class="progress-number">{appState.stage > 3 ? '✓' : '03'}</span>
          <span><strong>Stackings</strong><small>Review interactions</small></span>
        </div>
        <div class="progress-rail"></div>
        <div class:current={appState.stage === 4} class:complete={appState.stage > 4} class="progress-step">
          <span class="progress-number">{appState.stage > 4 ? '✓' : '04'}</span>
          <span><strong>Best fits</strong><small>Set thresholds</small></span>
        </div>
        <div class="progress-rail"></div>
        <div class:current={appState.stage === 5} class:complete={appState.stage > 5} class="progress-step">
          <span class="progress-number">{appState.stage > 5 ? '✓' : '05'}</span>
          <span><strong>Interactions</strong><small>Review energy matrix</small></span>
        </div>
        <div class="progress-rail"></div>
        <div class:current={appState.stage === 6} class="progress-step">
          <span class="progress-number">06</span>
          <span><strong>Gaussian</strong><small>Validate classes</small></span>
        </div>
      </div>

      {#if viewerError}
        <div class="viewer-error" role="alert"><strong>Viewer unavailable</strong><span>{viewerError}</span></div>
      {/if}
      {#if seamlessState.syncNotice}
        <div class="viewer-error" role="alert"><strong>Project synchronization</strong><span>{seamlessState.syncNotice}</span></div>
      {/if}
      {#each Object.entries(seamlessState.decisionErrors).filter(([key]) => decisionErrorIsVisible(key)) as [key, message] (`${key}:${message}`)}
        <div class="viewer-error" role="alert"><strong>Stored decision unavailable</strong><span>{message}</span></div>
      {/each}
      {#each seamlessState.stateGraphErrors as issue (`${issue.name}:${issue.state}:${issue.message}`)}
        <div class="viewer-error" role="alert"><strong>{issue.name} {issue.state}</strong><span>{issue.message}</span></div>
      {/each}

      {#if appState.stage > 1 && appState.stage <= 3 && (!appState.structureReady || appState.loading)}
        <div hidden><Stage1Structure {viewer} startupPdbId={null} autoContinue={false} onContinue={proceedToRanges} /></div>
      {/if}
      {#if appState.stage === 1}
        <Stage1Structure {viewer} startupPdbId={startupPdbIdForRestore} autoContinue={false} onContinue={proceedToRanges} />
      {:else if appState.stage === 2}
        <Stage2Ranges onBack={goBack} onComplete={complete} />
      {:else if appState.stage === 3}
        <Stage3Stackings onBack={backToRanges} onContinue={continueToStage4} />
      {:else if appState.stage === 4}
        <Stage4BestFits onBack={backToStackings} onContinue={continueToStage5} />
      {:else if appState.stage === 5}
        <Stage5Energetics onBack={backToStage4} onContinue={continueToStage6} />
      {:else}
        <section class="stage-content" aria-live="polite">
          <div class="eyebrow">Stage {appState.stage}</div>
          <h1>Shared project restored</h1>
          <p class="stage-description">This stage is waiting for its editor to load.</p>
        </section>
      {/if}

      <div class="sidebar-bottom"><span class="sidebar-bottom-mark">P</span><span>Stage {appState.stage} · Shared project</span><span class="sidebar-bottom-status">●</span></div>
    </aside>

    <section class="viewer-pane" aria-label="Molecular viewer workspace" inert={!seamlessState.ready}>
      <div class="viewer-toolbar">
        <div class="viewer-title"><span class="viewer-title-icon">◈</span><span>{appState.pdbId || 'Structure viewer'}</span>{#if appState.structureReady}<span class="loaded-badge">LOADED</span>{/if}</div>
        <div class="viewer-toolbar-right"><span class="viewer-mode"><span>●</span> {appState.stage >= 4 ? 'Prepared range view' : appState.stage === 1 ? 'Structure view' : appState.stage === 2 ? 'Pick mode' : 'Interaction review'}</span><span class="toolbar-divider"></span><span class="viewer-version">Mol*</span></div>
      </div>
      <div class="viewer-stage">
        <Viewer onReady={handleViewerReady} onError={(message) => viewerError = message} onPick={handlePick} />
        {#if seamlessState.ready && appState.stage >= 4 && seamlessState.computedFreshness.range_pdbs_aa !== true}
          <div class="viewer-empty"><strong>Preparing all-atom structure…</strong><span>Waiting for fresh prepared range PDBs.</span></div>
        {:else if seamlessState.ready && appState.stage < 4 && !appState.structureReady}
          <div class="viewer-empty" class:viewer-empty-loading={appState.loading}>
            <div class="empty-illustration">
              <div class="molecule-orbit orbit-one"></div><div class="molecule-orbit orbit-two"></div>
              <span class="molecule-node node-one"></span><span class="molecule-node node-two"></span><span class="molecule-node node-three"></span><span class="molecule-core"></span>
            </div>
            <strong>{appState.loading ? 'Preparing the structure' : 'Your structure will appear here'}</strong>
            <span>{appState.loading ? 'Downloading coordinates and assembling the selected unit.' : 'Enter a PDB code to begin exploring.'}</span>
          </div>
        {/if}
        <div class="viewer-corner-label"><span>{appState.stage >= 4 ? 'PREPARED PDBs' : 'RCSB PDB'}</span><span class="corner-separator">/</span><span>Mol*</span></div>
      </div>
      <div class="viewer-footer">
        <div class="viewer-footer-left"><span class="footer-live-dot"></span>{appState.structureReady ? `Ready${appState.pdbId ? ` · ${appState.pdbId}` : ''}` : 'Awaiting structure'}</div>
        <div class="axis-indicator"><span class="axis-y">Y</span><span class="axis-x">X</span><span class="axis-z">Z</span><i></i></div>
        <div class="viewer-footer-right">Molecular graphics <span>•</span> Phase 1</div>
      </div>
    </section>
    {#if !seamlessState.ready}
      <div class="sync-overlay" role="status" aria-live="polite">
        <div class="sync-overlay-card">
          <span class="sync-overlay-spinner" aria-hidden="true"></span>
          <strong>{seamlessState.status === 'disconnected' ? 'Connection lost' : 'Synchronizing shared project'}</strong>
          <span>{seamlessState.status === 'disconnected'
            ? 'Your last project view is preserved. Editing is paused until shared state is reread.'
            : 'Your project view is preserved while shared values are reread.'}</span>
        </div>
      </div>
    {/if}
  </div>
</main>
