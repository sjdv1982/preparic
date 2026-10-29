<script lang="ts">
  import Stage1Structure from './components/Stage1Structure.svelte';
  import Stage2Ranges from './components/Stage2Ranges.svelte';
  import Stage3Stackings from './components/Stage3Stackings.svelte';
  import Viewer from './components/Viewer.svelte';
  import { enterStage2, enterStage3, appState, invalidateStackingState, returnToStage1, returnToStage2 } from './lib/state.svelte';
  import { cancelSharedCifWrite, deactivateProjectDecisionAndWait, prepareProjectTransition, seamlessState, setSharedBiounitIndex, setSharedCif, setSharedPdbCode, setSharedRangeRows, setSharedStackingSelectionAndWait, waitForShareWrite } from './lib/seamless.svelte';
  import { findChain } from './lib/chains';
  import { cloneRangeRows, updateRow, validateRows } from './lib/ranges';
  import { mergeStackings, observeStage2RangeRows, resolveSelectedStackings } from './lib/stackings';
  import type { ChainKey, PickInfo, RangeRow, RangeSpec } from './lib/types';
  import type { MolViewer } from './lib/molstar/viewer';

  const startupQuery = new URLSearchParams(window.location.search);
  const startupPdbId = startupQuery.get('pdb');
  const startupAutoContinue = startupQuery.get('continue') === '1';
  let startupPdbIdForRestore = $state<string | null>(null);

  let seamlessStateInitialized = false;
  let startupStageRestoreHandled = false;
  let appliedRangeRevision = 0;
  let rangeObservation = { active: false, signature: null as string | null };

  $effect(() => {
    if (!seamlessState.ready) return;

    if (!seamlessStateInitialized) {
      seamlessStateInitialized = true;
      if (!seamlessState.pdbCode && seamlessState.projectStage === 1
        && startupPdbId && /^[0-9][A-Za-z0-9]{3}$/.test(startupPdbId)) {
        startupPdbIdForRestore = startupPdbId.toUpperCase();
        appState.pdbId = startupPdbIdForRestore;
        setSharedPdbCode(startupPdbIdForRestore);
      } else if (!seamlessState.pdbCodePendingWrite) {
        appState.pdbId = seamlessState.pdbCode;
      }
      if (!seamlessState.biounitIndexPendingWrite) {
        appState.biounitIndex = seamlessState.biounitIndex;
      }
      return;
    }

    if (seamlessState.pdbCode !== appState.pdbId && !seamlessState.pdbCodePendingWrite) {
      cancelSharedCifWrite();
      appState.pdbId = seamlessState.pdbCode;
      appState.assemblies = [];
      appState.assemblyId = '';
      appState.structureReady = false;
      appState.loading = false;
      appState.error = '';
      returnToStage1();
      return;
    }

    if (seamlessState.biounitIndex !== appState.biounitIndex) {
      appState.biounitIndex = seamlessState.biounitIndex;
      if (appState.stage === 2) returnToStage1();
    }
  });

  $effect(() => {
    if (!seamlessState.ready || !seamlessStateInitialized) return;
    setSharedBiounitIndex(appState.biounitIndex);
  });

  $effect(() => {
    if (appState.stage !== 2 || !seamlessState.rangesReady) return;
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

  $effect(() => {
    if (!viewer || !seamlessState.ready) return;
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
    <aside class="sidebar" inert={!seamlessState.ready}>
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
      {#each seamlessState.stateGraphErrors as issue (`${issue.name}:${issue.state}:${issue.message}`)}
        <div class="viewer-error" role="alert"><strong>{issue.name} {issue.state}</strong><span>{issue.message}</span></div>
      {/each}

      {#if appState.stage === 1 || !appState.structureReady}
        <Stage1Structure {viewer} startupPdbId={startupPdbIdForRestore} autoContinue={false} onContinue={proceedToRanges} />
      {:else if appState.stage === 2}
        <Stage2Ranges onBack={goBack} onComplete={complete} />
      {:else if appState.stage === 3}
        <Stage3Stackings onBack={backToRanges} />
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
        <div class="viewer-toolbar-right"><span class="viewer-mode"><span>●</span> {appState.stage === 1 ? 'Structure view' : appState.stage === 2 ? 'Pick mode' : 'Interaction review'}</span><span class="toolbar-divider"></span><span class="viewer-version">Mol*</span></div>
      </div>
      <div class="viewer-stage">
        <Viewer onReady={handleViewerReady} onError={(message) => viewerError = message} onPick={handlePick} />
        {#if !appState.structureReady}
          <div class="viewer-empty" class:viewer-empty-loading={appState.loading}>
            <div class="empty-illustration">
              <div class="molecule-orbit orbit-one"></div><div class="molecule-orbit orbit-two"></div>
              <span class="molecule-node node-one"></span><span class="molecule-node node-two"></span><span class="molecule-node node-three"></span><span class="molecule-core"></span>
            </div>
            <strong>{appState.loading ? 'Preparing the structure' : 'Your structure will appear here'}</strong>
            <span>{appState.loading ? 'Downloading coordinates and assembling the selected unit.' : 'Enter a PDB code to begin exploring.'}</span>
          </div>
        {/if}
        <div class="viewer-corner-label"><span>RCSB PDB</span><span class="corner-separator">/</span><span>Mol*</span></div>
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
