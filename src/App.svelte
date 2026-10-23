<script lang="ts">
  import Stage1Structure from './components/Stage1Structure.svelte';
  import Stage2Ranges from './components/Stage2Ranges.svelte';
  import Stage3Stackings from './components/Stage3Stackings.svelte';
  import Viewer from './components/Viewer.svelte';
  import { enterStage2, enterStage3, appState, invalidateStackingState, returnToStage1, returnToStage2 } from './lib/state.svelte';
  import { cancelSharedCifWrite, seamlessState, setSharedBiounitIndex, setSharedCif, setSharedPdbCode, setSharedRangeRows, setSharedRangeSpecs, setSharedStackingSelection } from './lib/seamless.svelte';
  import { findChain } from './lib/chains';
  import { cloneRangeRows, updateRow } from './lib/ranges';
  import { mergeStackings, observeStage2RangeRows, resolveSelectedStackings } from './lib/stackings';
  import type { ChainKey, PickInfo, RangeRow, RangeSpec } from './lib/types';
  import type { MolViewer } from './lib/molstar/viewer';

  const startupQuery = new URLSearchParams(window.location.search);
  const startupPdbId = startupQuery.get('pdb');
  const autoContinue = startupQuery.get('continue') === '1';
  if (startupPdbId && /^[0-9][A-Za-z0-9]{3}$/.test(startupPdbId)) {
    appState.pdbId = startupPdbId.toUpperCase();
  }

  let seamlessStateInitialized = false;
  let startupStageRestoreHandled = false;
  let startupStage2Restored = false;
  let appliedRangeRevision = 0;
  let rangeObservation = { active: false, signature: null as string | null };

  $effect(() => {
    if (!seamlessState.ready) return;

    if (!seamlessStateInitialized) {
      seamlessStateInitialized = true;
      if (!seamlessState.pdbCodePendingWrite && !appState.pdbId) {
        appState.pdbId = seamlessState.pdbCode;
      }
      if (!seamlessState.biounitIndexPendingWrite) {
        appState.biounitIndex = seamlessState.biounitIndex;
      }
      return;
    }

    if (seamlessState.pdbCode !== appState.pdbId) {
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
    setSharedPdbCode(appState.pdbId);
    setSharedBiounitIndex(appState.biounitIndex);
  });

  $effect(() => {
    if (appState.stage !== 2 || !seamlessState.rangesReady) return;
    if (seamlessState.rangeRevision > appliedRangeRevision) {
      appliedRangeRevision = seamlessState.rangeRevision;
      appState.rows = cloneRangeRows(seamlessState.rangeRows);
      setSharedRangeSpecs(appState.rows);
      return;
    }
    setSharedRangeRows(appState.rows);
    setSharedRangeSpecs(appState.rows);
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
    if (appState.stage !== 2) return;
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
      return;
    }
    if (cell.column === 'chain') {
      const chain = findChain(appState.chains, pick.labelAsymId, pick.instanceId);
      if (chain) appState.rows = updateRow(appState.rows, row.id, { chain: chain.name });
      return;
    }
    const focusChain = appState.chains.find(chain => chain.name === row.chain);
    if (!focusChain || focusChain.labelAsymId !== pick.labelAsymId || focusChain.instanceId !== pick.instanceId) return;
    appState.rows = updateRow(appState.rows, row.id, { [cell.column]: String(pick.labelSeqId) });
  }

  function proceedToRanges() {
    if (!viewer) return;
    const chains = viewer.listChains();
    if (chains.length === 0) {
      appState.error = 'This unit has no polymer chains with residue numbering.';
      appState.stage = 1;
      return;
    }
    appState.chains = chains;
    appliedRangeRevision = seamlessState.rangeRevision;
    enterStage2(seamlessState.rangesReady ? seamlessState.rangeRows : undefined);
  }

  $effect(() => {
    if (startupStageRestoreHandled
      || !seamlessState.ready
      || !seamlessStateInitialized
      || !appState.structureReady
      || appState.loading
      || !seamlessState.rangesReady) return;

    if (!startupStage2Restored) {
      if (appState.stage === 3) {
        startupStageRestoreHandled = true;
        return;
      }
      startupStage2Restored = true;
      if (seamlessState.rangeRows.length === 0) {
        startupStageRestoreHandled = true;
        return;
      }
      if (appState.stage === 1) proceedToRanges();
    }

    if (!seamlessState.stackingsReady || !seamlessState.stackingSelectionReady) return;
    startupStageRestoreHandled = true;
    if (appState.stage !== 2) return;

    const hasSavedStackingSelection = seamlessState.stackingSelectionNonEmpty
      && !seamlessState.stackingSelectionError
      && appState.stackings.length > 0
      && appState.stackingSelection.length === appState.stackings.length;
    if (hasSavedStackingSelection) enterStage3();
  });

  function goBack() {
    setSharedRangeRows([]);
    setSharedRangeSpecs([]);
    setSharedStackingSelection([]);
    returnToStage1();
    void viewer?.render({ focus: null });
  }

  function complete(): void {
    if (!seamlessState.stackingSelectionReady || !seamlessState.stackingSelectionNonEmpty) {
      const selection = mergeStackings([], [], appState.stackings).stackingSelection;
      appState.stackingSelection = selection;
      setSharedStackingSelection(selection);
    }
    enterStage3();
  }

  function backToRanges(): void {
    setSharedStackingSelection([]);
    appState.stackingSelection = appState.stackings.map(() => false);
    returnToStage2();
  }

  $effect(() => {
    if (!viewer) return;
    if (appState.stage === 3) {
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
    <aside class="sidebar">
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
        <div class:current={appState.stage === 3} class="progress-step">
          <span class="progress-number">03</span>
          <span><strong>Stackings</strong><small>Review interactions</small></span>
        </div>
      </div>

      {#if viewerError}
        <div class="viewer-error" role="alert"><strong>Viewer unavailable</strong><span>{viewerError}</span></div>
      {/if}

      {#if appState.stage === 1}
        <Stage1Structure {viewer} {startupPdbId} {autoContinue} onContinue={proceedToRanges} />
      {:else if appState.stage === 2}
        <Stage2Ranges onBack={goBack} onComplete={complete} />
      {:else}
        <Stage3Stackings onBack={backToRanges} />
      {/if}

      <div class="sidebar-bottom"><span class="sidebar-bottom-mark">P</span><span>{appState.stage === 3 ? 'Stage 3 · Stacking review' : 'Phase 1 · Local selection'}</span><span class="sidebar-bottom-status">●</span></div>
    </aside>

    <section class="viewer-pane" aria-label="Molecular viewer workspace">
      <div class="viewer-toolbar">
        <div class="viewer-title"><span class="viewer-title-icon">◈</span><span>{appState.pdbId || 'Structure viewer'}</span>{#if appState.structureReady}<span class="loaded-badge">LOADED</span>{/if}</div>
        <div class="viewer-toolbar-right"><span class="viewer-mode"><span>●</span> {appState.stage === 1 ? 'Structure view' : appState.stage === 2 ? 'Pick mode' : 'Interaction review'}</span><span class="toolbar-divider"></span><span class="viewer-version">Mol*</span></div>
      </div>
      {#if seamlessState.rangeTransformError}
        <div class="viewer-error" role="alert"><strong>Range extraction failed</strong><span>{seamlessState.rangeTransformError}</span></div>
      {/if}
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
  </div>
</main>
