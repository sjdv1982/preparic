<script lang="ts">
  import Stage1Structure from './components/Stage1Structure.svelte';
  import Stage2Ranges from './components/Stage2Ranges.svelte';
  import Viewer from './components/Viewer.svelte';
  import { enterStage2, appState, returnToStage1 } from './lib/state.svelte';
  import { findChain } from './lib/chains';
  import { updateRow } from './lib/ranges';
  import type { ChainKey, PickInfo, RangeSpec } from './lib/types';
  import type { MolViewer } from './lib/molstar/viewer';

  let viewer = $state<MolViewer | null>(null);
  let viewerError = $state('');

  function handleViewerReady(created: MolViewer) {
    viewer = created;
    viewerError = '';
  }

  function handlePick(pick: PickInfo) {
    if (appState.stage !== 2) return;
    const cell = appState.activeCell;
    if (!cell) return;
    const row = appState.rows.find(item => item.id === cell.rowId);
    if (!row) return;
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
    enterStage2();
  }

  function goBack() {
    returnToStage1();
    void viewer?.render({ focus: null });
  }

  function complete(): void {
    // Reserved for the following phase, which will consume this selection.
  }

  $effect(() => {
    if (!viewer) return;
    if (appState.stage !== 2) {
      viewer.setGranularity('chain');
      void viewer.render({ focus: null });
      return;
    }

    const cell = appState.activeCell;
    const row = cell ? appState.rows.find(item => item.id === cell.rowId) : undefined;
    const chain = row ? appState.chains.find(item => item.name === row.chain) : undefined;
    const isRangeCell = !!cell && (cell.column === 'start' || cell.column === 'end');
    if (!isRangeCell || !chain) {
      viewer.setGranularity('chain');
      void viewer.render({ focus: null });
      return;
    }

    const alwaysVisible: RangeSpec[] = [];
    for (const pinned of appState.rows) {
      if (!pinned.alwaysVisible) continue;
      const pinnedChain = appState.chains.find(item => item.name === pinned.chain);
      const start = Number(pinned.start);
      const end = Number(pinned.end);
      if (pinnedChain && Number.isInteger(start) && Number.isInteger(end) && pinnedChain.seqIds.has(start) && pinnedChain.seqIds.has(end)) {
        alwaysVisible.push({
          labelAsymId: pinnedChain.labelAsymId,
          instanceId: pinnedChain.instanceId,
          start: Math.min(start, end),
          end: Math.max(start, end),
        });
      }
    }
    const focus: ChainKey = { labelAsymId: chain.labelAsymId, instanceId: chain.instanceId };
    viewer.setGranularity('residue');
    void viewer.render({ focus, alwaysVisible, others: appState.otherChainsMode });
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
        <div class:current={appState.stage === 2} class="progress-step">
          <span class="progress-number">02</span>
          <span><strong>Ranges</strong><small>Select residues</small></span>
        </div>
      </div>

      {#if viewerError}
        <div class="viewer-error" role="alert"><strong>Viewer unavailable</strong><span>{viewerError}</span></div>
      {/if}

      {#if appState.stage === 1}
        <Stage1Structure {viewer} onContinue={proceedToRanges} />
      {:else}
        <Stage2Ranges onBack={goBack} onComplete={complete} />
      {/if}

      <div class="sidebar-bottom"><span class="sidebar-bottom-mark">P</span><span>Phase 1 · Local selection</span><span class="sidebar-bottom-status">●</span></div>
    </aside>

    <section class="viewer-pane" aria-label="Molecular viewer workspace">
      <div class="viewer-toolbar">
        <div class="viewer-title"><span class="viewer-title-icon">◈</span><span>{appState.pdbId || 'Structure viewer'}</span>{#if appState.structureReady}<span class="loaded-badge">LOADED</span>{/if}</div>
        <div class="viewer-toolbar-right"><span class="viewer-mode"><span>●</span> {appState.stage === 1 ? 'Structure view' : 'Pick mode'}</span><span class="toolbar-divider"></span><span class="viewer-version">Mol*</span></div>
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
  </div>
</main>
