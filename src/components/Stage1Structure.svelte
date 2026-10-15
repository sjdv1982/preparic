<script lang="ts">
  import { onDestroy } from 'svelte';
  import { appState } from '../lib/state.svelte';
  import { setSharedBiounitIndex, setSharedPdbCode } from '../lib/seamless.svelte';
  import type { MolViewer } from '../lib/molstar/viewer';

  let { viewer, onContinue, startupPdbId, autoContinue }: {
    viewer: MolViewer | null;
    onContinue: () => void;
    startupPdbId: string | null;
    autoContinue: boolean;
  } = $props();
  let debounce: ReturnType<typeof setTimeout> | undefined;
  let requestId = 0;
  let lastHandledPdbId: string | null = null;
  let autoContinueHandled = false;

  const pdbPattern = /^[0-9][A-Za-z0-9]{3}$/;

  $effect(() => {
    if (!viewer) return;
    const id = appState.pdbId.trim().toUpperCase();
    if (appState.assemblies.length > 0) {
      lastHandledPdbId = id;
      return;
    }
    if (id === lastHandledPdbId) return;
    queuePdbLoad(id, 0);
  });

  function handleInput(event: Event) {
    const value = (event.currentTarget as HTMLInputElement).value.trim();
    const id = value.toUpperCase();
    queuePdbLoad(id, 300);
    setSharedPdbCode(id);
  }

  function queuePdbLoad(value: string, delay: number) {
    const id = value.trim().toUpperCase();
    lastHandledPdbId = id;
    const request = ++requestId;
    if (debounce) clearTimeout(debounce);
    appState.pdbId = id;
    appState.assemblies = [];
    appState.assemblyId = '';
    appState.structureReady = false;
    appState.error = '';
    appState.loading = false;
    appState.stage = 1;
    appState.rows = [];
    appState.chains = [];
    appState.activeCell = null;
    if (viewer) void viewer.clear();
    if (!pdbPattern.test(id)) return;

    const continueAfterLoad = Boolean(
      autoContinue &&
      !autoContinueHandled &&
      startupPdbId &&
      startupPdbId.trim().toUpperCase() === id,
    );
    if (continueAfterLoad) autoContinueHandled = true;
    const load = () => {
      debounce = undefined;
      void loadEntry(id, request, continueAfterLoad);
    };
    if (delay > 0) debounce = setTimeout(load, delay);
    else load();
  }

  async function loadEntry(id: string, request: number, continueAfterLoad = false) {
    if (!viewer || request !== requestId) return;
    appState.loading = true;
    appState.error = '';
    let shouldContinue = false;
    try {
      const assemblies = await viewer.loadPdb(id);
      if (request !== requestId || !assemblies) return;
      appState.assemblies = assemblies;
      const requestedIndex = appState.biounitIndex;
      const selectedIndex = Number.isSafeInteger(requestedIndex) && requestedIndex >= 1 && requestedIndex <= assemblies.length
        ? requestedIndex
        : 1;
      const selectedAssembly = assemblies[selectedIndex - 1];
      if (!selectedAssembly) throw new Error('No biological units or polymer chains were found.');
      appState.biounitIndex = selectedIndex;
      appState.assemblyId = selectedAssembly.id;
      await viewer.showAssembly(selectedAssembly.id);
      if (request === requestId) {
        appState.structureReady = true;
        shouldContinue = continueAfterLoad;
      }
    } catch (error) {
      if (request === requestId) appState.error = error instanceof Error ? error.message : String(error);
    } finally {
      if (request === requestId) appState.loading = false;
    }
    if (shouldContinue && request === requestId) onContinue();
  }

  async function handleAssemblyChange(event: Event) {
    const id = (event.currentTarget as HTMLSelectElement).value;
    await changeAssembly(id);
  }

  async function changeAssembly(id: string) {
    const selectedIndex = appState.assemblies.findIndex(assembly => assembly.id === id) + 1;
    if (selectedIndex < 1) return;
    const request = ++requestId;
    if (debounce) clearTimeout(debounce);
    appState.assemblyId = id;
    appState.biounitIndex = selectedIndex;
    appState.structureReady = false;
    appState.loading = true;
    appState.error = '';
    setSharedBiounitIndex(selectedIndex);
    try {
      await viewer?.showAssembly(id);
      if (request === requestId) appState.structureReady = true;
    } catch (error) {
      if (request === requestId) appState.error = error instanceof Error ? error.message : String(error);
    } finally {
      if (request === requestId) appState.loading = false;
    }
  }

  $effect(() => {
    const index = appState.biounitIndex;
    const assemblies = appState.assemblies;
    if (!viewer || !appState.structureReady || appState.loading || assemblies.length === 0) return;
    if (!Number.isSafeInteger(index) || index < 1 || index > assemblies.length) {
      appState.biounitIndex = 1;
      return;
    }
    const selectedAssembly = assemblies[index - 1];
    if (selectedAssembly && selectedAssembly.id !== appState.assemblyId) {
      void changeAssembly(selectedAssembly.id);
    }
  });

  function assemblyLabel(index: number): string {
    const assembly = appState.assemblies[index];
    if (!assembly) return '';
    const names = assembly.participants.map(participant => `${participant.labelAsymId}${participant.count > 1 ? ` ×${participant.count}` : ''}`);
    const visibleNames = names.length > 6 ? `${names.slice(0, 6).join(', ')}, …` : names.join(', ');
    return `${index + 1} — ${visibleNames || 'no polymer chains'}`;
  }

  function assemblyDescription(index: number): string {
    const assembly = appState.assemblies[index];
    if (!assembly) return '';
    const names = assembly.participants.map(participant => `${participant.labelAsymId}${participant.count > 1 ? ` ×${participant.count}` : ''}`);
    return [names.join(', '), assembly.details].filter(Boolean).join(' · ');
  }

  onDestroy(() => {
    if (debounce) clearTimeout(debounce);
    requestId += 1;
  });
</script>

<section class="stage-content">
  <div class="eyebrow">Stage 01 <span>Structure</span></div>
  <h1>Choose a structure</h1>
  <p class="stage-description">Load a PDB entry, then choose the biological assembly to explore.</p>

  <label class="field-label" for="pdb-code">PDB code</label>
  <div class="input-with-prefix">
    <span class="input-prefix" aria-hidden="true">↗</span>
    <input id="pdb-code" class="text-input" type="text" maxlength="4" autocomplete="off" spellcheck="false" placeholder="e.g. 2JDQ" value={appState.pdbId} oninput={handleInput} disabled={!viewer} />
  </div>
  <div class="field-note">Four-character code from the Protein Data Bank</div>

  {#if appState.loading}
    <div class="status-card loading-card"><span class="spinner"></span><span>{appState.assemblyId ? 'Building selected unit…' : 'Fetching structure…'}</span></div>
  {:else if appState.error}
    <div class="status-card error-card" role="alert"><span class="status-icon">!</span><span>{appState.error}</span></div>
  {:else if appState.pdbId && !pdbPattern.test(appState.pdbId)}
    <div class="field-note validation-note">Enter a valid four-character PDB code.</div>
  {/if}

  {#if appState.assemblies.length > 0}
    <label class="field-label assembly-label" for="assembly-select">Biological unit</label>
    <select id="assembly-select" class="select-input" value={appState.assemblyId} onchange={handleAssemblyChange} disabled={appState.loading || !viewer}>
      {#each appState.assemblies as assembly, index (assembly.id)}
        <option value={assembly.id} title={assemblyDescription(index)}>{assemblyLabel(index)}</option>
      {/each}
    </select>
    <div class="field-note">{appState.assemblies.find(item => item.id === appState.assemblyId)?.details ?? 'Deposited biological assembly'}</div>
  {/if}

  <div class="stage-footer">
    <div class="footer-meta"><span class="ready-dot" class:ready={appState.structureReady}></span>{appState.structureReady ? 'Structure ready' : viewer ? 'Waiting for a PDB code' : 'Starting viewer…'}</div>
    <button class="button button-primary" onclick={onContinue} disabled={!appState.structureReady || appState.loading}>
      Continue <span aria-hidden="true">→</span>
    </button>
  </div>
</section>
