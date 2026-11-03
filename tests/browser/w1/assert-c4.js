// Invoke after C4 reconnect settles; fixture: accepted 1B7F at stage 2.
// Read-only assertion. It never clicks, writes shares, or repairs local state.
async () => {
  const { seamlessState: shared } = await import('/src/lib/seamless.svelte.ts');
  const { appState: app } = await import('/src/lib/state.svelte.ts');
  const failures = [];
  const check = (condition, message) => { if (!condition) failures.push(message); };
  check(shared.status === 'live' && shared.ready, 'reconnect did not finish');
  check(shared.projectStage === 2 && app.stage === 2, 'idle follower did not reach accepted stage 2');
  check(shared.pdbCode === '1B7F' && app.pdbId === '1B7F', 'abandoned 2JDQ intent survived in shared/local code');
  const input = document.querySelector('#pdb-code');
  // This component may unmount once the viewer has restored; absence is legal.
  if (input) check(input.value === '1B7F', 'mounted hidden PDB input retains 2JDQ');
  check(shared.structureIdReady && shared.structureId === '1B7F', 'stored CIF identity is not accepted 1B7F');
  check(app.structureReady && !app.loading, 'accepted viewer did not finish restoring');
  const title = document.querySelector('.viewer-title')?.textContent ?? '';
  const footer = document.querySelector('.sidebar-bottom')?.textContent ?? '';
  const panel = document.querySelector('.stage-content .eyebrow')?.textContent ?? '';
  check(title.includes('1B7F') && !title.includes('2JDQ'), 'viewer title retains abandoned structure');
  check(/Stage\s+2\s*·/.test(footer) && /Stage\s+02/.test(panel), 'stage-2 footer/panel mismatch');
  check(!shared.pendingInputs.includes('pdbcode') && !shared.pdbCodePendingWrite, 'abandoned code remains pending');
  return { pass: failures.length === 0, failures, code: shared.pdbCode,
    stage: shared.projectStage, inputValue: input?.value ?? null, title: title.trim(),
    notice: shared.syncNotice };
}
