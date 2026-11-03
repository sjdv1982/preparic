// Read-only observation in live Vite page, via Chrome DevTools MCP only.
async () => {
  const { seamlessState: s } = await import('/src/lib/seamless.svelte.ts');
  const { appState: a } = await import('/src/lib/state.svelte.ts');
  const visible = e => !!e && e.getClientRects().length > 0;
  const table = document.querySelector('table[aria-label="Editable fragment and pair RMSD thresholds"]');
  const copy = [...document.querySelectorAll('button')].find(e => e.textContent.trim() === 'Copy computed values');
  const failures = [];
  if (s.ready && a.stage === 4) {
    if (!table) failures.push('P4.3: stage 4 table absent');
    if (!copy || copy.disabled !== (s.computedFreshness.best_fit_thresholds !== true)) failures.push('P4.3: Copy availability differs from S9 freshness');
    if (table && [...table.querySelectorAll('input')].some(e => e.disabled)) failures.push('P4.3: threshold input disabled');
  }
  if (s.ready && a.stage >= 4 && !s.computedFreshness.range_pdbs_aa) {
    const pending = [...document.querySelectorAll('.viewer-empty')].filter(visible).map(e => e.textContent).join(' ');
    if (!/preparing all-atom structure/i.test(pending)) failures.push('P11: stale S23 lacks preparing wording');
  }
  const hash = async text => [...new Uint8Array(await crypto.subtle.digest('SHA-256', new TextEncoder().encode(text)))].map(n => n.toString(16).padStart(2, '0')).join('');
  const pdbs = s.rangePdbsAa?.pdbs ?? [];
  return {
    pass: failures.length === 0, failures,
    snapshot: {
      at: new Date().toISOString(), ready: s.ready, stage: a.stage, derivedStage: s.projectStage,
      freshness: { ...s.computedFreshness }, decisions: JSON.parse(JSON.stringify(s.projectSnapshot)),
      graphErrors: JSON.parse(JSON.stringify(s.stateGraphErrors)),
      computed: JSON.parse(JSON.stringify(s.computed)),
      copy: copy ? { disabled: copy.disabled, visible: visible(copy) } : null,
      inputs: table ? [...table.querySelectorAll('input')].map(e => ({ label: e.getAttribute('aria-label'), value: e.value, disabled: e.disabled })) : [],
      tableText: table?.textContent ?? null,
      overlayFragments: [...a.stage4OverlayFragments],
      alerts: [...document.querySelectorAll('[role="alert"], .stage4-error, .stage4-warning')].filter(visible).map(e => e.textContent.trim()),
      diffMarks: table?.querySelectorAll('.stage4-diff').length ?? 0,
      s23: await Promise.all(pdbs.map(async p => ({ range: p.range, name: p.name, bytes: new TextEncoder().encode(p.pdb).length, sha256: await hash(p.pdb) }))),
    },
  };
}
