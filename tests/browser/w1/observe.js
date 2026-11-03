// Chrome DevTools MCP evaluate_script sourcePath; read-only, no app setters.
async () => {
  const { seamlessState: shared } = await import('/src/lib/seamless.svelte.ts');
  const { appState: app } = await import('/src/lib/state.svelte.ts');
  const text = selector => document.querySelector(selector)?.textContent?.trim() ?? null;
  const visible = element => !!element && element.getClientRects().length > 0
    && getComputedStyle(element).visibility !== 'hidden';
  const input = document.querySelector('#pdb-code');
  const snapshot = {
    url: location.href, at: new Date().toISOString(),
    status: shared.status, ready: shared.ready,
    derivedStage: shared.projectStage, renderedStage: app.stage,
    footer: text('.sidebar-bottom'), panel: text('.stage-content .eyebrow'),
    sharedCode: shared.pdbCode, localCode: app.pdbId,
    structureId: shared.structureId, structureIdReady: shared.structureIdReady,
    structureReady: app.structureReady, loading: app.loading,
    input: input ? { value: input.value, disabled: input.disabled, visible: visible(input) } : null,
    viewerTitle: text('.viewer-title'), viewerFooter: text('.viewer-footer-left'),
    viewerPending: text('.viewer-empty'),
    overlay: visible(document.querySelector('.sync-overlay')),
    sidebarInert: document.querySelector('.sidebar')?.inert ?? null,
    decisions: JSON.parse(JSON.stringify(shared.projectSnapshot)),
    freshness: { ...shared.computedFreshness }, pendingInputs: [...shared.pendingInputs],
    hasInitialStateGraph: shared.hasInitialStateGraph,
    graphErrors: JSON.parse(JSON.stringify(shared.stateGraphErrors)),
    decisionErrors: { ...shared.decisionErrors },
    alerts: [...document.querySelectorAll('[role="alert"]')]
      .filter(visible).map(element => element.textContent.trim()),
    notice: shared.syncNotice,
    actions: [...document.querySelectorAll('.stage-content button')].map(button => ({
      text: button.textContent.trim(), disabled: button.disabled,
      blocked: button.disabled || !!button.closest('[inert]') || !visible(button),
    })),
  };
  const failures = [];
  const check = (condition, message) => { if (!condition) failures.push(message); };
  if (!shared.ready) {
    check(snapshot.overlay, 'P9: unread/disconnected project lacks read-only overlay');
    check(snapshot.sidebarInert === true, 'P9: unread/disconnected sidebar is interactive');
    check(!Object.values(snapshot.freshness).some(Boolean), 'P5/P9: result trusted before resynchronization');
  } else {
    check(app.stage === shared.projectStage, 'P1/P3: panel state differs from derived shared stage');
    check(snapshot.footer?.includes(`Stage ${shared.projectStage} ·`), 'P3: footer differs from shared stage');
    const panelStage = snapshot.panel?.match(/Stage\s+0?(\d)/i)?.[1];
    if (shared.projectStage <= 4) check(Number(panelStage) === shared.projectStage, 'P3: visible panel differs from footer');
    if (!app.loading && !shared.pdbCodePendingWrite && snapshot.pendingInputs.length === 0) {
      check(app.pdbId === shared.pdbCode, 'P1/P9: settled local PDB code differs from accepted shared code');
      if (input) check(input.value === shared.pdbCode, 'P9: hidden/visible PDB input retains abandoned intent');
    }
    if (app.stage >= 4 && !snapshot.freshness.range_pdbs_aa) {
      check(/preparing all-atom structure/i.test(snapshot.viewerPending ?? ''), 'P11: stale S23 lacks preparing wording');
    }
  }
  if (!snapshot.hasInitialStateGraph) check(snapshot.graphErrors.length === 0, 'P6: failure shown before first state graph');
  return { pass: failures.length === 0, failures, snapshot };
}
