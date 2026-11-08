import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import {
  connectToSeamless,
  deactivateProjectDecisionAndWait,
  seamlessState,
  setProjectDecision,
  setProjectDecisionDebounced,
  type SharedCell,
  type SharedValue,
  type ShareContext,
} from '../../src/lib/seamless.svelte';
import { appState } from '../../src/lib/state.svelte';
import { cloneThresholdMarker } from '../../src/lib/stage4';
import { cloneStage5Decision, stage5CanAdvance } from '../../src/lib/stage5';
import {
  classMarkerForContinue,
  copyComputedGaussianClasses,
  gaussianClassDiffCount,
  reconcileGaussianClassDraft,
  revertGaussianClassRow,
  switchGaussianClass,
  type GaussianClassMarker,
} from '../../src/lib/stage6';

const rna = { chain: 'A', start: 101, end: 102 };
const protein = { chain: 'B', start: 10, end: 30 };
const rangeRows = [rna, protein];
const rnaBasis = { pdbcode: '1B7F', rna };
const projectBasis = { pdbcode: '1B7F', biounit: 1, rows: rangeRows };
const key = '["A",101,2,30]';
const selectedKeys = [key];
const classMarker = {
  schema: 1,
  basis: selectedKeys,
  active: true,
  rows: [{ key, rrm: 90 as const, tilt30: true, disp45: false }],
};
const thresholdMarker = {
  schema: 1,
  basis: rnaBasis,
  active: true,
  pairs: [],
};
const clip = { schema: 1, basis: rnaBasis, clip1: 0, clip2: 0 };

class FakeCell implements SharedCell {
  readonly: boolean;
  marker: number;
  url?: string;
  onchange: (value: SharedValue) => void = () => undefined;
  writes: string[] = [];

  constructor(public value: SharedValue, marker = 1, readonly = false) {
    this.marker = marker;
    this.readonly = readonly;
  }

  set(value: string): boolean {
    if (this.readonly) return false;
    this.writes.push(value);
    return true;
  }

  emit(value: SharedValue, marker = this.marker): void {
    this.value = value;
    this.marker = marker;
    this.onchange(value);
  }
}

function json(value: unknown): string {
  return JSON.stringify(value);
}

function setupStage6(options: {
  configureContext?: (context: any) => void;
  shareList?: (keys: string[]) => string[];
  suppressInitialEmits?: string[];
} = {}) {
  const rangeMarker = { schema: 1, basis: { pdbcode: '1B7F', biounit: 1 }, active: true, rows: rangeRows };
  const selection = {
    schema: 1,
    basis: projectBasis,
    active: true,
    selected: { [key]: true },
  };
  const interaction = {
    schema: 1,
    basis: selectedKeys,
    active: true,
    rows: ['frag-1-AC'],
    ranges: ['range-2'],
    grid: [[true]],
  };
  const bestFits = { schema: 1, echo: { fit: '1b7f-current' }, rows: [] };
  const gaussianValidation = {
    schema: 1,
    echo: { classes: classMarker.rows, clip, selected_keys: selectedKeys, best_fits: bestFits.echo },
    rows: [],
  };
  const socket = { readyState: 1 };
  const context = {
    self: { sharelist: [] as string[], onsharelist: (_keys: string[]) => undefined, ws: socket as unknown as WebSocket },
    pdbcode: new FakeCell('1B7F'),
    biounit: new FakeCell('1'),
    ranges: new FakeCell(json(rangeMarker)),
    stacking_selection: new FakeCell(json(selection)),
    threshold_edits: new FakeCell(json(thresholdMarker)),
    best_fit_clip: new FakeCell(json(clip)),
    interaction_matrix: new FakeCell(json(interaction)),
    gaussian_classes: new FakeCell(json(classMarker)),
    energy_threshold_edits: new FakeCell(null),
    structure_id: new FakeCell('1B7F', 1, true),
    best_fits: new FakeCell(json(bestFits), 1, true),
    gaussian_validation: new FakeCell(json(gaussianValidation), 1, true),
  };
  options.configureContext?.(context);
  const connector = (_updateServer?: string | number | null, _restServer?: string | number | null, namespace?: string) => {
    if (namespace !== 'preparic') throw new Error(`Unexpected namespace: ${namespace}`);
    return context as unknown as ShareContext;
  };

  connectToSeamless(connector);
  const keys = options.shareList?.(Object.keys(context)) ?? Object.keys(context);
  context.self.onsharelist(keys);
  for (const [name, cell] of Object.entries(context)) {
    if (name === 'self' || options.suppressInitialEmits?.includes(name)) continue;
    (cell as FakeCell).emit((cell as FakeCell).value, 1);
  }
  return { context, classes: context.gaussian_classes as FakeCell, thresholds: context.threshold_edits as FakeCell };
}

const computedScores = {
  schema: 1,
  filter_set_version: 'model-1',
  filter_names: [
    'rrm50', 'rrm60', 'rrm70', 'rrm80',
    'rrm90', 'rrm90-tilt30', 'rrm90-disp45', 'rrm90-tilt30-disp45',
    'rrm95', 'rrm95-tilt30', 'rrm95-disp45', 'rrm95-tilt30-disp45',
  ],
  classes: [
    { key: 'anchor-A', rrm: 90, tilt30: true, disp45: false,
      auto_prefilters: { tilt30: true, disp45: false }, filter_name: 'rrm90-tilt30' },
    { key: 'anchor-B', rrm: 80, filter_name: 'rrm80' },
  ],
};

function stage5GateFixture() {
  const fragments = {
    schema: 1,
    echo: rnaBasis,
    fragments: [{ number: 1, key: 'frag-1-AC', sequence: 'AC', residues: [{ resid: 101 }, { resid: 102 }] }],
  };
  const rangeMarker = { schema: 1, basis: { pdbcode: '1B7F', biounit: 1 }, active: true, rows: rangeRows };
  const selection = { schema: 1, basis: projectBasis, active: true, selected: { [key]: true } };
  const thresholdEdits = {
    schema: 1, basis: rnaBasis, active: true, pairs: [],
  };
  const interactionMatrix = {
    schema: 1, basis: selectedKeys, active: true,
    rows: ['frag-1-AC'], ranges: ['range-2'], grid: [[true]],
  };
  const snapshot = {
    pdbcode: '1B7F', biounit: 1, ranges: rangeMarker, stacking_selection: selection,
    threshold_edits: thresholdEdits, best_fit_clip: clip,
    interaction_matrix: interactionMatrix, gaussian_classes: null,
  };
  const energyEdits = { schema: 1, basis: rnaBasis,
    ranges: [{ row: protein, values: { 'frag-1-AC': '1.00' } }] };
  const stackings = [{ rna_resname: 'A', rna_resid: 101, protein_range_pdb_id: 2, protein_resid: 30 }];
  const energetics = { schema: 1, ranges: [{ range: 'range-2', row: protein, rows: [{ key: 'frag-1-AC',
    reference: { energy: 1, error: null }, bestfit: { energy: 1, error: null } }] }] };
  const state = {
    ready: true, stage: 5, energeticsFresh: true, bestFitsFresh: true, geometryFresh: true,
    gaussianFresh: true, gaussianBackendFailed: false, validatorFailed: false,
    fragments, rangeRows, stackings, energetics, energyEdits,
  };
  return { snapshot, state };
}

describe('W7 class editor helpers', () => {
  it('drops flags at RRM-80 and reseeds current automatic flags on return to RRM-90/95', () => {
    const edited = { key: 'anchor-A', rrm: 90 as const, tilt30: false, disp45: true };
    const low = switchGaussianClass(edited, 80, { tilt30: true, disp45: false });
    expect(low).toEqual({ key: 'anchor-A', rrm: 80 });
    expect(switchGaussianClass(low, 90, { tilt30: true, disp45: false }))
      .toEqual({ key: 'anchor-A', rrm: 90, tilt30: true, disp45: false });

    const retained = { schema: 1, basis: ['anchor-A'], active: false, rows: [edited] };
    expect(classMarkerForContinue(retained, computedScores, ['anchor-A'])?.rows).toEqual([edited]);
    const reseeded = classMarkerForContinue(retained, computedScores, ['anchor-A', 'anchor-B']);
    expect(reseeded?.rows).toEqual([
      { key: 'anchor-A', rrm: 90, tilt30: true, disp45: false },
      { key: 'anchor-B', rrm: 80 },
    ]);
  });

  it('reverts only the chosen row and reports row-level differences', () => {
    const marker = {
      schema: 1, basis: ['anchor-A', 'anchor-B'], active: true,
      rows: [
        { key: 'anchor-A', rrm: 90, tilt30: false, disp45: true },
        { key: 'anchor-B', rrm: 50 },
      ],
    };
    expect(gaussianClassDiffCount(marker, computedScores)).toBe(2);
    const reverted = revertGaussianClassRow(marker, 'anchor-A', computedScores);
    expect(reverted?.rows).toEqual([
      { key: 'anchor-A', rrm: 90, tilt30: true, disp45: false },
      { key: 'anchor-B', rrm: 50 },
    ]);
    expect(gaussianClassDiffCount(reverted, computedScores)).toBe(1);
  });

  it('returns no copy for cancellation and prepares one full computed replacement after confirmation', () => {
    const marker = { schema: 1, basis: ['anchor-A', 'anchor-B'], active: true,
      rows: [{ key: 'anchor-A', rrm: 50 }, { key: 'anchor-B', rrm: 90, tilt30: false, disp45: true }] };
    expect(copyComputedGaussianClasses(marker, computedScores, false)).toBeNull();
    expect(copyComputedGaussianClasses(marker, computedScores, true)?.rows).toEqual([
      { key: 'anchor-A', rrm: 90, tilt30: true, disp45: false },
      { key: 'anchor-B', rrm: 80 },
    ]);
  });
});

describe('W7 first-entry readiness gate', () => {
  it('requires fresh S17/S25, complete true S24 cells, and a nonfailed T27 while allowing absent S18', () => {
    const { snapshot, state } = stage5GateFixture();
    expect(snapshot.gaussian_classes).toBeNull();
    expect(stage5CanAdvance(snapshot, state)).toEqual({ allowed: true, reasons: [] });
    expect(stage5CanAdvance(snapshot, { ...state, geometryFresh: false }).allowed).toBe(false);
    expect(stage5CanAdvance(snapshot, { ...state, gaussianFresh: false }).allowed).toBe(false);
    expect(stage5CanAdvance(snapshot, { ...state, validatorFailed: true }).allowed).toBe(false);
    expect(stage5CanAdvance(snapshot, { ...state, gaussianBackendFailed: true }).allowed).toBe(false);

    const incomplete = structuredClone(state) as any;
    incomplete.energyEdits.ranges[0].values['frag-1-AC'] = null;
    expect(stage5CanAdvance(snapshot, incomplete).allowed).toBe(false);
  });
});

describe('W7 S18 pending writes and validation freshness', () => {
  beforeEach(() => vi.useFakeTimers());
  afterEach(() => {
    vi.clearAllTimers();
    vi.useRealTimers();
    appState.stage = 1;
  });

  it('keeps the editor pending immediately, coalesces edits, and writes only at 500 ms', async () => {
    const { classes } = setupStage6();
    expect(seamlessState.projectStage).toBe(6);
    expect(seamlessState.computedFreshness.gaussian_validation).toBe(true);
    expect(seamlessState.computedFreshness).not.toHaveProperty('refe_validation');

    const first = { ...classMarker, rows: [{ key, rrm: 90 as const, tilt30: false, disp45: false }] };
    const latest = { ...classMarker, rows: [{ key, rrm: 95 as const, tilt30: true, disp45: true }] };
    const firstWrite = setProjectDecisionDebounced('gaussian_classes', first, 500);
    const latestWrite = setProjectDecisionDebounced('gaussian_classes', latest, 500);
    await expect(firstWrite).resolves.toBe(false);
    expect(classes.writes).toEqual([]);
    expect(seamlessState.pendingInputs).toContain('gaussian_classes');
    expect(seamlessState.computedFreshness.gaussian_validation).toBe(false);
    await vi.advanceTimersByTimeAsync(499);
    expect(classes.writes).toEqual([]);
    await vi.advanceTimersByTimeAsync(1);
    expect(classes.writes).toHaveLength(1);
    expect(JSON.parse(classes.writes[0])).toEqual(latest);
    expect(seamlessState.computedFreshness.gaussian_validation).toBe(false);

    classes.emit(json(latest), 2);
    await expect(latestWrite).resolves.toBe(true);
    expect(seamlessState.projectSnapshot.gaussian_classes).toEqual(latest);
    expect(seamlessState.pendingInputs).not.toContain('gaussian_classes');
  });

  it('test_s18_edit_during_inflight_put_persists_latest_marker', async () => {
    const { classes } = setupStage6();
    const first = { ...classMarker, rows: [{ key, rrm: 90 as const, tilt30: false, disp45: false }] };
    const latest = { ...classMarker, rows: [{ key, rrm: 95 as const, tilt30: true, disp45: true }] };
    const firstWrite = setProjectDecisionDebounced('gaussian_classes', first, 500);
    await vi.advanceTimersByTimeAsync(500);
    expect(classes.writes).toHaveLength(1);
    expect(JSON.parse(classes.writes[0])).toEqual(first);

    const latestWrite = setProjectDecisionDebounced('gaussian_classes', latest, 500);
    expect(latestWrite).toBe(firstWrite);
    await vi.advanceTimersByTimeAsync(500);
    expect(classes.writes).toHaveLength(1);

    classes.emit(json(first), 2);
    await Promise.resolve();
    expect(classes.writes).toHaveLength(2);
    expect(JSON.parse(classes.writes[1])).toEqual(latest);
    expect(seamlessState.pendingInputs).toContain('gaussian_classes');

    classes.emit(json(latest), 3);
    await expect(firstWrite).resolves.toBe(true);
    await expect(latestWrite).resolves.toBe(true);
    expect(seamlessState.projectSnapshot.gaussian_classes).toEqual(latest);
    expect(seamlessState.pendingInputs).not.toContain('gaussian_classes');
  });

  it('test_s18_prior_put_ack_preserves_newer_component_draft_and_subsequent_edit', () => {
    const component = readFileSync(resolve(process.cwd(), 'src/components/Stage6Gaussian.svelte'), 'utf8');
    expect(component).toContain('reconcileGaussianClassDraft(marker, untrack(() => classDraft), writePending)');
    expect(component).toContain("seamlessState.pendingInputs.includes('gaussian_classes')");

    const first: GaussianClassMarker = { ...classMarker, schema: 1, rows: [{ key, rrm: 90, tilt30: false, disp45: false }] };
    const newer: GaussianClassMarker = { ...classMarker, schema: 1, rows: [{ key, rrm: 95, tilt30: true, disp45: false }] };
    const afterFirstLocalEdit = newer;

    // At t=700, the first PUT is acknowledged while the t=600 edit is still
    // queued. The component effect sees the old ACK, but must retain the draft.
    const afterOldAck = reconcileGaussianClassDraft(first, afterFirstLocalEdit, true);
    expect(afterOldAck).toEqual(newer);

    // A third user edit at t=800 must start from the preserved t=600 value,
    // not the stale t=500 ACK that used to overwrite it.
    const third = {
      ...afterOldAck!,
      rows: [{ ...afterOldAck!.rows[0], tilt30: false }],
    };
    expect(third.rows).toEqual([{ key, rrm: 95, tilt30: false, disp45: false }]);
    expect(reconcileGaussianClassDraft(first, third, true)).toEqual(third);

    // Once the latest marker is acknowledged, the component can hydrate it.
    expect(reconcileGaussianClassDraft(third, first, false)).toEqual(third);
  });

  it('discards queued edits on a newer remote value and Back', async () => {
    const { classes } = setupStage6();
    const delayed = { ...classMarker, rows: [{ key, rrm: 90 as const, tilt30: false, disp45: true }] };
    const saved = setProjectDecisionDebounced('gaussian_classes', delayed, 500);
    const remote = { ...classMarker, rows: [{ key, rrm: 95 as const, tilt30: false, disp45: false }] };
    classes.emit(json(remote), 2);
    await expect(saved).resolves.toBe(false);
    await vi.advanceTimersByTimeAsync(500);
    expect(classes.writes).toEqual([]);
    expect(seamlessState.projectSnapshot.gaussian_classes).toEqual(remote);

    const queuedAgain = setProjectDecisionDebounced('gaussian_classes', delayed, 500);
    const back = deactivateProjectDecisionAndWait('gaussian_classes');
    const inactive = { ...remote, active: false };
    classes.emit(json(inactive), 3);
    await expect(queuedAgain).resolves.toBe(false);
    await expect(back).resolves.toBe(true);
    await vi.advanceTimersByTimeAsync(500);
    expect(classes.writes).toHaveLength(1);
    expect(JSON.parse(classes.writes[0])).toEqual(inactive);
    expect(seamlessState.projectStage).toBe(5);
  });

  it('adopts a newer conflicting class share without a reference validation cell', async () => {
    const { classes } = setupStage6();
    const delayed = { ...classMarker, rows: [{ key, rrm: 95 as const, tilt30: false, disp45: true }] };
    const saved = setProjectDecisionDebounced('gaussian_classes', delayed, 500);
    await vi.advanceTimersByTimeAsync(500);
    expect(classes.writes).toHaveLength(1);
    const remote = { ...classMarker, rows: [{ key, rrm: 90 as const, tilt30: false, disp45: true }] };
    classes.emit(json(remote), 2);
    await expect(saved).resolves.toBe(false);
    expect(seamlessState.projectSnapshot.gaussian_classes).toEqual(remote);
    expect(seamlessState.syncNotice).toMatch(/newer shared gaussian classes/i);
    await vi.advanceTimersByTimeAsync(500);
    expect(classes.writes).toHaveLength(1);

  });
});

describe('W7 S9 threshold decision snapshots', () => {
  it('copies reactive JSON markers without the structuredClone DataCloneError', () => {
    const reactive = (value: any): any => {
      if (value === null || typeof value !== 'object') return value;
      return new Proxy(value, {
        get(target, key, receiver) {
          return reactive(Reflect.get(target, key, receiver));
        },
      });
    };
    const marker = reactive(thresholdMarker);

    expect(() => structuredClone(marker)).toThrow();
    const copy = cloneThresholdMarker(marker);
    expect(copy).toEqual(thresholdMarker);
    expect(copy).not.toBe(marker);
    expect(copy.pairs).not.toBe(marker.pairs);
    expect(copy).not.toHaveProperty('bases');
  });
});

describe('W7 Stage 6 S24 reactive snapshots', () => {
  it('hydrates and edits the S24 table through the JSON proxy-safe clone helper', () => {
    const reactive = (value: any): any => {
      if (value === null || typeof value !== 'object') return value;
      return new Proxy(value, {
        get(target, key, receiver) {
          return reactive(Reflect.get(target, key, receiver));
        },
      });
    };
    const marker = reactive({
      schema: 1,
      basis: rnaBasis,
      ranges: [{ row: protein, values: { 'frag-1-AC': '-16.60' } }],
    });
    expect(() => structuredClone(marker)).toThrow();
    const copy = cloneStage5Decision(marker);
    expect(copy).toEqual({
      schema: 1,
      basis: rnaBasis,
      ranges: [{ row: protein, values: { 'frag-1-AC': '-16.60' } }],
    });
    expect(copy).not.toBe(marker);
    expect(copy.ranges).not.toBe(marker.ranges);

    const component = readFileSync(resolve(process.cwd(), 'src/components/Stage6Gaussian.svelte'), 'utf8');
    expect(component).toContain('cloneStage5Decision');
    expect(component).not.toContain('structuredClone(');
  });
});

describe('W7 initial computed-share reads', () => {
  afterEach(() => {
    vi.clearAllTimers();
    vi.unstubAllGlobals();
    appState.stage = 1;
  });

  it('reads S9 by its custom file route when the share-list key uses the route name', async () => {
    const fitThresholds = {
      schema: 1,
      echo: { pdb_code: '1b7f', rna, resids: [101, 102] },
      pairs: [],
      bases: [],
    };
    const fragments = {
      schema: 1,
      echo: { pdbcode: '1B7F', rna },
      residues: [101, 102],
      fragments: [],
    };
    const rangeSpecs = { schema: 1, pdbcode: '1B7F', biounit: 1, rows: rangeRows };
    const fetchMock = vi.fn(async (input: string | URL) => ({
      status: 200,
      ok: true,
      headers: { get: (name: string) => name.toLowerCase() === 'x-seamless-marker' ? '2' : null },
      text: async () => {
        expect(new URL(String(input)).pathname).toBe('/best-fit-thresholds.json');
        return json(fitThresholds);
      },
    }));
    vi.stubGlobal('fetch', fetchMock);

    setupStage6({
      configureContext: context => {
        context.self.server = 'http://preparic.test';
        context.range_specs = new FakeCell(json(rangeSpecs), 1, true);
        context.fragments = new FakeCell(json(fragments), 1, true);
        context.best_fit_thresholds = new FakeCell(json(fitThresholds), 1, true);
        context.best_fit_thresholds.url = '/best-fit-thresholds.json';
      },
      shareList: keys => keys.filter(key => key !== 'best_fit_thresholds'),
      suppressInitialEmits: ['best_fit_thresholds'],
    });

    await vi.waitFor(() => expect(seamlessState.computedFreshness.best_fit_thresholds).toBe(true));
    expect(fetchMock).toHaveBeenCalledTimes(1);
    expect(seamlessState.computed.best_fit_thresholds).toEqual(fitThresholds);
  });
});

it('W8b stage 6 exposes no reference-anchor editor while keeping Gaussian validation', () => {
  const source = readFileSync(resolve(process.cwd(), 'src/components/Stage6Gaussian.svelte'), 'utf8');
  const stage4 = readFileSync(resolve(process.cwd(), 'src/components/Stage4BestFits.svelte'), 'utf8');
  expect(source).not.toContain('refe_validation');
  expect(source).not.toContain('nuc1_ovrmsd');
  expect(source).not.toContain('nuc2_ovrmsd');
  expect(source).toContain('gaussian_validation');
  expect(source).toContain('<Stage4BestFits embedded={true} readOnly={true}');
  expect(stage4).toContain('readOnly = false');
  expect(stage4).toContain("editValue('pairs'");
  expect(stage4).toContain('<output aria-label="Overlap cRMSD threshold">');
  expect(stage4).toContain('<output aria-label="Overlap ovRMSD threshold">');
  expect(source).toContain('aria-label="Editable Gaussian classes by selected stacking"');
  expect(source).not.toContain('setEnergyThresholdEdits');
  expect(source).not.toContain('copyEnergyThresholds');
  expect(source).not.toContain('editEnergy');
  expect(source).not.toContain('Editable energy thresholds');
  const energyStart = source.indexOf('<section class="stage6-section" aria-label="Read-only energy threshold results">');
  const energyEnd = source.indexOf('</section>', energyStart);
  const energyResults = source.slice(energyStart, energyEnd);
  expect(energyResults).not.toMatch(/<input\b|<button\b/);
});
