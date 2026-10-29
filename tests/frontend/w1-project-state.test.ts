import { describe, expect, it } from 'vitest';
import {
  deriveProjectStage,
  markerIsWellFormed,
  markerMatchesUpstream,
  echoIsFresh,
  deferredWriteIsCurrent,
  visibleTransformerErrors,
} from '../../src/lib/project-state';

// These are schema-1 decision documents from the W0 canonical representations.
// Result values and state-graph status deliberately do not participate in stage derivation.
const rna = { chain: 'A', start: 101, end: 104 };
const protein = { chain: 'B', start: 10, end: 30 };
const rows = [rna, protein];
const keys = ['frag-1-AC', 'frag-2-CG', 'frag-3-GU'];
const rnaBasis = { pdbcode: '1B7F', rna };

function project() {
  return {
    pdbcode: '1B7F',
    biounit: 1,
    ranges: { schema: 1, basis: { pdbcode: '1B7F', biounit: 1 }, active: true, rows },
    stacking_selection: {
      schema: 1,
      basis: { pdbcode: '1B7F', biounit: 1, rows },
      active: true,
      selected: { 'stack-A': true, 'stack-B': false },
    },
    threshold_edits: {
      schema: 1, basis: rnaBasis, active: true,
      pairs: [
        { down: 1, up: 2, sequence: 'ACG', ovrmsd: null, crmsd: null },
        { down: 2, up: 3, sequence: 'CGU', ovrmsd: null, crmsd: null },
      ],
      bases: keys.map((fragment, i) => ({ fragment: i + 1, sequence: ['AC', 'CG', 'GU'][i], nuc1_ovrmsd: null, nuc2_ovrmsd: null })),
    },
    best_fit_clip: { schema: 1, basis: rnaBasis, clip1: 0, clip2: 0 },
    interaction_matrix: {
      schema: 1, basis: ['stack-A'], active: true,
      rows: [...keys], ranges: ['range-2'], grid: [[true], [true], [false]],
    },
    gaussian_classes: { schema: 1, basis: ['stack-A'], active: true, rows: [{ key: 'stack-A', rrm: 50 }] },
  };
}

describe('W1: the shared stage is derived solely from the valid edit-cell prefix (OC P3; ST 3)', () => {
  it('restores stage 6 from mounted decisions without computed results', () => {
    expect(deriveProjectStage(project())).toBe(6);
  });

  it('keeps a blank threshold skeleton at stage 4 while computation and filling are pending', () => {
    const state = project();
    state.interaction_matrix.active = false;
    state.gaussian_classes.active = false;
    expect(deriveProjectStage(state)).toBe(4);
  });

  it.each([
    ['ranges', 1], ['stacking_selection', 2], ['threshold_edits', 3],
    ['interaction_matrix', 4], ['gaussian_classes', 5],
  ] as const)('Back deactivating %s retains the decision and yields stage %s', (name, stage) => {
    const state = project();
    const before = structuredClone(state[name]);
    state[name].active = false;
    expect(deriveProjectStage(state)).toBe(stage);
    expect({ ...state[name], active: true }).toEqual(before);
  });

  it('distinguishes an undecided marker from an intentionally empty keyed selection', () => {
    const state: any = project();
    state.stacking_selection = null;
    expect(deriveProjectStage(state)).toBe(2);
    state.stacking_selection = { ...project().stacking_selection, selected: {} };
    state.threshold_edits.active = false;
    expect(deriveProjectStage(state)).toBe(3);
  });

  it('falls back on a delayed range edit and restores on exact reversion without a computed update', () => {
    const state = project();
    state.ranges.rows = [rna, { ...protein, end: 31 }];
    expect(deriveProjectStage(state)).toBe(2);
    state.ranges.rows = rows;
    expect(deriveProjectStage(state)).toBe(6);
  });

  it('does not reactivate an orphan beyond a missing or inactive prefix', () => {
    const state = project();
    state.threshold_edits.active = false;
    expect(state.interaction_matrix.active).toBe(true);
    expect(state.gaussian_classes.active).toBe(true);
    expect(deriveProjectStage(state)).toBe(3);
  });

  it.each(['pdbcode', 'biounit'] as const)('invalidates M2 when the project %s changes', (field) => {
    const state: any = project();
    state[field] = field === 'pdbcode' ? '2ABC' : 2;
    expect(deriveProjectStage(state)).toBe(1);
  });

  it('does not restore thresholds from another project that happens to have the same RNA addresses', () => {
    const state = project();
    state.threshold_edits.basis = { ...rnaBasis, pdbcode: '2ABC' };
    expect(deriveProjectStage(state)).toBe(3);
  });

  it('keeps thresholds independent of protein ranges after the selection has been reviewed again', () => {
    const state = project();
    const newRows = [rna, { ...protein, end: 31 }];
    state.ranges.rows = newRows;
    state.stacking_selection.basis.rows = newRows;
    state.interaction_matrix.active = false;
    expect(deriveProjectStage(state)).toBe(4);
  });

  it('compares selected keys, rather than all selection keys, for matrix and class restoration', () => {
    const state = project();
    state.stacking_selection.selected['stack-B'] = true;
    expect(deriveProjectStage(state)).toBe(4);
  });

  it('checks the addressed clip extent even when old and new clips have equal lengths', () => {
    const state = project();
    state.best_fit_clip.clip1 = 1;
    state.interaction_matrix.rows = keys.slice(0, 2);
    state.interaction_matrix.grid = [[true], [false]];
    expect(deriveProjectStage(state)).toBe(4);
    state.interaction_matrix.rows = keys.slice(1);
    expect(deriveProjectStage(state)).toBe(6);
  });

  it.each([
    { rows: ['frag-2-CG', 'frag-1-AC', 'frag-3-GU'] },
    { rows: ['frag-1-AC', 'frag-1-AC', 'frag-3-GU'] },
    { ranges: ['range-3'] },
    { grid: [[true], [true]] },
    { grid: [[true], [1], [false]] },
  ])('refuses a malformed or mismatched matrix without result-dependent lock checks: %j', patch => {
    const state: any = project();
    Object.assign(state.interaction_matrix, patch);
    expect(deriveProjectStage(state)).toBe(4);
  });
});

describe('W1: decision reads are validated before use (OC P7/P8; W0 schema 1)', () => {
  it('accepts the canonical stage markers, including permitted threshold blanks', () => {
    const state = project();
    for (const kind of ['ranges', 'stacking_selection', 'threshold_edits', 'interaction_matrix', 'gaussian_classes'] as const) {
      expect(markerIsWellFormed(kind, state[kind])).toBe(true);
      expect(markerMatchesUpstream(kind, state[kind], state)).toBe(true);
    }
  });

  it.each([
    { schema: 2 }, { active: 1 }, { basis: null },
  ])('does not partially consume an unsupported or malformed marker: %j', patch => {
    const state: any = project();
    Object.assign(state.ranges, patch);
    expect(deriveProjectStage(state)).toBe(1);
  });

  it('rejects legacy bare lists instead of silently assigning a basis', () => {
    expect(markerIsWellFormed('stacking_selection', [true, false])).toBe(false);
    expect(markerIsWellFormed('ranges', rows)).toBe(false);
  });
});

describe('W1: freshness uses the actual complete input echo (OC P5; ST 6)', () => {
  const fitEcho = { pdb_code: '1b7f', rna, resids: [{ chain: 'A', resid: 101, insertion_code: '', base: 'A' }] };
  const energyEcho = { best_fits: fitEcho, protein_ranges: [protein] };

  it('allows an exact reversion and ignores object property order', () => {
    expect(echoIsFresh(fitEcho, { resids: fitEcho.resids, rna, pdb_code: '1b7f' })).toBe(true);
    expect(echoIsFresh(energyEcho, structuredClone(energyEcho))).toBe(true);
  });

  it('invalidates an energy table immediately on a local protein edit, even when fragment keys are unchanged', () => {
    expect(echoIsFresh(energyEcho, { ...energyEcho, protein_ranges: [{ ...protein, end: 31 }] })).toBe(false);
  });

  it('does not treat a missing, partial, or differently ordered input lineage as fresh', () => {
    expect(echoIsFresh(undefined, fitEcho)).toBe(false);
    expect(echoIsFresh({ pdb_code: '1b7f' }, fitEcho)).toBe(false);
    expect(echoIsFresh({ ...energyEcho, protein_ranges: [protein, rna] }, { ...energyEcho, protein_ranges: [rna, protein] })).toBe(false);
  });
});

describe('W1: deferred writes belong to their original stage and revision (ST 4.1)', () => {
  const valid = { stageAtIssue: 4, currentStage: 4, live: true, revisionAtIssue: 7, currentRevision: 7 };
  it('admits an unchanged connected intent', () => expect(deferredWriteIsCurrent(valid)).toBe(true));
  it.each([
    { currentStage: 3 }, { live: false }, { currentRevision: 8 },
  ])('drops a deferred intent after %j', patch => {
    expect(deferredWriteIsCurrent({ ...valid, ...patch })).toBe(false);
  });
});


describe('W1: state-graph errors are scoped to their stage and pending inputs (OC P6; CB 8.4)', () => {
  const failed = (name: string) => ({ path: [name], state: 'failed', exception: `${name} failed` });
  const visible = (stage: number, nodes: any[], pendingInputs = new Set<string>(), hasInitialStateGraph = true) =>
    visibleTransformerErrors({ stage, hasInitialStateGraph, pendingInputs, nodes });

  it('hides every failure until the initial state graph has arrived', () => {
    expect(visible(6, [failed('best_fit')], new Set(), false)).toEqual([]);
  });

  it.each([
    ['structure_id', 1], ['parse_ranges', 2], ['ranges_to_pdb', 2],
    ['aareduce_ranges', 3], ['detect_stacking', 3], ['refe_fragments', 3],
    ['list_fragments', 3], ['share_range_pdbs_aa', 3], ['best_fit', 4],
    ['validate_threshold_edits', 4], ['energetics_ligands', 5],
    ['energetics_range_2', 5], ['merge_energetics', 5], ['interaction_handoff', 5],
    ['validate_energy_threshold_edits', 5], ['gaussian_check', 6], ['gaussian_scores', 6],
  ] as const)('shows %s beginning at stage %s and never before it', (name, startStage) => {
    expect(visible(startStage, [failed(name)]).map(error => error.name)).toContain(name);
    if (startStage > 1) expect(visible(startStage - 1, [failed(name)])).toEqual([]);
    expect(visible(6, [failed(name)]).map(error => error.name)).toContain(name);
  });

  it('hides an old failure while a relevant local range write is unconfirmed', () => {
    expect(visible(4, [failed('best_fit')], new Set(['ranges']))).toEqual([]);
    expect(visible(4, [failed('best_fit')], new Set(['gaussian_classes'])).map(error => error.name)).toContain('best_fit');
  });

  it('shows an input conversion failure as blocked even without an exception string', () => {
    const errors = visible(3, [{ path: ['detect_stacking'], state: 'blocked', exception: '' }]);
    expect(errors).toHaveLength(1);
    expect(errors[0].state).toBe('blocked');
    expect(errors[0].message).toMatch(/blocked/i);
  });
});
