/** W5 E1/E3/M5/M6 contract acceptance. Authored without executing tests. */
import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { describe, expect, it } from 'vitest';
import {
  adaptInteractionMatrix, candidateAnchors, cloneStage5Decision, makeInteractionMatrix, reconcileEnergyThresholdEdits,
  recomputeInteractionMatrix, toggleStage5MatrixCell, validateStage5InteractionMatrix,
  stage5CanAdvance, stage5RangeIds, stage5Rows, validEnergyEdit,
  validateStage5MatrixRead, setStage5UntilNextAnchor,
} from '../../src/lib/stage5';

const keys = Array.from({ length: 9 }, (_, i) => `frag-${i + 1}-UU`);
const rna = { chain: 'A', start: 101, end: 110 };
const rangeRows = [rna, { chain: 'B', start: 1, end: 20 }, { chain: 'C', start: 5, end: 25 }];
const basis = { pdbcode: '1B7F', rna };
const fragments = { schema: 1, echo: basis, fragments: keys.map(key => ({ key })) };
const clip = { schema: 1, basis, clip1: 0, clip2: 0 };
const selection = { schema: 1, active: true,
  basis: { pdbcode: '1B7F', biounit: 1, rows: rangeRows }, selected: { 'anchor-A': true } };
const ranges = ['range-2', 'range-3'];
const anchors = [{ resid: 104, range: 'range-2', stackings: ['anchor-A'], locks: [{ fragment: keys[3], nucleotide: 1 }] },
  { resid: 105, range: 'range-3', stackings: ['anchor-A'], locks: [{ fragment: keys[4], nucleotide: 1 }] }];

function matrix(rows = keys): any {
  return { schema: 1, basis: ['anchor-A'], active: true, rows: [...rows], ranges: [...ranges],
    grid: rows.map(() => [false, false]) };
}
function energetics(): any {
  return { schema: 1, ranges: ranges.map((range, index) => ({ range, row: rangeRows[index + 1],
    rows: keys.map(key => ({ key, reference: { energy: 3, error: null }, bestfit: { energy: -1, error: null } })) })) };
}
function gateFixture(): { snapshot: any; state: any } {
  const value = matrix();
  value.grid[3][0] = true;
  const edits = reconcileEnergyThresholdEdits(null, fragments, rangeRows)!;
  edits.ranges[0].values[keys[3]] = '-16.60';
  const snapshot = { pdbcode: '1B7F', biounit: 1,
    ranges: { schema: 1, active: true, basis: { pdbcode: '1B7F', biounit: 1 }, rows: rangeRows },
    stacking_selection: structuredClone(selection),
    threshold_edits: { schema: 1, basis, active: true,
      pairs: keys.slice(1).map((_, i) => ({ down: i + 1, up: i + 2, sequence: 'UUU', ovrmsd: '1.00', crmsd: '1.00' })) },
    best_fit_clip: clip, interaction_matrix: value, gaussian_classes: null };
  return { snapshot, state: { ready: true, stage: 5, energeticsFresh: true, bestFitsFresh: true,
    geometryFresh: true, gaussianFresh: true, gaussianClassesValid: true, validatorFailed: false,
    fragments, rangeRows, energetics: energetics(), energyEdits: edits } };
}

describe('W5 E1 addressed S24 reconciliation', () => {
  it('opens empty and preserves values by row when ranges change index', () => {
    const edits = reconcileEnergyThresholdEdits(null, fragments, rangeRows)!;
    expect(edits.ranges.map(entry => entry.values)).toEqual(rangeRows.slice(1).map(() => Object.fromEntries(keys.map(key => [key, null]))));
    edits.ranges[0].values[keys[0]] = '-16.60';
    edits.ranges[1].values[keys[1]] = '21.00';
    const before = structuredClone(edits);
    const added = { chain: 'D', start: 10, end: 30 };
    const actual = reconcileEnergyThresholdEdits(edits, fragments, [rna, rangeRows[2], rangeRows[1], added])!;
    expect(actual.ranges[0].values[keys[1]]).toBe('21.00');
    expect(actual.ranges[1].values[keys[0]]).toBe('-16.60');
    expect(Object.values(actual.ranges[2].values).every(value => value === null)).toBe(true);
    expect(edits).toEqual(before);
    const changed = reconcileEnergyThresholdEdits(edits, fragments, [rna, { ...rangeRows[1], end: 21 }, rangeRows[2]])!;
    expect(Object.values(changed.ranges[0].values).every(value => value === null)).toBe(true);
    expect(changed.ranges[1].values[keys[1]]).toBe('21.00');
  });
  it.each([{ ...basis, pdbcode: 'OTHER' }, { ...basis, rna: { ...rna, end: 111 } }])('resets the whole table on stale basis %j', echo => {
    const edits = reconcileEnergyThresholdEdits(null, fragments, rangeRows)!;
    edits.ranges[0].values[keys[0]] = '4.00';
    const result = reconcileEnergyThresholdEdits(edits, { ...fragments, echo }, rangeRows)!;
    expect(result.ranges.every(entry => Object.values(entry.values).every(value => value === null))).toBe(true);
  });
  it.each(['-16.60', '0', '0.00', '+21', '0.1'])('accepts finite ≤2-decimal edit %s', value => expect(validEnergyEdit(value)).toBe(true));
  it.each(['1.234', 'NaN', 'Infinity', '-Infinity', '', 'bad', '1e2'])('rejects rather than rounds %s', value => expect(validEnergyEdit(value)).toBe(false));
});

describe('W5 M5 structural read refusal', () => {
  it('copies reactive S14 and S24 JSON markers without structuredClone proxy errors', () => {
    const reactive = (value: any): any => {
      if (value === null || typeof value !== 'object') return value;
      return new Proxy(value, {
        get(target, key, receiver) {
          return reactive(Reflect.get(target, key, receiver));
        },
      });
    };
    const markers = [matrix(), reconcileEnergyThresholdEdits(null, fragments, rangeRows)!];
    for (const marker of markers) {
      const proxied = reactive(marker);
      expect(() => structuredClone(proxied)).toThrow();
      const copy = cloneStage5Decision(proxied as any);
      expect(copy).toEqual(marker);
      expect(copy).not.toBe(proxied);
      expect(copy.ranges).not.toBe(proxied.ranges);
    }
  });

  it('builds addressed rows and ranges from clip decisions', () => {
    expect(stage5Rows(fragments, { clip1: 1, clip2: 2 })).toEqual(keys.slice(1, -2));
    expect(stage5RangeIds(rangeRows)).toEqual(ranges);
    const stacks = [{ rna_resname: 'G', rna_resid: 104, protein_range_pdb_id: 2, protein_resid: 8 },
      { rna_resname: 'G', rna_resid: 105, protein_range_pdb_id: 3, protein_resid: 9 }];
    const selected = { ...selection, selected: Object.fromEntries(stacks.map(row =>
      [JSON.stringify([row.rna_resname, row.rna_resid, row.protein_range_pdb_id, row.protein_resid]), true])) };
    const mapped = { ...fragments, fragments: keys.map((key, index) => ({ key,
      residues: [{ resid: 101 + index }] })) };
    const expected = matrix();
    expected.basis = Object.keys(selected.selected).sort();
    expect(makeInteractionMatrix(mapped, clip, rangeRows, selected, stacks)).toEqual(expected);
  });
  it.each(['schema', 'row-count', 'ragged', 'zero', 'one', 'null', 'string', 'duplicate-row', 'order', 'range-order', 'duplicate-range', 'extent'])('refuses %s without mutating the read', kind => {
    const value = matrix();
    if (kind === 'schema') value.schema = 2;
    if (kind === 'row-count') value.grid.pop();
    if (kind === 'ragged') value.grid[0].pop();
    if (['zero', 'one', 'null', 'string'].includes(kind)) value.grid[0][0] = ({ zero: 0, one: 1, null: null, string: 'true' } as any)[kind];
    if (kind === 'duplicate-row') value.rows[0] = value.rows[1];
    if (kind === 'order') value.rows.reverse();
    if (kind === 'range-order') value.ranges.reverse();
    if (kind === 'duplicate-range') value.ranges[1] = value.ranges[0];
    if (kind === 'extent') value.rows.pop();
    const before = structuredClone(value);
    expect(validateStage5MatrixRead(value, keys, ranges)).toBeTypeOf('string');
    expect(value).toEqual(before);
  });
  it('refuses an opposite-end clip with the same row count', () => {
    expect(validateStage5MatrixRead(matrix(keys.slice(1)), keys.slice(0, -1), ranges)).not.toBeNull();
  });
});

describe('W5 M6 clip adaptation', () => {
  it('adapts a shifted extent when one old edge is removed and the opposite edge is uncovered', () => {
    const old = matrix(keys.slice(2, 7));
    old.grid.forEach((row: boolean[], index: number) => {
      row[0] = index < 4; // 3' cap stopped short of the old boundary
      row[1] = index > 0; // 5' cap stopped short; 3' end touched the old boundary
    });
    const energy = energetics();
    energy.ranges.forEach((range: any) => range.rows.forEach((row: any) => {
      row.reference.energy = 3;
      row.bestfit.energy = 3;
    }));
    const oldInner = structuredClone(old.grid.slice(1, 4));
    const shiftedRows = keys.slice(1, 6); // adds frag-2, drops frag-7
    const shifted = adaptInteractionMatrix(old, shiftedRows, anchors, energy)!;
    expect(shifted.rows).toEqual(shiftedRows);
    expect(shifted.grid.map((row: boolean[]) => row[0])).toEqual([true, true, true, true, true]);
    expect(shifted.grid.map((row: boolean[]) => row[1])).toEqual([false, false, true, true, true]);
    expect(shifted.grid.slice(2, 5)).toEqual(oldInner);
    expect(validateStage5InteractionMatrix(shifted, shiftedRows, anchors, energy)).toBeNull();
  });

  it.each(['high', 'low', 'error'])('grows only touched outer ends, with %s new rows', stop => {
    const old = matrix(keys.slice(2, 7));
    old.grid.forEach((row: boolean[]) => { row[0] = true; });
    old.grid.slice(1, -1).forEach((row: boolean[]) => { row[1] = true; });
    const energy = energetics();
    if (stop !== 'high') {
      energy.ranges[0].rows[1].reference = { energy: stop === 'low' ? 2.999 : null, error: stop === 'error' ? 'bad ligand' : null };
      energy.ranges[0].rows[1].bestfit = { energy: 2.999, error: null };
    }
    const before = structuredClone(old);
    const result = adaptInteractionMatrix(old, keys, anchors, energy)!;
    expect(result).not.toBeNull();
    expect(result.grid.map(row => row[0])).toEqual(keys.map((_, index) => stop === 'high' || index >= 2));
    expect(result.grid.map(row => row[1])).toEqual(keys.map((_, index) => index >= 3 && index <= 5));
    expect(result.grid.slice(2, 7)).toEqual(old.grid);
    expect(adaptInteractionMatrix(result, old.rows, anchors, energy)).toEqual(old);
    expect(adaptInteractionMatrix(old, old.rows, anchors, energy)).toEqual(old);
    expect(adaptInteractionMatrix(old, keys.slice(5), anchors, energy)).toBeNull();
    expect(old).toEqual(before);
  });
});

describe('W5 E3 stage 5→6 gates', () => {
  it('allows valid typed true cells without a computed-threshold freshness dependency', () => {
    const { snapshot, state } = gateFixture();
    expect(stage5CanAdvance(snapshot, state)).toEqual({ allowed: true, reasons: [] });
  });
  it.each(['ready', 'energeticsFresh', 'bestFitsFresh', 'geometryFresh', 'gaussianFresh'])('requires %s', key => {
    const { snapshot, state } = gateFixture();
    state[key] = false;
    expect(stage5CanAdvance(snapshot, state).allowed).toBe(false);
  });
  it('does not require an existing Gaussian class table for first entry to stage 6', () => {
    const { snapshot, state } = gateFixture();
    expect(snapshot.gaussian_classes).toBeNull();
    expect(stage5CanAdvance(snapshot, state)).toEqual({ allowed: true, reasons: [] });
  });
  it('closes first entry when best-fit stacking or Gaussian scoring failed', () => {
    const { snapshot, state } = gateFixture();
    expect(stage5CanAdvance(snapshot, { ...state, gaussianBackendFailed: true }).allowed).toBe(false);
  });
  it('blocks a failed T27 and the wrong stage', () => {
    const { snapshot, state } = gateFixture();
    expect(stage5CanAdvance(snapshot, { ...state, validatorFailed: true }).allowed).toBe(false);
    expect(stage5CanAdvance(snapshot, { ...state, stage: 4 }).allowed).toBe(false);
  });
  it.each([null, '1.234', 'NaN'])('requires a valid value for every true cell: %j', value => {
    const { snapshot, state } = gateFixture();
    state.energyEdits.ranges[0].values[keys[3]] = value;
    expect(stage5CanAdvance(snapshot, state).allowed).toBe(false);
  });
  it('ignores empty false cells but closes when a cap makes one true', () => {
    const { snapshot, state } = gateFixture();
    expect(stage5CanAdvance(snapshot, state).allowed).toBe(true);
    snapshot.interaction_matrix.grid[4][0] = true;
    expect(stage5CanAdvance(snapshot, state).allowed).toBe(false);
    state.energyEdits.ranges[0].values[keys[4]] = '0.00';
    expect(stage5CanAdvance(snapshot, state).allowed).toBe(true);
  });
  it('refuses a matrix that does not match the current clip', () => {
    const { snapshot, state } = gateFixture();
    snapshot.best_fit_clip = { ...clip, clip1: 1 };
    expect(stage5CanAdvance(snapshot, state).allowed).toBe(false);
  });
  it('refuses a true cell with a row-local scoring error', () => {
    const { snapshot, state } = gateFixture();
    state.energetics.ranges[0].rows[3].bestfit = { energy: null, error: 'bad ligand' };
    expect(stage5CanAdvance(snapshot, state).allowed).toBe(false);
  });
  it('refuses thresholds from another RNA/PDB basis', () => {
    const { snapshot, state } = gateFixture();
    state.energyEdits = { ...state.energyEdits, basis: { ...basis, pdbcode: 'OTHER' } };
    expect(stage5CanAdvance(snapshot, state).allowed).toBe(false);
  });
});


describe('W5 P5.3/P5.4 range columns, anchor locks, and Recompute', () => {
  const mapped = { ...fragments, fragments: keys.map((key, index) => ({ key,
    residues: [{ resid: 101 + 2 * index }, { resid: 102 + 2 * index }] })) };
  function anchorAt(row: number, range = 'range-2', locks = [row]): any {
    return { resid: 99 + 2 * row, range, stackings: [`a${row}-${range}`],
      locks: locks.map(n => ({ fragment: keys[n - 1], nucleotide: 1 })) };
  }
  function emptyFor(items: any[]): any {
    return { ...matrix(), basis: items.flatMap(item => item.stackings), ranges: [...ranges],
      grid: keys.map(() => ranges.map(() => false)) };
  }
  function vector(): any {
    const result = energetics();
    const values: Array<[number, number] | null> = [[4, 0], [2.999, -2.999], [3, 0], [0, 0], [0, -3], [25, 0], null, [-4, 0], [0, 0]];
    result.ranges[0].rows.forEach((row: any, index: number) => {
      row.reference = { energy: values[index]?.[0] ?? null, error: values[index] === null ? 'bad ligand' : null };
      row.bestfit = { energy: values[index]?.[1] ?? null, error: values[index] === null ? 'bad ligand' : null };
    });
    result.ranges[1].rows.forEach((row: any, index: number) => {
      row.reference = { energy: index === 4 ? 20 : 0, error: null };
      row.bestfit = { energy: 0, error: null };
    });
    return result;
  }
  function columnRows(value: any, column: number): string[] {
    return value.rows.filter((_: string, row: number) => value.grid[row][column]);
  }
  it('keeps unique range columns while preserving every selected stacking job and anchor lock', () => {
    const addresses = [[115, 2, 8], [107, 3, 9], [105, 2, 10], [105, 2, 11]];
    const stacks = addresses.map(([resid, range, protein]) => ({ rna_resname: 'G', rna_resid: resid,
      protein_range_pdb_id: range, protein_resid: protein }));
    const originalKeys = addresses.map(([resid, range, protein]) => JSON.stringify(['G', resid, range, protein]));
    const selected = { ...selection, selected: Object.fromEntries(originalKeys.map(key => [key, true])) };
    const before = structuredClone(selected);
    const actual = candidateAnchors(mapped, selected, stacks, keys)!;
    expect(actual.map(item => item.resid)).toEqual([105, 107, 115]);
    expect(actual.map(item => item.range)).toEqual(['range-2', 'range-3', 'range-2']);
    expect(actual.map(item => item.stackings.length)).toEqual([2, 1, 1]);
    expect(actual.flatMap(item => item.stackings).sort()).toEqual([...originalKeys].sort());
    const value = makeInteractionMatrix(mapped, clip, rangeRows, selected, stacks)!;
    expect(value.ranges).toEqual(['range-2', 'range-3']);
    expect(value.grid.every(row => row.length === 2)).toBe(true);
    expect(value.basis).toEqual([...originalKeys].sort());
    expect(selected).toEqual(before);
    expect(validateStage5MatrixRead(value, keys, ['range-2', 'range-3'])).toBeNull();
  });
  it('maps a shared residue to both fragment candidates and their nucleotide positions', () => {
    const shared = structuredClone(mapped);
    shared.fragments[2].residues[1].resid = 107;
    const stack = { rna_resname: 'G', rna_resid: 107, protein_range_pdb_id: 2, protein_resid: 8 };
    const selected = { ...selection, selected: { [JSON.stringify(['G', 107, 2, 8])]: true } };
    expect(candidateAnchors(shared, selected, [stack], keys)![0].locks).toEqual([
      { fragment: keys[2], nucleotide: 2 }, { fragment: keys[3], nucleotide: 1 },
    ]);
    expect(candidateAnchors(shared, selected, [stack], keys.slice(3))![0].locks).toEqual([
      { fragment: keys[3], nucleotide: 1 },
    ]);
    expect(candidateAnchors(shared, selected, [stack], keys.slice(0, 3))![0].locks).toEqual([
      { fragment: keys[2], nucleotide: 2 },
    ]);
    expect(candidateAnchors(shared, selected, [stack], keys.slice(4))).toBeNull();
  });
  it('W8a restores an uncovered low-energy lock and retains the selected stacking and interior cap', () => {
    const full = [anchorAt(4, 'range-2', [3, 4]), anchorAt(6, 'range-3')];
    const clipped = [{ ...full[0], locks: [full[0].locks[1]] }, full[1]];
    const energy = energetics();
    energy.ranges.forEach((column: any) => column.rows.forEach((row: any) => {
      row.reference = { energy: 0, error: null }; row.bestfit = { energy: 0, error: null };
    }));
    const old = { ...emptyFor(full), rows: keys.slice(3), grid: keys.slice(3).map(() => [false, false]) };
    const value = recomputeInteractionMatrix(old, clipped, energy)!;
    value.grid[1][1] = true;
    const before = structuredClone(value);
    const restored = adaptInteractionMatrix(value, keys, full, energy)!;
    expect(restored.basis).toEqual(value.basis);
    expect(restored.grid[2][0]).toBe(true);
    expect(restored.grid.slice(3)).toEqual(value.grid);
    expect(adaptInteractionMatrix(restored, value.rows, clipped, energy)).toEqual(value);
    expect(value).toEqual(before);
    energy.ranges[0].rows[2].bestfit = { energy: null, error: 'bad candidate' };
    const restoredWithError = adaptInteractionMatrix(value, keys, full, energy)!;
    expect(restoredWithError.grid[2][0]).toBe(true);
    expect(restoredWithError.grid.slice(3)).toEqual(value.grid);
    expect(validateStage5InteractionMatrix(restoredWithError, keys, full, energy, ranges)).toMatch(/107/);
    expect(validateStage5InteractionMatrix(restored, keys, full, energy, ranges)).toMatch(/107/);
    expect(validateStage5InteractionMatrix(value, value.rows, clipped, energy, ranges)).toBeNull();
  });
  it('M1 includes exact magnitude 3 and positive clash, excludes low/error stop rows and unreachable rows', () => {
    const items = [anchorAt(4)];
    const energy = vector();
    const actual = recomputeInteractionMatrix(emptyFor(items), items, energy)!;
    expect(actual.ranges).toEqual(ranges);
    expect(columnRows(actual, 0)).toEqual(keys.slice(2, 6));
    expect(validateStage5InteractionMatrix(actual, keys, items, energy)).toBeNull();
    const clipped = { ...emptyFor(items), rows: keys.slice(3, 6), grid: keys.slice(3, 6).map(() => ranges.map(() => false)) };
    expect(columnRows(recomputeInteractionMatrix(clipped, items, energy), 0)).toEqual(keys.slice(3, 6));
  });
  it('M2 allows independent ranges to overlap and keeps unanchored range columns empty', () => {
    const items = [anchorAt(4), anchorAt(5, 'range-3')];
    const actual = recomputeInteractionMatrix(emptyFor(items), items, vector())!;
    expect(columnRows(actual, 0)).toEqual(keys.slice(2, 6));
    expect(columnRows(actual, 1)).toEqual([keys[4]]);
    expect(actual.ranges).toEqual(ranges);
    expect(columnRows(actual, 2)).toEqual([]);
  });
  it.each([[[3, 4], [5, 6]], [[3, 4], [4, 5]]])('M3 joins touching/overlapping locks in their shared range column (%j)', (left, right) => {
    const items = [anchorAt(3, 'range-2', left), anchorAt(5, 'range-2', right)];
    const actual = recomputeInteractionMatrix(emptyFor(items), items, energetics())!;
    expect(actual.ranges).toEqual(ranges);
    expect(columnRows(actual, 0)).toEqual(keys);
    expect(columnRows(actual, 1)).toEqual([]);
    expect(validateStage5InteractionMatrix(actual, keys, items, energetics())).toBeNull();
    expect(toggleStage5MatrixCell(actual, 3, 0, items, energetics())).toBeNull(); // cannot remove any joined lock
  });
  it('refuses removing a lock, discontiguous growth, an error cell, and out-of-bounds edits without mutation', () => {
    const items = [anchorAt(4)];
    const energy = vector();
    const value = recomputeInteractionMatrix(emptyFor(items), items, energy)!;
    const before = structuredClone(value);
    expect(toggleStage5MatrixCell(value, 3, 0, items, energy)).toBeNull();
    expect(toggleStage5MatrixCell(value, 0, 0, items, energy)).toBeNull();
    expect(toggleStage5MatrixCell(value, 6, 0, items, energy)).toBeNull();
    expect(toggleStage5MatrixCell(value, -1, 0, items, energy)).toBeNull();
    expect(toggleStage5MatrixCell(value, 0, 1, items, energy)).toBeNull();
    expect(value).toEqual(before);
    const bad = structuredClone(value);
    bad.grid[0][0] = true;
    expect(validateStage5InteractionMatrix(bad, keys, items, energy)).not.toBeNull();
    bad.grid[0][0] = false;
    bad.grid[6][0] = true;
    expect(validateStage5InteractionMatrix(bad, keys, items, energy)).toMatch(/error/i);
    expect(recomputeInteractionMatrix(emptyFor([anchorAt(7)]), [anchorAt(7)], energy)).toBeNull();
  });
  it('refuses duplicate range columns rather than normalizing an old marker', () => {
    const items = [anchorAt(3), anchorAt(6)];
    const value = emptyFor(items);
    value.ranges = ['range-2', 'range-2'];
    value.grid = keys.map(() => [false, false]);
    expect(validateStage5InteractionMatrix(value, keys, items, energetics())).not.toBeNull();
    expect(adaptInteractionMatrix(value, keys, items, energetics())).toBeNull();
    expect(recomputeInteractionMatrix(value, items, energetics())).toBeNull();
  });

  it('Until next anchor joins and can unjoin disjoint same-range locks', () => {
    const items = [anchorAt(3), anchorAt(6)];
    const energy = energetics();
    for (const range of energy.ranges) {
      for (const row of range.rows) {
        row.reference.energy = 0.5;
        row.bestfit.energy = 0.5;
      }
    }
    const separated = recomputeInteractionMatrix(emptyFor(items), items, energy)!;
    separated.grid[1][0] = true; // retain an outer cap before the first lock
    separated.grid[6][0] = true; // retain an outer cap after the second lock
    const joined = setStage5UntilNextAnchor(separated, 0, true, items, energy)!;
    expect(joined.grid.map(row => row[0])).toEqual([false, true, true, true, true, true, true, false, false]);
    expect(joined.grid.map(row => row[1])).toEqual(Array(9).fill(false));
    const split = setStage5UntilNextAnchor(joined, 0, false, items, energy)!;
    expect(split.grid.map(row => row[0])).toEqual([false, true, true, false, false, true, true, false, false]);
    expect(split.grid.map(row => row[1])).toEqual(Array(9).fill(false));
    expect(validateStage5InteractionMatrix(split, keys, items, energy)).toBeNull();
  });

  it('keeps four anchors in two shared range columns and links successive same-range locks', () => {
    const items = [anchorAt(3, 'range-2'), anchorAt(4, 'range-3'),
      anchorAt(7, 'range-2'), anchorAt(8, 'range-3')];
    const energy = energetics();
    for (const range of energy.ranges) for (const row of range.rows) {
      row.reference.energy = 0.5;
      row.bestfit.energy = 0.5;
    }
    const separated = recomputeInteractionMatrix(emptyFor(items), items, energy)!;
    expect(separated.ranges).toEqual(ranges);
    expect(separated.grid).toHaveLength(keys.length);
    expect(separated.grid.every(row => row.length === 2)).toBe(true);
    expect(separated.grid[2][0]).toBe(true);
    expect(separated.grid[6][0]).toBe(true);
    expect(separated.grid[3][1]).toBe(true);
    expect(separated.grid[7][1]).toBe(true);

    const joined = setStage5UntilNextAnchor(separated, 0, true, items, energy)!;
    expect(joined.grid.slice(2, 7).every(row => row[0])).toBe(true);
    expect(joined.grid[4][1]).toBe(false);
    const split = setStage5UntilNextAnchor(joined, 0, false, items, energy)!;
    expect(split.grid[2][0]).toBe(true);
    expect(split.grid[6][0]).toBe(true);
    expect(split.grid.slice(3, 6).every(row => !row[0])).toBe(true);
    expect(validateStage5InteractionMatrix(split, keys, items, energy)).toBeNull();
  });

  it.each([[[3, 4], [5, 6]], [[3, 4], [4, 5]]])(
    'forces Until next anchor on when selected locks touch or overlap (%j)', (left, right) => {
      const items = [anchorAt(3, 'range-2', left), anchorAt(5, 'range-2', right)];
      const energy = energetics();
      const joined = recomputeInteractionMatrix(emptyFor(items), items, energy)!;
      expect(joined.grid.map(row => row[0])).toEqual(keys.map(() => true));
      expect(setStage5UntilNextAnchor(joined, 0, false, items, energy)).toBeNull();
    },
  );
});

describe('W4/W5 integrated contract wiring', () => {
  const source = (path: string) => readFileSync(resolve(process.cwd(), path), 'utf8');

  it('4→5 action requires fresh computed thresholds for typed pair thresholds', () => {
    const app = source('src/App.svelte');
    const start = app.indexOf('async function continueToStage5()');
    const end = app.indexOf('async function backToStage4()', start);
    const transition = app.slice(start, end);
    expect(transition).toContain('canContinueStage4');
    expect(transition).toContain('thresholdsFresh: seamlessState.computedFreshness.best_fit_thresholds');
  });

  it('renders one accessible range column with the selected anchors and locks in its header', () => {
    const stage5 = source('src/components/Stage5Energetics.svelte');
    expect(stage5).toContain('cloneStage5Decision');
    expect(stage5).not.toContain('structuredClone(');
    expect(stage5).toContain('aria-label="ATTRACT scores, thresholds, and interaction matrix by protein range"');
    expect(stage5).toContain('aria-label={`Show ${rangeId}, chain ${proteinRange.chain} ${proteinRange.start} to ${proteinRange.end}`}');
    expect(stage5).toContain('Range {item.index + 2} · chain {item.row.chain}');
    expect(stage5).toContain('Anchor {anchor.resid}:');
    expect(stage5).toContain('lock.fragment} (nt ${lock.nucleotide})');
    expect(stage5).toContain('Join anchor ${anchor.resid} to anchor ${nextAnchor.resid} in ${item.rangeId}');
    expect(stage5).not.toContain('aria-label="Per-anchor interaction selections"');
  });

  it('test_continue_4_to_5_rechecks_fresh_energetics_before_any_write', () => {
    const app = source('src/App.svelte');
    const start = app.indexOf('async function continueToStage5()');
    const end = app.indexOf('async function backToStage4()', start);
    const transition = app.slice(start, end);
    const gate = transition.indexOf('canContinueStage4');
    expect(gate).toBeGreaterThanOrEqual(0);
    expect(transition).toContain('energeticsFresh: seamlessState.computedFreshness.energetics');
    for (const write of ['setEnergyThresholdEditsAwaited', 'prepareProjectTransition', "setProjectDecisionAwaited('interaction_matrix'"]) {
      expect(gate).toBeLessThan(transition.indexOf(write));
    }
  });

  it('test_continue_4_to_5_rechecks_threshold_completeness_and_clip', () => {
    const app = source('src/App.svelte');
    const stage4 = source('src/lib/stage4.ts');
    const start = app.indexOf('async function continueToStage5()');
    const end = app.indexOf('async function backToStage4()', start);
    const transition = app.slice(start, end);
    expect(transition).toContain('canContinueStage4');
    expect(transition).toContain('selectedCandidateGroups');
    expect(stage4).toContain('thresholdsCompleteForClip');
    expect(stage4).toContain('validClip(clip, fragmentCount)');
  });

  it('Stage 5 keeps S24 editable and omits the S11 RMSD editor', () => {
    const stage5 = source('src/components/Stage5Energetics.svelte');
    const importedChildren = [...stage5.matchAll(/import\s+(\w+)\s+from\s+['"](.+\.svelte)['"]/g)]
      .filter(match => new RegExp(`<${match[1]}[\\s/>]`).test(stage5))
      .map(match => readFileSync(resolve(process.cwd(), 'src/components', match[2]), 'utf8'));
    const editor = [stage5, ...importedChildren].join('\n');
    expect(editor).toContain('energy threshold');
    expect(editor).toMatch(/<input[^>]*energy threshold/i);
    expect(editor).not.toMatch(/nuc1_ovrmsd|crmsd|ovrmsd|Overlap cRMSD threshold|Overlap ovRMSD threshold/i);
  });

  it('test_stage5_displays_computed_energy_threshold_beside_edit', () => {
    const stage5 = source('src/components/Stage5Energetics.svelte');
    const cell = stage5.match(/<td class="stage5-threshold-cell"[^>]*>[\s\S]*?<\/td>/)?.[0] ?? '';
    expect(cell).toMatch(/<input[^>]*energy threshold/i);
    expect(cell).toMatch(/<(?:span|small)[^>]*>[^<]*\{computed\}[^<]*<\/(?:span|small)>/s);
  });

  it('renders the Until next anchor control for stage 5', () => {
    const stage5 = source('src/components/Stage5Energetics.svelte');
    expect(stage5).toContain('Until next anchor');
    expect(stage5).toContain('setStage5UntilNextAnchor');
  });
});
