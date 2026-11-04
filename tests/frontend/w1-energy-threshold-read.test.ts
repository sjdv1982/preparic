import { afterEach, describe, expect, it, vi } from 'vitest';
import {
  connectToSeamless,
  seamlessState,
  type SharedCell,
  type SharedValue,
  type ShareContext,
} from '../../src/lib/seamless.svelte';
import { appState } from '../../src/lib/state.svelte';

const rna = { chain: 'A', start: 101, end: 104 };
const protein = { chain: 'B', start: 10, end: 30 };
const rows = [rna, protein];
const basis = { pdbcode: '1B7F', rna };
const energyThresholdEdits = {
  schema: 1,
  basis,
  ranges: [{ row: protein, values: { 'frag-1-AC': null } }],
};

class Cell implements SharedCell {
  readonly = false;
  marker = 1;
  writes: string[] = [];
  onchange: (value: SharedValue) => void = () => undefined;

  constructor(public value: SharedValue, public url: string) {}

  set(value: string): boolean {
    this.writes.push(value);
    return true;
  }

  emit(value: SharedValue, marker = this.marker): void {
    this.value = value;
    this.marker = marker;
    this.onchange(value);
  }
}

function setup(stage: 5 | 6) {
  const doc = (value: unknown) => JSON.stringify(value);
  const values: Record<string, SharedValue> = {
    pdbcode: '1B7F',
    biounit: '1',
    ranges: doc({ schema: 1, basis: { pdbcode: '1B7F', biounit: 1 }, active: true, rows }),
    stacking_selection: doc({ schema: 1, basis: { pdbcode: '1B7F', biounit: 1, rows }, active: true, selected: {} }),
    threshold_edits: doc({
      schema: 1, basis, active: true, pairs: [],
    }),
    best_fit_clip: doc({ schema: 1, basis, clip1: 0, clip2: 0 }),
    interaction_matrix: doc({
      schema: 1, basis: [], active: true, rows: ['frag-1-AC'], ranges: ['range-2'], grid: [[true]],
    }),
    energy_threshold_edits: doc(energyThresholdEdits),
    gaussian_classes: stage === 6 ? doc({ schema: 1, basis: [], active: true, rows: [] }) : null,
    structure_id: '1B7F',
  };
  const cells = Object.fromEntries(
    Object.entries(values).map(([key, value]) => [key, new Cell(value, `/preparic/${key}`)]),
  );
  const socket = { readyState: 1 };
  const context = {
    ...cells,
    self: {
      sharelist: [],
      onsharelist: (_keys: string[]) => undefined,
      server: 'http://w1.invalid',
      ws: socket as unknown as WebSocket,
    },
  } as unknown as ShareContext;

  connectToSeamless(() => context);
  context.self.onsharelist(Object.keys(cells));
  for (const cell of Object.values(cells)) cell.emit(cell.value);
  return { cells, context };
}

afterEach(() => {
  vi.clearAllTimers();
  vi.useRealTimers();
  appState.stage = 1;
});

describe('W1 S24 energy-threshold decision reads (OC G8/G9; P7/P8)', () => {
  it.each([5, 6] as const)('validates without changing derived stage %s', stage => {
    vi.useFakeTimers();
    const { cells } = setup(stage);
    const energyCell = cells.energy_threshold_edits as Cell;
    const invalidValues = [
      '{',
      JSON.stringify({ ...energyThresholdEdits, schema: 2 }),
      JSON.stringify({ schema: 1, ranges: energyThresholdEdits.ranges }),
    ];

    expect(seamlessState.ready).toBe(true);
    expect(seamlessState.projectStage).toBe(stage);
    expect(seamlessState.energyThresholdEdits).toEqual(energyThresholdEdits);

    for (const invalid of invalidValues) {
      energyCell.emit(invalid, energyCell.marker + 1);
      expect(seamlessState.ready).toBe(true);
      expect(seamlessState.projectStage).toBe(stage);
      expect(appState.stage).toBe(stage);
      expect(seamlessState.energyThresholdEdits).toBeNull();
      expect(seamlessState.decisionErrors.energy_threshold_edits).toMatch(/stored energy threshold edits are invalid/i);

      // Unrelated reads must not clear a refused decision's explanation.
      (cells.pdbcode as Cell).emit('1B7F', (cells.pdbcode as Cell).marker + 1);
      expect(seamlessState.decisionErrors.energy_threshold_edits).toMatch(/stored energy threshold edits are invalid/i);
    }

    energyCell.emit(JSON.stringify(energyThresholdEdits), energyCell.marker + 1);
    expect(seamlessState.projectStage).toBe(stage);
    expect(seamlessState.energyThresholdEdits).toEqual(energyThresholdEdits);
    expect(seamlessState.decisionErrors.energy_threshold_edits).toBeUndefined();
    expect(Object.values(cells).flatMap(cell => (cell as Cell).writes)).toEqual([]);
  });
});
