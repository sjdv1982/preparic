import { beforeEach, describe, expect, it } from 'vitest';
import { appState, clearStackingState, enterStage2, enterStage3, invalidateStackingState, returnToStage1, returnToStage2 } from '../../src/lib/state.svelte';
import {
  mergeStackings,
  normalizeStackings,
  observeStage2RangeRows,
  resolveSelectedStackings,
  stackingKey,
  uniqueStackingParticipantRanges,
  type Stage2RangeObservation,
} from '../../src/lib/stackings';
import { stackingDashSegments } from '../../src/lib/molstar/stacking-geometry';
import type { ChainInfo, RangeRow, Stacking, StackingViewerInput } from '../../src/lib/types';

function stacking(values: Partial<Stacking> = {}): Stacking {
  return {
    rna_resname: 'G',
    rna_resid: 4,
    protein_range_pdb_id: 2,
    protein_resid: 93,
    distance: 4.2,
    tilt_angle: 22,
    slip_angle: 19,
    rna_centroid: [1, 2, 3],
    protein_centroid: [4, 5, 6],
    ...values,
  };
}

function rangeRow(values: Partial<RangeRow> = {}): RangeRow {
  return {
    id: 'range-1',
    alwaysVisible: false,
    chain: 'A',
    start: '10',
    end: '20',
    ...values,
  };
}

describe('stacking boundary normalization', () => {
  it('parses JSON and returns normalized plain rows in deterministic order', () => {
    const later = stacking({ protein_range_pdb_id: 3, protein_resid: 7 });
    const earlier = stacking({ protein_resid: 5 });
    expect(normalizeStackings(JSON.stringify([later, earlier]))).toEqual([earlier, later]);
    expect(normalizeStackings([])).toEqual([]);
  });

  it.each([
    null,
    {},
    'not JSON',
    [null],
    [{ ...stacking(), protein_range_pdb_id: 1 }],
    [{ ...stacking(), distance: Number.NaN }],
    [{ ...stacking(), tilt_angle: 91 }],
    [{ ...stacking(), rna_centroid: [1, 2] }],
    [{ ...stacking(), rna_resname: ' ' }],
    [stacking(), stacking()],
  ])('rejects malformed table values (%#)', value => {
    expect(normalizeStackings(value)).toBeNull();
  });
});

describe('stacking identity and selection merge', () => {
  it('uses participant identity and ignores measurement changes', () => {
    const original = stacking();
    const measuredAgain = stacking({ distance: 5.5, tilt_angle: 40, slip_angle: 71 });
    expect(stackingKey(original)).toBe(stackingKey(measuredAgain));
    expect(stackingKey(original)).not.toBe(stackingKey(stacking({ protein_resid: 94 })));
  });

  it('preserves selection through backend reordering and measurement-only updates', () => {
    const first = stacking({ rna_resid: 4 });
    const second = stacking({ rna_resid: 8, distance: 5.4, tilt_angle: 41 });
    const next = mergeStackings(
      [first, second],
      [false, true],
      [
        stacking({ rna_resid: 8, distance: 4, tilt_angle: 20 }),
        stacking({ rna_resid: 4, distance: 5.8, tilt_angle: 42 }),
      ],
    );
    expect(next.stackings.map(row => row.rna_resid)).toEqual([8, 4]);
    expect(next.stackingSelection).toEqual([true, false]);
    expect(next.stackings[0].distance).toBe(4);
  });

  it('auto-selects only genuinely new rows at distance <= 5 Å and slip < 60° (strict)', () => {
    const previous = stacking({ rna_resid: 1 });
    const result = mergeStackings([previous], [false], [
      previous,
      stacking({ rna_resid: 2, distance: 5, slip_angle: 59.99 }),
      stacking({ rna_resid: 3, distance: 5, slip_angle: 60 }),
      stacking({ rna_resid: 4, distance: 5.01, slip_angle: 10 }),
      stacking({ rna_resid: 5, distance: 4, tilt_angle: 45, slip_angle: 10 }),
    ]);
    expect(result.stackingSelection).toEqual([false, true, false, false, true]);
  });

  it('does not preserve misaligned selection arrays', () => {
    const result = mergeStackings([stacking()], [], [stacking()]);
    expect(result.stackingSelection).toEqual([true]);
  });
});

describe('stage transitions and range-edit invalidation', () => {
  beforeEach(() => {
    appState.stage = 1;
    appState.rows = [];
    appState.staleStackingRevision = null;
    appState.stackings = [stacking()];
    appState.stackingSelection = [true];
  });

  it('enters stage 3 and returns to ranges without losing range or stacking state', () => {
    const rows = [rangeRow()];
    enterStage2(rows);
    enterStage3();
    expect(appState.stage).toBe(3);
    returnToStage2();
    expect(appState.stage).toBe(2);
    expect(appState.rows).toEqual(rows);
    expect(appState.stackings).toHaveLength(1);
    expect(appState.stackingSelection).toEqual([true]);
  });

  it('clears stacking rows and selection when returning to stage 1', () => {
    returnToStage1();
    expect(appState.stackings).toEqual([]);
    expect(appState.stackingSelection).toEqual([]);
  });

  it('treats the first Stage 2 snapshot as initialization and later content changes as edits', () => {
    let observation: Stage2RangeObservation = { active: false, signature: null };
    const initial = rangeRow();
    const entered = observeStage2RangeRows(observation, 2, [initial]);
    observation = entered.observation;
    expect(entered.changed).toBe(false);

    const idOnlyChange = observeStage2RangeRows(observation, 2, [rangeRow({ id: 'different-id' })]);
    observation = idOnlyChange.observation;
    expect(idOnlyChange.changed).toBe(false);

    const localEdit = observeStage2RangeRows(observation, 2, [rangeRow({ end: '25' })]);
    observation = localEdit.observation;
    expect(localEdit.changed).toBe(true);

    const sharedEdit = observeStage2RangeRows(observation, 2, [rangeRow({ end: '30' })]);
    observation = sharedEdit.observation;
    expect(sharedEdit.changed).toBe(true);

    const leftStage2 = observeStage2RangeRows(observation, 3, [rangeRow({ end: '30' })]);
    const returned = observeStage2RangeRows(leftStage2.observation, 2, [rangeRow({ end: '30' })]);
    expect(returned.changed).toBe(false);
  });

  it('clears the local table and records the stale backend revision on invalidation', () => {
    invalidateStackingState(9);
    expect(appState.stackings).toEqual([]);
    expect(appState.stackingSelection).toEqual([]);
    expect(appState.staleStackingRevision).toBe(9);
    clearStackingState();
    expect(appState.staleStackingRevision).toBeNull();
  });
});

describe('stacking viewer participant mapping', () => {
  const chains: ChainInfo[] = [
    { name: 'A', labelAsymId: 'A', instanceId: '1', seqIds: new Set([1, 4, 8]), seqMin: 1, seqMax: 8 },
    { name: 'B', labelAsymId: 'B', instanceId: '1', seqIds: new Set([20, 30, 40]), seqMin: 20, seqMax: 40 },
    { name: 'C', labelAsymId: 'C', instanceId: '2', seqIds: new Set([80, 83, 90]), seqMin: 80, seqMax: 90 },
  ];

  it('maps one-based range PDB IDs after empty editor rows are omitted', () => {
    const rows = [
      rangeRow({ chain: 'A', start: '1', end: '8' }),
      rangeRow({ id: 'empty', chain: '', start: '', end: '' }),
      rangeRow({ id: 'protein-b', chain: 'B', start: '20', end: '40' }),
      rangeRow({ id: 'protein-c', chain: 'C', start: '80', end: '90' }),
    ];
    const result = resolveSelectedStackings(
      [stacking({ protein_range_pdb_id: 3, protein_resid: 83 })],
      [true],
      rows,
      chains,
    );
    expect(result).toEqual([{
      rna: { labelAsymId: 'A', instanceId: '1', start: 4, end: 4 },
      protein: { labelAsymId: 'C', instanceId: '2', start: 83, end: 83 },
      rnaCentroid: [1, 2, 3],
      proteinCentroid: [4, 5, 6],
    }]);
  });

  it('sends only selected rows and skips participants that do not belong to their range', () => {
    const rows = [
      rangeRow({ chain: 'A', start: '1', end: '8' }),
      rangeRow({ id: 'protein-b', chain: 'B', start: '20', end: '40' }),
    ];
    const result = resolveSelectedStackings(
      [stacking(), stacking({ rna_resid: 4, protein_resid: 41 })],
      [false, true],
      rows,
      chains,
    );
    expect(result).toEqual([]);
    expect(resolveSelectedStackings([stacking()], [], rows, chains)).toEqual([]);
  });

  it('de-duplicates overlapping participant residue expressions in first-seen order', () => {
    const rna = { labelAsymId: 'A', instanceId: '1', start: 4, end: 4 };
    const protein = { labelAsymId: 'B', instanceId: '1', start: 30, end: 30 };
    const anotherProtein = { labelAsymId: 'C', instanceId: '2', start: 83, end: 83 };
    const inputs: StackingViewerInput[] = [
      { rna, protein, rnaCentroid: [0, 0, 0], proteinCentroid: [1, 0, 0] },
      { rna: { ...rna }, protein: { ...protein }, rnaCentroid: [0, 0, 0], proteinCentroid: [1, 0, 0] },
      { rna: { ...rna }, protein: anotherProtein, rnaCentroid: [0, 0, 0], proteinCentroid: [1, 0, 0] },
    ];
    expect(uniqueStackingParticipantRanges(inputs)).toEqual([rna, protein, anotherProtein]);
  });
});

describe('stacking connector dash geometry', () => {
  it('generates separated, evenly spaced cylinder segments along the centroid vector', () => {
    expect(stackingDashSegments([0, 0, 0], [8, 0, 0], 4, 0.5)).toEqual([
      { start: [0, 0, 0], end: [1, 0, 0] },
      { start: [2, 0, 0], end: [3, 0, 0] },
      { start: [4, 0, 0], end: [5, 0, 0] },
      { start: [6, 0, 0], end: [7, 0, 0] },
    ]);
  });

  it('omits degenerate and invalid connector segments', () => {
    expect(stackingDashSegments([1, 2, 3], [1, 2, 3])).toEqual([]);
    expect(stackingDashSegments([0, 0, 0], [1, 0, 0], 0)).toEqual([]);
    expect(stackingDashSegments([0, 0, 0], [1, 0, 0], 4, 0)).toEqual([]);
  });
});
