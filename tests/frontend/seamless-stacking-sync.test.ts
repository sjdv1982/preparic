import { beforeEach, describe, expect, it } from 'vitest';
import {
  connectToSeamless,
  seamlessState,
  setSharedStackingSelection,
  type SharedCell,
  type SharedValue,
  type ShareConnector,
  type ShareContext,
} from '../../src/lib/seamless.svelte';
import { appState, invalidateStackingState } from '../../src/lib/state.svelte';
import type { Stacking } from '../../src/lib/types';

class FakeCell implements SharedCell {
  value: SharedValue;
  readonly: boolean;
  marker = 0;
  onchange: (value: SharedValue) => void = () => undefined;
  writes: string[] = [];
  rejectWrites = false;

  constructor(value: SharedValue = undefined, readonly = false) {
    this.value = value;
    this.readonly = readonly;
  }

  set(value: string): boolean {
    if (this.readonly || this.rejectWrites) return false;
    this.value = value;
    this.writes.push(value);
    return true;
  }

  emit(value: SharedValue, marker: number): void {
    this.value = value;
    this.marker = marker;
    this.onchange(value);
  }
}

interface FakeConnector {
  context: ShareContext;
  stackings: FakeCell;
  selection: FakeCell;
}

function createFakeConnector(): FakeConnector {
  const stackings = new FakeCell(undefined, true);
  const selection = new FakeCell();
  const context: ShareContext = {
    self: {
      sharelist: [],
      onsharelist: () => undefined,
    },
    pdbcode: new FakeCell('1B7F'),
    biounit: new FakeCell('1'),
    ranges: new FakeCell(''),
    stackings,
    stacking_selection: selection,
  };
  const fakeConnector: ShareConnector = (_updateServer, _restServer, namespace) => {
    if (namespace !== 'preparic') throw new Error(`Unexpected namespace: ${namespace}`);
    return context;
  };
  connectToSeamless(fakeConnector);
  context.self.onsharelist(['pdbcode', 'biounit', 'ranges', 'stackings', 'stacking_selection']);
  return { context, stackings, selection };
}

function stacking(values: Partial<Stacking> = {}): Stacking {
  return {
    rna_resname: 'G',
    rna_resid: 4,
    protein_range_pdb_id: 2,
    protein_resid: 93,
    distance: 4.2,
    tilt_angle: 22,
    slip_angle: 71,
    rna_centroid: [1, 2, 3],
    protein_centroid: [4, 5, 6],
    ...values,
  };
}

function emitTable(connector: FakeConnector, value: unknown, marker: number): void {
  connector.stackings.emit(JSON.stringify(value), marker);
}

describe('Seamless stacking share synchronization', () => {
  beforeEach(() => {
    appState.stackings = [];
    appState.stackingSelection = [];
    appState.staleStackingRevision = null;
    seamlessState.stackingRevision = 0;
    seamlessState.stackingSelectionRevision = 0;
    seamlessState.stackingsReady = false;
    seamlessState.stackingsError = '';
    seamlessState.stackingSelectionReady = false;
    seamlessState.stackingSelectionError = '';
    seamlessState.stackingSelectionPendingWrite = false;
    seamlessState.stackingTransformError = '';
  });

  it('attaches the read-only table and writable selection shares with separate revisions', () => {
    const connector = createFakeConnector();
    expect(connector.stackings.readonly).toBe(true);
    expect(connector.selection.readonly).toBe(false);

    emitTable(connector, [stacking()], 11);
    connector.selection.emit('[false]', 5);

    expect(seamlessState.stackingRevision).toBe(11);
    expect(seamlessState.stackingSelectionRevision).toBe(5);
    expect(appState.stackingSelection).toEqual([false]);
  });

  it('auto-selects new rows and pushes the complete list to the backend', () => {
    const connector = createFakeConnector();
    const selectedByDefault = stacking();
    const notSelected = stacking({ rna_resid: 8, distance: 5.2 });

    emitTable(connector, [selectedByDefault, notSelected], 1);

    expect(appState.stackingSelection).toEqual([true, false]);
    expect(connector.selection.writes).toEqual(['[true,false]']);
  });

  it('applies an exact-length remote boolean list without writing it back', () => {
    const connector = createFakeConnector();
    emitTable(connector, [stacking(), stacking({ rna_resid: 8 })], 1);
    const writesBeforeRemoteUpdate = [...connector.selection.writes];

    connector.selection.emit('[false,true]', 2);

    expect(appState.stackingSelection).toEqual([false, true]);
    expect(connector.selection.writes).toEqual(writesBeforeRemoteUpdate);
    expect(seamlessState.stackingSelectionError).toBe('');
  });

  it('queues a length mismatch through an empty table and applies it to the next non-empty table', () => {
    const connector = createFakeConnector();
    emitTable(connector, [stacking()], 1);
    const writesBeforeQueue = [...connector.selection.writes];

    connector.selection.emit('[false,true]', 2);
    expect(appState.stackingSelection).toEqual([true]);
    expect(seamlessState.stackingSelectionError).toBe('');

    emitTable(connector, [], 2);
    expect(appState.stackings).toEqual([]);
    expect(connector.selection.writes).toEqual(writesBeforeQueue);

    emitTable(connector, [stacking(), stacking({ rna_resid: 8 })], 3);
    expect(appState.stackingSelection).toEqual([false, true]);
    expect(connector.selection.writes).toEqual(writesBeforeQueue);
  });

  it('rejects invalid selection types without replacing the current selection', () => {
    const connector = createFakeConnector();
    emitTable(connector, [stacking()], 1);
    const currentSelection = [...appState.stackingSelection];

    connector.selection.emit('[true,"false"]', 2);

    expect(appState.stackingSelection).toEqual(currentSelection);
    expect(seamlessState.stackingSelectionReady).toBe(false);
    expect(seamlessState.stackingSelectionError).toMatch(/only booleans/);
  });

  it('pushes checkbox changes and suppresses the matching backend echo', () => {
    const connector = createFakeConnector();
    emitTable(connector, [stacking()], 1);
    appState.stackingSelection = [false];

    expect(setSharedStackingSelection(appState.stackingSelection)).toBe(true);
    expect(connector.selection.writes).toEqual(['[true]', '[false]']);
    expect(seamlessState.stackingSelectionPendingWrite).toBe(true);

    connector.selection.emit('[false]', 2);

    expect(connector.selection.writes).toEqual(['[true]', '[false]']);
    expect(seamlessState.stackingSelectionPendingWrite).toBe(false);
    expect(appState.stackingSelection).toEqual([false]);
  });

  it('ignores an older backend read that races with an optimistic selection write', () => {
    const connector = createFakeConnector();
    emitTable(connector, [stacking()], 1);
    expect(seamlessState.stackingSelectionPendingWrite).toBe(true);

    connector.selection.emit('[]', 0);
    expect(appState.stackingSelection).toEqual([true]);
    expect(seamlessState.stackingSelectionPendingWrite).toBe(true);

    connector.selection.emit('[true]', 1);
    expect(seamlessState.stackingSelectionPendingWrite).toBe(false);
    expect(connector.selection.writes).toEqual(['[true]']);
  });

  it('surfaces rejected writes without discarding the optimistic checkbox state', () => {
    const connector = createFakeConnector();
    emitTable(connector, [stacking()], 1);
    connector.selection.rejectWrites = true;
    appState.stackingSelection = [false];

    expect(setSharedStackingSelection(appState.stackingSelection)).toBe(false);
    expect(appState.stackingSelection).toEqual([false]);
    expect(seamlessState.stackingSelectionPendingWrite).toBe(false);
    expect(seamlessState.stackingSelectionError).toMatch(/rejected/);
  });

  it('ignores the stale table marker after a range edit, then accepts the next table revision', () => {
    const connector = createFakeConnector();
    const oldRows = [stacking()];
    emitTable(connector, oldRows, 10);

    invalidateStackingState(seamlessState.stackingRevision);
    emitTable(connector, oldRows, 10);
    expect(appState.stackings).toEqual([]);
    expect(appState.staleStackingRevision).toBe(10);

    emitTable(connector, oldRows, 11);
    expect(appState.stackings).toEqual(oldRows);
    expect(appState.staleStackingRevision).toBeNull();
  });

  it('reports malformed table JSON and detector failures separately', () => {
    const connector = createFakeConnector();
    emitTable(connector, [stacking()], 1);
    const validRows = [...appState.stackings];

    connector.stackings.emit('{not json', 2);
    connector.context.self.onstategraph?.({
      nodes: [
        { path: ['ranges_to_pdb'], exception: 'range error' },
        { path: ['detect_stacking'], exception: 'detector error' },
      ],
    });

    expect(appState.stackings).toEqual(validRows);
    expect(seamlessState.stackingsError).toMatch(/valid JSON/);
    expect(seamlessState.rangeTransformError).toBe('range error');
    expect(seamlessState.stackingTransformError).toBe('detector error');
  });
});
