import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import {
  connectToSeamless,
  resetStage4ActiveNucleotide,
  seamlessState,
  setStage4ActiveNucleotide,
  type SharedCell,
  type SharedValue,
  type ShareConnector,
  type ShareContext,
} from '../../src/lib/seamless.svelte';
import { appState } from '../../src/lib/state.svelte';
import { makeThresholdSkeleton, type FragmentIdentity } from '../../src/lib/stage4';

const rows = [{ chain: 'A', start: 101, end: 104 }];
const rangeBasis = { pdbcode: '1B7F', biounit: 1 };
const rna = { chain: 'A', start: 101, end: 104 };
const projectBasis = { pdbcode: '1B7F', rna };

class FakeCell implements SharedCell {
  value: SharedValue;
  readonly: boolean;
  marker: number;
  onchange: (value: SharedValue) => void = () => undefined;
  writes: string[] = [];

  constructor(value: SharedValue, marker = 1, readonly = false) {
    this.value = value;
    this.marker = marker;
    this.readonly = readonly;
  }

  set(value: string): boolean {
    if (this.readonly) return false;
    this.writes.push(value);
    return true;
  }

  emit(value: SharedValue, marker: number): void {
    this.value = value;
    this.marker = marker;
    this.onchange(value);
  }
}

function makeMarker(value: unknown): string {
  return JSON.stringify(value);
}

function setupStage4() {
  const inventory = {
    schema: 1,
    echo: { pdbcode: '1B7F', rna },
    residues: [],
    fragments: [{ number: 1, key: 'frag-1-AC', sequence: 'AC', residues: [] }],
    adjacent_pairs: [],
  };
  const thresholds = makeThresholdSkeleton(inventory);
  if (!thresholds) throw new Error('failed to build test threshold marker');
  const selection = {
    schema: 1,
    basis: { ...rangeBasis, rows },
    active: true,
    selected: { '["C",102,2,30]': true },
  };
  const initialClip = { schema: 1, basis: projectBasis, clip1: 0, clip2: 0 };
  const socket = { readyState: 1 };
  const context = {
    self: { sharelist: [] as string[], onsharelist: (_keys: string[]) => undefined, ws: socket as unknown as WebSocket },
    pdbcode: new FakeCell('1B7F'),
    biounit: new FakeCell('1'),
    ranges: new FakeCell(makeMarker({ schema: 1, basis: rangeBasis, active: true, rows })),
    stacking_selection: new FakeCell(makeMarker(selection)),
    threshold_edits: new FakeCell(makeMarker(thresholds)),
    best_fit_clip: new FakeCell(makeMarker(initialClip)),
    interaction_matrix: new FakeCell(null),
    energy_threshold_edits: new FakeCell(null),
    gaussian_classes: new FakeCell(null),
    structure_id: new FakeCell('1B7F', 1, true),
    stage4_active_nucleotide: new FakeCell(null),
  };
  const connector: ShareConnector = (_updateServer, _restServer, namespace) => {
    if (namespace !== 'preparic') throw new Error(`Unexpected namespace: ${namespace}`);
    return context as unknown as ShareContext;
  };

  connectToSeamless(connector);
  const keys = [
    'pdbcode', 'biounit', 'ranges', 'stacking_selection', 'threshold_edits', 'best_fit_clip',
    'interaction_matrix', 'energy_threshold_edits', 'gaussian_classes', 'structure_id',
    'stage4_active_nucleotide',
  ];
  context.self.onsharelist(keys);
  context.pdbcode!.emit('1B7F', 1);
  context.biounit!.emit('1', 1);
  context.ranges!.emit(makeMarker({ schema: 1, basis: rangeBasis, active: true, rows }), 1);
  context.stacking_selection!.emit(makeMarker(selection), 1);
  context.threshold_edits!.emit(makeMarker(thresholds), 1);
  context.best_fit_clip!.emit(makeMarker(initialClip), 1);
  context.interaction_matrix!.emit(null, 1);
  context.energy_threshold_edits!.emit(null, 1);
  context.gaussian_classes!.emit(null, 1);
  context.structure_id!.emit('1B7F', 1);
  context.stage4_active_nucleotide!.emit(null, 1);
  return {
    context, socket, ranges: context.ranges as FakeCell,
    active: context.stage4_active_nucleotide as FakeCell,
  };
}

const n102: FragmentIdentity = { chain: 'A', resid: 102, insertion_code: '', base: 'C' };
const n103: FragmentIdentity = { chain: 'A', resid: 103, insertion_code: '', base: 'G' };

describe('S26 stage4_active_nucleotide (shared view state)', () => {
  beforeEach(() => vi.useFakeTimers());
  afterEach(() => {
    vi.clearAllTimers();
    vi.useRealTimers();
    appState.stage = 1;
  });

  it('shows a local choice at once and writes it after the debounce, without a pending input', async () => {
    const { active } = setupStage4();
    expect(seamlessState.projectStage).toBe(4);
    expect(setStage4ActiveNucleotide(n102)).toBe(true);
    expect(seamlessState.stage4ActiveNucleotide).toEqual(n102);
    expect(seamlessState.pendingInputs).not.toContain('stage4_active_nucleotide');
    await vi.advanceTimersByTimeAsync(299);
    expect(active.writes).toEqual([]);
    await vi.advanceTimersByTimeAsync(1);
    expect(active.writes.map(value => JSON.parse(value))).toEqual([n102]);
    active.emit(JSON.stringify(n102), 2);
    expect(seamlessState.stage4ActiveNucleotide).toEqual(n102);
  });

  it('coalesces slider movement into one write of the last value', async () => {
    const { active } = setupStage4();
    setStage4ActiveNucleotide(n102);
    await vi.advanceTimersByTimeAsync(200);
    setStage4ActiveNucleotide(n103);
    await vi.advanceTimersByTimeAsync(300);
    expect(active.writes.map(value => JSON.parse(value))).toEqual([n103]);
  });

  it('adopts a remote choice and drops the queued local one', async () => {
    const { active } = setupStage4();
    setStage4ActiveNucleotide(n102);
    active.emit(JSON.stringify(n103), 2);
    expect(seamlessState.stage4ActiveNucleotide).toEqual(n103);
    await vi.advanceTimersByTimeAsync(600);
    expect(active.writes).toEqual([]);
    expect(seamlessState.syncNotice).toBe('');
  });

  it('keeps a newer local choice when the echo of an earlier write arrives', async () => {
    const { active } = setupStage4();
    setStage4ActiveNucleotide(n102);
    await vi.advanceTimersByTimeAsync(300);
    setStage4ActiveNucleotide(n103);
    active.emit(JSON.stringify(n102), 2);
    expect(seamlessState.stage4ActiveNucleotide).toEqual(n103);
    await vi.advanceTimersByTimeAsync(300);
    expect(active.writes.map(value => JSON.parse(value))).toEqual([n102, n103]);
  });

  it('drops a queued choice when the stage changes', async () => {
    const { active, ranges } = setupStage4();
    setStage4ActiveNucleotide(n102);
    ranges.emit(makeMarker({
      schema: 1, basis: rangeBasis, active: true, rows: [{ chain: 'A', start: 102, end: 104 }],
    }), 2);
    expect(seamlessState.projectStage).not.toBe(4);
    await vi.advanceTimersByTimeAsync(600);
    expect(active.writes).toEqual([]);
    expect(setStage4ActiveNucleotide(n103)).toBe(false);
  });

  it('treats a malformed value as null (the default is shown)', () => {
    const { active } = setupStage4();
    active.emit(JSON.stringify({ resid: 'x' }), 2);
    expect(seamlessState.stage4ActiveNucleotide).toBeNull();
    expect(active.writes).toEqual([]);
  });

  it('reset writes null at once, cancels a queued choice and skips a redundant write', async () => {
    const { active } = setupStage4();
    resetStage4ActiveNucleotide();
    expect(active.writes).toEqual([]); // already null on the server

    active.emit(JSON.stringify(n102), 2);
    setStage4ActiveNucleotide(n103);
    resetStage4ActiveNucleotide();
    expect(seamlessState.stage4ActiveNucleotide).toBeNull();
    expect(active.writes).toEqual(['null']);
    await vi.advanceTimersByTimeAsync(600);
    expect(active.writes).toEqual(['null']);
    active.emit(null, 3);
    expect(seamlessState.stage4ActiveNucleotide).toBeNull();
  });
});
