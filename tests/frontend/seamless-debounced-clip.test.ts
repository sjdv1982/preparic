import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import {
  connectToSeamless,
  flushDebouncedProjectDecision,
  seamlessState,
  setProjectDecisionDebounced,
  type SharedCell,
  type SharedValue,
  type ShareConnector,
  type ShareContext,
} from '../../src/lib/seamless.svelte';
import { appState } from '../../src/lib/state.svelte';
import { makeThresholdSkeleton } from '../../src/lib/stage4';

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
  };
  const connector: ShareConnector = (_updateServer, _restServer, namespace) => {
    if (namespace !== 'preparic') throw new Error(`Unexpected namespace: ${namespace}`);
    return context as unknown as ShareContext;
  };

  connectToSeamless(connector);
  const keys = [
    'pdbcode', 'biounit', 'ranges', 'stacking_selection', 'threshold_edits', 'best_fit_clip',
    'interaction_matrix', 'energy_threshold_edits', 'gaussian_classes', 'structure_id',
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
  return { context, socket, clip: context.best_fit_clip as FakeCell, ranges: context.ranges as FakeCell };
}

describe('Stage 4 debounced clip writes', () => {
  beforeEach(() => vi.useFakeTimers());
  afterEach(() => {
    vi.clearAllTimers();
    vi.useRealTimers();
    appState.stage = 1;
  });

  it('waits for the debounce, then writes and waits for acknowledgment', async () => {
    const { clip } = setupStage4();
    expect(seamlessState.projectStage).toBe(4);
    const next = { schema: 1, basis: projectBasis, clip1: 1, clip2: 0 };

    const saved = setProjectDecisionDebounced('best_fit_clip', next, 500);
    expect(clip.writes).toEqual([]);
    expect(seamlessState.pendingInputs).toContain('best_fit_clip');
    await vi.advanceTimersByTimeAsync(499);
    expect(clip.writes).toEqual([]);
    await vi.advanceTimersByTimeAsync(1);
    expect(clip.writes).toHaveLength(1);
    expect(JSON.parse(clip.writes[0])).toEqual(next);

    clip.emit(makeMarker(next), 2);
    await expect(saved).resolves.toBe(true);
    expect(seamlessState.pendingInputs).not.toContain('best_fit_clip');
  });

  it.each(['stage', 'connection', 'revision'])('drops the queued edit after a %s change', async change => {
    const { context, socket, clip, ranges } = setupStage4();
    const next = { schema: 1, basis: projectBasis, clip1: 1, clip2: 0 };
    const saved = setProjectDecisionDebounced('best_fit_clip', next, 500);

    if (change === 'stage') {
      ranges.emit(makeMarker({
        schema: 1, basis: rangeBasis, active: true, rows: [{ chain: 'A', start: 102, end: 104 }],
      }), 2);
    } else if (change === 'connection') {
      await vi.advanceTimersByTimeAsync(250); // observe the live websocket
      socket.readyState = 3;
      await vi.advanceTimersByTimeAsync(250);
      expect(seamlessState.status).toBe('disconnected');
    } else {
      clip.emit(makeMarker({ schema: 1, basis: projectBasis, clip1: 0, clip2: 0 }), 2);
    }

    await expect(saved).resolves.toBe(false);
    await vi.advanceTimersByTimeAsync(600);
    expect(clip.writes).toEqual([]);
    expect(seamlessState.pendingInputs).not.toContain('best_fit_clip');
    expect(context.best_fit_clip).toBe(clip);
  });

  it('flushes a queued clip before a dependent transition', async () => {
    const { clip } = setupStage4();
    const next = { schema: 1, basis: projectBasis, clip1: 0, clip2: 1 };
    const queued = setProjectDecisionDebounced('best_fit_clip', next, 500);

    const flush = flushDebouncedProjectDecision('best_fit_clip');
    expect(clip.writes).toHaveLength(1);
    clip.emit(makeMarker(next), 2);

    await expect(queued).resolves.toBe(true);
    await expect(flush).resolves.toBe(true);
  });
});
