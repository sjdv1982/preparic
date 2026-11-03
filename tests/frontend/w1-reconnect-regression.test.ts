import { afterEach, describe, expect, it, vi } from 'vitest';
import {
  connectToSeamless, seamlessState,
  type SharedCell, type SharedValue, type ShareContext,
} from '../../src/lib/seamless.svelte';

const rna = { chain: 'A', start: 101, end: 104 };
const rows = [rna];

class Cell implements SharedCell {
  readonly = false;
  marker = 1;
  writes: string[] = [];
  onchange: (value: SharedValue) => void = () => undefined;
  constructor(public value: SharedValue, public url: string) {}
  set(value: string) { this.writes.push(value); return true; }
  emit(value: SharedValue) { this.value = value; this.onchange(value); }
}

function setup() {
  const doc = (value: unknown) => JSON.stringify(value);
  const values: Record<string, SharedValue> = {
    pdbcode: '1B7F', biounit: '1',
    ranges: doc({ schema: 1, basis: { pdbcode: '1B7F', biounit: 1 }, active: true, rows }),
    stacking_selection: doc({ schema: 1, basis: { pdbcode: '1B7F', biounit: 1, rows }, active: true, selected: {} }),
    threshold_edits: null,
    best_fit_clip: doc({ schema: 1, basis: { pdbcode: '1B7F', rna }, clip1: 0, clip2: 0 }),
    interaction_matrix: null, energy_threshold_edits: null, gaussian_classes: null,
    structure_id: '1B7F',
    fragments: doc({ schema: 1, echo: { pdbcode: '1B7F', rna }, residues: [], fragments: [], adjacent_pairs: [] }),
  };
  const cells = Object.fromEntries(Object.entries(values).map(([key, value]) => [key, new Cell(value, `/preparic/${key}`)]));
  const socket = { readyState: 1 };
  const context = {
    ...cells,
    self: { sharelist: [], onsharelist: (_keys: string[]) => undefined,
      server: 'http://w1.invalid', ws: socket as unknown as WebSocket },
  } as unknown as ShareContext;
  connectToSeamless(() => context);
  context.self.onsharelist(Object.keys(cells));
  for (const cell of Object.values(cells)) cell.emit(cell.value);
  return { context, cells, socket };
}

afterEach(() => { vi.clearAllTimers(); vi.useRealTimers(); vi.unstubAllGlobals(); });

describe('W1 P5/P9 reconnect requires a current-connection read of each result', () => {
  it.each(['omitted', 'failed'] as const)('keeps retained matching payload stale when its reread is %s', async mode => {
    vi.useFakeTimers();
    const { context, cells, socket } = setup();
    expect(seamlessState.ready).toBe(true);
    expect(seamlessState.computedFreshness.fragments).toBe(true);
    const oldPayload = JSON.parse(JSON.stringify(seamlessState.computed.fragments));
    await vi.advanceTimersByTimeAsync(250); // observe the live socket before disconnect
    socket.readyState = 3;
    await vi.advanceTimersByTimeAsync(250);
    expect(seamlessState.status).toBe('disconnected');
    expect(seamlessState.ready).toBe(false);
    expect(seamlessState.computed.fragments).toEqual(oldPayload); // retained behind overlay
    expect(seamlessState.computedFreshness.fragments).not.toBe(true);

    const requested: string[] = [];
    vi.stubGlobal('fetch', vi.fn(async (url: URL) => {
      const key = url.pathname.split('/').pop()!;
      requested.push(key);
      if (key === 'fragments' && mode === 'failed') return new Response('unavailable', { status: 503 });
      const value = cells[key].value;
      return new Response(value === null ? null : String(value), {
        status: value === null ? 204 : 200, headers: { 'x-seamless-marker': '1' },
      });
    }));
    socket.readyState = 1;
    context.self.onsharelist(Object.keys(cells).filter(key => mode !== 'omitted' || key !== 'fragments'));
    await vi.advanceTimersByTimeAsync(0);
    // Decisions may unlock editing after the barrier, but cannot certify this result.
    expect(seamlessState.ready).toBe(true);
    expect(seamlessState.projectStage).toBe(3);
    expect(seamlessState.computed.fragments).toEqual(oldPayload);
    expect(seamlessState.computedFreshness.fragments).not.toBe(true);
    expect(requested.includes('fragments')).toBe(mode === 'failed');
    expect(Object.values(cells).flatMap(cell => cell.writes)).toEqual([]);

    // An actual new delivery on this connection restores freshness without editing decisions.
    if (mode === 'failed') {
      cells.fragments.emit(cells.fragments.value);
      expect(seamlessState.computedFreshness.fragments).toBe(true);
    }
  });
});
