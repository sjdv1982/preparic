import { rangeDefinitionObject, rangeDefinitionSignature, parseRangeDefinition, serializeRangeDefinition } from './range-definition';
import type { RangeRow } from './types';

type SharedValue = string | Blob | null | undefined;
type ShareKey = 'pdbcode' | 'biounit' | 'ranges';

interface SharedCell {
  value: SharedValue;
  readonly: boolean;
  onchange: (value: SharedValue) => void;
  set: (value: string) => boolean;
}

interface StateGraphNode {
  path: Array<string | number>;
  exception?: string | null;
}

interface StateGraph {
  nodes?: StateGraphNode[];
}

interface ShareContext {
  self: {
    sharelist: string[];
    onsharelist: (keys: string[]) => void;
    onstategraph?: (graph: StateGraph) => void;
  };
  pdbcode?: SharedCell;
  biounit?: SharedCell;
  ranges?: SharedCell;
  range_specs?: SharedCell;
  'cif.cif'?: SharedCell;
}

type ShareAddress = string | number | null;
type ShareConnector = (
  updateServer?: ShareAddress,
  restServer?: ShareAddress,
  namespace?: string,
) => ShareContext;

type IdleDeadlineLike = { didTimeout: boolean; timeRemaining: () => number };
type IdleWindow = Window & {
  requestIdleCallback?: (callback: (deadline: IdleDeadlineLike) => void, options?: { timeout: number }) => number;
  cancelIdleCallback?: (handle: number) => void;
};

declare global {
  interface Window {
    connect_seamless?: ShareConnector;
  }
}

export const seamlessState = $state({
  status: 'connecting' as 'connecting' | 'connected' | 'unavailable' | 'waiting-for-shares',
  ready: false,
  pdbCode: '',
  biounitIndex: 1,
  pdbCodePendingWrite: false,
  biounitIndexPendingWrite: false,
  rangesReady: false,
  rangeRows: [] as RangeRow[],
  rangeRevision: 0,
  rangeError: '',
  rangeTransformError: '',
});

let client: ShareContext | undefined;
let pendingPdbCode: string | undefined;
let pendingBiounitIndex: number | undefined;
let pendingCif: string | undefined;
let cifWriteCallbackId: number | undefined;
let cifWriteCallbackKind: 'idle' | 'timeout' | undefined;
let cifWriteGeneration = 0;
let rangeSignature: string | undefined;
let rangeSpecsSignature: string | undefined;
const received = new Set<ShareKey>();

function normalizePdbCode(value: SharedValue): string {
  return typeof value === 'string' ? value.trim().toUpperCase() : '';
}

function parseBiounitIndex(value: SharedValue): number {
  const raw = typeof value === 'string' ? value.trim() : '';
  const index = raw === '' ? Number.NaN : Number(raw);
  return Number.isSafeInteger(index) && index >= 1 ? index : 1;
}

function updateReady(): void {
  if (!received.has('pdbcode') || !received.has('biounit')) return;
  seamlessState.ready = true;
  seamlessState.status = 'connected';
}

function receiveStateGraph(graph: StateGraph): void {
  const failedNode = graph.nodes?.find(node =>
    node.path.includes('ranges_to_pdb')
    && typeof node.exception === 'string'
    && node.exception.length > 0
  );
  seamlessState.rangeTransformError = failedNode?.exception ?? '';
}

function receivePdbCode(value: SharedValue): void {
  const code = normalizePdbCode(value);
  received.add('pdbcode');
  if (pendingPdbCode !== undefined && pendingPdbCode !== code) {
    seamlessState.pdbCode = pendingPdbCode;
  } else {
    pendingPdbCode = undefined;
    seamlessState.pdbCodePendingWrite = false;
    seamlessState.pdbCode = code;
  }
  updateReady();
}

function receiveBiounitIndex(value: SharedValue): void {
  const index = parseBiounitIndex(value);
  received.add('biounit');
  if (pendingBiounitIndex !== undefined && pendingBiounitIndex !== index) {
    seamlessState.biounitIndex = pendingBiounitIndex;
  } else {
    pendingBiounitIndex = undefined;
    seamlessState.biounitIndexPendingWrite = false;
    seamlessState.biounitIndex = index;
  }
  updateReady();
}

function receiveRanges(value: SharedValue): void {
  if (value !== undefined && value !== null && typeof value !== 'string') {
    seamlessState.rangesReady = false;
    seamlessState.rangeError = 'The shared range definition is not YAML text.';
    return;
  }
  try {
    const rows = parseRangeDefinition(typeof value === 'string' ? value : '');
    const signature = rangeDefinitionSignature(rows);
    seamlessState.rangeError = '';
    seamlessState.rangesReady = true;
    if (signature === rangeSignature) return;
    rangeSignature = signature;
    seamlessState.rangeRows = rows;
    seamlessState.rangeRevision += 1;
  } catch (error) {
    seamlessState.rangesReady = false;
    seamlessState.rangeError = error instanceof Error ? error.message : 'Could not read the YAML range definition.';
  }
}

function attachShare(key: ShareKey, keys: string[]): void {
  if (!keys.includes(key)) return;
  const cell = client?.[key];
  if (!cell) return;
  cell.onchange = key === 'pdbcode'
    ? receivePdbCode
    : key === 'biounit'
      ? receiveBiounitIndex
      : receiveRanges;
  if (cell.value !== undefined) cell.onchange(cell.value);
}

function flushCifWrite(generation: number): void {
  cifWriteCallbackId = undefined;
  cifWriteCallbackKind = undefined;
  if (generation !== cifWriteGeneration) return;
  const cell = client?.['cif.cif'];
  const cif = pendingCif;
  if (cif === undefined || !cell || cell.readonly) return;
  let accepted = false;
  try {
    accepted = cell.set(cif);
  } catch {
    // Keep the latest CIF queued; export errors must not affect the viewer.
  }
  if (accepted) {
    if (pendingCif === cif) pendingCif = undefined;
    return;
  }
  // A share entry can be rejected while the browser client is replacing a stale entry.
  // Retry until the current writable entry is ready; accepted writes handle transport
  // retries in the browser client.
  if (pendingCif !== cif || typeof window === 'undefined') return;
  cifWriteCallbackKind = 'timeout';
  cifWriteCallbackId = window.setTimeout(() => flushCifWrite(generation), 250);
}

function scheduleCifWrite(): void {
  if (pendingCif === undefined || cifWriteCallbackId !== undefined) return;
  const cell = client?.['cif.cif'];
  if (!cell || cell.readonly || typeof window === 'undefined') return;
  // Prefer idle time, but cap the delay so an active viewer cannot starve the upload.
  const generation = cifWriteGeneration;
  const idleWindow = window as IdleWindow;
  if (idleWindow.requestIdleCallback) {
    cifWriteCallbackKind = 'idle';
    cifWriteCallbackId = idleWindow.requestIdleCallback(
      () => flushCifWrite(generation),
      { timeout: 1000 }
    );
  } else {
    cifWriteCallbackKind = 'timeout';
    cifWriteCallbackId = window.setTimeout(() => flushCifWrite(generation), 100);
  }
}

function attachCifShare(keys: string[]): void {
  if (!keys.includes('cif.cif')) return;
  scheduleCifWrite();
}

if (typeof window === 'undefined' || !window.connect_seamless) {
  seamlessState.status = 'unavailable';
} else {
  try {
    client = window.connect_seamless(null, null, 'preparic');
    client.self.onstategraph = receiveStateGraph;
    client.self.onsharelist = keys => {
      attachShare('pdbcode', keys);
      attachShare('biounit', keys);
      attachShare('ranges', keys);
      attachCifShare(keys);
      if (!keys.includes('pdbcode') || !keys.includes('biounit')) {
        seamlessState.status = 'waiting-for-shares';
      }
    };
  } catch {
    seamlessState.status = 'unavailable';
  }
}

export function setSharedCif(cif: string): void {
  pendingCif = cif;
  scheduleCifWrite();
}

export function cancelSharedCifWrite(): void {
  cifWriteGeneration += 1;
  const idleWindow = typeof window === 'undefined' ? undefined : window as IdleWindow;
  if (cifWriteCallbackId !== undefined) {
    if (cifWriteCallbackKind === 'idle' && idleWindow?.cancelIdleCallback) {
      idleWindow.cancelIdleCallback(cifWriteCallbackId);
    } else if (cifWriteCallbackKind === 'timeout') {
      window.clearTimeout(cifWriteCallbackId);
    }
  }
  cifWriteCallbackId = undefined;
  cifWriteCallbackKind = undefined;
  pendingCif = undefined;
}

export function setSharedPdbCode(value: string): boolean {
  const code = value.trim().toUpperCase();
  const cell = client?.pdbcode;
  if (!client || !seamlessState.ready || !cell) {
    pendingPdbCode = code;
    seamlessState.pdbCodePendingWrite = true;
    seamlessState.pdbCode = code;
    return true;
  }
  if (cell.readonly) return false;
  if (normalizePdbCode(cell.value) === code) {
    pendingPdbCode = undefined;
    seamlessState.pdbCodePendingWrite = false;
    seamlessState.pdbCode = code;
    return true;
  }
  pendingPdbCode = code;
  seamlessState.pdbCodePendingWrite = true;
  const accepted = cell.set(code);
  if (accepted) seamlessState.pdbCode = code;
  return accepted;
}

export function setSharedBiounitIndex(value: number): boolean {
  const index = Number.isSafeInteger(value) && value >= 1 ? value : 1;
  const cell = client?.biounit;
  if (!client || !seamlessState.ready || !cell) {
    pendingBiounitIndex = index;
    seamlessState.biounitIndexPendingWrite = true;
    seamlessState.biounitIndex = index;
    return true;
  }
  if (cell.readonly) return false;
  if (parseBiounitIndex(cell.value) === index) {
    pendingBiounitIndex = undefined;
    seamlessState.biounitIndexPendingWrite = false;
    seamlessState.biounitIndex = index;
    return true;
  }
  pendingBiounitIndex = index;
  seamlessState.biounitIndexPendingWrite = true;
  const accepted = cell.set(String(index));
  if (accepted) seamlessState.biounitIndex = index;
  return accepted;
}

export function setSharedRangeRows(rows: RangeRow[]): boolean {
  const cell = client?.ranges;
  if (!client || !seamlessState.rangesReady || !cell || cell.readonly) return false;
  const signature = rangeDefinitionSignature(rows);
  if (signature === rangeSignature) return true;
  try {
    const accepted = cell.set(serializeRangeDefinition(rows));
    if (!accepted) return false;
    rangeSignature = signature;
    seamlessState.rangeRows = rows.map(row => ({ ...row }));
    return true;
  } catch {
    return false;
  }
}

export function setSharedRangeSpecs(rows: RangeRow[]): boolean {
  const cell = client?.range_specs;
  if (!client || !seamlessState.rangesReady || !cell || cell.readonly) return false;
  const value = JSON.stringify(rangeDefinitionObject(rows));
  if (value === rangeSpecsSignature) return true;
  try {
    const accepted = cell.set(value);
    if (!accepted) return false;
    rangeSpecsSignature = value;
    return true;
  } catch {
    return false;
  }
}
