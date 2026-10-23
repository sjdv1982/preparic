import { rangeDefinitionObject, rangeDefinitionSignature, parseRangeDefinition, serializeRangeDefinition } from './range-definition';
import { mergeStackings, normalizeStackings } from './stackings';
import { appState } from './state.svelte';
import type { RangeRow, Stacking } from './types';

export type SharedValue = string | Blob | null | undefined;
type ShareKey = 'pdbcode' | 'biounit' | 'ranges' | 'stackings' | 'stacking_selection';

export interface SharedCell {
  value: SharedValue;
  readonly: boolean;
  onchange: (value: SharedValue) => void;
  set: (value: string) => boolean;
  marker?: number;
}

interface StateGraphNode {
  path: Array<string | number>;
  exception?: string | null;
}

interface StateGraph {
  nodes?: StateGraphNode[];
}

export interface ShareContext {
  self: {
    sharelist: string[];
    onsharelist: (keys: string[]) => void;
    onstategraph?: (graph: StateGraph) => void;
  };
  pdbcode?: SharedCell;
  biounit?: SharedCell;
  ranges?: SharedCell;
  range_specs?: SharedCell;
  stackings?: SharedCell;
  stacking_selection?: SharedCell;
  'cif.cif'?: SharedCell;
}

type ShareAddress = string | number | null;
export type ShareConnector = (
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
  stackingsReady: false,
  stackingsError: '',
  stackingRevision: 0,
  stackingSelectionReady: false,
  stackingSelectionRevision: 0,
  stackingSelectionError: '',
  stackingSelectionPendingWrite: false,
  stackingTransformError: '',
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
let stackingsSignature: string | undefined;
let backendStackingSelectionSignature: string | undefined;
let pendingStackingSelectionWriteSignature: string | undefined;
let pendingStackingSelectionWriteMarker: number | undefined;
let pendingBackendStackingSelection: boolean[] | undefined;
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
  const failedRangeNode = graph.nodes?.find(node =>
    node.path.includes('ranges_to_pdb')
    && typeof node.exception === 'string'
    && node.exception.length > 0
  );
  const failedStackingNode = graph.nodes?.find(node =>
    node.path.some(part => String(part).includes('detect_stacking'))
    && typeof node.exception === 'string'
    && node.exception.length > 0
  );
  seamlessState.rangeTransformError = failedRangeNode?.exception ?? '';
  seamlessState.stackingTransformError = failedStackingNode?.exception ?? '';
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

function updateBackendRevision(current: number, marker?: number): number {
  // Seamless markers are per shared cell, which keeps table and selection
  // revisions independent and lets the range-edit barrier reject old values.
  if (Number.isSafeInteger(marker) && (marker as number) >= 0) return marker as number;
  return current + 1;
}

function sharedValueSignature(value: SharedValue): string {
  if (typeof value === 'string') return `text:${value}`;
  if (value === null) return 'null';
  if (value === undefined) return 'undefined';
  try {
    return `other:${JSON.stringify(value)}`;
  } catch {
    return `other:${String(value)}`;
  }
}

function parseBooleanSelection(value: SharedValue): boolean[] | null {
  if (typeof value !== 'string') return null;
  let parsed: unknown;
  try {
    parsed = JSON.parse(value) as unknown;
  } catch {
    return null;
  }
  if (!Array.isArray(parsed) || !parsed.every(item => typeof item === 'boolean')) return null;
  return [...parsed];
}

function writeStackingSelection(selection: readonly boolean[]): boolean {
  if (!Array.isArray(selection) || !selection.every(item => typeof item === 'boolean')) {
    seamlessState.stackingSelectionError = 'The stacking selection must contain only booleans.';
    return false;
  }
  if (selection.length !== appState.stackings.length) {
    seamlessState.stackingSelectionError = 'The stacking selection does not match the current table length.';
    return false;
  }

  const cell = client?.stacking_selection;
  if (!cell) {
    seamlessState.stackingSelectionPendingWrite = false;
    seamlessState.stackingSelectionError = 'The shared stacking selection is not available.';
    return false;
  }
  if (cell.readonly) {
    seamlessState.stackingSelectionPendingWrite = false;
    seamlessState.stackingSelectionError = 'The shared stacking selection is read-only.';
    return false;
  }

  const serialized = JSON.stringify([...selection]);
  if (pendingStackingSelectionWriteSignature === serialized) return true;
  const currentSelection = parseBooleanSelection(cell.value);
  if (pendingStackingSelectionWriteSignature === undefined
    && currentSelection !== null
    && JSON.stringify(currentSelection) === serialized) {
    backendStackingSelectionSignature = serialized;
    seamlessState.stackingSelectionError = '';
    return true;
  }
  if (pendingStackingSelectionWriteSignature === undefined
    && backendStackingSelectionSignature === serialized) {
    seamlessState.stackingSelectionError = '';
    return true;
  }

  pendingStackingSelectionWriteSignature = serialized;
  pendingStackingSelectionWriteMarker = Number.isSafeInteger(cell.marker) ? cell.marker : undefined;
  seamlessState.stackingSelectionPendingWrite = true;
  let accepted = false;
  try {
    accepted = cell.set(serialized);
  } catch {
    accepted = false;
  }
  if (!accepted) {
    if (pendingStackingSelectionWriteSignature === serialized) {
      pendingStackingSelectionWriteSignature = undefined;
      pendingStackingSelectionWriteMarker = undefined;
      seamlessState.stackingSelectionPendingWrite = false;
    }
    seamlessState.stackingSelectionError = 'The shared stacking selection rejected the update.';
    return false;
  }
  seamlessState.stackingSelectionError = '';
  return true;
}

function receiveStackings(value: SharedValue, marker?: number): void {
  seamlessState.stackingRevision = updateBackendRevision(seamlessState.stackingRevision, marker);
  const stackings = normalizeStackings(value);
  if (!stackings) {
    seamlessState.stackingsReady = false;
    seamlessState.stackingsError = 'The shared stacking table is not a valid JSON list of stacking rows.';
    return;
  }

  seamlessState.stackingsReady = true;
  seamlessState.stackingsError = '';
  if (appState.staleStackingRevision !== null
    && seamlessState.stackingRevision <= appState.staleStackingRevision) return;
  if (appState.staleStackingRevision !== null) appState.staleStackingRevision = null;

  const signature = JSON.stringify(stackings);
  const sameAsLastBackendTable = signature === stackingsSignature;
  const sameAsLocalTable = JSON.stringify(appState.stackings) === signature;
  const hasQueuedSelection = pendingBackendStackingSelection !== undefined;
  if (sameAsLastBackendTable && sameAsLocalTable
    && appState.stackingSelection.length === appState.stackings.length
    && !hasQueuedSelection) return;
  stackingsSignature = signature;

  const merged = mergeStackings(appState.stackings, appState.stackingSelection, stackings);
  appState.stackings = merged.stackings;
  appState.stackingSelection = merged.stackingSelection;

  let shouldPushSelection = true;
  if (pendingBackendStackingSelection !== undefined) {
    if (stackings.length === 0) {
      // A queued backend selection is associated with a later non-empty table.
      shouldPushSelection = false;
    } else {
      const queued = pendingBackendStackingSelection;
      pendingBackendStackingSelection = undefined;
      if (queued.length === stackings.length) {
        appState.stackingSelection = [...queued];
        backendStackingSelectionSignature = JSON.stringify(queued);
        seamlessState.stackingSelectionReady = true;
        seamlessState.stackingSelectionError = '';
        shouldPushSelection = false;
      } else {
        seamlessState.stackingSelectionError = `The queued backend selection has ${queued.length} values, but the next non-empty table has ${stackings.length} rows.`;
      }
    }
  }

  if (shouldPushSelection) writeStackingSelection(appState.stackingSelection);
}

function receiveStackingSelection(value: SharedValue, marker?: number): void {
  seamlessState.stackingSelectionRevision = updateBackendRevision(
    seamlessState.stackingSelectionRevision,
    marker,
  );
  const selection = parseBooleanSelection(value);
  const signature = selection === null ? sharedValueSignature(value) : JSON.stringify(selection);
  seamlessState.stackingSelectionReady = selection !== null;
  let newerThanPendingWrite = false;

  if (pendingStackingSelectionWriteSignature !== undefined) {
    if (signature === pendingStackingSelectionWriteSignature) {
      backendStackingSelectionSignature = signature;
      pendingStackingSelectionWriteSignature = undefined;
      pendingStackingSelectionWriteMarker = undefined;
      seamlessState.stackingSelectionPendingWrite = false;
      seamlessState.stackingSelectionError = '';
      return;
    }
    if (Number.isSafeInteger(marker)
      && pendingStackingSelectionWriteMarker !== undefined
      && (marker as number) <= pendingStackingSelectionWriteMarker) return;
    // A distinct backend value arrived while our write was pending; it is newer
    // than the last value we observed and therefore wins. Equal/older markers
    // above are stale reads that raced with the optimistic local write.
    newerThanPendingWrite = true;
    pendingStackingSelectionWriteSignature = undefined;
    pendingStackingSelectionWriteMarker = undefined;
    seamlessState.stackingSelectionPendingWrite = false;
  }
  if (!newerThanPendingWrite && signature === backendStackingSelectionSignature) return;
  backendStackingSelectionSignature = signature;

  if (selection === null) {
    pendingBackendStackingSelection = undefined;
    seamlessState.stackingSelectionError = 'The shared stacking selection must be a JSON list containing only booleans.';
    return;
  }
  if (selection.length === appState.stackings.length) {
    pendingBackendStackingSelection = undefined;
    appState.stackingSelection = selection;
    seamlessState.stackingSelectionError = '';
    return;
  }

  pendingBackendStackingSelection = selection;
  seamlessState.stackingSelectionError = '';
}

function attachShare(key: ShareKey, keys: string[]): void {
  if (!keys.includes(key)) return;
  const cell = client?.[key];
  if (!cell) return;
  switch (key) {
    case 'pdbcode':
      cell.onchange = receivePdbCode;
      break;
    case 'biounit':
      cell.onchange = receiveBiounitIndex;
      break;
    case 'ranges':
      cell.onchange = receiveRanges;
      break;
    case 'stackings':
      cell.onchange = value => receiveStackings(value, cell.marker);
      break;
    case 'stacking_selection':
      cell.onchange = value => receiveStackingSelection(value, cell.marker);
      break;
  }
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

function attachShareContext(context: ShareContext): void {
  client = context;
  stackingsSignature = undefined;
  backendStackingSelectionSignature = undefined;
  pendingStackingSelectionWriteSignature = undefined;
  pendingStackingSelectionWriteMarker = undefined;
  pendingBackendStackingSelection = undefined;
  seamlessState.stackingSelectionPendingWrite = false;
  context.self.onstategraph = receiveStateGraph;
  context.self.onsharelist = keys => {
    attachShare('pdbcode', keys);
    attachShare('biounit', keys);
    attachShare('ranges', keys);
    // Selection must be observed before stackings so a persisted exact-length
    // value can supersede defaults for the table being attached.
    attachShare('stacking_selection', keys);
    attachShare('stackings', keys);
    attachCifShare(keys);
    if (!keys.includes('pdbcode') || !keys.includes('biounit')) {
      seamlessState.status = 'waiting-for-shares';
    }
  };
}

export function connectToSeamless(connect: ShareConnector): void {
  attachShareContext(connect(null, null, 'preparic'));
}

if (typeof window === 'undefined' || !window.connect_seamless) {
  seamlessState.status = 'unavailable';
} else {
  try {
    connectToSeamless(window.connect_seamless);
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

export function setSharedStackingSelection(selection: readonly boolean[]): boolean {
  // An explicit checkbox edit supersedes a delayed backend list that was
  // waiting for a table of a different length.
  pendingBackendStackingSelection = undefined;
  return writeStackingSelection(selection);
}
