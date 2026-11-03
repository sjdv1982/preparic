import { parse, stringify } from 'yaml';
import { mergeStackings, normalizeStackings, stackingKey } from './stackings';
import { appState } from './state.svelte';
import {
  deferredWriteIsCurrent,
  deriveProjectStage,
  echoIsFresh,
  markerIsWellFormed,
  markerMatchesUpstream,
  visibleTransformerErrors,
  type ProjectSnapshot,
  type ProjectStage,
} from './project-state';
import type { RangeRow, Stacking } from './types';

export type SharedValue = string | Blob | null | undefined;
type DecisionKey = 'ranges' | 'stacking_selection' | 'threshold_edits' | 'best_fit_clip' | 'interaction_matrix' | 'gaussian_classes';
type ShareKey = 'pdbcode' | 'biounit' | 'ranges' | 'range_specs' | 'stackings' | 'stacking_selection'
  | 'threshold_edits' | 'best_fit_clip' | 'interaction_matrix' | 'gaussian_classes' | 'energy_threshold_edits'
  | 'structure_id' | 'range_pdbs_aa' | 'fragments' | 'best_fit_thresholds' | 'best_fits' | 'energetics'
  | 'reference_fragment_pdbs'
  | 'energy_thresholds' | 'interactions' | 'bestfit_stackings' | 'gaussian_validation' | 'refe_validation'
  | 'gaussian_scores';

export interface SharedCell {
  value: SharedValue;
  readonly: boolean;
  onchange: (value: SharedValue) => void;
  set: (value: string) => boolean;
  marker?: number;
  url?: string;
  binary?: boolean;
}

interface StateGraphNode {
  path: Array<string | number>;
  state?: string;
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
    ws?: WebSocket | null;
    server?: string;
  };
  pdbcode?: SharedCell;
  biounit?: SharedCell;
  ranges?: SharedCell;
  range_specs?: SharedCell;
  stackings?: SharedCell;
  stacking_selection?: SharedCell;
  threshold_edits?: SharedCell;
  best_fit_clip?: SharedCell;
  interaction_matrix?: SharedCell;
  gaussian_classes?: SharedCell;
  energy_threshold_edits?: SharedCell;
  structure_id?: SharedCell;
  range_pdbs_aa?: SharedCell;
  fragments?: SharedCell;
  reference_fragment_pdbs?: SharedCell;
  best_fit_thresholds?: SharedCell;
  best_fits?: SharedCell;
  energetics?: SharedCell;
  energy_thresholds?: SharedCell;
  interactions?: SharedCell;
  bestfit_stackings?: SharedCell;
  gaussian_validation?: SharedCell;
  refe_validation?: SharedCell;
  gaussian_scores?: SharedCell;
  'cif.cif'?: SharedCell;
}

type ShareAddress = string | number | null;
interface StructureCommit {
  code: string;
  cif: string;
  codeAcknowledged: boolean;
  cifAcknowledged: boolean;
  cifConflict: boolean;
  retryingCif: boolean;
}
export interface SharedStructureTarget {
  stageAtIssue: ProjectStage;
  liveAtIssue: boolean;
  connectionRevision: number;
  biounitIndex: number;
  pdbCodeRevision: number;
  cifRevision: number;
  pdbCodeMarker: number | undefined;
  cifMarker: number | undefined;
}
interface PendingShareWrite {
  expected: unknown;
  serialized: string;
  baseMarker: number | undefined;
  stageAtIssue: ProjectStage;
  revision: number;
  promise: Promise<boolean>;
  resolve: (accepted: boolean) => void;
}
export type ShareConnector = (
  updateServer?: ShareAddress,
  restServer?: ShareAddress,
  namespace?: string,
) => ShareContext;

declare global {
  interface Window {
    connect_seamless?: ShareConnector;
  }
}

export const seamlessState = $state({
  status: 'connecting' as 'connecting' | 'live' | 'disconnected' | 'unavailable' | 'waiting-for-shares',
  ready: false,
  projectStage: 1 as ProjectStage,
  projectSnapshot: {
    pdbcode: '',
    biounit: 1,
    ranges: null,
    stacking_selection: null,
    threshold_edits: null,
    best_fit_clip: { schema: 1, basis: null, clip1: 0, clip2: 0 },
    interaction_matrix: null,
    gaussian_classes: null,
  } as ProjectSnapshot,
  decisionErrors: {} as Partial<Record<DecisionKey, string>>,
  stateGraphErrors: [] as Array<{ name: string; state: string; message: string }>,
  hasInitialStateGraph: false,
  pendingInputs: [] as string[],
  syncNotice: '',
  structureId: '',
  structureIdReady: false,
  rangeSpecs: null as unknown,
  rangePdbsAa: null as unknown,
  computed: {} as Partial<Record<ShareKey, unknown>>,
  computedFreshness: {} as Partial<Record<ShareKey, boolean>>,
  pdbCode: '',
  biounitIndex: 1,
  pdbCodePendingWrite: false,
  biounitIndexPendingWrite: false,
  rangesReady: false,
  rangeRows: [] as RangeRow[],
  rangeRevision: 0,
  stackingsReady: false,
  hasStackingsSnapshot: false,
  stackingsError: '',
  stackingRevision: 0,
  stackingSelectionReady: false,
  stackingSelectionNonEmpty: false,
  stackingSelectionRevision: 0,
  stackingSelectionError: '',
  stackingSelectionPendingWrite: false,
  stackingTransformError: '',
  rangeError: '',
  rangeTransformError: '',
});

let client: ShareContext | undefined;
let pendingBiounitIndex: number | undefined;
let pendingCif: string | undefined;
let pendingPdbCode: string | undefined;
let pendingStructureCommit: StructureCommit | undefined;
const pendingShareWrites = new Map<string, PendingShareWrite>();
const latestWriteRevisions = new Map<string, number>();
const shareReceivers = new Map<ShareKey, (value: SharedValue) => void>();
let nextWriteRevision = 0;
let resyncInProgress = false;
let applyingAuthoritativeRead = false;
let rangeSignature: string | undefined;
let stackingsSignature: string | undefined;
let backendStackingSelectionSignature: string | undefined;
let pendingStackingSelectionWriteSignature: string | undefined;
let pendingStackingSelectionWriteMarker: number | undefined;
let pendingBackendStackingSelection: Record<string, boolean> | undefined;
const received = new Set<ShareKey>();
const requiredEditShares: readonly ShareKey[] = [
  'pdbcode', 'biounit', 'ranges', 'stacking_selection', 'threshold_edits',
  'best_fit_clip', 'interaction_matrix', 'energy_threshold_edits', 'gaussian_classes',
];
const requiredComputedShares: readonly ShareKey[] = ['structure_id'];
let websocketPoll: ReturnType<typeof setInterval> | undefined;
let hasConnectedBefore = false;
let websocketWasOpen = false;
let structureConnectionRevision = 0;
let structureShareRevisions = { pdbcode: 0, cif: 0 };

function normalizePdbCode(value: SharedValue): string {
  return typeof value === 'string' ? value.trim().toUpperCase() : '';
}

function parseBiounitIndex(value: SharedValue): number {
  const raw = typeof value === 'string' ? value.trim() : '';
  const index = raw === '' ? Number.NaN : Number(raw);
  return Number.isSafeInteger(index) && index >= 1 ? index : 1;
}

function parseJsonValue(value: SharedValue): unknown | undefined {
  if (value === null || value === undefined) return value;
  if (typeof value !== 'string') return typeof Blob !== 'undefined' && value instanceof Blob ? undefined : value;
  try {
    return JSON.parse(value) as unknown;
  } catch {
    return value;
  }
}

function parseMarkerValue(value: SharedValue, kind: DecisionKey): unknown | undefined {
  if (value === null || value === undefined) return null;
  if (typeof value !== 'string') return undefined;
  try {
    if (kind === 'ranges') return parse(value) as unknown;
    return JSON.parse(value) as unknown;
  } catch {
    return undefined;
  }
}

function decisionCell(key: DecisionKey): SharedCell | undefined {
  return client?.[key];
}

function writeDecision(key: DecisionKey, marker: unknown): boolean {
  const cell = decisionCell(key);
  if (!cell || cell.readonly || !markerIsWellFormed(key, marker)) return false;
  let serialized: string;
  try {
    serialized = key === 'ranges' ? stringify(marker) : JSON.stringify(marker);
  } catch {
    return false;
  }
  const existing = pendingShareWrites.get(key);
  if (existing && stableJson(existing.expected) === stableJson(marker)) return true;
  if (!pendingShareWrites.has(key)
    && stableJson(canonicalShareValue(key, cell.value)) === stableJson(marker)) return true;
  const pending = beginPendingShareWrite(key, serialized, marker);
  let accepted = false;
  try { accepted = cell.set(serialized); } catch { accepted = false; }
  if (!accepted) {
    resolvePendingShareWrite(key, false);
    return false;
  }
  pending.serialized = serialized;
  return true;
}

function currentSnapshot(): ProjectSnapshot {
  return seamlessState.projectSnapshot;
}

function refreshProjectStage(): void {
  seamlessState.projectStage = deriveProjectStage(currentSnapshot());
  appState.stage = seamlessState.projectStage;
  refreshVisibleErrors();
}

function recordDecision(key: DecisionKey, value: unknown): void {
  const snapshot = seamlessState.projectSnapshot as unknown as Record<DecisionKey, unknown>;
  snapshot[key] = value;
  delete seamlessState.decisionErrors[key];
  refreshProjectStage();
  refreshComputedFreshness();
}

function markPendingInput(key: string): void {
  if (!seamlessState.pendingInputs.includes(key)) seamlessState.pendingInputs = [...seamlessState.pendingInputs, key];
  refreshVisibleErrors();
}

function clearPendingInput(key: string): void {
  seamlessState.pendingInputs = seamlessState.pendingInputs.filter(item => item !== key);
  refreshVisibleErrors();
}

function canonicalShareValue(key: string, value: SharedValue): unknown {
  if (key === 'pdbcode') return normalizePdbCode(value);
  if (key === 'biounit') return parseBiounitIndex(value);
  if (key === 'cif.cif') return value;
  if (key === 'ranges' || key === 'stacking_selection' || key === 'threshold_edits'
    || key === 'best_fit_clip' || key === 'interaction_matrix' || key === 'gaussian_classes') {
    return parseMarkerValue(value, key);
  }
  return parseJsonValue(value);
}

function serializeFreshShareValue(key: string, value: SharedValue): string {
  if (key === 'cif.cif') return typeof value === 'string' ? value : '';
  if (key === 'pdbcode') return normalizePdbCode(value);
  if (key === 'biounit') return String(parseBiounitIndex(value));
  if (key === 'ranges') return typeof value === 'string' ? value : stringify(value ?? null);
  return typeof value === 'string' ? value : JSON.stringify(value ?? null);
}

function resolvePendingShareWrite(key: string, accepted: boolean, conflict = false): void {
  const pending = pendingShareWrites.get(key);
  if (!pending) return;
  pendingShareWrites.delete(key);
  clearPendingInput(key);
  pending.resolve(accepted);
  if (conflict) {
    seamlessState.syncNotice = `A newer shared ${key.replaceAll('_', ' ')} value replaced this edit. The accepted project value is now shown.`;
  }
}

function beginPendingShareWrite(key: string, serialized: string, expected: unknown): PendingShareWrite {
  const existing = pendingShareWrites.get(key);
  if (existing && stableJson(existing.expected) === stableJson(expected)) return existing;
  if (existing) resolvePendingShareWrite(key, false, true);
  const cell = key === 'cif.cif' ? client?.['cif.cif'] : client?.[key as ShareKey];
  let resolve!: (accepted: boolean) => void;
  const promise = new Promise<boolean>(done => { resolve = done; });
  const pending: PendingShareWrite = {
    expected,
    serialized,
    baseMarker: Number.isSafeInteger(cell?.marker) ? cell?.marker : undefined,
    stageAtIssue: seamlessState.projectStage,
    revision: ++nextWriteRevision,
    promise,
    resolve,
  };
  pendingShareWrites.set(key, pending);
  latestWriteRevisions.set(key, pending.revision);
  markPendingInput(key);
  return pending;
}

function observePendingWrite(key: string, value: SharedValue, marker?: number, authoritative = false): boolean {
  const pending = pendingShareWrites.get(key);
  if (!pending) return false;
  const receivedValue = canonicalShareValue(key, value);
  if (stableJson(receivedValue) === stableJson(pending.expected)) {
    resolvePendingShareWrite(key, true);
    return false;
  }
  if (authoritative && pending.baseMarker !== undefined && Number.isSafeInteger(marker)
    && (marker as number) <= pending.baseMarker) return true;
  // A fresh authoritative mismatch must stay pending until resynchronizeShares
  // can queue the freshly read server value after the client's retained retry.
  if (authoritative) return false;
  if (!authoritative && pending.baseMarker !== undefined && Number.isSafeInteger(marker)
    && (marker as number) <= pending.baseMarker) {
    return true;
  }
  resolvePendingShareWrite(key, false, true);
  return false;
}

export function waitForShareWrite(key: string): Promise<boolean> {
  return pendingShareWrites.get(key)?.promise ?? Promise.resolve(true);
}

function refreshVisibleErrors(): void {
  seamlessState.stateGraphErrors = visibleTransformerErrors({
    stage: seamlessState.projectStage,
    hasInitialStateGraph: seamlessState.hasInitialStateGraph,
    pendingInputs: new Set(seamlessState.pendingInputs),
    nodes: lastStateGraphNodes,
  });
  seamlessState.rangeTransformError = seamlessState.stateGraphErrors
    .find(issue => issue.name === 'ranges_to_pdb' || issue.name === 'parse_ranges')?.message ?? '';
  seamlessState.stackingTransformError = seamlessState.stateGraphErrors
    .find(issue => issue.name === 'detect_stacking')?.message ?? '';
}

let lastStateGraphNodes: StateGraphNode[] = [];
function updateReady(): void {
  if (resyncInProgress) {
    seamlessState.ready = false;
    seamlessState.status = 'waiting-for-shares';
    return;
  }
  // An empty stage 1 must stay recoverable even when structure_id cannot be computed yet.
  const structureIdReadRequired = !(seamlessState.projectStage === 1 && !seamlessState.pdbCode)
    && !requiredComputedShares.every(key => received.has(key));
  if (!requiredEditShares.every(key => received.has(key)) || structureIdReadRequired) {
    seamlessState.ready = false;
    seamlessState.status = 'waiting-for-shares';
    return;
  }
  seamlessState.ready = true;
  seamlessState.status = 'live';
  refreshProjectStage();
  refreshComputedFreshness();
  if (pendingBiounitIndex !== undefined) setSharedBiounitIndex(pendingBiounitIndex);
}

function suspendForDisconnect(): void {
  if (!hasConnectedBefore || seamlessState.status === 'disconnected') return;
  if (!pendingStructureCommit) {
    pendingCif = undefined;
    pendingPdbCode = undefined;
    seamlessState.pdbCodePendingWrite = false;
    clearPendingInput('pdbcode');
  }
  structureConnectionRevision += 1;
  seamlessState.ready = false;
  seamlessState.status = 'disconnected';
  seamlessState.computedFreshness = {};
  seamlessState.structureIdReady = false;
  seamlessState.hasInitialStateGraph = false;
  seamlessState.stateGraphErrors = [];
  lastStateGraphNodes = [];
  seamlessState.stackingsReady = false;
  seamlessState.rangesReady = false;
  seamlessState.stackingSelectionReady = false;
  stackingsSignature = undefined;
  backendStackingSelectionSignature = undefined;
  pendingBackendStackingSelection = undefined;
  received.clear();
  refreshVisibleErrors();
}

function watchShareConnection(): void {
  const open = client?.self.ws?.readyState === 1;
  if (open) {
    websocketWasOpen = true;
    return;
  }
  if (websocketWasOpen) suspendForDisconnect();
}

function computedValue(key: ShareKey): unknown {
  return seamlessState.computed[key];
}

function valueEcho(key: ShareKey): unknown {
  const value = computedValue(key);
  return value && typeof value === 'object' && 'echo' in value
    ? (value as { echo: unknown }).echo
    : undefined;
}

function primaryRange(): unknown {
  const marker = seamlessState.projectSnapshot.ranges;
  return marker && typeof marker === 'object' && 'rows' in marker && Array.isArray(marker.rows)
    ? marker.rows[0] ?? null
    : null;
}

function currentRangeSpecs(): unknown {
  const snapshot = seamlessState.projectSnapshot;
  const marker = snapshot.ranges;
  if (!markerMatchesUpstream('ranges', marker, snapshot)) return null;
  return {
    schema: 1,
    pdbcode: seamlessState.pdbCode,
    biounit: seamlessState.biounitIndex,
    rows: marker && typeof marker === 'object' && 'rows' in marker ? marker.rows : [],
  };
}

function selectionKeys(): string[] {
  const marker = seamlessState.projectSnapshot.stacking_selection;
  if (!marker || typeof marker !== 'object' || !('selected' in marker) || !marker.selected || typeof marker.selected !== 'object') return [];
  return Object.entries(marker.selected as Record<string, unknown>)
    .filter(([, selected]) => selected === true).map(([key]) => key).sort();
}

function currentComputedEcho(key: ShareKey): unknown | null {
  const snapshot = seamlessState.projectSnapshot;
  const firstRange = primaryRange();
  const rangeSpecs = currentRangeSpecs();
  const fragments = computedValue('fragments');
  const fragmentResidues = fragments && typeof fragments === 'object' && 'residues' in fragments
    ? (fragments as { residues: unknown }).residues
    : null;
  const bestFitsEcho = valueEcho('best_fits');
  const rangeRows = snapshot.ranges && typeof snapshot.ranges === 'object' && 'rows' in snapshot.ranges
    ? (snapshot.ranges as { rows: Array<{ chain: string; start: number; end: number }> }).rows
    : [];
  const proteinRanges = rangeRows.slice(1);
  switch (key) {
    case 'structure_id': return seamlessState.pdbCode || null;
    case 'stackings': {
      const structureId = computedValue('structure_id');
      if (typeof structureId !== 'string' || !rangeSpecs || typeof rangeSpecs !== 'object') return null;
      const rs = rangeSpecs as { rows?: unknown; biounit?: unknown };
      if (!Array.isArray(rs.rows)) return null;
      return { structure_id: structureId, pdbcode: seamlessState.pdbCode, biounit: rs.biounit, rows: rs.rows };
    }
    case 'fragments': return firstRange ? { pdbcode: seamlessState.pdbCode, rna: firstRange } : null;
    case 'reference_fragment_pdbs': {
      const fragmentMap = computedValue('fragments');
      return fragmentMap && typeof fragmentMap === 'object' && 'echo' in fragmentMap
        ? (fragmentMap as { echo: unknown }).echo : null;
    }
    case 'best_fit_thresholds':
    case 'best_fits':
      return fragmentResidues === null || !firstRange ? null
        : { pdb_code: seamlessState.pdbCode.toLowerCase(), rna: firstRange, resids: fragmentResidues };
    case 'energetics':
    case 'energy_thresholds':
    case 'bestfit_stackings':
      return bestFitsEcho ? { best_fits: bestFitsEcho, protein_ranges: proteinRanges } : null;
    case 'gaussian_scores': {
      const classesEcho = valueEcho('bestfit_stackings');
      const clip = snapshot.best_fit_clip;
      return classesEcho && clip ? { bestfit_stackings: classesEcho, selected_keys: selectionKeys(), clip } : null;
    }
    case 'gaussian_validation':
      return snapshot.gaussian_classes && snapshot.best_fit_clip
        ? { classes: (snapshot.gaussian_classes as { rows?: unknown }).rows, clip: snapshot.best_fit_clip,
            selected_keys: selectionKeys(), best_fits: bestFitsEcho }
        : null;
    case 'refe_validation':
      return snapshot.threshold_edits && snapshot.best_fit_clip
        ? { bases: (snapshot.threshold_edits as { bases?: unknown }).bases, clip: snapshot.best_fit_clip,
            selected_keys: selectionKeys(), best_fits: bestFitsEcho }
        : null;
    default: return null;
  }
}

function refreshComputedFreshness(): void {
  const keys: ShareKey[] = [
    'stackings', 'fragments', 'reference_fragment_pdbs', 'best_fit_thresholds', 'best_fits', 'energetics', 'energy_thresholds',
    'bestfit_stackings', 'gaussian_scores', 'gaussian_validation', 'refe_validation', 'range_pdbs_aa',
  ];
  const freshness: Partial<Record<ShareKey, boolean>> = {};
  freshness.structure_id = seamlessState.ready
    && computedValue('structure_id') === seamlessState.pdbCode && !!seamlessState.pdbCode;
  const currentRanges = currentRangeSpecs();
  freshness.range_specs = seamlessState.ready && stableJson(computedValue('range_specs')) === stableJson(currentRanges);
  for (const key of keys) {
    const expected = currentComputedEcho(key);
    freshness[key] = seamlessState.ready && echoIsFresh(valueEcho(key), expected);
  }
  const reduced = computedValue('range_pdbs_aa');
  freshness.range_pdbs_aa = !!reduced && typeof reduced === 'object'
    && echoIsFresh((reduced as { range_specs?: unknown }).range_specs, currentRanges)
    && (reduced as { pdb_code?: unknown }).pdb_code === seamlessState.pdbCode.toLowerCase();
  seamlessState.computedFreshness = freshness;
}

function receiveStateGraph(graph: StateGraph): void {
  lastStateGraphNodes = graph.nodes ?? [];
  seamlessState.hasInitialStateGraph = true;
  refreshVisibleErrors();
}

function finishStructureCommitIfAcknowledged(): void {
  const commit = pendingStructureCommit;
  if (!commit?.codeAcknowledged || !commit.cifAcknowledged) return;
  pendingStructureCommit = undefined;
  if (!pendingStructureCommit) clearPendingInput('pdbcode');
  clearPendingInput('cif.cif');
  if (pendingPdbCode === commit.code) {
    pendingPdbCode = undefined;
    seamlessState.pdbCodePendingWrite = false;
  }
}

function retryCifWhileCommitOwnsCode(): void {
  const commit = pendingStructureCommit;
  const cell = client?.['cif.cif'];
  if (!commit?.codeAcknowledged || commit.cifAcknowledged || !commit.cifConflict || commit.retryingCif
    || seamlessState.projectStage !== 1 || normalizePdbCode(appState.pdbId) !== commit.code
    || !cell || cell.readonly) return;
  commit.retryingCif = true;
  if (!pendingShareWrites.has('cif.cif')) beginPendingShareWrite('cif.cif', commit.cif, commit.cif);
  let accepted = false;
  try { accepted = cell.set(commit.cif); } catch { accepted = false; }
  commit.retryingCif = false;
  if (!accepted) {
    pendingStructureCommit = undefined;
    resolvePendingShareWrite('pdbcode', false);
    resolvePendingShareWrite('cif.cif', false);
  }
}

function receiveCifCommit(value: SharedValue): void {
  structureShareRevisions.cif += 1;
  if (observePendingWrite('cif.cif', value, client?.['cif.cif']?.marker, applyingAuthoritativeRead)) return;
  const commit = pendingStructureCommit;
  if (!commit) return;
  if (typeof value === 'string' && value === commit.cif) {
    commit.cifAcknowledged = true;
    finishStructureCommitIfAcknowledged();
  } else {
    commit.cifConflict = true;
    retryCifWhileCommitOwnsCode();
  }
}

function receivePdbCode(value: SharedValue): void {
  structureShareRevisions.pdbcode += 1;
  if (observePendingWrite('pdbcode', value, client?.pdbcode?.marker, applyingAuthoritativeRead)) {
    if (applyingAuthoritativeRead) { received.add('pdbcode'); updateReady(); }
    return;
  }
  const code = normalizePdbCode(value);
  received.add('pdbcode');
  const commit = pendingStructureCommit;
  if (commit && code === commit.code) {
    commit.codeAcknowledged = true;
    if (pendingPdbCode === commit.code) seamlessState.pdbCodePendingWrite = false;
    retryCifWhileCommitOwnsCode();
    finishStructureCommitIfAcknowledged();
  } else if (commit && code !== commit.code) {
    pendingStructureCommit = undefined;
    clearPendingInput('pdbcode');
    clearPendingInput('cif.cif');
    resolvePendingShareWrite('cif.cif', false, true);
    if (pendingPdbCode === commit.code) {
      pendingPdbCode = undefined;
      seamlessState.pdbCodePendingWrite = false;
    }
  } else if (pendingPdbCode === undefined || pendingPdbCode === code) {
    pendingPdbCode = undefined;
    seamlessState.pdbCodePendingWrite = false;
  }
  seamlessState.pdbCode = code;
  seamlessState.projectSnapshot.pdbcode = code;
  if (!pendingStructureCommit) clearPendingInput('pdbcode');
  updateReady();
}

function receiveBiounitIndex(value: SharedValue): void {
  if (observePendingWrite('biounit', value, client?.biounit?.marker, applyingAuthoritativeRead)) {
    if (applyingAuthoritativeRead) { received.add('biounit'); updateReady(); }
    return;
  }
  const index = parseBiounitIndex(value);
  received.add('biounit');
  pendingBiounitIndex = undefined;
  seamlessState.biounitIndexPendingWrite = false;
  seamlessState.biounitIndex = index;
  seamlessState.projectSnapshot.biounit = index;
  clearPendingInput('biounit');
  updateReady();
}

function receiveRanges(value: SharedValue): void {
  if (observePendingWrite('ranges', value, client?.ranges?.marker, applyingAuthoritativeRead)) {
    if (applyingAuthoritativeRead) { received.add('ranges'); seamlessState.rangesReady = true; updateReady(); }
    return;
  }
  const parsed = parseMarkerValue(value, 'ranges');
  const valid = parsed !== undefined && markerIsWellFormed('ranges', parsed);
  const marker = valid ? parsed : null;
  const signature = stableJson(marker);
  seamlessState.rangeError = valid ? '' : 'The stored ranges decision is invalid; review the project ranges.';
  seamlessState.rangesReady = true;
  received.add('ranges');
  if (valid) clearPendingInput('ranges');
  else seamlessState.decisionErrors.ranges = seamlessState.rangeError;
  recordDecision('ranges', marker);
  const rows = marker && typeof marker === 'object' && 'rows' in marker && Array.isArray(marker.rows)
    ? marker.rows.map((row, index) => ({
        id: `shared-${index + 1}`,
        alwaysVisible: false,
        chain: row.chain,
        start: String(row.start),
        end: String(row.end),
      }))
    : [];
  if (signature !== rangeSignature) {
    rangeSignature = signature;
    seamlessState.rangeRows = rows;
    seamlessState.rangeRevision += 1;
    if (seamlessState.ready && appState.stage === 2) appState.rows = rows.map(row => ({ ...row }));
  }
  updateReady();
}

function receiveDecision(key: Exclude<DecisionKey, 'ranges'>, value: SharedValue): void {
  if (observePendingWrite(key, value, client?.[key]?.marker, applyingAuthoritativeRead)) {
    if (applyingAuthoritativeRead) { received.add(key); updateReady(); }
    return;
  }
  const parsed = parseMarkerValue(value, key);
  const valid = parsed !== undefined && markerIsWellFormed(key, parsed);
  received.add(key);
  if (valid) clearPendingInput(key);
  else seamlessState.decisionErrors[key] = `The stored ${key.replaceAll('_', ' ')} decision is invalid; review it before continuing.`;
  recordDecision(key, valid ? parsed : null);
  refreshComputedFreshness();
  updateReady();
}

function receiveComputed(key: ShareKey, value: SharedValue): void {
  const parsed = parseJsonValue(value);
  seamlessState.computed[key] = parsed;
  received.add(key);
  if (key === 'range_specs') seamlessState.rangeSpecs = parsed;
  if (key === 'range_pdbs_aa') seamlessState.rangePdbsAa = parsed;
  if (key === 'structure_id') {
    seamlessState.structureId = typeof parsed === 'string' ? parsed : '';
    seamlessState.structureIdReady = true;
  }
  refreshComputedFreshness();
  updateReady();
}

function updateBackendRevision(current: number, marker?: number): number {
  // Seamless markers are per shared cell, which keeps table and selection
  // revisions independent and lets the range-edit barrier reject old values.
  if (Number.isSafeInteger(marker) && (marker as number) >= 0) return marker as number;
  return current + 1;
}

function stableJson(value: unknown): string {
  if (Array.isArray(value)) return `[${value.map(stableJson).join(',')}]`;
  if (value && typeof value === 'object') {
    return `{${Object.keys(value).sort().map(key => `${JSON.stringify(key)}:${stableJson((value as Record<string, unknown>)[key])}`).join(',')}}`;
  }
  return JSON.stringify(value) ?? 'undefined';
}

function parseStackingMarker(value: SharedValue): unknown | undefined {
  if (value === null || value === undefined) return null;
  return parseJsonValue(value);
}

function markerSelection(marker: unknown): Record<string, boolean> | null {
  if (!markerIsWellFormed('stacking_selection', marker) || !marker || typeof marker !== 'object'
    || !('selected' in marker) || !marker.selected || typeof marker.selected !== 'object') return null;
  return marker.selected as Record<string, boolean>;
}

function rangeRowsBasis(): Array<{ chain: string; start: number; end: number }> | null {
  const marker = seamlessState.projectSnapshot.ranges;
  if (!marker || typeof marker !== 'object' || !('rows' in marker) || !Array.isArray(marker.rows)) return null;
  return marker.rows as Array<{ chain: string; start: number; end: number }>;
}

function makeStackingMarker(selection: readonly boolean[], active: boolean): unknown | null {
  const rows = rangeRowsBasis();
  if (!rows || selection.length !== appState.stackings.length
    || !selection.every(item => typeof item === 'boolean')) return null;
  const selected: Record<string, boolean> = {};
  appState.stackings.forEach((stacking, index) => { selected[stackingKey(stacking)] = selection[index]; });
  return {
    schema: 1,
    basis: { pdbcode: seamlessState.pdbCode, biounit: seamlessState.biounitIndex, rows },
    active,
    selected,
  };
}

function writeStackingSelection(selection: readonly boolean[], active = true): boolean {
  const marker = makeStackingMarker(selection, active);
  const cell = client?.stacking_selection;
  if (!marker || !cell) {
    seamlessState.stackingSelectionPendingWrite = false;
    seamlessState.stackingSelectionError = 'The stacking selection does not match the current table or range decision.';
    return false;
  }
  if (cell.readonly) {
    seamlessState.stackingSelectionPendingWrite = false;
    seamlessState.stackingSelectionError = 'The shared stacking selection is read-only.';
    return false;
  }
  const serialized = JSON.stringify(marker);
  const signature = stableJson(marker);
  if (pendingStackingSelectionWriteSignature === signature) return true;
  const current = parseStackingMarker(cell.value);
  if (pendingStackingSelectionWriteSignature === undefined && stableJson(current) === signature) {
    recordDecision('stacking_selection', marker);
    backendStackingSelectionSignature = signature;
    seamlessState.stackingSelectionError = '';
    return true;
  }
  pendingStackingSelectionWriteSignature = signature;
  pendingStackingSelectionWriteMarker = Number.isSafeInteger(cell.marker) ? cell.marker : undefined;
  seamlessState.stackingSelectionPendingWrite = true;
  beginPendingShareWrite('stacking_selection', serialized, marker);
  let accepted = false;
  try { accepted = cell.set(serialized); } catch { accepted = false; }
  if (!accepted) {
    pendingStackingSelectionWriteSignature = undefined;
    pendingStackingSelectionWriteMarker = undefined;
    seamlessState.stackingSelectionPendingWrite = false;
    resolvePendingShareWrite('stacking_selection', false);
    seamlessState.stackingSelectionError = 'The shared stacking selection rejected the update.';
    return false;
  }
  recordDecision('stacking_selection', marker);
  seamlessState.stackingSelectionError = '';
  return true;
}

function receiveStackings(value: SharedValue, marker?: number): void {
  seamlessState.stackingRevision = updateBackendRevision(seamlessState.stackingRevision, marker);
  const payload = parseJsonValue(value);
  seamlessState.computed.stackings = payload;
  const rows = payload && typeof payload === 'object' && !Array.isArray(payload) && 'rows' in payload
    ? (payload as { rows: unknown }).rows
    : payload;
  const stackings = normalizeStackings(rows);
  received.add('stackings');
  refreshComputedFreshness();
  if (!stackings) {
    seamlessState.stackingsReady = false;
    seamlessState.stackingsError = 'The shared stacking table is not a valid JSON list of stacking rows.';
    return;
  }

  const selectionCell = client?.stacking_selection;
  if (!seamlessState.stackingSelectionReady && selectionCell?.value !== undefined) {
    receiveStackingSelection(selectionCell.value, selectionCell.marker);
  }

  seamlessState.stackingsReady = true;
  seamlessState.hasStackingsSnapshot = true;
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
  const storedSelection = markerSelection(seamlessState.projectSnapshot.stacking_selection);
  const selectionApplies = storedSelection !== null
    && markerMatchesUpstream('stacking_selection', seamlessState.projectSnapshot.stacking_selection, currentSnapshot());
  appState.stackingSelection = stackings.map((row, index) => {
    const key = stackingKey(row);
    if (selectionApplies && storedSelection && Object.hasOwn(storedSelection, key)) return storedSelection[key];
    return merged.stackingSelection[index];
  });

  if (pendingBackendStackingSelection !== undefined && stackings.length > 0) {
    const queued = pendingBackendStackingSelection;
    pendingBackendStackingSelection = undefined;
    appState.stackingSelection = stackings.map(stacking => queued[stackingKey(stacking)] ?? false);
    backendStackingSelectionSignature = stableJson(queued);
    seamlessState.stackingSelectionReady = true;
    seamlessState.stackingSelectionError = '';
  }
}

function receiveStackingSelection(value: SharedValue, marker?: number): void {
  if (observePendingWrite('stacking_selection', value, marker, applyingAuthoritativeRead)) {
    if (applyingAuthoritativeRead) { received.add('stacking_selection'); updateReady(); }
    return;
  }
  seamlessState.stackingSelectionRevision = updateBackendRevision(
    seamlessState.stackingSelectionRevision,
    marker,
  );
  const parsed = parseStackingMarker(value);
  const selection = markerSelection(parsed);
  const signature = stableJson(parsed ?? null);
  let newerThanPendingWrite = false;

  if (pendingStackingSelectionWriteSignature !== undefined) {
    if (signature === pendingStackingSelectionWriteSignature) {
      seamlessState.stackingSelectionReady = selection !== null;
      seamlessState.stackingSelectionNonEmpty = selection !== null && Object.values(selection).some(Boolean);
      backendStackingSelectionSignature = signature;
      if (selection !== null) {
        clearPendingInput('stacking_selection');
        recordDecision('stacking_selection', parsed);
      }
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
  seamlessState.stackingSelectionReady = selection !== null;
  seamlessState.stackingSelectionNonEmpty = selection !== null && Object.values(selection).some(Boolean);
  if (!newerThanPendingWrite && signature === backendStackingSelectionSignature) return;
  backendStackingSelectionSignature = signature;
  received.add('stacking_selection');
  if (markerSelection(parsed) !== null) clearPendingInput('stacking_selection');
  recordDecision('stacking_selection', selection === null ? null : parsed);

  if (selection === null) {
    pendingBackendStackingSelection = undefined;
    seamlessState.stackingSelectionError = parsed === null ? '' : 'The shared stacking selection decision is invalid.';
    return;
  }
  const appliesToCurrentRanges = selection !== null
    && markerMatchesUpstream('stacking_selection', parsed, currentSnapshot());
  if (appState.stackings.length > 0) {
    pendingBackendStackingSelection = undefined;
    appState.stackingSelection = appState.stackings.map(stacking => {
      const stored = appliesToCurrentRanges ? selection[stackingKey(stacking)] : undefined;
      return stored ?? (stacking.distance <= 5 && stacking.slip_angle < 60);
    });
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
  const install = (handler: (value: SharedValue) => void) => {
    shareReceivers.set(key, handler);
    cell.onchange = value => {
      if (!resyncInProgress) handler(value);
    };
  };
  switch (key) {
    case 'pdbcode':
      install(receivePdbCode);
      break;
    case 'biounit':
      install(receiveBiounitIndex);
      break;
    case 'ranges':
      install(receiveRanges);
      break;
    case 'stackings':
      install(value => receiveStackings(value, cell.marker));
      break;
    case 'stacking_selection':
      install(value => receiveStackingSelection(value, cell.marker));
      break;
    case 'threshold_edits':
    case 'best_fit_clip':
    case 'interaction_matrix':
    case 'gaussian_classes':
      install(value => receiveDecision(key, value));
      break;
    case 'energy_threshold_edits':
      install(value => {
        seamlessState.computed.energy_threshold_edits = parseJsonValue(value) ?? null;
        received.add(key);
        updateReady();
      });
      break;
    case 'range_specs':
    case 'structure_id':
    case 'range_pdbs_aa':
    case 'fragments':
    case 'reference_fragment_pdbs':
    case 'best_fit_thresholds':
    case 'best_fits':
    case 'energetics':
    case 'energy_thresholds':
    case 'interactions':
    case 'bestfit_stackings':
    case 'gaussian_validation':
    case 'refe_validation':
    case 'gaussian_scores':
      install(value => receiveComputed(key, value));
      break;
  }
}

async function readShareFresh(key: ShareKey): Promise<void> {
  const cell = client?.[key];
  const receiver = shareReceivers.get(key);
  if (!cell || !cell.url || !receiver || typeof fetch !== 'function') {
    throw new Error(`Shared ${key} cannot be reread.`);
  }
  const base = client?.self.server ?? (typeof window === 'undefined' ? '' : window.location.origin);
  const url = new URL(cell.url, `${base.replace(/\/$/, '')}/`);
  let lastError: unknown;
  for (let attempt = 0; attempt < 3; attempt += 1) {
    try {
      const response = await fetch(url, { cache: 'no-store' });
      let value: SharedValue;
      if (response.status === 204 || response.status === 404) value = null;
      else if (!response.ok) throw new Error(`HTTP ${response.status}`);
      else value = cell.binary ? await response.blob() : await response.text();
      const markerHeader = response.headers.get('x-seamless-marker');
      const serverMarker = markerHeader === null ? undefined : Number(markerHeader);
      if (Number.isSafeInteger(serverMarker) && Number.isSafeInteger(cell.marker)
        && (serverMarker as number) < (cell.marker as number)) {
        throw new Error('A newer share update arrived during the reread.');
      }
      if (Number.isSafeInteger(serverMarker)) cell.marker = serverMarker as number;
      cell.value = value;
      applyingAuthoritativeRead = true;
      receiver(value);
      applyingAuthoritativeRead = false;
      return;
    } catch (error) {
      applyingAuthoritativeRead = false;
      lastError = error;
    }
  }
  throw lastError instanceof Error ? lastError : new Error(`Shared ${key} could not be reread.`);
}

async function readPendingCifFresh(): Promise<void> {
  const cell = client?.['cif.cif'];
  if (!cell?.url) throw new Error('The shared CIF cannot be reread.');
  const base = client?.self.server ?? (typeof window === 'undefined' ? '' : window.location.origin);
  const response = await fetch(new URL(cell.url, `${base.replace(/\/$/, '')}/`), { cache: 'no-store' });
  if (!response.ok) throw new Error(`The shared CIF could not be reread (HTTP ${response.status}).`);
  const value = await response.text();
  const markerHeader = response.headers.get('x-seamless-marker');
  const serverMarker = markerHeader === null ? undefined : Number(markerHeader);
  if (Number.isSafeInteger(serverMarker)) cell.marker = serverMarker as number;
  cell.value = value;
  applyingAuthoritativeRead = true;
  receiveCifCommit(value);
  applyingAuthoritativeRead = false;
}

async function resynchronizeShares(keys: string[]): Promise<void> {
  const ordered: ShareKey[] = [
    'pdbcode', 'biounit', 'ranges', 'stacking_selection', 'threshold_edits',
    'best_fit_clip', 'interaction_matrix', 'gaussian_classes', 'energy_threshold_edits',
    'range_specs', 'structure_id', 'range_pdbs_aa', 'stackings', 'fragments', 'reference_fragment_pdbs',
    'best_fit_thresholds', 'best_fits', 'energetics', 'energy_thresholds', 'interactions',
    'bestfit_stackings', 'gaussian_validation', 'refe_validation', 'gaussian_scores',
  ];
  const tasks = ordered.filter(key => keys.includes(key)).map(readShareFresh);
  if (pendingStructureCommit && keys.includes('cif.cif')) tasks.push(readPendingCifFresh());
  const results = await Promise.allSettled(tasks);
  const failed = results.filter(result => result.status === 'rejected').length;
  resyncInProgress = false;
  if (failed > 0) {
    seamlessState.syncNotice = `Could not reread ${failed} shared value${failed === 1 ? '' : 's'} after reconnect. Editing stays paused until they can be read.`;
  } else {
    seamlessState.syncNotice = '';
  }
  updateReady();

  if (!seamlessState.ready) return;
  for (const [key, pending] of [...pendingShareWrites]) {
    const cell = key === 'cif.cif' ? client?.['cif.cif'] : client?.[key as ShareKey];
    if (!cell) continue;
    const freshServerValue = cell.value;
    const currentValue = canonicalShareValue(key, freshServerValue);
    if (stableJson(currentValue) === stableJson(pending.expected)) {
      resolvePendingShareWrite(key, true);
      continue;
    }
    const currentRevision = latestWriteRevisions.get(key) ?? pending.revision;
    const admissible = deferredWriteIsCurrent({
      stageAtIssue: pending.stageAtIssue,
      currentStage: seamlessState.projectStage,
      live: seamlessState.ready,
      revisionAtIssue: pending.revision,
      currentRevision,
    });
    if (admissible && (pending.baseMarker === undefined || cell.marker === pending.baseMarker)) {
      let accepted = false;
      try { accepted = cell.set(pending.serialized); } catch { accepted = false; }
      if (accepted) continue;
    }
    let correctionQueued = false;
    try { correctionQueued = cell.set(serializeFreshShareValue(key, freshServerValue)); } catch { correctionQueued = false; }
    cell.value = freshServerValue;
    resolvePendingShareWrite(key, false, true);
    if (key !== 'cif.cif') shareReceivers.get(key as ShareKey)?.(freshServerValue);
    else receiveCifCommit(freshServerValue);
    if (!correctionQueued) {
      seamlessState.syncNotice = `A stale ${key.replaceAll('_', ' ')} retry could not be corrected to the freshly read shared value.`;
    }
  }
}

function attachShareContext(context: ShareContext): void {
  client = context;
  shareReceivers.clear();
  stackingsSignature = undefined;
  backendStackingSelectionSignature = undefined;
  pendingStackingSelectionWriteSignature = undefined;
  pendingStackingSelectionWriteMarker = undefined;
  pendingBackendStackingSelection = undefined;
  seamlessState.stackingSelectionPendingWrite = false;
  context.self.onstategraph = receiveStateGraph;
  context.self.onsharelist = keys => {
    const reconnect = hasConnectedBefore;
    hasConnectedBefore = true;
    if (reconnect) structureConnectionRevision += 1;
    resyncInProgress = reconnect;
    if (reconnect) {
      seamlessState.ready = false;
      seamlessState.computedFreshness = {};
      seamlessState.structureIdReady = false;
      seamlessState.hasInitialStateGraph = false;
      seamlessState.stateGraphErrors = [];
      lastStateGraphNodes = [];
      received.clear();
      seamlessState.status = 'waiting-for-shares';
      seamlessState.stackingsReady = false;
      seamlessState.rangesReady = false;
      stackingsSignature = undefined;
      backendStackingSelectionSignature = undefined;
      pendingBackendStackingSelection = undefined;
      seamlessState.stackingSelectionReady = false;
    }
    attachShare('pdbcode', keys);
    attachShare('biounit', keys);
    attachShare('ranges', keys);
    // Selection must be observed before stackings so a persisted exact-length
    // value can supersede defaults for the table being attached.
    attachShare('stacking_selection', keys);
    attachShare('stackings', keys);
    for (const key of [
      'threshold_edits', 'best_fit_clip', 'interaction_matrix', 'gaussian_classes', 'energy_threshold_edits',
      'range_specs', 'structure_id', 'range_pdbs_aa', 'fragments', 'reference_fragment_pdbs', 'best_fit_thresholds', 'best_fits',
      'energetics', 'energy_thresholds', 'interactions', 'bestfit_stackings', 'gaussian_validation',
      'refe_validation', 'gaussian_scores',
    ] as const) attachShare(key, keys);
    const cifCell = context['cif.cif'];
    if (cifCell) cifCell.onchange = value => { if (!resyncInProgress) receiveCifCommit(value); };
    if (!requiredEditShares.every(key => keys.includes(key))
      || !requiredComputedShares.every(key => keys.includes(key))) {
      seamlessState.status = 'waiting-for-shares';
    }
    if (reconnect && context.self.ws) void resynchronizeShares(keys);
    else resyncInProgress = false;
  };
}

export function connectToSeamless(connect: ShareConnector): void {
  structureConnectionRevision += 1;
  structureShareRevisions = { pdbcode: 0, cif: 0 };
  for (const key of [...pendingShareWrites.keys()]) resolvePendingShareWrite(key, false);
  pendingShareWrites.clear();
  latestWriteRevisions.clear();
  seamlessState.pendingInputs = [];
  seamlessState.ready = false;
  seamlessState.status = 'connecting';
  seamlessState.projectSnapshot = {
    pdbcode: '',
    biounit: 1,
    ranges: null,
    stacking_selection: null,
    threshold_edits: null,
    best_fit_clip: { schema: 1, basis: null, clip1: 0, clip2: 0 },
    interaction_matrix: null,
    gaussian_classes: null,
  };
  seamlessState.computed = {};
  seamlessState.computedFreshness = {};
  seamlessState.decisionErrors = {};
  seamlessState.rangeError = '';
  seamlessState.stackingSelectionError = '';
  seamlessState.stackingSelectionReady = false;
  seamlessState.stackingSelectionNonEmpty = false;
  seamlessState.stackingsReady = false;
  seamlessState.hasStackingsSnapshot = false;
  seamlessState.rangesReady = false;
  seamlessState.structureIdReady = false;
  seamlessState.hasInitialStateGraph = false;
  received.clear();
  appState.stackings = [];
  appState.stackingSelection = [];
  appState.staleStackingRevision = null;
  resyncInProgress = false;
  hasConnectedBefore = false;
  websocketWasOpen = false;
  refreshProjectStage();
  attachShareContext(connect(null, null, 'preparic'));
  if (websocketPoll !== undefined) clearInterval(websocketPoll);
  websocketPoll = setInterval(watchShareConnection, 250);
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
}

/** Explicitly read the file-authoritative CIF; dotted shares are not auto-read. */
export async function readSharedCif(): Promise<string> {
  const cell = client?.['cif.cif'];
  if (!cell) throw new Error('The shared CIF cell is not available.');
  if (!cell.url || typeof fetch !== 'function') throw new Error('The shared CIF could not be read.');
  const base = client?.self.server ?? (typeof window === 'undefined' ? '' : window.location.origin);
  const response = await fetch(new URL(cell.url, `${base.replace(/\/$/, '')}/`), { cache: 'no-store' });
  if (!response.ok) throw new Error(`The shared CIF could not be read (HTTP ${response.status}).`);
  const cif = await response.text();
  if (!cif.trim()) throw new Error('The stored shared CIF is empty.');
  return cif;
}

/** Capture the live project target before starting a CIF download. */
export function captureSharedStructureTarget(): SharedStructureTarget {
  return {
    stageAtIssue: seamlessState.projectStage,
    liveAtIssue: seamlessState.ready && seamlessState.status === 'live',
    connectionRevision: structureConnectionRevision,
    biounitIndex: seamlessState.biounitIndex,
    pdbCodeRevision: structureShareRevisions.pdbcode,
    cifRevision: structureShareRevisions.cif,
    pdbCodeMarker: client?.pdbcode?.marker,
    cifMarker: client?.['cif.cif']?.marker,
  };
}

export function sharedStructureTargetIsCurrent(target: SharedStructureTarget): boolean {
  const codeCell = client?.pdbcode;
  const cifCell = client?.['cif.cif'];
  return target.liveAtIssue
    && seamlessState.ready
    && seamlessState.status === 'live'
    && seamlessState.projectStage === target.stageAtIssue
    && target.connectionRevision === structureConnectionRevision
    && target.biounitIndex === seamlessState.biounitIndex
    && target.pdbCodeRevision === structureShareRevisions.pdbcode
    && target.cifRevision === structureShareRevisions.cif
    && !!codeCell && !!cifCell
    && target.pdbCodeMarker === codeCell.marker
    && target.cifMarker === cifCell.marker;
}

/** Commit the downloaded structure together with its matching project code. */
export function commitSharedStructure(value: string, target: SharedStructureTarget): boolean {
  const code = normalizePdbCode(value);
  const cifCell = client?.['cif.cif'];
  const codeCell = client?.pdbcode;
  if (!sharedStructureTargetIsCurrent(target)) {
    pendingCif = undefined;
    pendingPdbCode = undefined;
    seamlessState.pdbCodePendingWrite = false;
    clearPendingInput('pdbcode');
    seamlessState.syncNotice = 'The shared project changed while this structure was loading; the downloaded CIF was discarded.';
    return false;
  }
  if (!code || code !== normalizePdbCode(appState.pdbId) || !codeCell) return false;
  if (normalizePdbCode(codeCell.value) === code && !pendingStructureCommit) {
    pendingCif = undefined;
    pendingPdbCode = undefined;
    pendingStructureCommit = undefined;
    seamlessState.pdbCode = code;
    seamlessState.projectSnapshot.pdbcode = code;
    seamlessState.pdbCodePendingWrite = false;
    clearPendingInput('pdbcode');
    return true;
  }
  if (seamlessState.projectStage !== 1 || !pendingCif || !cifCell || cifCell.readonly || codeCell.readonly) return false;
  const cif = pendingCif;
  pendingPdbCode = code;
  pendingStructureCommit = {
    code,
    cif,
    codeAcknowledged: false,
    cifAcknowledged: false,
    cifConflict: false,
    retryingCif: false,
  };
  cifCell.onchange = value => { if (!resyncInProgress) receiveCifCommit(value); };
  markPendingInput('pdbcode');
  markPendingInput('cif.cif');
  beginPendingShareWrite('cif.cif', cif, cif);
  beginPendingShareWrite('pdbcode', code, code);
  let cifAccepted = false;
  try { cifAccepted = cifCell.set(cif); } catch { cifAccepted = false; }
  if (!cifAccepted) {
    pendingStructureCommit = undefined;
    resolvePendingShareWrite('pdbcode', false);
    resolvePendingShareWrite('cif.cif', false);
    return false;
  }
  let codeAccepted = false;
  try { codeAccepted = codeCell.set(code); } catch { codeAccepted = false; }
  if (!codeAccepted) {
    pendingStructureCommit = undefined;
    resolvePendingShareWrite('pdbcode', false);
    resolvePendingShareWrite('cif.cif', false);
    return false;
  }
  pendingCif = undefined;
  seamlessState.pdbCode = code;
  seamlessState.projectSnapshot.pdbcode = code;
  seamlessState.pdbCodePendingWrite = true;
  refreshProjectStage();
  refreshComputedFreshness();
  return true;
}

export function cancelSharedCifWrite(): void {
  pendingCif = undefined;
  if (!pendingStructureCommit) {
    pendingPdbCode = undefined;
    seamlessState.pdbCodePendingWrite = false;
    clearPendingInput('pdbcode');
  }
}

export function setSharedPdbCode(value: string): boolean {
  if (!seamlessState.ready || seamlessState.status !== 'live' || seamlessState.projectStage !== 1) return false;
  const code = value.trim().toUpperCase();
  pendingPdbCode = code;
  seamlessState.pdbCodePendingWrite = true;
  return true;
}

export function setSharedBiounitIndex(value: number): boolean {
  const index = Number.isSafeInteger(value) && value >= 1 ? value : 1;
  const cell = client?.biounit;
  if (!client || !seamlessState.ready || !cell) {
    pendingBiounitIndex = index;
    seamlessState.biounitIndexPendingWrite = true;
    seamlessState.biounitIndex = index;
    seamlessState.projectSnapshot.biounit = index;
    refreshProjectStage();
    refreshComputedFreshness();
    return true;
  }
  if (cell.readonly) return false;
  const existingWrite = pendingShareWrites.get('biounit');
  if (existingWrite && existingWrite.expected === index) return true;
  if (!pendingShareWrites.has('biounit') && parseBiounitIndex(cell.value) === index) {
    pendingBiounitIndex = undefined;
    seamlessState.biounitIndexPendingWrite = false;
    seamlessState.biounitIndex = index;
    return true;
  }
  pendingBiounitIndex = index;
  seamlessState.biounitIndexPendingWrite = true;
  markPendingInput('biounit');
  seamlessState.projectSnapshot.biounit = index;
  refreshProjectStage();
  refreshComputedFreshness();
  beginPendingShareWrite('biounit', String(index), index);
  let accepted = false;
  try { accepted = cell.set(String(index)); } catch { accepted = false; }
  if (accepted) seamlessState.biounitIndex = index;
  else {
    resolvePendingShareWrite('biounit', false);
    clearPendingInput('biounit');
  }
  return accepted;
}

export function setSharedRangeRows(rows: RangeRow[]): boolean {
  const cell = client?.ranges;
  if (!client || !seamlessState.rangesReady || !cell || cell.readonly) return false;
  const completeRows: Array<{ chain: string; start: number; end: number }> = [];
  let invalidRow = false;
  for (const row of rows) {
    const chain = appState.chains.find(item => item.name === row.chain);
    const start = Number(row.start.trim());
    const end = Number(row.end.trim());
    const observed = chain && Number.isSafeInteger(start) && Number.isSafeInteger(end) && start <= end
      ? [...chain.seqIds].filter(resid => resid >= start && resid <= end).length
      : 0;
    if (!chain || !row.start.trim() || !row.end.trim() || !Number.isSafeInteger(start) || !Number.isSafeInteger(end)
      || !chain.seqIds.has(start) || !chain.seqIds.has(end) || start > end || observed !== end - start + 1) {
      invalidRow = true;
      continue;
    }
    completeRows.push({ chain: row.chain, start, end });
  }
  if (rows.length === 0 || invalidRow || completeRows.length !== rows.length) {
    seamlessState.rangeError = 'Complete every range with an observed chain and ordered residue endpoints before saving.';
    return false;
  }
  const marker = {
    schema: 1,
    basis: { pdbcode: seamlessState.pdbCode, biounit: seamlessState.biounitIndex },
    active: true,
    rows: completeRows,
  };
  if (!markerIsWellFormed('ranges', marker)) {
    seamlessState.rangeError = 'Completed ranges must be unique and use valid residue endpoints.';
    return false;
  }
  const signature = stableJson(marker);
  if (signature === rangeSignature) return true;
  if (!writeDecision('ranges', marker)) {
    return false;
  }
  seamlessState.rangeError = '';
  rangeSignature = signature;
  seamlessState.rangeRows = completeRows.map((row, index) => ({ id: `shared-${index + 1}`, alwaysVisible: false, chain: row.chain, start: String(row.start), end: String(row.end) }));
  seamlessState.rangeRevision += 1;
  recordDecision('ranges', marker);
  return true;
}

export function setProjectDecision(key: DecisionKey, marker: unknown): boolean {
  if (!markerIsWellFormed(key, marker)) return false;
  const current = seamlessState.projectSnapshot[key];
  if (stableJson(current) === stableJson(marker)) return true;
  if (!writeDecision(key, marker)) {
    return false;
  }
  recordDecision(key, marker);
  return true;
}

export async function setProjectDecisionAwaited(key: DecisionKey, marker: unknown): Promise<boolean> {
  if (!setProjectDecision(key, marker)) return false;
  return waitForShareWrite(key);
}

export function deactivateProjectDecision(key: DecisionKey): boolean {
  const current = seamlessState.projectSnapshot[key];
  if (!current || typeof current !== 'object' || !('active' in current)) return true;
  return setProjectDecision(key, { ...current, active: false });
}

export async function deactivateProjectDecisionAndWait(key: DecisionKey): Promise<boolean> {
  const current = seamlessState.projectSnapshot[key];
  if (!current || typeof current !== 'object' || !('active' in current) || current.active !== true) return true;
  return setProjectDecisionAwaited(key, { ...current, active: false });
}

export async function prepareProjectTransition(targetStage: 2 | 3 | 4 | 5 | 6): Promise<boolean> {
  const stageMarkers: Array<[DecisionKey, number]> = [
    ['ranges', 2], ['stacking_selection', 3], ['threshold_edits', 4],
    ['interaction_matrix', 5], ['gaussian_classes', 6],
  ];
  for (const [key, stage] of stageMarkers) {
    if (stage <= targetStage) continue;
    if (!await deactivateProjectDecisionAndWait(key)) return false;
  }
  return true;
}

export function setSharedStackingSelection(selection: readonly boolean[]): boolean {
  // An explicit checkbox edit supersedes a delayed backend list that was
  // waiting for a table of a different length.
  pendingBackendStackingSelection = undefined;
  if (selection.length === 0) backendStackingSelectionSignature = undefined;
  const accepted = writeStackingSelection(selection);
  if (accepted) {
    seamlessState.stackingSelectionReady = true;
    seamlessState.stackingSelectionNonEmpty = selection.length > 0;
  }
  return accepted;
}

export async function setSharedStackingSelectionAndWait(selection: readonly boolean[]): Promise<boolean> {
  if (!setSharedStackingSelection(selection)) return false;
  return waitForShareWrite('stacking_selection');
}
