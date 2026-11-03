/** Pure state and validation rules for the shared Preparic project.
 *
 * The stage is derived from the persisted decision markers alone.  Computed
 * shares and state-graph status are deliberately absent from this module's
 * stage calculation.
 */

export type ProjectStage = 1 | 2 | 3 | 4 | 5 | 6;
export type RangeAddress = { chain: string; start: number; end: number };

export interface ProjectSnapshot {
  pdbcode: string;
  biounit: number;
  ranges: unknown;
  stacking_selection: unknown;
  threshold_edits: unknown;
  best_fit_clip: unknown;
  interaction_matrix: unknown;
  gaussian_classes: unknown;
}

type AnyRecord = Record<string, unknown>;
type MarkerKind = 'ranges' | 'stacking_selection' | 'threshold_edits' | 'best_fit_clip' | 'interaction_matrix' | 'gaussian_classes';

const TRANSFORMER_STAGE: Record<string, ProjectStage> = {
  structure_id: 1,
  parse_ranges: 2,
  ranges_to_pdb: 2,
  pdbs_to_folder: 2,
  detect_stacking: 3,
  aareduce_ranges: 3,
  refe_fragments: 3,
  list_fragments: 3,
  share_range_pdbs_aa: 3,
  best_fit: 4,
  best_fit_tables: 4,
  best_fit_thresholds: 4,
  parse_best_fits: 4,
  validate_threshold_edits: 4,
  clip_fragments: 4,
  energetics_ligands: 5,
  energetics: 5,
  energetics_range: 5,
  merge_energetics: 5,
  interaction_handoff: 5,
  energy_thresholds: 5,
  validate_energy_threshold_edits: 5,
  detect_stacking_bestfit: 5,
  validate_gaussian_classes: 6,
  gaussian_check_jobs: 6,
  gaussian_check_job_json: 6,
  gaussian_semantic_identity: 6,
  gaussian_check_transformer: 6,
  gaussian_check: 6,
  gaussian_scores_transformer: 5,
  merge_gaussian_validation: 6,
  validate_refe_anchors: 6,
  gaussian_scores: 6,
};

const TRANSFORMER_INPUTS: Record<string, ReadonlySet<string>> = {
  structure_id: new Set(['cif.cif']),
  parse_ranges: new Set(['ranges', 'biounit']),
  ranges_to_pdb: new Set(['cif.cif', 'range_specs']),
  pdbs_to_folder: new Set(['range_pdbs']),
  detect_stacking: new Set(['range_pdbs', 'pdbcode', 'biounit']),
  aareduce_ranges: new Set(['range_pdb_files', 'aareduce_files']),
  refe_fragments: new Set(['range-1-aa.pdb']),
  list_fragments: new Set(['pdbcode', 'range-1-aa.pdb']),
  share_range_pdbs_aa: new Set(['range_specs', 'exclude_pdb_code', 'range_pdb_files_aa']),
  best_fit: new Set(['pdbcode', 'ranges', 'stacking_selection', 'range-1-aa.pdb']),
  best_fit_tables: new Set(['best_fit', 'range-1-aa.pdb', 'pdbcode']),
  best_fit_thresholds: new Set(['best_fit_tables']),
  parse_best_fits: new Set(['best_fit', 'range-1-aa.pdb', 'pdbcode']),
  validate_threshold_edits: new Set(['threshold_edits', 'fragments']),
  clip_fragments: new Set(['threshold_edits', 'best_fit_clip', 'parse_best_fits']),
  energetics_ligands: new Set(['refe_fragments', 'best_fit', 'fragments']),
  energetics_range: new Set(['receptor', 'energetics_ligands']),
  energetics: new Set(['energetics_range_results', 'best_fits', 'range_specs', 'fragments']),
  merge_energetics: new Set(['energetics_range_results', 'best_fits', 'range_specs', 'fragments']),
  interaction_handoff: new Set(['interaction_matrix', 'stacking_selection', 'stackings', 'best_fit_clip', 'fragments', 'energetics']),
  energy_thresholds: new Set(['energetics', 'best_fit']),
  validate_energy_threshold_edits: new Set(['energy_threshold_edits', 'fragments', 'protein_ranges']),
  detect_stacking_bestfit: new Set(['best_fits', 'proteins', 'fragments', 'range_specs']),
  validate_gaussian_classes: new Set(['gaussian_classes', 'gaussian_model_json']),
  gaussian_check_jobs: new Set(['best_fits', 'fragments', 'stackings', 'stacking_selection', 'best_fit_clip']),
  gaussian_check_job_json: new Set(['gaussian_check_jobs']),
  gaussian_semantic_identity: new Set(['gaussian_model_json']),
  gaussian_check_transformer: new Set(['gaussian_check_job_json', 'best_fit', 'proteins', 'exclude_pdb_code', 'fraglib_identities', 'gaussian_semantic_identity']),
  gaussian_check: new Set(['gaussian_check_transformer']),
  gaussian_scores_transformer: new Set(['bestfit_stackings', 'stackings', 'stacking_selection', 'best_fit_clip', 'fragments', 'gaussian_model_json']),
  merge_gaussian_validation: new Set(['gaussian_check', 'gaussian_check_jobs', 'gaussian_classes', 'stacking_selection', 'best_fit_clip', 'best_fits', 'gaussian_model_json', 'validate_gaussian_classes']),
  validate_refe_anchors: new Set(['threshold_edits', 'best_fit_thresholds', 'fragments', 'stackings', 'stacking_selection', 'best_fit_clip', 'best_fits']),
  gaussian_scores: new Set(['bestfit_stackings', 'stackings', 'stacking_selection', 'best_fit_clip', 'fragments', 'gaussian_model_json']),
};

function isRecord(value: unknown): value is AnyRecord {
  return value !== null && typeof value === 'object' && !Array.isArray(value);
}

function hasOnlyKeys(value: AnyRecord, keys: readonly string[]): boolean {
  return Object.keys(value).length === keys.length && keys.every(key => Object.hasOwn(value, key));
}

function isInteger(value: unknown): value is number {
  return typeof value === 'number' && Number.isSafeInteger(value);
}

function isUnique<T>(items: readonly T[]): boolean {
  return new Set(items).size === items.length;
}

function isRangeAddress(value: unknown): value is RangeAddress {
  if (!isRecord(value) || !hasOnlyKeys(value, ['chain', 'start', 'end'])) return false;
  return typeof value.chain === 'string' && value.chain.length > 0
    && isInteger(value.start) && isInteger(value.end) && value.end >= value.start;
}

function sameValue(left: unknown, right: unknown): boolean {
  return stableStringify(left) === stableStringify(right);
}

function stableStringify(value: unknown): string {
  if (Array.isArray(value)) return `[${value.map(stableStringify).join(',')}]`;
  if (isRecord(value)) {
    return `{${Object.keys(value).sort().map(key => `${JSON.stringify(key)}:${stableStringify(value[key])}`).join(',')}}`;
  }
  return JSON.stringify(value) ?? 'undefined';
}

function validBasisProject(value: unknown): value is { pdbcode: string; biounit: number; rows?: RangeAddress[] } {
  if (!isRecord(value)) return false;
  const hasRows = Object.hasOwn(value, 'rows');
  if (!hasOnlyKeys(value, hasRows ? ['pdbcode', 'biounit', 'rows'] : ['pdbcode', 'biounit'])) return false;
  return typeof value.pdbcode === 'string' && value.pdbcode.length > 0
    && isInteger(value.biounit) && value.biounit >= 1
    && (!hasRows || (Array.isArray(value.rows) && value.rows.length > 0
      && value.rows.every(isRangeAddress)
      && isUnique(value.rows.map(row => `${row.chain}\u0000${row.start}\u0000${row.end}`))));
}

function validRnaBasis(value: unknown): value is { pdbcode: string; rna: RangeAddress } {
  return isRecord(value) && hasOnlyKeys(value, ['pdbcode', 'rna'])
    && typeof value.pdbcode === 'string' && value.pdbcode.length > 0
    && isRangeAddress(value.rna);
}

function validRoundedPositive(value: unknown): boolean {
  if (typeof value !== 'string' || !/^(?:[1-9]\d*(?:\.\d{1,2})?|0\.(?:[1-9]\d?|0[1-9]))$/.test(value)) return false;
  const numeric = Number(value);
  return Number.isFinite(numeric) && numeric > 0;
}

function validFragmentKey(value: unknown): value is string {
  return typeof value === 'string' && /^frag-[1-9]\d*-[ACGU]{2}$/.test(value);
}

function validSelectionBasis(value: unknown): value is string[] {
  return Array.isArray(value) && value.every(item => typeof item === 'string' && item.length > 0) && isUnique(value);
}

function selectedKeys(marker: unknown): string[] | null {
  if (!isRecord(marker) || !isRecord(marker.selected)) return null;
  if (!Object.values(marker.selected).every(value => typeof value === 'boolean')) return null;
  const selected = marker.selected;
  return Object.keys(selected).filter(key => selected[key] === true).sort();
}

function selectedBasisMatches(basis: unknown, selection: unknown): boolean {
  const keys = selectedKeys(selection);
  return keys !== null && validSelectionBasis(basis)
    && sameValue([...basis].sort(), keys);
}

function isRangeMarker(value: unknown): value is AnyRecord & { active: boolean; rows: RangeAddress[] } {
  return isRecord(value) && hasOnlyKeys(value, ['schema', 'basis', 'active', 'rows'])
    && value.schema === 1 && validBasisProject(value.basis)
    && typeof value.active === 'boolean'
    && Array.isArray(value.rows) && value.rows.length > 0
    && value.rows.every(isRangeAddress)
    && isUnique(value.rows.map(row => `${row.chain}\u0000${row.start}\u0000${row.end}`));
}

function isSelectionMarker(value: unknown): value is AnyRecord & { active: boolean } {
  return isRecord(value) && hasOnlyKeys(value, ['schema', 'basis', 'active', 'selected'])
    && value.schema === 1 && validBasisProject(value.basis)
    && typeof value.active === 'boolean'
    && isRecord(value.selected)
    && Object.keys(value.selected).every(key => key.length > 0)
    && Object.values(value.selected).every(item => typeof item === 'boolean');
}

function validThresholdMarker(value: unknown): boolean {
  if (!isRecord(value) || !hasOnlyKeys(value, ['schema', 'basis', 'active', 'pairs', 'bases']) || value.schema !== 1
    || !validRnaBasis(value.basis) || typeof value.active !== 'boolean'
    || !Array.isArray(value.pairs) || !Array.isArray(value.bases)) return false;

  const pairIds: string[] = [];
  for (const pair of value.pairs) {
    if (!isRecord(pair) || !hasOnlyKeys(pair, ['down', 'up', 'sequence', 'ovrmsd', 'crmsd'])
      || !isInteger(pair.down) || pair.down < 1 || !isInteger(pair.up) || pair.up !== pair.down + 1
      || typeof pair.sequence !== 'string' || !/^[ACGU]{3}$/.test(pair.sequence)
      || !(pair.ovrmsd === null || validRoundedPositive(pair.ovrmsd))
      || !(pair.crmsd === null || validRoundedPositive(pair.crmsd))) return false;
    pairIds.push(`${pair.down}:${pair.up}`);
  }
  if (!isUnique(pairIds)) return false;

  const baseIds: number[] = [];
  for (const base of value.bases) {
    if (!isRecord(base) || !hasOnlyKeys(base, ['fragment', 'sequence', 'nuc1_ovrmsd', 'nuc2_ovrmsd'])
      || !isInteger(base.fragment) || base.fragment < 1
      || typeof base.sequence !== 'string' || !/^[ACGU]{2}$/.test(base.sequence)
      || !(base.nuc1_ovrmsd === null || validRoundedPositive(base.nuc1_ovrmsd))
      || !(base.nuc2_ovrmsd === null || validRoundedPositive(base.nuc2_ovrmsd))) return false;
    baseIds.push(base.fragment);
  }
  return isUnique(baseIds);
}

function validClip(value: unknown): boolean {
  if (!isRecord(value) || !hasOnlyKeys(value, ['schema', 'basis', 'clip1', 'clip2']) || value.schema !== 1
    || !isInteger(value.clip1) || value.clip1 < 0 || !isInteger(value.clip2) || value.clip2 < 0) return false;
  if (value.basis === null) return value.clip1 === 0 && value.clip2 === 0;
  return validRnaBasis(value.basis);
}

function validMatrix(value: unknown): boolean {
  if (!isRecord(value) || !hasOnlyKeys(value, ['schema', 'basis', 'active', 'rows', 'ranges', 'grid'])
    || value.schema !== 1 || !validSelectionBasis(value.basis) || typeof value.active !== 'boolean'
    || !Array.isArray(value.rows) || !value.rows.every(validFragmentKey) || !isUnique(value.rows)
    || !Array.isArray(value.ranges) || !value.ranges.every(item => typeof item === 'string' && /^range-[1-9]\d*$/.test(item))
    || !Array.isArray(value.grid)) return false;
  const rangeCount = value.ranges.length;
  return value.grid.every((row: unknown) => Array.isArray(row) && row.length === rangeCount
    && row.every((cell: unknown) => typeof cell === 'boolean'));
}

function validClasses(value: unknown): boolean {
  if (!isRecord(value) || !hasOnlyKeys(value, ['schema', 'basis', 'active', 'rows'])
    || value.schema !== 1 || !validSelectionBasis(value.basis) || typeof value.active !== 'boolean'
    || !Array.isArray(value.rows)) return false;
  const keys: string[] = [];
  for (const row of value.rows) {
    if (!isRecord(row) || typeof row.key !== 'string' || row.key.length === 0
      || ![50, 60, 70, 80, 90, 95].includes(row.rrm as number)) return false;
    keys.push(row.key);
    if (row.rrm === 90 || row.rrm === 95) {
      if (!hasOnlyKeys(row, ['key', 'rrm', 'tilt30', 'disp45'])
        || typeof row.tilt30 !== 'boolean' || typeof row.disp45 !== 'boolean') return false;
    } else if (!hasOnlyKeys(row, ['key', 'rrm'])) return false;
  }
  return isUnique(keys) && sameValue([...keys].sort(), [...value.basis].sort());
}

/** Strict read validator for the stage-5 energy-threshold decision (S24). */
export function energyThresholdEditsIsWellFormed(value: unknown): boolean {
  if (value === null) return true;
  if (!isRecord(value) || !hasOnlyKeys(value, ['schema', 'basis', 'ranges']) || value.schema !== 1
    || !validRnaBasis(value.basis) || !Array.isArray(value.ranges)) return false;

  const rangeIds: string[] = [];
  for (const entry of value.ranges) {
    if (!isRecord(entry) || !hasOnlyKeys(entry, ['row', 'values'])
      || !isRangeAddress(entry.row) || !isRecord(entry.values)) return false;
    const row = entry.row;
    rangeIds.push(`${row.chain}\u0000${row.start}\u0000${row.end}`);
    for (const [fragment, threshold] of Object.entries(entry.values)) {
      if (!validFragmentKey(fragment)) return false;
      if (threshold === null) continue;
      if (typeof threshold !== 'string' || !/^[+-]?(?:0|[1-9]\d*)(?:\.\d{1,2})?$/.test(threshold)
        || !Number.isFinite(Number(threshold))) return false;
    }
  }
  return isUnique(rangeIds);
}

/** Strict read validator for the edit markers consumed by the stage graph. */
export function markerIsWellFormed(kind: MarkerKind, value: unknown): boolean {
  if (kind !== 'best_fit_clip' && value === null) return true;
  switch (kind) {
    case 'ranges': return isRangeMarker(value);
    case 'stacking_selection': return isSelectionMarker(value);
    case 'threshold_edits': return validThresholdMarker(value);
    case 'best_fit_clip': return validClip(value);
    case 'interaction_matrix': return validMatrix(value);
    case 'gaussian_classes': return validClasses(value);
  }
}

function currentRnaRange(snapshot: ProjectSnapshot): RangeAddress | null {
  if (!isRangeMarker(snapshot.ranges) || !Array.isArray(snapshot.ranges.rows)) return null;
  const row = snapshot.ranges.rows[0];
  return isRangeAddress(row) ? row : null;
}

function matrixRowsMatchClip(matrix: AnyRecord, snapshot: ProjectSnapshot): boolean {
  if (!Array.isArray(matrix.rows) || !isRecord(snapshot.threshold_edits)
    || !Array.isArray(snapshot.threshold_edits.bases)
    || !markerIsWellFormed('best_fit_clip', snapshot.best_fit_clip)) return false;
  const clip = snapshot.best_fit_clip as AnyRecord & { clip1: number; clip2: number };
  const rna = currentRnaRange(snapshot);
  if (!rna || !validRnaBasis(clip.basis)
    || !sameValue(clip.basis, { pdbcode: snapshot.pdbcode, rna })) return false;
  const bases = snapshot.threshold_edits.bases as AnyRecord[];
  const allKeys = bases.map(base => `frag-${base.fragment}-${base.sequence}`);
  if (new Set(allKeys).size !== allKeys.length || clip.clip1 + clip.clip2 >= allKeys.length) return false;
  const expectedRows = allKeys.slice(clip.clip1, allKeys.length - clip.clip2);
  const expectedRanges = isRangeMarker(snapshot.ranges)
    ? snapshot.ranges.rows.slice(1).map((_, index) => `range-${index + 2}`)
    : [];
  return sameValue(matrix.rows, expectedRows) && sameValue(matrix.ranges, expectedRanges)
    && Array.isArray(matrix.grid) && matrix.grid.length === expectedRows.length
    && matrix.grid.every((row: unknown) => Array.isArray(row) && row.length === expectedRanges.length);
}

/** Check one marker against only the persisted upstream decisions it records. */
export function markerMatchesUpstream(kind: MarkerKind, value: unknown, snapshot: ProjectSnapshot): boolean {
  if (!markerIsWellFormed(kind, value) || !isRecord(value)) return false;
  if (kind !== 'best_fit_clip' && value.active !== true) return false;
  if (kind === 'ranges') {
    return sameValue(value.basis, { pdbcode: snapshot.pdbcode, biounit: snapshot.biounit });
  }
  if (!isRangeMarker(snapshot.ranges) || snapshot.ranges.active !== true) return false;
  if (kind === 'stacking_selection') {
    const expected = {
      pdbcode: snapshot.pdbcode,
      biounit: snapshot.biounit,
      rows: snapshot.ranges.rows,
    };
    return sameValue(value.basis, expected);
  }
  const rna = currentRnaRange(snapshot);
  if (!rna) return false;
  if (kind === 'threshold_edits') {
    return sameValue(value.basis, { pdbcode: snapshot.pdbcode, rna });
  }
  if (kind === 'best_fit_clip') {
    return validRnaBasis(value.basis)
      && sameValue(value.basis, { pdbcode: snapshot.pdbcode, rna });
  }
  if (!isSelectionMarker(snapshot.stacking_selection)
    || snapshot.stacking_selection.active !== true) return false;
  if (kind === 'interaction_matrix') {
    return selectedBasisMatches(value.basis, snapshot.stacking_selection)
      && markerMatchesUpstream('threshold_edits', snapshot.threshold_edits, snapshot)
      && matrixRowsMatchClip(value, snapshot);
  }
  if (kind === 'gaussian_classes') {
    const keys = selectedKeys(snapshot.stacking_selection);
    return keys !== null && selectedBasisMatches(value.basis, snapshot.stacking_selection)
      && Array.isArray(value.rows)
      && sameValue(value.rows.map((row: AnyRecord) => row.key).sort(), keys);
  }
  return false;
}

/** Highest stage supported by a well-formed, active prefix of decision markers. */
export function deriveProjectStage(snapshot: ProjectSnapshot): ProjectStage {
  let stage: ProjectStage = 1;
  if (!markerMatchesUpstream('ranges', snapshot.ranges, snapshot)) return stage;
  stage = 2;
  if (!markerMatchesUpstream('stacking_selection', snapshot.stacking_selection, snapshot)) return stage;
  stage = 3;
  if (!markerMatchesUpstream('threshold_edits', snapshot.threshold_edits, snapshot)) return stage;
  stage = 4;
  if (!markerMatchesUpstream('interaction_matrix', snapshot.interaction_matrix, snapshot)) return stage;
  stage = 5;
  if (!markerMatchesUpstream('gaussian_classes', snapshot.gaussian_classes, snapshot)) return stage;
  return 6;
}

/** Values, not hashes: a result is fresh exactly when its echoed inputs match. */
export function echoIsFresh(echo: unknown, currentEcho: unknown): boolean {
  return echo !== null && echo !== undefined && sameValue(echo, currentEcho);
}

export interface DeferredWriteGuard {
  stageAtIssue: ProjectStage;
  currentStage: ProjectStage;
  live: boolean;
  revisionAtIssue: number;
  currentRevision: number;
  connectionAtIssue?: number;
  currentConnection?: number;
}

/** Deferred edits belong to their original stage and cell revision. */
export function deferredWriteIsCurrent(guard: DeferredWriteGuard): boolean {
  const connectionMatches = guard.connectionAtIssue === undefined && guard.currentConnection === undefined
    || Number.isSafeInteger(guard.connectionAtIssue)
      && guard.connectionAtIssue === guard.currentConnection;
  return guard.live && guard.stageAtIssue === guard.currentStage
    && Number.isSafeInteger(guard.revisionAtIssue)
    && guard.revisionAtIssue === guard.currentRevision
    && connectionMatches;
}

export interface StateGraphNode {
  path: Array<string | number>;
  state?: string;
  exception?: string | null;
}

export interface VisibleTransformerError {
  name: string;
  state: string;
  message: string;
}

export interface TransformerErrorOptions {
  stage: ProjectStage;
  hasInitialStateGraph: boolean;
  pendingInputs: ReadonlySet<string>;
  nodes: readonly StateGraphNode[];
}

/**
 * Return in-scope failures only.  The state graph never changes the stage; it
 * only supplies banners once the first graph has arrived.
 */
export function visibleTransformerErrors(options: TransformerErrorOptions): VisibleTransformerError[] {
  if (!options.hasInitialStateGraph) return [];
  const errors: VisibleTransformerError[] = [];
  for (const node of options.nodes) {
    const name = node.path.map(String).find(part => Object.hasOwn(TRANSFORMER_STAGE, part)
      || /^energetics_range_[1-9]\d*$/.test(part));
    if (!name) continue;
    const templateName = name.startsWith('energetics_range_') ? 'energetics_range' : name;
    if (TRANSFORMER_STAGE[templateName] > options.stage) continue;
    const failed = node.state === 'failed' || (node.state === undefined && !!node.exception);
    const blocked = node.state === 'blocked';
    if (!failed && !blocked) continue;
    const inputs = TRANSFORMER_INPUTS[templateName] ?? new Set<string>();
    if ([...inputs].some(input => options.pendingInputs.has(input))) continue;
    errors.push({
      name,
      state: blocked ? 'blocked' : 'failed',
      message: blocked ? 'blocked' : node.exception || 'Transformation failed.',
    });
  }
  return errors;
}

// Keep the list visible to runtime callers that need to validate names without
// accepting aliases from historical plans (T17 and S16 were removed).
export const ACTIVE_TRANSFORMER_NAMES = Object.freeze(Object.keys(TRANSFORMER_STAGE));
