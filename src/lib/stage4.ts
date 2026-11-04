import { stackingKey } from './stackings';
import type { Stacking } from './types';

export interface FragmentIdentity {
  chain: string;
  resid: number;
  insertion_code: string;
  base: 'A' | 'C' | 'G' | 'U';
}

export interface FragmentRow {
  number: number;
  key: string;
  sequence: string;
  residues: [FragmentIdentity, FragmentIdentity];
}

export interface FragmentInventory {
  schema: 1;
  echo: { pdbcode: string; rna: { chain: string; start: number; end: number } };
  residues: FragmentIdentity[];
  fragments: FragmentRow[];
  adjacent_pairs: Array<{ down: number; up: number; sequence: string }>;
}

export interface ThresholdPair {
  down: number;
  up: number;
  sequence: string;
  ovrmsd: string | null;
  crmsd: string | null;
}

export interface ThresholdMarker {
  schema: 1;
  basis: FragmentInventory['echo'];
  active: boolean;
  pairs: ThresholdPair[];
}

/** Copy a JSON decision marker while removing reactive proxies from its object graph. */
export function cloneThresholdMarker(marker: ThresholdMarker): ThresholdMarker {
  return JSON.parse(JSON.stringify(marker)) as ThresholdMarker;
}

export interface ClipMarker {
  schema: 1;
  basis: FragmentInventory['echo'] | null;
  clip1: number;
  clip2: number;
}

type UnknownRecord = Record<string, unknown>;

function isRecord(value: unknown): value is UnknownRecord {
  return value !== null && typeof value === 'object' && !Array.isArray(value);
}

function sameAddress(left: unknown, right: unknown): boolean {
  if (!isRecord(left) || !isRecord(right)) return false;
  return left.chain === right.chain && left.start === right.start && left.end === right.end;
}

/** Build an exact stage-3 selection map, retaining choices only for live keys. */
export function defaultStackingSelection(
  stackings: readonly Stacking[],
  retained: Readonly<Record<string, boolean>> = {},
): Record<string, boolean> {
  const selection: Record<string, boolean> = {};
  for (const stacking of stackings) {
    const key = stackingKey(stacking);
    const previous = retained[key];
    selection[key] = typeof previous === 'boolean'
      ? previous
      : stacking.distance <= 5 && stacking.slip_angle < 60;
  }
  return selection;
}

/** Make stage 4's blank, keyed table from T24 before the fit finishes. */
export function makeThresholdSkeleton(
  inventory: unknown,
  pdbCode?: string,
  active = true,
): ThresholdMarker | null {
  if (!isRecord(inventory) || inventory.schema !== 1 || !isRecord(inventory.echo)
    || !isRecord(inventory.echo.rna) || !Array.isArray(inventory.adjacent_pairs)
    || !Array.isArray(inventory.fragments) || inventory.fragments.length === 0) return null;
  const echo = inventory.echo as FragmentInventory['echo'];
  if (typeof echo.pdbcode !== 'string' || !echo.pdbcode
    || typeof echo.rna.chain !== 'string' || typeof echo.rna.start !== 'number'
    || !Number.isSafeInteger(echo.rna.start) || typeof echo.rna.end !== 'number'
    || !Number.isSafeInteger(echo.rna.end) || echo.rna.end < echo.rna.start
    || (pdbCode && echo.pdbcode !== pdbCode.toUpperCase())) return null;
  const pairs: ThresholdPair[] = [];
  for (const pair of inventory.adjacent_pairs) {
    if (!isRecord(pair) || typeof pair.down !== 'number' || !Number.isSafeInteger(pair.down)
      || typeof pair.up !== 'number' || !Number.isSafeInteger(pair.up)
      || pair.up !== Number(pair.down) + 1 || typeof pair.sequence !== 'string') return null;
    pairs.push({ down: Number(pair.down), up: Number(pair.up), sequence: pair.sequence, ovrmsd: null, crmsd: null });
  }
  return { schema: 1, basis: { pdbcode: echo.pdbcode, rna: { ...echo.rna } }, active, pairs };
}

/** RMSD values are finite positive decimals with no more than two decimals. */
export function validateThresholdValue(value: unknown): boolean {
  if (value === null) return true;
  if (typeof value !== 'string' || !/^(?:0|[1-9]\d*)(?:\.\d{1,2})?$/.test(value)) return false;
  const numeric = Number(value);
  return Number.isFinite(numeric) && numeric > 0;
}

/** Match every edit row to T24; null cells are valid while the user fills the table. */
export function validateThresholdEdits(edits: unknown, inventory: unknown): boolean {
  if (!isRecord(edits) || edits.schema !== 1 || typeof edits.active !== 'boolean'
    || Object.keys(edits).sort().join(',') !== 'active,basis,pairs,schema'
    || !isRecord(edits.basis) || !isRecord(inventory) || inventory.schema !== 1
    || !isRecord(inventory.echo) || !sameAddress(edits.basis.rna, inventory.echo.rna)
    || edits.basis.pdbcode !== inventory.echo.pdbcode
    || !Array.isArray(edits.pairs) || !Array.isArray(inventory.adjacent_pairs)) return false;
  if (edits.pairs.length !== inventory.adjacent_pairs.length) return false;
  for (let index = 0; index < edits.pairs.length; index += 1) {
    const row = edits.pairs[index];
    const source = inventory.adjacent_pairs[index];
    if (!isRecord(row) || !isRecord(source)
      || row.down !== source.down || row.up !== source.up || row.sequence !== source.sequence
      || !(row.ovrmsd === null || (typeof row.ovrmsd === 'string' && validateThresholdValue(row.ovrmsd)))
      || !(row.crmsd === null || (typeof row.crmsd === 'string' && validateThresholdValue(row.crmsd)))) return false;
  }
  return true;
}

type CandidateGroupSet = readonly (readonly number[])[];
type CandidateInput = readonly number[] | CandidateGroupSet;

function asCandidateGroups(value: CandidateInput): number[][] | null {
  if (!Array.isArray(value)) return null;
  if (value.every(candidate => typeof candidate === 'number')) {
    return value.map(candidate => [candidate as number]);
  }
  if (value.every(group => Array.isArray(group))) {
    return value.map(group => [...group as readonly number[]]);
  }
  return null;
}

function keepsEveryCandidateGroup(clip: ClipMarker, fragmentCount: number, value: CandidateInput): boolean {
  const groups = asCandidateGroups(value);
  if (groups === null) return false;
  return groups.every(group => group.length > 0
    && group.every(number => Number.isSafeInteger(number) && number >= 1 && number <= fragmentCount)
    && group.some(number => number > clip.clip1 && number <= fragmentCount - clip.clip2));
}

function validClip(
  clip: unknown,
  fragmentCount: number,
  candidates?: CandidateInput,
): clip is ClipMarker {
  if (!isRecord(clip) || clip.schema !== 1
    || !Number.isSafeInteger(clip.clip1) || Number(clip.clip1) < 0
    || !Number.isSafeInteger(clip.clip2) || Number(clip.clip2) < 0
    || Number(clip.clip1) + Number(clip.clip2) >= fragmentCount) return false;
  return candidates === undefined || keepsEveryCandidateGroup(clip as unknown as ClipMarker, fragmentCount, candidates);
}

/** Validate the stage-4 extent and require one kept candidate for every selected stacking. */
export function validClipForCandidates(
  clip: unknown,
  fragmentCount: number,
  candidates: CandidateInput,
): clip is ClipMarker {
  return validClip(clip, fragmentCount, candidates);
}

/** Return one candidate-fragment group per selected stacking key. */
export function selectedCandidateGroups(
  inventory: unknown,
  selectedKeys: readonly string[],
  stackings: unknown,
): number[][] {
  if (!isRecord(inventory) || !Array.isArray(inventory.fragments)) return [];
  const rows = Array.isArray(stackings) ? stackings
    : isRecord(stackings) && Array.isArray(stackings.rows) ? stackings.rows : null;
  if (!rows) return [];
  const candidatesByKey = new Map(selectedKeys.map(key => [key, new Set<number>()]));
  for (const row of rows) {
    if (!isRecord(row)) continue;
    const key = stackingKey(row as unknown as Stacking);
    const candidates = candidatesByKey.get(key);
    if (!candidates) continue;
    const resid = row.rna_resid;
    let base = typeof row.rna_resname === 'string' ? row.rna_resname.trim().toUpperCase() : '';
    if (base.length === 2 && (base[0] === 'R' || base[0] === 'D')) base = base[1];
    if (!Number.isSafeInteger(resid) || !/^[ACGU]$/.test(base)) continue;
    for (const fragment of inventory.fragments) {
      if (!isRecord(fragment) || !Array.isArray(fragment.residues)) continue;
      if (fragment.residues.some(identity => isRecord(identity)
        && identity.resid === resid && identity.base === base)
        && typeof fragment.number === 'number' && Number.isSafeInteger(fragment.number)) {
        candidates.add(Number(fragment.number));
      }
    }
  }
  return [...candidatesByKey.values()].map(group => [...group].sort((left, right) => left - right));
}

export function selectedCandidateFragments(
  inventory: unknown,
  selectedKeys: readonly string[],
  stackings: unknown,
): number[] {
  return [...new Set(selectedCandidateGroups(inventory, selectedKeys, stackings).flat())]
    .sort((left, right) => left - right);
}

/** Apply the sole automatic stage-3→4 correction: clamp a clip off anchors. */
export function clipForSelectedAnchors(
  clip: ClipMarker,
  inventory: unknown,
  selectedKeys: readonly string[],
  stackings: unknown,
): ClipMarker {
  if (!isRecord(inventory) || !Array.isArray(inventory.fragments)) return clip;
  const groups = selectedCandidateGroups(inventory, selectedKeys, stackings);
  if (groups.length === 0 || groups.some(group => group.length === 0)) return clip;
  const count = inventory.fragments.length;
  const clip1Limit = Math.min(...groups.map(group => Math.max(...group) - 1));
  const clip2Limit = Math.min(...groups.map(group => count - Math.min(...group)));
  return {
    ...clip,
    clip1: Math.min(clip.clip1, clip1Limit),
    clip2: Math.min(clip.clip2, clip2Limit),
  };
}

/** True when the fresh range, selection, fragment list and viewer can enter stage 4. */
export function canEnterStage4(options: {
  selection: unknown;
  stackings: unknown;
  stackingsFresh: boolean;
  fragmentsFresh: boolean;
  rangePdbsFresh: boolean;
}): boolean {
  const { selection, stackings, stackingsFresh, fragmentsFresh, rangePdbsFresh } = options;
  if (!stackingsFresh || !fragmentsFresh || !rangePdbsFresh || !isRecord(selection)
    || selection.schema !== 1 || selection.active !== true || !isRecord(selection.selected)) return false;
  const rows = Array.isArray(stackings) ? stackings
    : isRecord(stackings) && Array.isArray(stackings.rows) ? stackings.rows : null;
  if (!rows || rows.length === 0 || Object.values(selection.selected).some(value => typeof value !== 'boolean')) return false;
  const keys = rows.map(row => isRecord(row) ? stackingKey(row as unknown as Stacking) : '');
  if (keys.some(key => key === '') || new Set(keys).size !== keys.length) return false;
  const selectedKeys = Object.keys(selection.selected);
  return selectedKeys.length === keys.length && keys.every(key => Object.hasOwn(selection.selected as object, key))
    && Object.values(selection.selected).some(value => value === true);
}

/** Gate 4→5 on typed values, valid clip, T10 success, and fresh energetics. */
export function canContinueStage4(options: {
  edits: unknown;
  fragmentCount: number;
  clip: unknown;
  selectedCandidates: CandidateInput;
  thresholdsFresh: boolean;
  energeticsFresh: boolean;
  thresholdValidatorFailed: boolean;
}): boolean {
  const { edits, fragmentCount, clip, selectedCandidates, thresholdsFresh, energeticsFresh, thresholdValidatorFailed } = options;
  if (!thresholdsFresh || !energeticsFresh || thresholdValidatorFailed
    || !validClip(clip, fragmentCount, selectedCandidates)
    || !isRecord(edits) || edits.active !== true || !Array.isArray(edits.pairs)) return false;
  const currentClip = clip as ClipMarker;
  if (isRecord(edits.basis) && currentClip.basis && !sameAddress(edits.basis.rna, currentClip.basis.rna)) return false;
  const kept = new Set(Array.from({ length: fragmentCount - currentClip.clip1 - currentClip.clip2 },
    (_, index) => index + currentClip.clip1 + 1));
  for (const row of edits.pairs) {
    if (!isRecord(row) || !Number.isSafeInteger(row.down) || !Number.isSafeInteger(row.up)) return false;
    if (kept.has(Number(row.down)) && kept.has(Number(row.up))
      && (typeof row.ovrmsd !== 'string' || typeof row.crmsd !== 'string'
        || !validateThresholdValue(row.ovrmsd) || !validateThresholdValue(row.crmsd))) return false;
  }
  return true;
}

/** Copy fresh computed values only after the caller's overwrite confirmation. */
export function copyComputedThresholds(
  edits: ThresholdMarker,
  computed: unknown,
  confirmed: boolean,
): ThresholdMarker {
  if (!confirmed || !isRecord(computed) || computed.schema !== 1
    || !Array.isArray(computed.pairs)) return edits;
  const computedPairs = new Map<string, UnknownRecord>();
  for (const value of computed.pairs) {
    if (isRecord(value) && Number.isSafeInteger(value.down) && Number.isSafeInteger(value.up)) {
      computedPairs.set(`${value.down}:${value.up}`, value);
    }
  }
  return {
    ...edits,
    pairs: edits.pairs.map(row => {
      const source = computedPairs.get(`${row.down}:${row.up}`);
      if (!source || typeof source.ovrmsd !== 'string' || typeof source.crmsd !== 'string'
        || !validateThresholdValue(source.ovrmsd) || !validateThresholdValue(source.crmsd)) return row;
      return { ...row, ovrmsd: source.ovrmsd, crmsd: source.crmsd };
    }),
  };
}

/** The stage-4 gate is independent of computed fit values; only Copy uses S9. */
function isThresholdPair(value: unknown): value is ThresholdPair {
  return isRecord(value)
    && Number.isSafeInteger(value.down) && Number.isSafeInteger(value.up)
    && typeof value.sequence === 'string'
    && (value.ovrmsd === null || typeof value.ovrmsd === 'string')
    && (value.crmsd === null || typeof value.crmsd === 'string');
}

function hasThresholdRows(value: unknown): value is Pick<ThresholdMarker, 'schema' | 'active' | 'pairs'> {
  return isRecord(value) && value.schema === 1 && typeof value.active === 'boolean'
    && Array.isArray(value.pairs) && value.pairs.every(isThresholdPair);
}

export function thresholdsCompleteForClip(
  edits: unknown,
  fragmentCount: number,
  clip: unknown,
): boolean {
  if (!hasThresholdRows(edits) || edits.active !== true
    || !Number.isSafeInteger(fragmentCount) || fragmentCount < 1
    || !validClip(clip, fragmentCount)) return false;
  const marker = edits;
  const currentClip = clip as ClipMarker;
  const kept = new Set(Array.from({ length: fragmentCount - currentClip.clip1 - currentClip.clip2 },
    (_, index) => index + currentClip.clip1 + 1));
  const filled = (value: unknown): boolean => typeof value === 'string' && validateThresholdValue(value);
  return marker.pairs.every(row => !kept.has(row.down) || !kept.has(row.up)
    || (filled(row.ovrmsd) && filled(row.crmsd)));
}
