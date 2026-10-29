import { definedRangeRows } from './range-definition';
import type { ChainInfo, RangeRow, RangeSpec, Stage, Stacking, StackingViewerInput } from './types';

export interface StackingTableState {
  stackings: Stacking[];
  stackingSelection: boolean[];
}

export interface Stage2RangeObservation {
  active: boolean;
  signature: string | null;
}

export interface Stage2RangeObservationResult {
  observation: Stage2RangeObservation;
  changed: boolean;
}

const initialRangeObservation: Stage2RangeObservation = { active: false, signature: null };

function isRecord(value: unknown): value is Record<string, unknown> {
  return typeof value === 'object' && value !== null && !Array.isArray(value);
}

function isFiniteNumber(value: unknown): value is number {
  return typeof value === 'number' && Number.isFinite(value);
}

function isInteger(value: unknown): value is number {
  return typeof value === 'number' && Number.isSafeInteger(value);
}

function normalizeCentroid(value: unknown): [number, number, number] | null {
  if (!Array.isArray(value) || value.length !== 3 || !value.every(isFiniteNumber)) return null;
  return [value[0], value[1], value[2]];
}

function normalizeStacking(value: unknown): Stacking | null {
  if (!isRecord(value)) return null;
  const {
    rna_resname,
    rna_resid,
    protein_range_pdb_id,
    protein_resid,
    distance,
    tilt_angle,
    slip_angle,
    rna_centroid,
    protein_centroid,
  } = value;
  const rnaCentroid = normalizeCentroid(rna_centroid);
  const proteinCentroid = normalizeCentroid(protein_centroid);
  if (typeof rna_resname !== 'string' || rna_resname.trim() === '') return null;
  if (!isInteger(rna_resid) || !isInteger(protein_resid)) return null;
  if (!isInteger(protein_range_pdb_id) || protein_range_pdb_id < 2) return null;
  if (!isFiniteNumber(distance) || distance < 0) return null;
  if (!isFiniteNumber(tilt_angle) || tilt_angle < 0 || tilt_angle > 90) return null;
  if (!isFiniteNumber(slip_angle) || slip_angle < 0 || slip_angle > 90) return null;
  if (!rnaCentroid || !proteinCentroid) return null;
  return {
    rna_resname: rna_resname.trim(),
    rna_resid,
    protein_range_pdb_id,
    protein_resid,
    distance,
    tilt_angle,
    slip_angle,
    rna_centroid: rnaCentroid,
    protein_centroid: proteinCentroid,
  };
}

/**
 * Parse and validate the JSON-compatible backend boundary. One malformed row
 * rejects the complete table so row and selection indices can never drift.
 */
export function normalizeStackings(value: unknown): Stacking[] | null {
  let parsed = value;
  if (typeof parsed === 'string') {
    try {
      parsed = JSON.parse(parsed) as unknown;
    } catch {
      return null;
    }
  }
  if (!Array.isArray(parsed)) return null;

  const normalized: Stacking[] = [];
  const seen = new Set<string>();
  for (const value of parsed) {
    const stacking = normalizeStacking(value);
    if (!stacking) return null;
    const key = stackingKey(stacking);
    if (seen.has(key)) return null;
    seen.add(key);
    normalized.push(stacking);
  }
  return normalized.sort(compareStackings);
}

/** Row identity follows the participant tuple, not measured geometry. */
export function stackingKey(stacking: Stacking): string {
  return JSON.stringify([
    stacking.rna_resname,
    stacking.rna_resid,
    stacking.protein_range_pdb_id,
    stacking.protein_resid,
  ]);
}

function compareStackings(left: Stacking, right: Stacking): number {
  return left.protein_range_pdb_id - right.protein_range_pdb_id
    || left.protein_resid - right.protein_resid
    || left.rna_resid - right.rna_resid
    || (left.rna_resname < right.rna_resname ? -1 : left.rna_resname > right.rna_resname ? 1 : 0);
}

function copyStacking(stacking: Stacking): Stacking {
  return {
    ...stacking,
    rna_centroid: [stacking.rna_centroid[0], stacking.rna_centroid[1], stacking.rna_centroid[2]],
    protein_centroid: [stacking.protein_centroid[0], stacking.protein_centroid[1], stacking.protein_centroid[2]],
  };
}

/** Merge a backend table while keeping choices for participants that remain. */
export function mergeStackings(
  previousStackings: readonly Stacking[],
  previousSelection: readonly boolean[],
  nextStackings: readonly Stacking[],
): StackingTableState {
  const previousByKey = new Map<string, boolean>();
  if (previousSelection.length === previousStackings.length) {
    previousStackings.forEach((stacking, index) => {
      previousByKey.set(stackingKey(stacking), previousSelection[index]);
    });
  }

  const stackings = nextStackings.map(copyStacking);
  const stackingSelection = stackings.map(stacking => {
    const previous = previousByKey.get(stackingKey(stacking));
    return previous ?? (stacking.distance <= 5.0 && stacking.slip_angle < 60.0);
  });
  return { stackings, stackingSelection };
}

/** Display-only numeric formatting; detector values remain unrounded. */
export function formatStackingNumber(value: number, fractionDigits = 2): string {
  return new Intl.NumberFormat('en-US', {
    minimumFractionDigits: fractionDigits,
    maximumFractionDigits: fractionDigits,
  }).format(value);
}

/**
 * Hash only range content that affects the Seamless detector. The local row ID
 * is UI identity and does not alter the extracted range PDBs.
 */
export function rangeRowsContentSignature(rows: readonly RangeRow[]): string {
  return JSON.stringify(rows.map(row => [row.alwaysVisible, row.chain, row.start, row.end]));
}

/** Treat the first Stage 2 snapshot as initialization; later edits invalidate. */
export function observeStage2RangeRows(
  previous: Stage2RangeObservation,
  stage: Stage,
  rows: readonly RangeRow[],
): Stage2RangeObservationResult {
  if (stage !== 2) {
    return { observation: initialRangeObservation, changed: false };
  }
  const signature = rangeRowsContentSignature(rows);
  if (!previous.active) {
    return { observation: { active: true, signature }, changed: false };
  }
  return {
    observation: { active: true, signature },
    changed: signature !== previous.signature,
  };
}

function participantRange(
  row: RangeRow | undefined,
  resid: number,
  chains: readonly ChainInfo[],
): RangeSpec | null {
  if (!row) return null;
  const chain = chains.find(item => item.name === row.chain);
  const rawStart = row.start.trim();
  const rawEnd = row.end.trim();
  const start = Number(rawStart);
  const end = Number(rawEnd);
  if (!chain || !rawStart || !rawEnd || !Number.isSafeInteger(start) || !Number.isSafeInteger(end)) return null;
  if (resid < Math.min(start, end) || resid > Math.max(start, end) || !chain.seqIds.has(resid)) return null;
  return {
    labelAsymId: chain.labelAsymId,
    instanceId: chain.instanceId,
    start: resid,
    end: resid,
  };
}

/** Resolve selected detector rows against the same ordered range list sent to Seamless. */
export function resolveSelectedStackings(
  stackings: readonly Stacking[],
  selection: readonly boolean[],
  rangeRows: readonly RangeRow[],
  chains: readonly ChainInfo[],
): StackingViewerInput[] {
  if (selection.length !== stackings.length) return [];
  const configuredRows = definedRangeRows(rangeRows);
  const rnaRow = configuredRows[0];
  if (!rnaRow) return [];

  const resolved: StackingViewerInput[] = [];
  stackings.forEach((stacking, index) => {
    if (selection[index] !== true) return;
    const proteinRow = configuredRows[stacking.protein_range_pdb_id - 1];
    const rna = participantRange(rnaRow, stacking.rna_resid, chains);
    const protein = participantRange(proteinRow, stacking.protein_resid, chains);
    if (!rna || !protein) return;
    resolved.push({
      rna,
      protein,
      rnaCentroid: [...stacking.rna_centroid],
      proteinCentroid: [...stacking.protein_centroid],
    });
  });
  return resolved;
}

/** Return distinct residue selections in first-seen order for one sticks representation. */
export function uniqueStackingParticipantRanges(
  stackings: readonly StackingViewerInput[],
): RangeSpec[] {
  const ranges = new Map<string, RangeSpec>();
  for (const stacking of stackings) {
    for (const range of [stacking.rna, stacking.protein]) {
      const key = JSON.stringify([range.labelAsymId, range.instanceId, range.start, range.end]);
      if (!ranges.has(key)) ranges.set(key, { ...range });
    }
  }
  return [...ranges.values()];
}
