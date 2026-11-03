import { markerIsWellFormed } from './project-state';

export type GaussianRrm = 50 | 60 | 70 | 80 | 90 | 95;
export type GaussianClassRow =
  | { key: string; rrm: 50 | 60 | 70 | 80 }
  | { key: string; rrm: 90 | 95; tilt30: boolean; disp45: boolean };

export interface GaussianClassMarker {
  schema: 1;
  basis: string[];
  active: boolean;
  rows: GaussianClassRow[];
}

export interface AutomaticGaussianPrefilters {
  tilt30: boolean;
  disp45: boolean;
}

type RecordValue = Record<string, any>;
const RRMS: readonly GaussianRrm[] = [50, 60, 70, 80, 90, 95];

function isRecord(value: unknown): value is RecordValue {
  return value !== null && typeof value === 'object' && !Array.isArray(value);
}

export function selectedGaussianKeys(selection: unknown): string[] {
  if (!isRecord(selection) || selection.active !== true || !isRecord(selection.selected)
    || !Object.values(selection.selected).every(value => typeof value === 'boolean')) return [];
  return Object.keys(selection.selected).filter(key => selection.selected[key] === true).sort();
}

export function editableGaussianRow(value: unknown): GaussianClassRow | null {
  if (!isRecord(value) || typeof value.key !== 'string' || !RRMS.includes(value.rrm)) return null;
  if (value.rrm === 90 || value.rrm === 95) {
    if (typeof value.tilt30 !== 'boolean' || typeof value.disp45 !== 'boolean') return null;
    return { key: value.key, rrm: value.rrm, tilt30: value.tilt30, disp45: value.disp45 };
  }
  return { key: value.key, rrm: value.rrm };
}

export function editableGaussianRows(values: unknown): GaussianClassRow[] | null {
  if (!Array.isArray(values)) return null;
  const rows: GaussianClassRow[] = [];
  for (const value of values) {
    const row = editableGaussianRow(value);
    if (!row) return null;
    rows.push(row);
  }
  return rows;
}

/**
 * Reconcile the S18 editor with the authoritative share. While a local write
 * is pending, the shared marker may echo an older in-flight PUT; keep the
 * current local draft so subsequent edits build on the latest user value.
 */
export function reconcileGaussianClassDraft(
  shared: unknown,
  draft: GaussianClassMarker | null,
  writePending: boolean,
): GaussianClassMarker | null {
  if (writePending || !markerIsWellFormed('gaussian_classes', shared)
    || !isRecord(shared) || shared.active !== true || !Array.isArray(shared.basis)) return draft;
  const rows = editableGaussianRows(shared.rows);
  return rows ? { schema: 1, basis: [...shared.basis], active: true, rows } : draft;
}

export function gaussianFilterName(row: GaussianClassRow): string {
  let name = `rrm${row.rrm}`;
  if (row.rrm === 90 || row.rrm === 95) {
    if (row.tilt30) name += '-tilt30';
    if (row.disp45) name += '-disp45';
  }
  return name;
}

/** Expose base classes from the backend's delivered legal filter names. */
export function gaussianRrmOptions(filterNames: readonly string[]): GaussianRrm[] {
  const values = new Set<GaussianRrm>();
  for (const name of filterNames) {
    const match = /^rrm(50|60|70|80|90|95)$/.exec(name);
    if (match) values.add(Number(match[1]) as GaussianRrm);
  }
  return [...values].sort((left, right) => left - right);
}

export function computedGaussianRows(scores: unknown, selectedKeys: readonly string[]): GaussianClassRow[] | null {
  if (!isRecord(scores) || scores.schema !== 1 || !Array.isArray(scores.classes)
    || !Array.isArray(scores.filter_names)) return null;
  const names = new Set(scores.filter_names.filter((name: unknown): name is string => typeof name === 'string'));
  const byKey = new Map<string, unknown>();
  for (const item of scores.classes) {
    if (!isRecord(item) || typeof item.key !== 'string' || byKey.has(item.key)) return null;
    byKey.set(item.key, item);
  }
  if (byKey.size !== selectedKeys.length || selectedKeys.some(key => !byKey.has(key))) return null;
  const rows: GaussianClassRow[] = [];
  for (const key of selectedKeys) {
    const computed = byKey.get(key);
    let row = editableGaussianRow(computed);
    if (row && (row.rrm === 90 || row.rrm === 95)) {
      const automatic = automaticGaussianPrefilters(computed);
      if (!automatic) return null;
      row = { key: row.key, rrm: row.rrm, ...automatic };
    }
    if (!row || !names.has(gaussianFilterName(row))) return null;
    rows.push(row);
  }
  return rows;
}

export function classMarkerForContinue(
  retained: unknown,
  scores: unknown,
  selectedKeys: readonly string[],
): GaussianClassMarker | null {
  if (markerIsWellFormed('gaussian_classes', retained) && isRecord(retained)
    && JSON.stringify([...retained.basis].sort()) === JSON.stringify([...selectedKeys].sort())) {
    const rows = editableGaussianRows(retained.rows);
    if (rows) return { schema: 1, basis: [...selectedKeys], active: true, rows };
  }
  const rows = computedGaussianRows(scores, selectedKeys);
  if (!rows) return null;
  return { schema: 1, basis: [...selectedKeys], active: true, rows };
}

/** Restore computed values for one row without disturbing any other edits. */
export function revertGaussianClassRow(
  marker: unknown,
  key: string,
  scores: unknown,
): GaussianClassMarker | null {
  if (!markerIsWellFormed('gaussian_classes', marker) || !isRecord(marker)
    || marker.active !== true || !Array.isArray(marker.basis) || !Array.isArray(marker.rows)) return null;
  const computed = computedGaussianRows(scores, marker.basis);
  const row = computed?.find(item => item.key === key);
  if (!row || !marker.rows.some((item: unknown) => isRecord(item) && item.key === key)) return null;
  return {
    schema: 1,
    basis: [...marker.basis],
    active: true,
    rows: marker.rows.map((item: GaussianClassRow) => item.key === key ? row : item),
  };
}

/** Produce the confirmed full-table copy; a cancelled action has no marker to write. */
export function copyComputedGaussianClasses(
  marker: unknown,
  scores: unknown,
  confirmed: boolean,
): GaussianClassMarker | null {
  if (!confirmed || !markerIsWellFormed('gaussian_classes', marker) || !isRecord(marker)
    || !Array.isArray(marker.basis)) return null;
  const rows = computedGaussianRows(scores, marker.basis);
  if (!rows) return null;
  return { schema: 1, basis: [...marker.basis], active: true, rows };
}

export function switchGaussianClass(
  current: GaussianClassRow,
  nextRrm: GaussianRrm,
  automatic: AutomaticGaussianPrefilters,
): GaussianClassRow {
  if (nextRrm < 90) return { key: current.key, rrm: nextRrm as 50 | 60 | 70 | 80 };
  if (current.rrm === nextRrm) return current;
  return {
    key: current.key,
    rrm: nextRrm as 90 | 95,
    tilt30: automatic.tilt30,
    disp45: automatic.disp45,
  };
}

export function automaticGaussianPrefilters(value: unknown): AutomaticGaussianPrefilters | null {
  if (!isRecord(value)) return null;
  const automatic = isRecord(value.auto_prefilters) ? value.auto_prefilters : value;
  return typeof automatic.tilt30 === 'boolean' && typeof automatic.disp45 === 'boolean'
    ? { tilt30: automatic.tilt30, disp45: automatic.disp45 } : null;
}

export function setGaussianFlag(
  row: GaussianClassRow,
  flag: 'tilt30' | 'disp45',
  value: boolean,
): GaussianClassRow {
  if (row.rrm < 90) return row;
  return { ...row, [flag]: value } as GaussianClassRow;
}

export function gaussianRowsDiffer(left: unknown, right: unknown): boolean {
  const a = editableGaussianRow(left);
  const b = editableGaussianRow(right);
  return !a || !b || JSON.stringify(a) !== JSON.stringify(b);
}

export function gaussianClassDiffCount(marker: unknown, computed: unknown): number {
  if (!isRecord(marker) || !Array.isArray(marker.rows)) return 0;
  const rows = computedGaussianRows(computed, marker.basis);
  if (!rows) return 0;
  const byKey = new Map(rows.map(row => [row.key, row]));
  return marker.rows.reduce((count: number, row: unknown) => {
    if (!isRecord(row) || typeof row.key !== 'string') return count;
    return count + (gaussianRowsDiffer(row, byKey.get(row.key)) ? 1 : 0);
  }, 0);
}

export function gaussianCandidateRows(scores: unknown, selectedKeys: readonly string[]): RecordValue[] {
  if (!isRecord(scores) || !Array.isArray(scores.candidate_scores)) return [];
  const selected = new Set(selectedKeys);
  return scores.candidate_scores.filter((row: unknown): row is RecordValue =>
    isRecord(row) && selected.has(row.key) && row.kept === true);
}

export function gaussianReferenceRows(scores: unknown, selectedKeys: readonly string[]): RecordValue[] {
  if (!isRecord(scores) || !Array.isArray(scores.reference_scores)) return [];
  const selected = new Set(selectedKeys);
  return scores.reference_scores.filter((row: unknown): row is RecordValue =>
    isRecord(row) && selected.has(row.key));
}

export function referenceAnchorState(value: unknown): 'green' | 'amber' | 'red' | 'blocked' | 'missing' | 'error' {
  if (!isRecord(value) || typeof value.status !== 'string') return 'error';
  if (['green', 'amber', 'red', 'blocked', 'missing'].includes(value.status)) return value.status as any;
  return 'error';
}
