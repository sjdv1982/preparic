import { energyThresholdEditsIsWellFormed, markerIsWellFormed, markerMatchesUpstream, type ProjectSnapshot } from './project-state';

export interface Stage5Range { chain: string; start: number; end: number }
export interface Stage5S14 {
  schema: 1;
  basis: string[];
  active: boolean;
  rows: string[];
  ranges: string[];
  grid: boolean[][];
}
export interface Stage5S24 {
  schema: 1;
  basis: { pdbcode: string; rna: Stage5Range };
  ranges: Array<{ row: Stage5Range; values: Record<string, string | null> }>;
}
export interface Stage5Anchor {
  resid: number;
  range: string;
  stackings: string[];
  locks: Array<{ fragment: string; nucleotide: number }>;
}

/** Copy a JSON decision marker while removing reactive proxies from its object graph. */
export function cloneStage5Decision<T extends Stage5S14 | Stage5S24>(marker: T): T {
  return JSON.parse(JSON.stringify(marker)) as T;
}

type RecordValue = Record<string, any>;

function isRecord(value: unknown): value is RecordValue {
  return value !== null && typeof value === 'object' && !Array.isArray(value);
}

function addressKey(row: Stage5Range): string {
  return `${row.chain}\u0000${row.start}\u0000${row.end}`;
}

function selectedKeys(selection: unknown): string[] | null {
  if (!isRecord(selection) || !isRecord(selection.selected)
    || !Object.values(selection.selected).every(value => typeof value === 'boolean')) return null;
  const selected = selection.selected;
  return Object.keys(selected).filter(key => selected[key] === true).sort();
}

export function stage5Rows(fragments: unknown, clip: unknown): string[] | null {
  if (!isRecord(fragments) || !Array.isArray(fragments.fragments) || !isRecord(clip)
    || !Number.isSafeInteger(clip.clip1) || !Number.isSafeInteger(clip.clip2)) return null;
  const rows = fragments.fragments.map((row: any) => row?.key);
  if (!rows.length || rows.some((key: unknown) => typeof key !== 'string')
    || new Set(rows).size !== rows.length || clip.clip1 < 0 || clip.clip2 < 0
    || clip.clip1 + clip.clip2 >= rows.length) return null;
  return rows.slice(clip.clip1, rows.length - clip.clip2);
}

export function stage5RangeIds(rangeRows: unknown): string[] | null {
  if (!Array.isArray(rangeRows) || rangeRows.length < 2) return null;
  return rangeRows.slice(1).map((_, index) => `range-${index + 2}`);
}

export function candidateAnchors(
  fragments: unknown,
  selection: unknown,
  stackings: unknown,
  matrixRows: readonly string[],
): Stage5Anchor[] | null {
  if (!isRecord(fragments) || !Array.isArray(fragments.fragments) || !isRecord(selection)
    || selection.active !== true || !isRecord(selection.selected)
    || !Object.values(selection.selected).every(value => typeof value === 'boolean')) return null;
  const stackingRows = Array.isArray(stackings) ? stackings
    : isRecord(stackings) && Array.isArray(stackings.rows) ? stackings.rows : null;
  if (!stackingRows) return null;
  const selected = new Set(Object.keys(selection.selected).filter(key => selection.selected[key] === true));
  const byResidue = new Map<number, Array<{ fragment: string; nucleotide: number }>>();
  for (const fragment of fragments.fragments) {
    if (!isRecord(fragment) || typeof fragment.key !== 'string' || !Array.isArray(fragment.residues)) return null;
    fragment.residues.forEach((residue: unknown, index: number) => {
      if (!isRecord(residue) || !Number.isSafeInteger(residue.resid)) return;
      const locks = byResidue.get(residue.resid as number) ?? [];
      locks.push({ fragment: fragment.key as string, nucleotide: index + 1 });
      byResidue.set(residue.resid as number, locks);
    });
  }
  const anchors = new Map<string, Stage5Anchor>();
  for (const raw of stackingRows) {
    if (!isRecord(raw) || typeof raw.rna_resname !== 'string' || !Number.isSafeInteger(raw.rna_resid)
      || !Number.isSafeInteger(raw.protein_range_pdb_id) || !Number.isSafeInteger(raw.protein_resid)) return null;
    const key = JSON.stringify([raw.rna_resname.trim(), raw.rna_resid,
      raw.protein_range_pdb_id, raw.protein_resid]);
    if (!selected.has(key)) continue;
    const range = `range-${raw.protein_range_pdb_id}`;
    const address = `${raw.rna_resid}\u0000${range}`;
    let anchor = anchors.get(address);
    if (!anchor) {
      anchor = { resid: raw.rna_resid as number, range, stackings: [], locks: [] };
      anchors.set(address, anchor);
    }
    if (!anchor.stackings.includes(key)) anchor.stackings.push(key);
    for (const lock of byResidue.get(raw.rna_resid as number) ?? []) {
      if (!anchor.locks.some(item => item.fragment === lock.fragment && item.nucleotide === lock.nucleotide)) {
        anchor.locks.push({ ...lock });
      }
    }
  }
  const rowSet = new Set(matrixRows);
  const result = [...anchors.values()].sort((left, right) => left.resid - right.resid
    || Number(left.range.slice(6)) - Number(right.range.slice(6)));
  for (const anchor of result) {
    anchor.locks = anchor.locks.filter(lock => rowSet.has(lock.fragment));
    if (anchor.locks.length === 0) return null;
  }
  return result;
}

export function makeInteractionMatrix(
  fragments: unknown,
  clip: unknown,
  rangeRows: unknown,
  selection: unknown,
  stackings: unknown,
): Stage5S14 | null {
  const rows = stage5Rows(fragments, clip);
  const proteinRanges = stage5RangeIds(rangeRows);
  const basis = selectedKeys(selection);
  const anchors = rows && basis ? candidateAnchors(fragments, selection, stackings, rows) : null;
  if (!rows || !proteinRanges || !basis || !anchors) return null;
  if (anchors.some(anchor => !proteinRanges.includes(anchor.range))) return null;
  const ranges = proteinRanges;
  return {
    schema: 1,
    basis,
    active: true,
    rows,
    ranges,
    grid: rows.map(() => ranges.map(() => false)),
  };
}

export function reconcileEnergyThresholdEdits(
  existing: unknown,
  fragments: unknown,
  rangeRows: unknown,
): Stage5S24 | null {
  if (!isRecord(fragments) || !isRecord(fragments.echo) || !Array.isArray(fragments.fragments)
    || !Array.isArray(rangeRows) || rangeRows.length < 2) return null;
  const basis = fragments.echo as Stage5S24['basis'];
  const keys = fragments.fragments.map((row: any) => row?.key);
  if (!keys.length || keys.some((key: unknown) => typeof key !== 'string')) return null;
  const previous = new Map<string, RecordValue>();
  if (isRecord(existing) && existing.schema === 1 && sameAddressBasis(existing.basis, basis)
    && Array.isArray(existing.ranges)) {
    for (const entry of existing.ranges) {
      if (isRecord(entry) && isRecord(entry.row) && isRecord(entry.values)) {
        previous.set(addressKey(entry.row as Stage5Range), entry.values);
      }
    }
  }
  const ranges = rangeRows.slice(1).map((row: any) => ({
    row: { chain: row.chain, start: row.start, end: row.end } as Stage5Range,
    values: Object.fromEntries(keys.map((key: string) => {
      const old = previous.get(addressKey(row as Stage5Range))?.[key];
      const accepted = old === null || (typeof old === 'string' && validEnergyEdit(old));
      return [key, accepted ? old : null];
    })),
  }));
  return { schema: 1, basis, ranges };
}

function sameAddressBasis(left: unknown, right: unknown): boolean {
  return isRecord(left) && isRecord(right) && left.pdbcode === right.pdbcode
    && isRecord(left.rna) && isRecord(right.rna)
    && addressKey(left.rna as Stage5Range) === addressKey(right.rna as Stage5Range);
}

export function validEnergyEdit(value: string): boolean {
  return /^[+-]?(?:0|[1-9]\d*)(?:\.\d{1,2})?$/.test(value)
    && Number.isFinite(Number(value));
}

export function validateStage5MatrixRead(
  value: unknown,
  expectedRows: readonly string[],
  expectedRanges: readonly string[],
): string | null {
  if (new Set(expectedRanges).size !== expectedRanges.length) return 'protein range order contains duplicates';
  if (!isRecord(value) || Object.keys(value).length !== 6
    || !['schema', 'basis', 'active', 'rows', 'ranges', 'grid'].every(key => Object.hasOwn(value, key))) {
    return 'schema check failed';
  }
  if (value.schema !== 1 || typeof value.active !== 'boolean') return 'schema or active-state check failed';
  if (!Array.isArray(value.basis) || value.basis.some(key => typeof key !== 'string')
    || new Set(value.basis).size !== value.basis.length) return 'selection basis check failed';
  if (!Array.isArray(value.rows) || value.rows.length !== expectedRows.length
    || value.rows.some((row, index) => row !== expectedRows[index])) return 'addressed rows or clip extent check failed';
  if (!Array.isArray(value.ranges) || value.ranges.length !== expectedRanges.length
    || value.ranges.some((range, index) => range !== expectedRanges[index])) return 'ordered protein ranges check failed';
  if (!Array.isArray(value.grid) || value.grid.length !== expectedRows.length) return 'row count check failed';
  for (const row of value.grid) {
    if (!Array.isArray(row) || row.length !== expectedRanges.length) return 'rectangular grid check failed';
    if (row.some(cell => typeof cell !== 'boolean')) return 'boolean cell check failed';
  }
  return null;
}

export function validateStage5InteractionMatrix(
  matrix: unknown,
  expectedRows: readonly string[],
  anchors: readonly Stage5Anchor[],
  energetics: unknown,
  expectedRanges?: readonly string[],
): string | null {
  const ranges = expectedRanges ?? (isRecord(matrix) && Array.isArray(matrix.ranges)
    ? matrix.ranges as string[] : [...new Set(anchors.map(anchor => anchor.range))]);
  const shapeError = validateStage5MatrixRead(matrix, expectedRows, ranges);
  if (shapeError) return shapeError;
  if (!isRecord(matrix) || !Array.isArray(matrix.grid)) return 'matrix structure is invalid';
  const rows = matrix.rows as string[];
  const matrixRanges = matrix.ranges as string[];
  const grid = matrix.grid as boolean[][];
  const rowIndex = new Map(rows.map((key, index) => [key, index]));
  const anchorsByRange = new Map<string, Stage5Anchor[]>();
  for (const anchor of anchors) {
    const column = matrixRanges.indexOf(anchor.range);
    if (column < 0) return `selected anchor ${anchor.resid} refers to an unknown protein range`;
    const lockPositions: number[] = [];
    for (const lock of anchor.locks) {
      const position = rowIndex.get(lock.fragment);
      if (position === undefined) return `selected anchor ${anchor.resid} has a clipped lock`;
      if (!grid[position][column]) return `selected anchor ${anchor.resid} is not locked in its column`;
      lockPositions.push(position);
    }
    if (!lockPositions.length) return `selected anchor ${anchor.resid} has no candidate locks`;
    const rangeAnchors = anchorsByRange.get(anchor.range) ?? [];
    rangeAnchors.push(anchor);
    anchorsByRange.set(anchor.range, rangeAnchors);
  }
  for (const [column, range] of matrixRanges.entries()) {
    const locks = (anchorsByRange.get(range) ?? []).flatMap(anchor => anchor.locks
      .map(lock => rowIndex.get(lock.fragment)).filter((position): position is number => position !== undefined));
    let index = 0;
    while (index < rows.length) {
      if (!grid[index][column]) { index += 1; continue; }
      const first = index;
      while (index + 1 < rows.length && grid[index + 1][column]) index += 1;
      const last = index;
      if (!locks.some(position => position >= first && position <= last)) {
        return `true run in ${range} has no selected anchor`;
      }
      for (let row = first; row <= last; row += 1) {
        if (!hasValidEnergies(energetics, range, rows[row])) {
          const lockOwner = (anchorsByRange.get(range) ?? []).find(anchor =>
            anchor.locks.some(lock => lock.fragment === rows[row]));
          return lockOwner
            ? `selected anchor ${lockOwner.resid} has a scoring error at ${rows[row]} / ${range}`
            : `scoring error at ${rows[row]} / ${range}`;
        }
      }
      index += 1;
    }
  }
  return null;
}

export function recomputeInteractionMatrix(
  matrix: unknown,
  anchors: readonly Stage5Anchor[],
  energetics: unknown,
): Stage5S14 | null {
  if (!isRecord(matrix) || !Array.isArray(matrix.rows) || !Array.isArray(matrix.ranges)
    || matrix.ranges.some((range: unknown) => typeof range !== 'string')
    || new Set(matrix.ranges).size !== matrix.ranges.length
    || anchors.some(anchor => !matrix.ranges.includes(anchor.range))) return null;
  const rows = matrix.rows as string[];
  const positions = new Map(rows.map((key, index) => [key, index]));
  const ranges = matrix.ranges as string[];
  const grid = rows.map(() => ranges.map(() => false));
  const locksByRange = new Map<string, Set<number>>();
  for (const anchor of anchors) {
    const locks = locksByRange.get(anchor.range) ?? new Set<number>();
    for (const lock of anchor.locks) {
      const position = positions.get(lock.fragment);
      if (position === undefined) return null;
      locks.add(position);
    }
    locksByRange.set(anchor.range, locks);
  }
  for (const anchor of anchors) {
    const column = ranges.indexOf(anchor.range);
    const locks = anchor.locks.map(lock => positions.get(lock.fragment));
    if (!locks.length || locks.some(position => position === undefined)) return null;
    for (const position of locks as number[]) {
      if (!hasValidEnergies(energetics, anchor.range, rows[position])) return null;
    }
    let first = Math.min(...locks as number[]);
    let last = Math.max(...locks as number[]);
    const otherLocks = new Set(locksByRange.get(anchor.range));
    for (const position of locks) if (position !== undefined) otherLocks.delete(position);
    while (first > 0 && !otherLocks.has(first - 1) && grows(energetics, anchor.range, rows[first - 1])) first -= 1;
    while (last + 1 < rows.length && !otherLocks.has(last + 1)
      && grows(energetics, anchor.range, rows[last + 1])) last += 1;
    for (let row = first; row <= last; row += 1) grid[row][column] = true;
  }
  return { ...matrix, grid } as Stage5S14;
}

export function toggleStage5MatrixCell(
  matrix: Stage5S14,
  row: number,
  column: number,
  anchors: readonly Stage5Anchor[],
  energetics: unknown,
): Stage5S14 | null {
  if (row < 0 || row >= matrix.rows.length || column < 0 || column >= matrix.ranges.length) return null;
  const next = { ...matrix, grid: matrix.grid.map(values => [...values]) };
  const value = !next.grid[row][column];
  next.grid[row][column] = value;
  return validateStage5InteractionMatrix(next, matrix.rows, anchors, energetics, matrix.ranges) === null ? next : null;
}

export interface UntilNextAnchorState { nextColumn: number | null; allowed: boolean; joined: boolean; forced: boolean }

function lockInterval(anchor: Stage5Anchor, rows: readonly string[]): { first: number; last: number } | null {
  const positions = anchor.locks.map(lock => rows.indexOf(lock.fragment));
  if (!positions.length || positions.some(position => position < 0)) return null;
  return { first: Math.min(...positions), last: Math.max(...positions) };
}

/** Describe the optional link to the next downstream anchor on the same range. */
export function untilNextAnchorState(
  matrix: Stage5S14,
  anchorIndex: number,
  anchors: readonly Stage5Anchor[],
): UntilNextAnchorState {
  const current = anchors[anchorIndex];
  if (!current) return { nextColumn: null, allowed: false, joined: false, forced: false };
  let nextColumn: number | null = null;
  for (let index = anchorIndex + 1; index < anchors.length; index += 1) {
    if (anchors[index].range === current.range) { nextColumn = index; break; }
  }
  if (nextColumn === null) return { nextColumn: null, allowed: false, joined: false, forced: false };
  const leftLocks = lockInterval(current, matrix.rows);
  const rightLocks = lockInterval(anchors[nextColumn], matrix.rows);
  if (!leftLocks || !rightLocks) return { nextColumn, allowed: false, joined: false, forced: false };
  const forced = rightLocks.first <= leftLocks.last + 1;
  const rangeColumn = matrix.ranges.indexOf(current.range);
  if (rangeColumn < 0) return { nextColumn, allowed: false, joined: false, forced: false };
  const first = Math.min(leftLocks.first, rightLocks.first);
  const last = Math.max(leftLocks.last, rightLocks.last);
  const joined = matrix.grid.slice(first, last + 1).every(row => row[rangeColumn] === true);
  return { nextColumn, allowed: true, joined, forced };
}

/** Join or split successive same-range anchors inside their shared range column. */
export function setStage5UntilNextAnchor(
  matrix: Stage5S14,
  anchorIndex: number,
  enabled: boolean,
  anchors: readonly Stage5Anchor[],
  energetics: unknown,
): Stage5S14 | null {
  const state = untilNextAnchorState(matrix, anchorIndex, anchors);
  if (!state.allowed || state.nextColumn === null || (!enabled && state.forced)) return null;
  if (enabled === state.joined) return matrix;
  const leftLocks = lockInterval(anchors[anchorIndex], matrix.rows);
  const rightLocks = lockInterval(anchors[state.nextColumn], matrix.rows);
  if (!leftLocks || !rightLocks) return null;
  const rangeColumn = matrix.ranges.indexOf(anchors[anchorIndex].range);
  if (rangeColumn < 0) return null;
  const next = { ...matrix, grid: matrix.grid.map(row => [...row]) };

  if (enabled) {
    const first = Math.min(leftLocks.first, rightLocks.first);
    const last = Math.max(leftLocks.last, rightLocks.last);
    for (let row = first; row <= last; row += 1) next.grid[row][rangeColumn] = true;
  } else {
    if (rightLocks.first <= leftLocks.last + 1) return null;
    for (let row = leftLocks.last + 1; row < rightLocks.first; row += 1) next.grid[row][rangeColumn] = false;
  }
  return validateStage5InteractionMatrix(next, matrix.rows, anchors, energetics, matrix.ranges) === null ? next : null;
}

function energyCell(energetics: unknown, range: string, fragment: string): any | null {
  if (!isRecord(energetics) || !Array.isArray(energetics.ranges)) return null;
  const rangeRow = energetics.ranges.find((item: any) => item?.range === range);
  return Array.isArray(rangeRow?.rows) ? rangeRow.rows.find((item: any) => item?.key === fragment) ?? null : null;
}

function hasValidEnergies(energetics: unknown, range: string, fragment: string): boolean {
  const row = energyCell(energetics, range, fragment);
  return !!row && ['reference', 'bestfit'].every(key => isRecord(row[key])
    && row[key].error == null && typeof row[key].energy === 'number' && Number.isFinite(row[key].energy));
}

function grows(energetics: unknown, range: string, fragment: string): boolean {
  const row = energyCell(energetics, range, fragment);
  if (!hasValidEnergies(energetics, range, fragment)) return false;
  return Math.abs(row.reference.energy) >= 1 || Math.abs(row.bestfit.energy) >= 1;
}

export function adaptInteractionMatrix(
  matrix: unknown,
  newRows: readonly string[],
  anchors: readonly RecordValue[],
  energetics: unknown,
): Stage5S14 | null {
  if (!isRecord(matrix) || !Array.isArray(matrix.rows) || !Array.isArray(matrix.ranges)
    || matrix.ranges.some((range: unknown) => typeof range !== 'string')
    || !Array.isArray(matrix.grid) || matrix.grid.length !== matrix.rows.length
    || matrix.grid.some((row: unknown) => !Array.isArray(row) || row.length !== matrix.ranges.length
      || row.some((cell: unknown) => typeof cell !== 'boolean'))
    || !newRows.length || new Set(newRows).size !== newRows.length) return null;
  const oldRows = matrix.rows as string[];
  const positions = new Map(newRows.map((key, index) => [key, index]));
  const mapped = oldRows.filter(key => positions.has(key)).map(key => positions.get(key)!);
  if (!mapped.length || mapped.some((position, index) => index > 0 && position !== mapped[index - 1] + 1)) return null;
  if (new Set(matrix.ranges).size !== matrix.ranges.length
    || anchors.some(anchor => !matrix.ranges.includes(anchor.range)
    || !Array.isArray(anchor.locks) || anchor.locks.some((lock: any) => !positions.has(lock.fragment)))) return null;
  const grid = newRows.map(() => matrix.ranges.map(() => false));
  const oldIndex = new Map(oldRows.map((key, index) => [key, index]));
  for (const [newIndex, key] of newRows.entries()) {
    const previous = oldIndex.get(key);
    if (previous !== undefined) grid[newIndex] = [...matrix.grid[previous]];
  }
  // Newly uncovered candidate rows are mandatory locks, even below the
  // matrix's normal growth threshold.
  for (const anchor of anchors as any[]) {
    const column = (matrix.ranges as string[]).indexOf(anchor.range);
    for (const lock of anchor.locks as Array<{ fragment: string }>) {
      const row = positions.get(lock.fragment);
      if (column < 0 || row === undefined) return null;
      // Preserve the mandatory lock even when its score is bad. The matrix
      // read validator reports the anchor-specific scoring error; dropping
      // the adapted matrix here would make the caller replace it with an
      // empty skeleton and lose addressed edits.
      grid[row][column] = true;
    }
  }
  if (mapped.length) {
    const first = mapped[0];
    const last = mapped[mapped.length - 1];
    matrix.ranges.forEach((range: string, column: number) => {
      if (matrix.grid[0][column]) {
        for (let index = first - 1; index >= 0 && grows(energetics, range, newRows[index]); index -= 1) grid[index][column] = true;
      }
      if (matrix.grid[matrix.grid.length - 1][column]) {
        for (let index = last + 1; index < newRows.length && grows(energetics, range, newRows[index]); index += 1) grid[index][column] = true;
      }
    });
  }
  return { ...matrix, rows: [...newRows], grid } as Stage5S14;
}

export function stage5CanAdvance(
  snapshot: ProjectSnapshot,
  state: {
    ready: boolean;
    stage: number;
    energeticsFresh: boolean;
    bestFitsFresh: boolean;
    geometryFresh: boolean;
    gaussianFresh: boolean;
    gaussianBackendFailed?: boolean;
    validatorFailed: boolean;
    fragments: unknown;
    rangeRows: unknown;
    stackings?: unknown;
    anchors?: readonly Stage5Anchor[];
    energetics: unknown;
    energyEdits: unknown;
  },
): { allowed: boolean; reasons: string[] } {
  const reasons: string[] = [];
  if (!state.ready || state.stage !== 5) reasons.push('The shared project is not live in stage 5.');
  if (!state.energeticsFresh) reasons.push('Energetics do not match the current best fits and protein ranges.');
  if (!state.bestFitsFresh || !state.geometryFresh) reasons.push('Best-fit geometry is stale.');
  if (!state.gaussianFresh) reasons.push('Fresh Gaussian scores are required.');
  if (state.gaussianBackendFailed) reasons.push('Best-fit stacking geometry or Gaussian scoring failed.');
  if (state.validatorFailed) reasons.push('Energy-threshold validation failed.');
  if (!energyThresholdEditsIsWellFormed(state.energyEdits)) reasons.push('The energy-threshold edits are malformed.');
  if (!isRecord(state.fragments) || !isRecord(state.fragments.echo)
    || !isRecord(state.energyEdits) || !sameAddressBasis(state.energyEdits.basis, state.fragments.echo)) {
    reasons.push('Energy-threshold edits do not match the current PDB and RNA range.');
  }
  const rows = stage5Rows(state.fragments, snapshot.best_fit_clip);
  const anchors = state.anchors ?? (rows && state.stackings !== undefined
    ? candidateAnchors(state.fragments, snapshot.stacking_selection, state.stackings, rows) : null);
  const ranges = stage5RangeIds(state.rangeRows);
  const matrix = snapshot.interaction_matrix;
  let matrixError = rows && ranges ? validateStage5MatrixRead(matrix, rows, ranges) : 'matrix extent is unavailable';
  if (!matrixError && anchors && state.energetics) {
    matrixError = validateStage5InteractionMatrix(matrix, rows!, anchors, state.energetics, ranges!);
  }
  if (state.stackings !== undefined && !anchors) matrixError = 'selected anchors are unavailable';
  if (matrixError || !markerIsWellFormed('interaction_matrix', matrix)
    || !markerMatchesUpstream('interaction_matrix', matrix, snapshot)) {
    reasons.push(`The interaction matrix is invalid${matrixError ? `: ${matrixError}` : '.'}`);
  } else if (state.energetics && isRecord(matrix)) {
    const edits = state.energyEdits;
    const currentRangeRows = Array.isArray(state.rangeRows) ? state.rangeRows.slice(1) : [];
    const editRows = isRecord(edits) && Array.isArray(edits.ranges) ? edits.ranges : [];
    for (let row = 0; row < matrix.rows.length; row += 1) {
      for (let column = 0; column < matrix.ranges.length; column += 1) {
        if (!matrix.grid[row][column]) continue;
        const rangeId = matrix.ranges[column];
        const rangeIndex = Number(String(rangeId).slice('range-'.length)) - 2;
        const rangeRow = currentRangeRows[rangeIndex];
        if (!hasValidEnergies(state.energetics, rangeId, matrix.rows[row])) {
          reasons.push(`The scoring result for ${matrix.rows[row]} / ${rangeId} is missing or contains an error.`);
          continue;
        }
        const editEntry = editRows.find((item: any) => isRecord(item) && isRecord(item.row)
          && isRecord(rangeRow) && addressKey(item.row as Stage5Range) === addressKey(rangeRow as Stage5Range));
        const value = editEntry?.values?.[matrix.rows[row]];
        if (value === null || value === undefined || typeof value !== 'string' || !validEnergyEdit(value)) {
          reasons.push(`Enter a valid energy threshold for ${matrix.rows[row]} / ${rangeId}.`);
        }
      }
    }
  }
  return { allowed: reasons.length === 0, reasons: [...new Set(reasons)] };
}
