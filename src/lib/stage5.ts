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
  const ranges = anchors.map(anchor => anchor.range);
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
): string | null {
  const expectedRanges = anchors.map(anchor => anchor.range);
  const shapeError = validateStage5MatrixRead(matrix, expectedRows, expectedRanges);
  if (shapeError) return shapeError;
  if (!isRecord(matrix) || !Array.isArray(matrix.grid)) return 'matrix structure is invalid';
  const rows = matrix.rows as string[];
  const grid = matrix.grid as boolean[][];
  const rowIndex = new Map(rows.map((key, index) => [key, index]));
  const intervals = new Map<string, Array<{ column: number; first: number; last: number }>>();
  for (const [column, anchor] of anchors.entries()) {
    const lockPositions: number[] = [];
    for (const lock of anchor.locks) {
      const position = rowIndex.get(lock.fragment);
      if (position === undefined) return `selected anchor ${anchor.resid} has a clipped lock`;
      if (!grid[position][column]) return `selected anchor ${anchor.resid} is not locked in its column`;
      lockPositions.push(position);
    }
    if (!lockPositions.length) return `selected anchor ${anchor.resid} has no candidate locks`;
    let index = 0;
    const columnIntervals: Array<{ first: number; last: number }> = [];
    while (index < rows.length) {
      if (!grid[index][column]) { index += 1; continue; }
      const first = index;
      while (index + 1 < rows.length && grid[index + 1][column]) index += 1;
      const last = index;
      if (!lockPositions.some(position => position >= first && position <= last)) {
        return `true run in ${anchor.range} has no selected anchor`;
      }
      for (let row = first; row <= last; row += 1) {
        if (!hasValidEnergies(energetics, anchor.range, rows[row])) {
          return `scoring error at ${rows[row]} / ${anchor.range}`;
        }
      }
      columnIntervals.push({ first, last });
      index += 1;
    }
    if (columnIntervals.length > 1) return `anchor column ${column} has more than one run`;
    const firstInterval = columnIntervals[0];
    if (firstInterval) {
      const rangeIntervals = intervals.get(anchor.range) ?? [];
      rangeIntervals.push({ column, ...firstInterval });
      intervals.set(anchor.range, rangeIntervals);
    }
  }
  for (const [range, rangeIntervals] of intervals) {
    const ordered = [...rangeIntervals].sort((left, right) => left.first - right.first
      || left.last - right.last || left.column - right.column);
    for (let index = 1; index < ordered.length; index += 1) {
      const left = ordered[index - 1];
      const right = ordered[index];
      if (right.first <= left.last + 1
        && (right.first !== left.first || right.last !== left.last)) {
        return `touching or overlapping anchors in ${range} must share a joined run`;
      }
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
    || matrix.ranges.length !== anchors.length
    || matrix.ranges.some((range: unknown, index: number) => range !== anchors[index].range)) return null;
  const rows = matrix.rows as string[];
  const positions = new Map(rows.map((key, index) => [key, index]));
  const grid = rows.map(() => anchors.map(() => false));
  const byRange = new Map<string, Array<{ column: number; first: number; last: number }>>();
  for (const [column, anchor] of anchors.entries()) {
    const locks = anchor.locks.map(lock => positions.get(lock.fragment));
    if (!locks.length || locks.some(position => position === undefined)) return null;
    for (const position of locks as number[]) {
      if (!hasValidEnergies(energetics, anchor.range, rows[position])) return null;
    }
    let first = Math.min(...locks as number[]);
    let last = Math.max(...locks as number[]);
    const otherLocks = new Set(anchors.flatMap((other, otherColumn) => otherColumn !== column
      && other.range === anchor.range ? other.locks.map(lock => positions.get(lock.fragment)!) : []));
    while (first > 0 && !otherLocks.has(first - 1) && grows(energetics, anchor.range, rows[first - 1])) first -= 1;
    while (last + 1 < rows.length && !otherLocks.has(last + 1)
      && grows(energetics, anchor.range, rows[last + 1])) last += 1;
    const sameRange = byRange.get(anchor.range) ?? [];
    sameRange.push({ column, first, last });
    byRange.set(anchor.range, sameRange);
  }
  for (const rangeIntervals of byRange.values()) {
    const ordered = [...rangeIntervals].sort((left, right) => left.first - right.first
      || left.last - right.last || left.column - right.column);
    const components: Array<typeof ordered> = [];
    for (const item of ordered) {
      const tail = components.at(-1);
      if (tail && item.first <= Math.max(...tail.map(value => value.last)) + 1) tail.push(item);
      else components.push([item]);
    }
    for (const component of components) {
      const first = Math.min(...component.map(item => item.first));
      const last = Math.max(...component.map(item => item.last));
      for (const item of component) for (let row = first; row <= last; row += 1) grid[row][item.column] = true;
    }
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
  const range = matrix.ranges[column];
  const sameRunColumns = matrix.ranges.map((item, index) => item === range
    && matrix.grid.every((line, rowIndex) => line[index] === matrix.grid[rowIndex][column]) ? index : -1)
    .filter(index => index >= 0);
  for (const target of sameRunColumns) next.grid[row][target] = value;
  return validateStage5InteractionMatrix(next, matrix.rows, anchors, energetics) === null ? next : null;
}

export interface UntilNextAnchorState { nextColumn: number | null; allowed: boolean; joined: boolean; forced: boolean }

function columnInterval(matrix: Stage5S14, column: number): { first: number; last: number } | null {
  const occupied = matrix.grid.map((row, index) => row[column] ? index : -1).filter(index => index >= 0);
  return occupied.length ? { first: occupied[0], last: occupied.at(-1)! } : null;
}

function lockInterval(anchor: Stage5Anchor, rows: readonly string[]): { first: number; last: number } | null {
  const positions = anchor.locks.map(lock => rows.indexOf(lock.fragment));
  if (!positions.length || positions.some(position => position < 0)) return null;
  return { first: Math.min(...positions), last: Math.max(...positions) };
}

/** Describe the optional link to the next downstream anchor on the same range. */
export function untilNextAnchorState(
  matrix: Stage5S14,
  column: number,
  anchors: readonly Stage5Anchor[],
): UntilNextAnchorState {
  const current = anchors[column];
  if (!current) return { nextColumn: null, allowed: false, joined: false, forced: false };
  let nextColumn: number | null = null;
  for (let index = column + 1; index < anchors.length; index += 1) {
    if (anchors[index].range === current.range) { nextColumn = index; break; }
  }
  if (nextColumn === null) return { nextColumn: null, allowed: false, joined: false, forced: false };
  const leftLocks = lockInterval(current, matrix.rows);
  const rightLocks = lockInterval(anchors[nextColumn], matrix.rows);
  if (!leftLocks || !rightLocks) return { nextColumn, allowed: false, joined: false, forced: false };
  const forced = rightLocks.first <= leftLocks.last + 1;
  const left = columnInterval(matrix, column);
  const right = columnInterval(matrix, nextColumn);
  const joined = !!left && !!right && left.first === right.first && left.last === right.last;
  return { nextColumn, allowed: true, joined, forced };
}

/** Join or split the current column and its next same-range anchor. */
export function setStage5UntilNextAnchor(
  matrix: Stage5S14,
  column: number,
  enabled: boolean,
  anchors: readonly Stage5Anchor[],
  energetics: unknown,
): Stage5S14 | null {
  const state = untilNextAnchorState(matrix, column, anchors);
  if (!state.allowed || state.nextColumn === null || (!enabled && state.forced)) return null;
  if (enabled === state.joined) return matrix;
  const range = anchors[column].range;
  const sameRangeColumns = anchors.map((anchor, index) => anchor.range === range ? index : -1)
    .filter(index => index >= 0);
  const leftPosition = sameRangeColumns.indexOf(column);
  const rightPosition = sameRangeColumns.indexOf(state.nextColumn);
  if (leftPosition < 0 || rightPosition !== leftPosition + 1) return null;

  const links = sameRangeColumns.slice(0, -1).map((leftColumn, index) => {
    const left = columnInterval(matrix, leftColumn);
    const right = columnInterval(matrix, sameRangeColumns[index + 1]);
    return !!left && !!right && left.first === right.first && left.last === right.last;
  });
  let componentStart = leftPosition;
  let componentEnd = rightPosition;
  while (componentStart > 0 && links[componentStart - 1]) componentStart -= 1;
  while (componentEnd < links.length && links[componentEnd]) componentEnd += 1;
  const componentColumns = sameRangeColumns.slice(componentStart, componentEnd + 1);
  const componentIntervals = componentColumns.map(target => columnInterval(matrix, target));
  if (componentIntervals.some(interval => interval === null)) return null;
  const first = Math.min(...componentIntervals.map(interval => interval!.first));
  const last = Math.max(...componentIntervals.map(interval => interval!.last));
  const next = { ...matrix, grid: matrix.grid.map(row => [...row]) };

  if (enabled) {
    for (const target of componentColumns) {
      for (let row = first; row <= last; row += 1) next.grid[row][target] = true;
    }
  } else {
    const leftColumns = sameRangeColumns.slice(componentStart, leftPosition + 1);
    const rightColumns = sameRangeColumns.slice(rightPosition, componentEnd + 1);
    const leftLocks = leftColumns.flatMap(target => anchors[target].locks.map(lock => matrix.rows.indexOf(lock.fragment)));
    const rightLocks = rightColumns.flatMap(target => anchors[target].locks.map(lock => matrix.rows.indexOf(lock.fragment)));
    if (!leftLocks.length || !rightLocks.length || leftLocks.some(index => index < 0) || rightLocks.some(index => index < 0)) return null;
    const leftLast = Math.max(...leftLocks);
    const rightFirst = Math.min(...rightLocks);
    if (rightFirst <= leftLast + 1) return null;
    for (const target of leftColumns) {
      for (let row = 0; row < matrix.rows.length; row += 1) {
        next.grid[row][target] = row >= first && row <= Math.min(last, leftLast);
      }
    }
    for (const target of rightColumns) {
      for (let row = 0; row < matrix.rows.length; row += 1) {
        next.grid[row][target] = row >= Math.max(first, rightFirst) && row <= last;
      }
    }
  }
  return validateStage5InteractionMatrix(next, matrix.rows, anchors, energetics) === null ? next : null;
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
    || !Array.isArray(matrix.grid) || matrix.grid.length !== matrix.rows.length
    || matrix.grid.some((row: unknown) => !Array.isArray(row) || row.length !== matrix.ranges.length
      || row.some((cell: unknown) => typeof cell !== 'boolean'))
    || !newRows.length || new Set(newRows).size !== newRows.length) return null;
  const oldRows = matrix.rows as string[];
  const positions = new Map(newRows.map((key, index) => [key, index]));
  const mapped = oldRows.filter(key => positions.has(key)).map(key => positions.get(key)!);
  if (!mapped.length || mapped.some((position, index) => index > 0 && position !== mapped[index - 1] + 1)) return null;
  if (anchors.some(anchor => !matrix.ranges.includes(anchor.range)
    || !Array.isArray(anchor.locks) || anchor.locks.some((lock: any) => !positions.has(lock.fragment)))) return null;
  if (anchors.length !== matrix.ranges.length
    || anchors.some((anchor, index) => anchor.range !== matrix.ranges[index])) return null;
  const grid = newRows.map(() => matrix.ranges.map(() => false));
  const oldIndex = new Map(oldRows.map((key, index) => [key, index]));
  for (const [newIndex, key] of newRows.entries()) {
    const previous = oldIndex.get(key);
    if (previous !== undefined) grid[newIndex] = [...matrix.grid[previous]];
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
  const columnsByRange = new Map<string, number[]>();
  matrix.ranges.forEach((range: string, column: number) => {
    const columns = columnsByRange.get(range) ?? [];
    columns.push(column);
    columnsByRange.set(range, columns);
  });
  for (const columns of columnsByRange.values()) {
    const intervals = columns.map(column => {
      const occupied = grid.map((row, index) => row[column] ? index : -1).filter(index => index >= 0);
      return occupied.length ? { columns: [column], first: occupied[0], last: occupied.at(-1)! } : null;
    }).filter((item): item is { columns: number[]; first: number; last: number } => item !== null)
      .sort((left, right) => left.first - right.first || left.last - right.last);
    const components: Array<{ columns: number[]; first: number; last: number }> = [];
    for (const interval of intervals) {
      const tail = components.at(-1);
      if (tail && interval.first <= tail.last + 1) {
        tail.last = Math.max(tail.last, interval.last);
        tail.columns.push(...interval.columns);
      } else components.push({ ...interval });
    }
    for (const component of components.filter(item => item.columns.length > 1)) {
      for (const column of component.columns) {
        for (let index = component.first; index <= component.last; index += 1) grid[index][column] = true;
      }
    }
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
    gaussianClassesValid: boolean;
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
  if (!state.gaussianFresh || !state.gaussianClassesValid) reasons.push('Fresh Gaussian scores and valid classes are required.');
  if (state.validatorFailed) reasons.push('Energy-threshold validation failed.');
  if (!energyThresholdEditsIsWellFormed(state.energyEdits)) reasons.push('The energy-threshold edits are malformed.');
  if (!isRecord(state.fragments) || !isRecord(state.fragments.echo)
    || !isRecord(state.energyEdits) || !sameAddressBasis(state.energyEdits.basis, state.fragments.echo)) {
    reasons.push('Energy-threshold edits do not match the current PDB and RNA range.');
  }
  const rows = stage5Rows(state.fragments, snapshot.best_fit_clip);
  const anchors = state.anchors ?? (rows && state.stackings !== undefined
    ? candidateAnchors(state.fragments, snapshot.stacking_selection, state.stackings, rows) : null);
  const ranges = anchors ? anchors.map(anchor => anchor.range) : stage5RangeIds(state.rangeRows);
  const matrix = snapshot.interaction_matrix;
  let matrixError = rows && ranges ? validateStage5MatrixRead(matrix, rows, ranges) : 'matrix extent is unavailable';
  if (!matrixError && anchors && state.energetics) {
    matrixError = validateStage5InteractionMatrix(matrix, rows!, anchors, state.energetics);
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
