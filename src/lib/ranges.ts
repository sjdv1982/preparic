import type { ActiveCell, ChainInfo, RangeColumn, RangeRow } from './types';

let rowSequence = 0;

export function createRangeRow(chain = ''): RangeRow {
  rowSequence += 1;
  return { id: `range-${rowSequence}`, alwaysVisible: false, chain, start: '', end: '' };
}

export function normalizeRow(row: RangeRow): RangeRow {
  const start = Number(row.start);
  const end = Number(row.end);
  if (row.start !== '' && row.end !== '' && Number.isFinite(start) && Number.isFinite(end) && start > end) {
    return { ...row, start: String(end), end: String(start) };
  }
  return row;
}

export function normalizeLeavingRow(rows: RangeRow[], previous: ActiveCell | null, next: ActiveCell | null): RangeRow[] {
  if (!previous || previous.rowId === next?.rowId) return rows;
  return rows.map(row => row.id === previous.rowId ? normalizeRow(row) : row);
}

export function setActiveCell(rows: RangeRow[], previous: ActiveCell | null, next: ActiveCell | null): RangeRow[] {
  return normalizeLeavingRow(rows, previous, next);
}

export function updateRow(rows: RangeRow[], rowId: string, patch: Partial<RangeRow>): RangeRow[] {
  return rows.map(row => row.id === rowId ? { ...row, ...patch } : row);
}

export function addRow(rows: RangeRow[], activeRowId?: string | null): RangeRow[] {
  const index = activeRowId ? rows.findIndex(row => row.id === activeRowId) : -1;
  const chain = index >= 0 ? rows[index].chain : '';
  const next = [...rows];
  next.splice(index >= 0 ? index + 1 : rows.length, 0, createRangeRow(chain));
  return next;
}

export function removeRow(rows: RangeRow[], rowId: string): RangeRow[] {
  return rows.filter(row => row.id !== rowId);
}

export function moveRow(rows: RangeRow[], rowId: string, direction: -1 | 1): RangeRow[] {
  const index = rows.findIndex(row => row.id === rowId);
  const target = index + direction;
  if (index < 0 || target < 0 || target >= rows.length) return rows;
  const next = [...rows];
  [next[index], next[target]] = [next[target], next[index]];
  return next;
}

export interface ValidationIssue {
  rowId: string;
  column: RangeColumn;
  message: string;
}

export function validateRows(rows: RangeRow[], chains: ChainInfo[], includeEmpty = true): ValidationIssue[] {
  const issues: ValidationIssue[] = [];
  const byName = new Map(chains.map(chain => [chain.name, chain]));
  for (const row of rows) {
    const chain = byName.get(row.chain);
    if (!chain && (includeEmpty || row.chain !== '')) {
      issues.push({ rowId: row.id, column: 'chain', message: row.chain ? 'Choose a valid chain.' : 'Choose a chain.' });
    }
    for (const column of ['start', 'end'] as const) {
      const raw = row[column].trim();
      const value = Number(raw);
      if (raw === '' && !includeEmpty) continue;
      if (raw === '' || !Number.isInteger(value)) {
        issues.push({ rowId: row.id, column, message: 'Enter an integer residue number.' });
      } else if (chain && !chain.seqIds.has(value)) {
        issues.push({ rowId: row.id, column, message: `Residue ${value} is not in chain ${chain.name}.` });
      }
    }
  }
  return issues;
}

export function isCellInvalid(issues: ValidationIssue[], rowId: string, column: RangeColumn): boolean {
  return issues.some(issue => issue.rowId === rowId && issue.column === column);
}

export function moveCell(rows: RangeRow[], cell: ActiveCell, key: 'ArrowUp' | 'ArrowDown' | 'ArrowLeft' | 'ArrowRight' | 'Tab', shift = false): ActiveCell | null {
  const columns: RangeColumn[] = ['chain', 'start', 'end'];
  const rowIndex = rows.findIndex(row => row.id === cell.rowId);
  const columnIndex = columns.indexOf(cell.column);
  if (rowIndex < 0) return null;
  if (key === 'ArrowUp') return { rowId: rows[Math.max(0, rowIndex - 1)].id, column: cell.column };
  if (key === 'ArrowDown') return { rowId: rows[Math.min(rows.length - 1, rowIndex + 1)].id, column: cell.column };
  if (key === 'ArrowLeft') return { rowId: rows[rowIndex].id, column: columns[Math.max(0, columnIndex - 1)] };
  if (key === 'ArrowRight') return { rowId: rows[rowIndex].id, column: columns[Math.min(columns.length - 1, columnIndex + 1)] };
  const step = shift ? -1 : 1;
  let nextColumn = columnIndex + step;
  let nextRow = rowIndex;
  if (nextColumn < 0 || nextColumn >= columns.length) {
    nextRow = Math.max(0, Math.min(rows.length - 1, rowIndex + step));
    nextColumn = step > 0 ? 0 : columns.length - 1;
  }
  return { rowId: rows[nextRow].id, column: columns[nextColumn] };
}
