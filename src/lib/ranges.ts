import type { ActiveCell, ChainInfo, RangeColumn, RangeRow } from './types';

let rowSequence = 0;

export function createRangeRow(chain = ''): RangeRow {
  rowSequence += 1;
  return { id: `range-${rowSequence}`, alwaysVisible: false, chain, start: '', end: '' };
}

export function cloneRangeRows(rows: RangeRow[]): RangeRow[] {
  return rows.map(row => ({ ...row, id: createRangeRow().id }));
}

export function normalizeRow(row: RangeRow): RangeRow {
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
  for (const [rowIndex, row] of rows.entries()) {
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
    const start = Number(row.start.trim());
    const end = Number(row.end.trim());
    if (row.start.trim() !== '' && row.end.trim() !== ''
      && Number.isInteger(start) && Number.isInteger(end) && start > end) {
      issues.push({ rowId: row.id, column: 'end', message: 'End residue must not come before the start residue.' });
    }
    // Range 1 is RNA. Its residue labels must cover the whole requested
    // interval; a gap inside the interval is just as invalid as a missing
    // endpoint. Protein ranges intentionally retain their existing behavior.
    if (rowIndex === 0 && chain && row.start.trim() !== '' && row.end.trim() !== ''
      && Number.isInteger(start) && Number.isInteger(end) && start <= end) {
      const missing = [];
      for (let resid = start; resid <= end; resid += 1) {
        if (!chain.seqIds.has(resid)) missing.push(resid);
      }
      if (missing.length > 0) {
        issues.push({
          rowId: row.id,
          column: 'end',
          message: `RNA range contains unobserved residue${missing.length === 1 ? '' : 's'}: ${missing.join(', ')}.`,
        });
      }
      const breaks = [...(chain.rnaBreakAfter ?? [])]
        .filter(resid => resid >= start && resid < end)
        .sort((left, right) => left - right);
      if (breaks.length > 0) {
        issues.push({
          rowId: row.id,
          column: 'end',
          message: `RNA range has a chemical chain break between residues ${breaks.map(resid => `${resid} and ${resid + 1}`).join(', ')}.`,
        });
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
