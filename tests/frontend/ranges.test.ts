import { describe, expect, it } from 'vitest';
import { addRow, createRangeRow, moveRow, normalizeLeavingRow, normalizeRow, removeRow, updateRow, validateRows } from '../../src/lib/ranges';
import type { ChainInfo, RangeRow } from '../../src/lib/types';

const chains: ChainInfo[] = [{
  name: 'A', labelAsymId: 'A', instanceId: '1_555', seqIds: new Set([1, 2, 10, 20]), seqMin: 1, seqMax: 20,
}];

function row(values: Partial<RangeRow> = {}): RangeRow {
  return { ...createRangeRow(), alwaysVisible: false, chain: 'A', start: '10', end: '2', ...values };
}

describe('range row normalization', () => {
  it('swaps start and end when leaving a row, but not when changing cells within it', () => {
    const current = row();
    const rows = [current];
    const active = { rowId: current.id, column: 'start' as const };
    expect(normalizeLeavingRow(rows, active, { rowId: current.id, column: 'end' })).toBe(rows);
    expect(normalizeLeavingRow(rows, active, null)[0]).toMatchObject({ start: '2', end: '10' });
    expect(normalizeRow(current)).toMatchObject({ start: '2', end: '10' });
    expect(updateRow(rows, current.id, { end: '20' })[0]).toMatchObject({ start: '10', end: '20' });
  });
});

describe('validateRows', () => {
  it('accepts label sequence ids that exist in the selected chain', () => {
    expect(validateRows([row({ start: '2', end: '10' })], chains)).toEqual([]);
  });

  it('reports unknown chains, missing cells on submit, and residues outside a chain', () => {
    const result = validateRows([row({ chain: 'Z', start: '', end: '3' })], chains);
    expect(result.map(issue => issue.column)).toEqual(['chain', 'start']);
    expect(validateRows([row({ start: '99', end: '2' })], chains).map(issue => issue.column)).toEqual(['start']);
    expect(validateRows([row({ chain: '', start: '', end: '' })], chains, false)).toEqual([]);
  });
});

describe('range row operations', () => {
  it('adds after the active row and copies its chain', () => {
    const first = row({ start: '', end: '' });
    const next = addRow([first], first.id);
    expect(next).toHaveLength(2);
    expect(next[1]).toMatchObject({ chain: 'A', start: '', end: '' });
  });

  it('removes and moves rows without mutating the input', () => {
    const first = row();
    const second = row({ chain: 'B' });
    const rows = [first, second];
    expect(moveRow(rows, second.id, -1).map(item => item.id)).toEqual([second.id, first.id]);
    expect(rows.map(item => item.id)).toEqual([first.id, second.id]);
    expect(removeRow(rows, first.id)).toEqual([second]);
    expect(updateRow(rows, first.id, { start: '1' })[0].start).toBe('1');
  });
});
