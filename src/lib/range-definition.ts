import { parse, stringify } from 'yaml';
import { createRangeRow } from './ranges';
import type { RangeRow } from './types';

function asText(value: unknown): string {
  if (typeof value === 'string') return value;
  if (typeof value === 'number' && Number.isFinite(value)) return String(value);
  return '';
}

function yamlValue(value: string): string | number {
  const trimmed = value.trim();
  if (trimmed !== '') {
    const number = Number(trimmed);
    if (Number.isSafeInteger(number)) return number;
  }
  return value;
}

function isDefinedRange(row: RangeRow): boolean {
  return row.alwaysVisible || row.chain !== '' || row.start !== '' || row.end !== '';
}

export function rangeDefinitionSignature(rows: RangeRow[]): string {
  return JSON.stringify(rows.filter(isDefinedRange).map(row => ({
    chain: row.chain,
    start: yamlValue(row.start),
    end: yamlValue(row.end),
    alwaysVisible: row.alwaysVisible,
  })));
}

export function serializeRangeDefinition(rows: RangeRow[]): string {
  return stringify({
    ranges: rows.filter(isDefinedRange).map(row => ({
      chain: row.chain,
      start: yamlValue(row.start),
      end: yamlValue(row.end),
      always_visible: row.alwaysVisible,
    })),
  }, { defaultStringType: 'QUOTE_DOUBLE' });
}

export function parseRangeDefinition(source: string): RangeRow[] {
  if (source.trim() === '') return [];
  const value: unknown = parse(source);
  if (!value || typeof value !== 'object' || Array.isArray(value)) {
    throw new Error('The YAML range definition must be a mapping with a ranges list.');
  }
  const ranges = (value as { ranges?: unknown }).ranges;
  if (!Array.isArray(ranges)) {
    throw new Error('The YAML range definition must contain a ranges list.');
  }

  return ranges.map((entry: unknown) => {
    if (!entry || typeof entry !== 'object' || Array.isArray(entry)) {
      throw new Error('Each YAML range must be a mapping.');
    }
    const row = entry as {
      chain?: unknown;
      start?: unknown;
      end?: unknown;
      always_visible?: unknown;
      alwaysVisible?: unknown;
    };
    const restored = createRangeRow(asText(row.chain));
    restored.start = asText(row.start);
    restored.end = asText(row.end);
    const alwaysVisible = row.always_visible ?? row.alwaysVisible;
    restored.alwaysVisible = alwaysVisible === true || (typeof alwaysVisible === 'string' && alwaysVisible.toLowerCase() === 'true');
    return restored;
  });
}

