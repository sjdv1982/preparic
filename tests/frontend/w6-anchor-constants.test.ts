import { readFileSync } from 'node:fs';
import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import { expect, test } from 'vitest';

import {
  BOUNDARY_RECOMPUTE_TOLERANCE,
  L_OV,
  MIN_CELL_FRACTION,
} from '../../src/lib/anchor-constants';

const pythonConstants = resolve(
  dirname(fileURLToPath(import.meta.url)),
  '../../../alaric/alaric/anchor_constants.py',
);

function pythonNumericConstant(name: string): number {
  const source = readFileSync(pythonConstants, 'utf8');
  const match = source.match(new RegExp(`^${name}\\s*=\\s*(.+?)\\s*(?:#.*)?$`, 'm'));
  if (!match) throw new Error(`Alaric anchor_constants.py does not define ${name}`);
  const expression = match[1].trim();
  if (/^(?:math\.)?sqrt\(3(?:\.0)?\)$/.test(expression)) return Math.sqrt(3);
  if (/^1\s*\/\s*3$/.test(expression)) return 1 / 3;
  const value = Number(expression);
  if (!Number.isFinite(value)) {
    throw new Error(`Unsupported Python expression for ${name}: ${expression}`);
  }
  return value;
}

test('Preparic mirrors Alaric A11 anchor constants', () => {
  expect(L_OV).toBe(pythonNumericConstant('L_OV'));
  expect(BOUNDARY_RECOMPUTE_TOLERANCE).toBe(
    pythonNumericConstant('BOUNDARY_RECOMPUTE_TOLERANCE'),
  );
  expect(MIN_CELL_FRACTION).toBe(pythonNumericConstant('MIN_CELL_FRACTION'));
  expect(L_OV).toBe(Math.sqrt(3));
  expect(BOUNDARY_RECOMPUTE_TOLERANCE).toBe(1e-5);
  expect(MIN_CELL_FRACTION).toBe(1 / 3);
});
