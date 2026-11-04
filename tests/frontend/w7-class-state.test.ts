// W7 initial acceptance vectors against existing state interfaces; not yet run.
import { expect, test } from 'vitest';
import {
  deferredWriteIsCurrent, echoIsFresh, markerIsWellFormed,
} from '../../src/lib/project-state';

function classes(rrm: number = 90): any {
  return { schema: 1, basis: ['anchor-A'], active: true, rows: [
    { key: 'anchor-A', rrm, ...(rrm >= 90 ? { tilt30: true, disp45: false } : {}) },
  ] };
}

test.each([50, 60, 70, 80, 90, 95])('W7 accepts complete class %s and inactive retained edits', rrm => {
  const value = classes(rrm);
  expect(markerIsWellFormed('gaussian_classes', value)).toBe(true);
  value.active = false;
  expect(markerIsWellFormed('gaussian_classes', value)).toBe(true);
});

test('W7 flags are present iff 90/95 and are strict booleans', () => {
  for (const rrm of [90, 95]) {
    for (const flag of ['tilt30', 'disp45']) {
      const absent = classes(rrm);
      delete absent.rows[0][flag];
      expect(markerIsWellFormed('gaussian_classes', absent)).toBe(false);
      for (const bad of [null, 0, 1, 'false']) {
        const value = classes(rrm);
        value.rows[0][flag] = bad;
        expect(markerIsWellFormed('gaussian_classes', value)).toBe(false);
      }
    }
  }
  for (const rrm of [50, 60, 70, 80]) {
    const value = classes(rrm);
    value.rows[0].tilt30 = false;
    expect(markerIsWellFormed('gaussian_classes', value)).toBe(false);
  }
});

test('W7 requires unique class keys and the exact selected basis', () => {
  const duplicate = classes();
  duplicate.rows.push({ ...duplicate.rows[0] });
  expect(markerIsWellFormed('gaussian_classes', duplicate)).toBe(false);
  const missing = classes();
  missing.rows = [];
  expect(markerIsWellFormed('gaussian_classes', missing)).toBe(false);
  const wrong = classes();
  wrong.basis = ['anchor-B'];
  expect(markerIsWellFormed('gaussian_classes', wrong)).toBe(false);
});

test('W7 old validation cannot match optimistic class or clip changes', () => {
  const oldEcho = { rows: classes().rows, clip: [0, 0] };
  expect(echoIsFresh(structuredClone(oldEcho), oldEcho)).toBe(true);
  const classEdit = structuredClone(oldEcho);
  classEdit.rows[0].disp45 = true;
  expect(echoIsFresh(oldEcho, classEdit)).toBe(false);
  const clipEdit = structuredClone(oldEcho);
  clipEdit.clip[0] = 1;
  expect(echoIsFresh(oldEcho, clipEdit)).toBe(false);
  expect(echoIsFresh(null, classEdit)).toBe(false);
});

test('W7 deferred class write is invalid after Back, remote revision or reconnect', () => {
  const guard = { stageAtIssue: 6 as const, currentStage: 6 as const, live: true,
    revisionAtIssue: 3, currentRevision: 3, connectionAtIssue: 1, currentConnection: 1 };
  expect(deferredWriteIsCurrent(guard)).toBe(true);
  expect(deferredWriteIsCurrent({ ...guard, currentStage: 5 })).toBe(false);
  expect(deferredWriteIsCurrent({ ...guard, currentRevision: 4 })).toBe(false);
  expect(deferredWriteIsCurrent({ ...guard, currentConnection: 2 })).toBe(false);
  expect(deferredWriteIsCurrent({ ...guard, live: false })).toBe(false);
});
