import { describe, expect, it } from 'vitest';
import {
  exportNeedsAcknowledgement,
  exportIsStale,
  exportReadiness,
  exportWarningAcknowledged,
} from '../../src/lib/export';

const digest = 'a'.repeat(64);
const shareChecksum = 'b'.repeat(64);
const warning = { kind: 'gaussian', address: 'stack-A', message: 'The chosen filter failed.' };

function admission(overrides: Record<string, unknown> = {}) {
  return {
    schema: 1,
    ready: true,
    blocks: [],
    warnings: [],
    digest,
    ...overrides,
  };
}

function readinessInput(overrides: Record<string, unknown> = {}) {
  return {
    stage: 6,
    connected: true,
    writePending: false,
    stage6EchoesFresh: true,
    admission: admission(),
    admissionState: 'complete',
    treeState: 'complete',
    admissionNodeChecksum: shareChecksum,
    admissionValueChecksum: shareChecksum,
    ...overrides,
  };
}

describe('E6 export readiness and X5 acknowledgment', () => {
  it('enables export only for the current complete, ready stage-6 admission and tree', () => {
    expect(exportReadiness(readinessInput()).enabled).toBe(true);
  });

  it.each([
    ['outside stage 6', { stage: 5 }],
    ['while disconnected', { connected: false }],
    ['while a write is pending or debounced', { writePending: true }],
    ['when stage-6 result echoes are stale', { stage6EchoesFresh: false }],
    ['before admission is complete', { admissionState: 'waiting' }],
    ['when the tree is not complete', { treeState: 'computing' }],
    ['when the published admission checksum is stale', { admissionNodeChecksum: 'c'.repeat(64) }],
    ['when the admission value has no checksum', { admissionValueChecksum: null }],
    ['when admission is not ready', { admission: admission({ ready: false }) }],
    ['when a hard block is present', { admission: admission({ blocks: [{ code: 'blocked' }] }) }],
    ['when admission has no digest', { admission: admission({ digest: null }) }],
  ])('keeps Export disabled %s', (_reason, changes) => {
    expect(exportReadiness(readinessInput(changes)).enabled).toBe(false);
  });

  it('passes backend X5 warning records through for display without deriving new warnings', () => {
    const warnings = [
      warning,
      { kind: 'pair-boundary', address: 'pair-1-2', message: 'Pair threshold is at its boundary.' },
      { kind: 'energy', address: 'frag-2/range-2', message: 'Best-fit energy does not pass.' },
    ];
    const result = exportReadiness(readinessInput({ admission: admission({ warnings }) }));
    expect(result.warnings).toEqual(warnings);
    expect(result.requiresAcknowledgement).toBe(true);
  });

  it('requires one explicit acknowledgment when backend warnings exist', () => {
    const withWarning = admission({ warnings: [warning] });
    expect(exportNeedsAcknowledgement(withWarning)).toBe(true);
    expect(exportWarningAcknowledged(withWarning, null)).toBe(false);
    expect(exportWarningAcknowledged(withWarning, digest)).toBe(true);
    expect(exportNeedsAcknowledgement(admission())).toBe(false);

    const unacknowledged = exportReadiness(readinessInput({ admission: withWarning }));
    expect(unacknowledged.enabled).toBe(true);
    expect(unacknowledged.requiresAcknowledgement).toBe(true);
  });

  it('binds acknowledgment to the digest while preserving it across change and revert', () => {
    const withWarning = admission({ warnings: [warning] });
    const changedDigest = 'd'.repeat(64);
    expect(exportWarningAcknowledged(withWarning, digest)).toBe(true);
    expect(exportWarningAcknowledged({ ...withWarning, digest: changedDigest }, digest)).toBe(false);
    // A revert restores the exact evidence identity and therefore its local acknowledgment.
    expect(exportWarningAcknowledged({ ...withWarning, digest }, digest)).toBe(true);
  });

  it('marks an existing downloaded digest stale only when the current admission digest differs', () => {
    expect(exportIsStale(digest, digest)).toBe(false);
    expect(exportIsStale(digest, 'd'.repeat(64))).toBe(true);
    expect(exportIsStale(null, 'd'.repeat(64))).toBe(false);
  });

  it('does not require acknowledgment for a snapshot with no X5 warnings', () => {
    const sentinelOnlyAdmission = admission({ warnings: [] });
    expect(exportNeedsAcknowledgement(sentinelOnlyAdmission)).toBe(false);
    expect(exportWarningAcknowledged(sentinelOnlyAdmission, null)).toBe(true);
  });
});
