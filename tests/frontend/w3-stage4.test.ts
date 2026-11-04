import { describe, expect, it } from 'vitest';
import {
  canContinueStage4, canEnterStage4, clipForSelectedAnchors, copyComputedThresholds,
  defaultStackingSelection, makeThresholdSkeleton, selectedCandidateFragments, selectedCandidateGroups, thresholdsCompleteForClip,
  validateThresholdEdits, validateThresholdValue,
} from '../../src/lib/stage4';
import { stackingKey } from '../../src/lib/stackings';
import { createRangeRow, validateRows } from '../../src/lib/ranges';
import { nameChains } from '../../src/lib/chains';
import { deferredWriteIsCurrent } from '../../src/lib/project-state';

const rna = { chain: 'A', start: 101, end: 104 };
const basis = { pdbcode: '1B7F', rna };
const residues = [...'ACGU'].map((base, i) => ({ chain: 'A', resid: 101 + i, insertion_code: '', base }));
const sequences = ['AC', 'CG', 'GU'];
const fragments: any = {
  schema: 1, echo: basis, residues,
  fragments: sequences.map((sequence, i) => ({ number: i + 1, key: `frag-${i + 1}-${sequence}`, sequence, residues: residues.slice(i, i + 2) })),
  adjacent_pairs: [{ down: 1, up: 2, sequence: 'ACG' }, { down: 2, up: 3, sequence: 'CGU' }],
};
const clip = { schema: 1 as const, basis, clip1: 0, clip2: 0 };

describe('W3 deferred clip writes', () => {
  const current = {
    stageAtIssue: 4 as const,
    currentStage: 4 as const,
    live: true,
    revisionAtIssue: 8,
    currentRevision: 8,
    connectionAtIssue: 3,
    currentConnection: 3,
  };

  it('admits an unchanged live clip intent', () => {
    expect(deferredWriteIsCurrent(current)).toBe(true);
  });

  it.each([
    { currentStage: 3 as const },
    { live: false },
    { currentRevision: 9 },
    { currentConnection: 4 },
  ])('drops a clip intent after %j', change => {
    expect(deferredWriteIsCurrent({ ...current, ...change })).toBe(false);
  });
});

function stacking(resid: number, distance = 5, slip = 59.999, tilt = 89) {
  return { rna_resname: 'C', rna_resid: resid, protein_range_pdb_id: 2, protein_resid: 30,
    distance, slip_angle: slip, tilt_angle: tilt, rna_centroid: [0, 0, distance], protein_centroid: [0, 0, 0] } as any;
}

function filled() {
  const table = makeThresholdSkeleton(fragments, '1B7F');
  if (!table) throw new Error('Expected a threshold skeleton for the valid inventory.');
  table.pairs.forEach((row: any) => { row.ovrmsd = '0.25'; row.crmsd = '0.50'; });
  table.bases.forEach((row: any) => { row.nuc1_ovrmsd = '0.75'; row.nuc2_ovrmsd = '1.00'; });
  return table;
}

describe('W3 stage 3 selection (OC P3.1)', () => {
  it('uses inclusive distance 5 and strict displacement 60 without a tilt condition', () => {
    const rows = [stacking(101), stacking(102, 5.001), stacking(103, 5, 60), stacking(104, 5, 0, 90)];
    expect(defaultStackingSelection(rows, {})).toEqual({
      [stackingKey(rows[0])]: true, [stackingKey(rows[1])]: false,
      [stackingKey(rows[2])]: false, [stackingKey(rows[3])]: true,
    });
  });

  it('merges retained decisions by key and drops absent keys', () => {
    const rows = [stacking(101), stacking(102)];
    expect(defaultStackingSelection(rows, { [stackingKey(rows[0])]: false, orphan: true })).toEqual({
      [stackingKey(rows[0])]: false, [stackingKey(rows[1])]: true,
    });
  });
});

describe('W3 RNA range validity agrees with extraction (OC P2.1/D7/D15)', () => {
  it('rejects an unobserved internal RNA residue even when both endpoints are observed', () => {
    const chains = [{ name: 'A', labelAsymId: 'A', instanceId: '1', seqIds: new Set([101, 102, 104, 105]), seqMin: 101, seqMax: 105 }];
    const row = { ...createRangeRow(), chain: 'A', start: '101', end: '105' };
    const issues = validateRows([row], chains);
    expect(issues.length).toBeGreaterThan(0);
    expect(JSON.stringify(issues)).toContain('103');
  });

  it('rejects a chemical O3′–P break recorded from the selected assembly', () => {
    const chains = [{ name: 'A', labelAsymId: 'A', instanceId: '1',
      seqIds: new Set([101, 102, 103, 104]), seqMin: 101, seqMax: 104,
      rnaBreakAfter: new Set([102]) }];
    const row = { ...createRangeRow(), chain: 'A', start: '101', end: '104' };
    expect(validateRows([row], chains).map(issue => issue.message)).toContain(
      'RNA range has a chemical chain break between residues 102 and 103.',
    );
  });

  it('accepts the selected assembly copy name and does not apply the RNA gap rule to protein rows', () => {
    const chains = nameChains([
      { labelAsymId: 'A', instanceId: '1', isIdentity: true, seqIds: [101, 102, 103, 104] },
      { labelAsymId: 'A', instanceId: '2', isIdentity: false, seqIds: [101, 102, 103, 104] },
      { labelAsymId: 'A', instanceId: '3', isIdentity: false, seqIds: [101, 102, 103, 104] },
      { labelAsymId: 'C', instanceId: '1', isIdentity: true, seqIds: [10, 11, 13, 14] },
    ]);
    expect(chains.map(row => row.name)).toEqual(['A', 'A-1', 'A-2', 'C']);
    expect(validateRows([
      { ...createRangeRow(), chain: 'A-2', start: '101', end: '104' },
      { ...createRangeRow(), chain: 'C', start: '10', end: '14' },
    ], chains)).toEqual([]);
  });
});

describe('W3 stage 4 decisions (OC P4.3, D4/D5)', () => {
  it('creates a keyed empty skeleton before any fit/computed values exist', () => {
    const table = makeThresholdSkeleton(fragments, '1B7F');
    expect(table).toEqual({ schema: 1, basis, active: true,
      pairs: fragments.adjacent_pairs.map((row: any) => ({ ...row, ovrmsd: null, crmsd: null })),
      bases: sequences.map((sequence, i) => ({ fragment: i + 1, sequence, nuc1_ovrmsd: null, nuc2_ovrmsd: null })),
    });
    expect(validateThresholdEdits(table, fragments)).toBe(true);
    expect(thresholdsCompleteForClip(table, 3, clip)).toBe(false);
  });

  it.each([null, '0.01', '0.1', '1', '1.00', '10.01'])('accepts a blank or finite positive two-decimal value %s', value => {
    expect(validateThresholdValue(value)).toBe(true);
  });

  it.each(['0', '-0.1', '0.001', '1e-2', 'NaN', 'Infinity', 0.25, true])('rejects %s without rounding or coercion', value => {
    expect(validateThresholdValue(value)).toBe(false);
  });

  it('requires only kept base/pair cells and closes when an empty row is uncovered', () => {
    const table = filled();
    table.bases[0].nuc1_ovrmsd = null;
    table.bases[0].nuc2_ovrmsd = null;
    table.pairs[0].ovrmsd = null;
    table.pairs[0].crmsd = null;
    expect(thresholdsCompleteForClip(table, 3, { ...clip, clip1: 1 })).toBe(true);
    expect(thresholdsCompleteForClip(table, 3, clip)).toBe(false);
  });

  it('refuses schema/key errors without partially accepting the table', () => {
    expect(validateThresholdEdits({ ...filled(), schema: 2 }, fragments)).toBe(false);
    expect(validateThresholdEdits({ ...filled(), schema: true }, fragments)).toBe(false);
    const table = filled();
    table.bases.push({ ...table.bases[0] });
    expect(validateThresholdEdits(table, fragments)).toBe(false);
  });

  it('copies only explicit rounded values after confirmation, preserving decision basis', () => {
    const table = filled();
    table.bases[0].nuc1_ovrmsd = '1.25';
    const computed: any = { schema: 1, echo: { pdb_code: '1b7f', rna, resids: residues },
      pairs: table.pairs.map((row: any) => ({ ...row, ovrmsd: '0.26', crmsd: '0.51', ovrmsd_raw: 0.25, crmsd_raw: 0.5 })),
      bases: table.bases.map((row: any) => ({ ...row, nuc1_ovrmsd: '0.76', nuc2_ovrmsd: '1.01', nuc1_ovrmsd_raw: 0.75, nuc2_ovrmsd_raw: 1 })),
    };
    const before = structuredClone(table);
    expect(copyComputedThresholds(table, computed, false)).toEqual(before);
    expect(copyComputedThresholds(table, null, true)).toEqual(before);
    const copied = copyComputedThresholds(table, computed, true);
    expect(copied.basis).toEqual(basis);
    expect(copied.bases[0].nuc1_ovrmsd).toBe('0.76');
    expect(copied.pairs[0].ovrmsd).toBe('0.26');
    expect(copied.bases[0]).not.toHaveProperty('nuc1_ovrmsd_raw');
    expect(table).toEqual(before);
  });

  it('retains one candidate at either end and minimally restores an excluded anchor', () => {
    const row = stacking(102);
    const key = stackingKey(row);
    expect(selectedCandidateFragments(fragments, [key], [row])).toEqual([1, 2]);
    expect(selectedCandidateGroups(fragments, [key], [row])).toEqual([[1, 2]]);
    for (const saved of [{ ...clip, clip1: 1 }, { ...clip, clip2: 2 }]) {
      expect(clipForSelectedAnchors(saved, fragments, [key], [row])).toEqual(saved);
    }
    const old = { ...clip, clip1: 2 };
    expect(clipForSelectedAnchors(old, fragments, [key], [row])).toEqual({ ...old, clip1: 1 });
    expect(old.clip1).toBe(2);
    const right = { ...stacking(103), rna_resname: 'G' };
    expect(selectedCandidateGroups(fragments, [stackingKey(right)], [right])).toEqual([[2, 3]]);
    expect(clipForSelectedAnchors({ ...clip, clip2: 1 }, fragments, [stackingKey(right)], [right]))
      .toEqual({ ...clip, clip2: 1 });
    expect(clipForSelectedAnchors({ ...clip, clip2: 2 }, fragments, [stackingKey(right)], [right]))
      .toEqual({ ...clip, clip2: 1 });
  });

  it('protects terminal singleton candidates and simultaneous selected anchors', () => {
    const rows = [{ ...stacking(101), rna_resname: 'A' }, { ...stacking(104), rna_resname: 'U' }];
    expect(selectedCandidateGroups(fragments, rows.map(stackingKey), rows)).toEqual([[1], [3]]);
    expect(clipForSelectedAnchors({ ...clip, clip1: 1, clip2: 1 }, fragments,
      [stackingKey(rows[0])], [rows[0]])).toEqual({ ...clip, clip2: 1 });
    expect(clipForSelectedAnchors({ ...clip, clip1: 1, clip2: 1 }, fragments,
      [stackingKey(rows[1])], [rows[1]])).toEqual({ ...clip, clip1: 1 });
    expect(clipForSelectedAnchors({ ...clip, clip1: 1, clip2: 1 }, fragments,
      rows.map(stackingKey), rows)).toEqual(clip);
  });

});

describe('W3 forward gates use current data without waiting for the fit', () => {
  const rows = [stacking(102)];
  const selection = { schema: 1, basis: { pdbcode: '1B7F', biounit: 1, rows: [rna, { chain: 'C', start: 30, end: 40 }] }, active: true, selected: { [stackingKey(rows[0])]: true } };
  const enter = { selection, stackings: rows, stackingsFresh: true, fragmentsFresh: true, rangePdbsFresh: true };

  it('opens 3→4 on fresh reference fragments and held fresh all-atom PDBs with no fit input', () => {
    expect(canEnterStage4(enter)).toBe(true);
  });

  it.each(['stackingsFresh', 'fragmentsFresh', 'rangePdbsFresh'] as const)('closes 3→4 while %s is false', field => {
    expect(canEnterStage4({ ...enter, [field]: false })).toBe(false);
  });

  it('requires exact selection keys and a nonempty selected subset', () => {
    expect(canEnterStage4({ ...enter, selection: { ...selection, selected: {} } })).toBe(false);
    expect(canEnterStage4({ ...enter, selection: { ...selection, selected: { [stackingKey(rows[0])]: false } } })).toBe(false);
    expect(canEnterStage4({ ...enter, selection: { ...selection, selected: { ...selection.selected, orphan: true } } })).toBe(false);
  });

  it('test_continue_4_to_5_accepts_complete_typed_thresholds_while_copy_source_stale', () => {
    const args = { edits: filled(), fragmentCount: 3, clip, selectedCandidates: [[1, 2]], energeticsFresh: true, thresholdValidatorFailed: false };
    expect(canContinueStage4(args)).toBe(true);
    expect(canContinueStage4({ ...args, energeticsFresh: false })).toBe(false);
    expect(canContinueStage4({ ...args, thresholdValidatorFailed: true })).toBe(false);
    expect(canContinueStage4({ ...args, clip: { ...clip, clip1: 1 } })).toBe(true);
    expect(canContinueStage4({ ...args, clip: { ...clip, clip1: 2 } })).toBe(false);
    expect(canContinueStage4({ ...args, edits: makeThresholdSkeleton(fragments, '1B7F') })).toBe(false);
  });
});
