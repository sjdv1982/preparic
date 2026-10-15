import { MolScriptBuilder as MS } from 'molstar/lib/mol-script/language/builder.js';
import type { Expression } from 'molstar/lib/mol-script/language/expression.js';
import { NucleicBackboneAtoms, ProteinBackboneAtoms } from 'molstar/lib/mol-model/structure/model/types.js';
import type { ChainKey, RangeSpec } from '../types';

const labelAsymId = MS.struct.atomProperty.macromolecular.label_asym_id();
const instanceId = MS.struct.atomProperty.core.instanceId();
const labelSeqId = MS.struct.atomProperty.macromolecular.label_seq_id();
const labelAtomId = MS.struct.atomProperty.macromolecular.label_atom_id();
// Keep CA and C1' (or C1*) as bond anchors; hide the rest of the backbone.
const excludedBackboneAtoms = [
  ...ProteinBackboneAtoms,
  ...NucleicBackboneAtoms,
].filter(atomName => atomName !== 'CA' && atomName !== "C1'" && atomName !== 'C1*');

export function chainExpression(chain: ChainKey): Expression {
  return MS.struct.generator.atomGroups({
    'chain-test': MS.core.logic.and([
      MS.core.rel.eq([labelAsymId, chain.labelAsymId]),
      MS.core.rel.eq([instanceId, chain.instanceId]),
    ]),
  });
}

export function rangeExpression(range: RangeSpec): Expression {
  return MS.struct.generator.atomGroups({
    'chain-test': MS.core.logic.and([
      MS.core.rel.eq([labelAsymId, range.labelAsymId]),
      MS.core.rel.eq([instanceId, range.instanceId]),
    ]),
    'residue-test': MS.core.rel.inRange([labelSeqId, range.start, range.end]),
  });
}

export function sidechainRangeExpression(range: RangeSpec): Expression {
  return MS.struct.modifier.intersectBy({
    0: rangeExpression(range),
    by: MS.struct.generator.atomGroups({
      'atom-test': MS.core.logic.not([
        MS.core.set.has([MS.set(...excludedBackboneAtoms), labelAtomId]),
      ]),
    }),
  });
}

export function unionExpression(expressions: Expression[]): Expression {
  const nonEmpty = expressions.length ? expressions : [MS.struct.generator.empty()];
  return MS.struct.modifier.union([
    nonEmpty.length === 1
      ? nonEmpty[0]
      : MS.struct.combinator.merge(nonEmpty.map(expression => MS.struct.modifier.union([expression]))),
  ]);
}

export function complementExpression(expression: Expression): Expression {
  return MS.struct.modifier.exceptBy({ 0: MS.struct.generator.all(), by: expression });
}
