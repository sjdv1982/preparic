import type { ChainInfo, ChainKey } from './types';

export interface ChainUnitInput extends ChainKey {
  isIdentity: boolean;
  seqIds: Iterable<number>;
}

/** Name copies of an asym id deterministically, retaining instanceId as identity. */
export function nameChains(units: ChainUnitInput[]): ChainInfo[] {
  const groups = new Map<string, ChainUnitInput[]>();
  for (const unit of units) {
    const group = groups.get(unit.labelAsymId) ?? [];
    group.push(unit);
    groups.set(unit.labelAsymId, group);
  }

  const result: ChainInfo[] = [];
  for (const [labelAsymId, copies] of groups) {
    const identityIndex = copies.findIndex(copy => copy.isIdentity);
    copies.forEach((copy, index) => {
      const name = copies.length === 1
        ? labelAsymId
        : copy.isIdentity
          ? labelAsymId
          : `${labelAsymId}-${index + (identityIndex >= 0 && index > identityIndex ? 0 : 1)}`;
      const seqIds = new Set([...copy.seqIds].filter(Number.isFinite));
      const sorted = [...seqIds].sort((a, b) => a - b);
      result.push({
        name,
        labelAsymId,
        instanceId: copy.instanceId,
        seqIds,
        seqMin: sorted[0] ?? 0,
        seqMax: sorted.at(-1) ?? 0,
      });
    });
  }
  return result;
}

export function findChain(chains: ChainInfo[], labelAsymId: string, instanceId: string): ChainInfo | undefined {
  return chains.find(chain => chain.labelAsymId === labelAsymId && chain.instanceId === instanceId);
}
