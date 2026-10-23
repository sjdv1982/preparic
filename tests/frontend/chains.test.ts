import { describe, expect, it } from 'vitest';
import { nameChains, type ChainUnitInput } from '../../src/lib/chains';

function copy(labelAsymId: string, instanceId: string, isIdentity: boolean, seqIds = [1, 2, 3]): ChainUnitInput {
  return { labelAsymId, instanceId, isIdentity, seqIds };
}

describe('nameChains', () => {
  it('keeps a single copy under its asym id', () => {
    expect(nameChains([copy('A', '1_555', true)])).toMatchObject([{ name: 'A', labelAsymId: 'A', instanceId: '1_555' }]);
  });

  it('names the identity copy A and later copies in order', () => {
    const chains = nameChains([
      copy('A', 'ASM-1', true),
      copy('A', 'ASM-2', false),
      copy('A', 'ASM-3', false),
    ]);
    expect(chains.map(chain => chain.name)).toEqual(['A', 'A-1', 'A-2']);
  });

  it('names copies without an identity operator from A-1', () => {
    const chains = nameChains([copy('A', 'ASM-4', false), copy('A', 'ASM-5', false)]);
    expect(chains.map(chain => chain.name)).toEqual(['A-1', 'A-2']);
  });

  it('keeps mixed chains and an identity copy that is not first', () => {
    const chains = nameChains([
      copy('A', 'ASM-2', false, [3, 5]),
      copy('B', 'ASM-1', true, [4, 6]),
      copy('A', 'ASM-1', true, [1, 2]),
      copy('A', 'ASM-3', false, [7]),
    ]);
    expect(chains.map(chain => chain.name)).toEqual(['A-1', 'A', 'A-2', 'B']);
    expect([...chains[0].seqIds]).toEqual([3, 5]);
    expect(chains[0].seqMin).toBe(3);
    expect(chains[0].seqMax).toBe(5);
  });
});
