import { describe, expect, it } from 'vitest';

describe('Mol* viewer module bootstrap', () => {
  it('defers custom shape-provider lookup until the viewer is initialized', async () => {
    const { MolViewer } = await import('../../src/lib/molstar/viewer');
    expect(typeof MolViewer.create).toBe('function');
  }, 20000);
});
