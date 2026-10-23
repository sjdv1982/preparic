# 1B7F backend fixtures

`structure.cif` is the 1B7F mmCIF saved by an earlier local deployment session
(`/tmp/preparic-cif-probe`). `deployment-ranges.yaml` is a copy of the selections
in `deployment-mounts/ranges.yaml` when this suite was created. Both are retained
here so the tests work after those ephemeral files disappear.

`ranges.yaml` shortens those selections for fast real aareduce runs:

- RNA chain A, residues 1–4.
- Protein chain C, residues 2–8.
- Protein chain C, residues 90–96.

`pdbs/range-{1,2,3}.pdb` are fixed extraction snapshots of those selections,
preserving label chain IDs, residue numbers and coordinates from the CIF. They
serve as expected extraction output and independent Bash input. The extraction
tests also assert residue identities and boundaries explicitly. The native Bash
test deletes an internal protein backbone O from a temporary copy and verifies
that the real reducer rebuilds it. Reduced coordinates are not golden files:
the tests check numbering, chains, heavy atoms and protein backbone completeness.
