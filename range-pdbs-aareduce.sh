#!/usr/bin/env bash
set -euo pipefail

mkdir -p RESULT
shopt -s nullglob
pdb_files=(pdbs/range-*.pdb)

if ((${#pdb_files[@]} == 0)); then
    printf 'No range PDB files were provided.\n' > RESULT/no-ranges.txt
    exit 0
fi

for pdb in "${pdb_files[@]}"; do
    name=${pdb##*/}
    output="RESULT/${name%.pdb}-aa.pdb"

    if [[ $name == range-1.pdb ]]; then
        python aareduce/aareduce.py "$pdb" "$output" \
            --rna --heavy --nalib --no-renumbering
    else
        python aareduce/aareduce.py "$pdb" "$output" \
            --heavy --pdb2pqr --no-renumbering
    fi
done
