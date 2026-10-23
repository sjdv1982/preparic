#!/usr/bin/env bash
# Native Bash assertions; run all cases or select empty/dispatch/failure/real.
set -euo pipefail

test_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
project=$(cd -- "$test_dir/../.." && pwd)
script="$project/range-pdbs-aareduce.sh"
fixtures="$test_dir/data/1b7f/pdbs"
scratch=$(mktemp -d)
trap 'rm -rf -- "$scratch"' EXIT

fail() { printf '%s\n' "$*" >&2; exit 1; }

# The fake python command isolates the Bash script’s behavior: 
#  it records arguments to verify the RNA/protein flags 
#  and can deliberately fail to check that the script stops and preserves the exit code.
mock_python() {
    mkdir -p bin
    cat > bin/python <<'MOCK'
#!/usr/bin/env bash
set -euo pipefail
printf '%s\n' "$*" >> calls
if [[ ${FAIL_ON:-} == "$2" ]]; then exit 23; fi
cp -- "$2" "$3"
MOCK
    chmod +x bin/python
    export PATH="$PWD/bin:$PATH"
}

run_case() (
    local case_name=$1
    mkdir -p "$scratch/$case_name/pdbs"
    cd -- "$scratch/$case_name"
    case "$case_name" in
        empty)
            # Non-range files must not be passed to aareduce.
            printf 'ignored\n' > pdbs/other.pdb
            mock_python
            bash "$script"
            [[ ! -e calls ]] || fail 'aareduce ran for empty input'
            [[ $(cat RESULT/no-ranges.txt) == 'No range PDB files were provided.' ]] || fail 'missing empty-input marker'
            [[ $(find RESULT -type f | wc -l) -eq 1 ]] || fail 'unexpected empty-input output'
            ;;
        dispatch)
            mock_python
            for n in 1 2 10; do printf 'PDB %s\n' "$n" > "pdbs/range-$n.pdb"; done
            printf 'ignored\n' > pdbs/other.pdb
            bash "$script"
            cat > expected-calls <<'CALLS'
aareduce/aareduce.py pdbs/range-1.pdb RESULT/range-1-aa.pdb --rna --heavy --nalib --no-renumbering
aareduce/aareduce.py pdbs/range-10.pdb RESULT/range-10-aa.pdb --heavy --pdb2pqr --no-renumbering
aareduce/aareduce.py pdbs/range-2.pdb RESULT/range-2-aa.pdb --heavy --pdb2pqr --no-renumbering
CALLS
            diff -u expected-calls calls
            for n in 1 2 10; do cmp "pdbs/range-$n.pdb" "RESULT/range-$n-aa.pdb"; done
            [[ $(find RESULT -type f | wc -l) -eq 3 ]] || fail 'unexpected output files'
            # The RNA rule is based on the filename, not iteration position.
            rm pdbs/range-1.pdb calls
            bash "$script"
            [[ $(head -n 1 calls) == *'--heavy --pdb2pqr --no-renumbering' ]] || fail 'protein-only input used RNA flags'
            ;;
        failure)
            mock_python
            for n in 1 2 3; do printf 'PDB\n' > "pdbs/range-$n.pdb"; done
            export FAIL_ON=pdbs/range-2.pdb
            local status=0
            bash "$script" || status=$?
            [[ $status -eq 23 ]] || fail "expected reducer exit 23, got $status"
            [[ $(wc -l < calls) -eq 2 ]] || fail 'continued after reduction failed'
            [[ ! -e RESULT/range-3-aa.pdb ]] || fail 'reduced a later range after failure'
            ;;
        real)
            cp "$fixtures"/*.pdb pdbs/
            ln -s "$project/aareduce" aareduce
            # Require real PDB2PQR to rebuild a missing internal backbone O.
            awk '!(/^ATOM/ && substr($0,13,4)==" O  " && substr($0,23,4)+0==4)' \
                pdbs/range-2.pdb > missing.pdb
            mv missing.pdb pdbs/range-2.pdb
            bash "$script"
            [[ $(find RESULT -type f | wc -l) -eq 3 ]] || fail 'expected three reduced PDBs'
            for n in 1 2 3; do
                local input="pdbs/range-$n.pdb" output="RESULT/range-$n-aa.pdb"
                [[ -s $output ]] || fail "missing $output"
                # Compare chain, residue number and insertion code in file order.
                awk '/^ATOM/{id=substr($0,22,6);if(id!=prev)print id;prev=id}' "$input" > input-ids
                awk '/^ATOM/{id=substr($0,22,6);if(id!=prev)print id;prev=id}' "$output" > output-ids
                diff -u input-ids output-ids
                awk '/^ATOM/{n++;name=substr($0,13,4);gsub(/ /,"",name);if(name ~ /^[0-9]*[HD]/)exit 1} END{if(!n)exit 1}' "$output" \
                    || fail "expected heavy atoms in $output"
            done
            awk '/^ATOM/ && substr($0,13,4)==" O  " && substr($0,23,4)+0==4{found=1} END{exit !found}' \
                RESULT/range-2-aa.pdb || fail 'PDB2PQR did not rebuild the missing O'
            ;;
        *) fail "Unknown test case: $case_name" ;;
    esac
    printf 'PASS %s\n' "$case_name"
)

if (($#)); then
    run_case "$1"
else
    for case_name in empty dispatch failure real; do run_case "$case_name"; done
fi
