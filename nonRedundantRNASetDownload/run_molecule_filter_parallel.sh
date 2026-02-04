#!/bin/bash
set -euo pipefail

INPUT_DIR="PDB_files_raw"
OUT_DIR="PDB_files"
FAIL_FILE="failed_molecule-filter.txt"
JOBS=4

mkdir -p "$OUT_DIR"
touch "$FAIL_FILE"
FILE_LIST=$(find "$INPUT_DIR" -name "*.cif")

printf '%s\n' "$FILE_LIST" | parallel --jobs "$JOBS" '
  id=$(basename {} .cif)
  out='"$OUT_DIR"'/$id.cif

  if grep -qx "$id" '"$FAIL_FILE"' 2>/dev/null; then
    echo "Skipping $id (previous failure)"

  elif [ ! -f "$out" ]; then
    echo -n "molecule-filter: processing $id... "
    if molecule-filter PDB_files_raw/$id.cif \
         --filter-by-poly-types polyribonucleotide \
         > "$out"
    then
      echo "OK"
    else
      echo "FAILED for $id"
      echo "$id" >> '"$FAIL_FILE"'
    fi
  fi
'
