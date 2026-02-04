#!/bin/bash
set -euo pipefail

INPUT_DIR="PDB_files"
OUT_DIR="dotbracket_files"
FAIL_FILE="logs_from_last_update/failed_RNApdbee.txt"
JOBS=4

mkdir -p "$OUT_DIR"
touch "$FAIL_FILE"

FILE_LIST=$(find "$INPUT_DIR" -name "*.cif")

printf '%s\n' "$FILE_LIST" | parallel --jobs "$JOBS" '
  infile={}
  id=$(basename "$infile" .cif)
  out='"$OUT_DIR"'/$id-2D-dotbracket.dbn

  if grep -qx "$id" '"$FAIL_FILE"' 2>/dev/null; then
    echo "Skipping $id (previous failure)"

  elif [ ! -f "$out" ]; then
    echo -e "RNApdbee: processing $id... "

    tmpdir=$(mktemp -d /tmp/rnapdbee.XXXXXX)

    if /opt/rnapdbee-standalone-old/rnapdbee \
         -i "$infile" \
         -o "$tmpdir" \
         -a DSSR
    then
      if mv "$tmpdir/0/strands.dbn" "$out"; then
        echo "OK"
      else
        echo "FAILED (no output file)"
        echo "$id" >> '"$FAIL_FILE"'
      fi
    else
      echo "FAILED (RNApdbee error)"
      echo "$id" >> '"$FAIL_FILE"'
    fi

    rm -rf "$tmpdir"
  fi
'
