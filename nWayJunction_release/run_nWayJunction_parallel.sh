#!/bin/bash
set -euo pipefail

INPUT_DIR="dotbracket_files"
OUT_DIR="output/single_records"
FAIL_FILE="logs_from_last_update/failed_nWayJunction.txt"
JOBS=4

touch "$FAIL_FILE"
FILE_LIST=$(find "$INPUT_DIR" -name "*.dbn")

printf '%s\n' "$FILE_LIST" | parallel --jobs "$JOBS" '
  id=$(basename {} -2D-dotbracket.dbn)
  out='"$OUT_DIR"'/$id.xml

  if grep -qx "$id" '"$FAIL_FILE"' 2>/dev/null; then
    echo "Skipping $id (previous failure)"

  elif [ ! -f "$out" ]; then
    echo -n "nWayJunction: processing $id... "

    if python3 ./nWayJunction_release/main.py SINGLE '"$INPUT_DIR"'/$id-2D-dotbracket.dbn
    then
      echo "OK"
    else
      echo "FAILED"
      echo "$id" >> '"$FAIL_FILE"'
    fi
  fi
'
