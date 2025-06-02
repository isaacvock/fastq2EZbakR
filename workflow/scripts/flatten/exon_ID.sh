#!/bin/bash
set -euo pipefail

input=$1
output=$2

awk -v FS="\t" '
  $3 != "exonic_part" {print $0}
  $3 == "exonic_part" {
    split($9, atr, "\"")
    $0 = $0 "; exon_id \"E" atr[2] atr[6] "\";"
    print $0
  }
' "$input" > "$output"