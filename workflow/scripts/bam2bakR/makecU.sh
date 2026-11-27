#!/bin/bash
set -euo pipefail

cpus=$1
pos_cutoff=$2
mutposout=$3
mutposfilter=$4
high_cutoff=$5

shift 5
cu_files=("$@")

if [ "${#cu_files[@]}" -eq 0 ]; then
    echo "Error: no _cU.csv.gz input files supplied to makecU.sh" >&2
    exit 1
fi

{
    echo "sample,rname,gloc,tp,trials,n"
    for f in "${cu_files[@]}"; do
        sname=$(basename "$f" _cU.csv.gz)
        pigz -d -c "$f" | awk -v s="$sname" 'NR > 1 { print s "," $0 }'
    done
} | pigz -p "$cpus" > "$mutposout"

pigz -d -c "$mutposout" \
    | awk -F "," \
        -v cutoff="$pos_cutoff" \
        -v upper="$high_cutoff" \
        'NR == 1 || ($(NF-1) >= cutoff && $(NF-1) <= upper) { print }' \
    | pigz -p "$cpus" > "$mutposfilter"

echo "** site-specific mutation file created: $mutposout"
