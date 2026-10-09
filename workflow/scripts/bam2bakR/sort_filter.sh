#!/bin/bash

set -euo pipefail

# Source the paths and variables:
cpus=$1
sample=$2
input=$3
output=$4 # Main bam file that actually matters
layout=${5:-PE}
sc=${6:-False}

# Make sure bam file is query sorted rather than coordinate sorted
    # This is important for mutation counting and other pipeline steps
    # that are parallelized via splitting up the bam file into chunks.
    # Query sorting is required by the custom awk script that splits up 
    # the bam file because the read pairs need to be kept together in the
    # bam file chunks.
# Mate repair and mate-unmapped filtering apply only to biological pairs.
if [[ "$layout" == "PE" && "$sc" != "True" ]]; then
    samtools sort -@ "$cpus" -n "$input" | \
    samtools fixmate -@ "$cpus" - - | \
    samtools view -@ "$cpus" -b -h -q 2 -F 0x4 -F 0x8 -F 0x100 -F 0x200 -F 0x400 -F 0x800 - > "$output"
else
    # Duplicate-flagged reads remain observations in single-cell mode.
    exclude=0xF04
    if [[ "$sc" == "True" ]]; then
        exclude=0xB04
    fi
    samtools sort -@ "$cpus" -n "$input" | \
    samtools view -@ "$cpus" -b -h -q 2 -F "$exclude" - > "$output"
fi

if [[ "$sc" == "True" ]]; then
    # Name sorting makes collisions adjacent. Fail before mutation counting,
    # where repeated qnames would otherwise be interpreted as biological mates.
    samtools view -@ "$cpus" "$output" | awk '
        $2 % 2 == 1 {print "Single-cell mode requires unpaired cDNA alignments: " $1 > "/dev/stderr"; exit 1}
        NR > 1 && $1 == previous {print "Repeated primary single-cell read identity: " $1 > "/dev/stderr"; exit 1}
        {previous=$1}'
fi

echo "* Reads filtered for sample $sample"
