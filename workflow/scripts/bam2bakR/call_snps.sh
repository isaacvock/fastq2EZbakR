#!/bin/bash
set -euo pipefail

cpus=$1
nsamps=$2
output_txt=$3
output_vcf=$4
mpileup_options="$5"
call_options="$6"
genome_fasta=$7

if [[ "$nsamps" -gt 0 ]]; then
    shift 8
    control_samples=("$@")
    mkdir -p ./results/snps
    bam_list=$(mktemp ./results/snps/bam_list.XXXXXX)
    trap 'rm -f "$bam_list"' EXIT

    printf "%s\n" "${control_samples[@]}" > "$bam_list"

    # shellcheck disable=SC2086
    bcftools mpileup --threads "$cpus" \
                     -f "$genome_fasta" \
                     -b "$bam_list" \
                     $mpileup_options \
                     -Ou \
    | bcftools call --threads "$cpus" $call_options -mv -Ov -o "$output_vcf"

    awk '$1 !~ /^#/ && length($4) == 1 {
            if (length($5) == 1) {
                print $4":"$5":"$1":"$2
            } else if (length($5) == 3 && $5 ~ /,/) {
                split($5, mut, ",")
                print $4":"mut[1]":"$1":"$2
                print $4":"mut[2]":"$1":"$2
            }
        }' "$output_vcf" \
        | LC_COLLATE=C sort -u > "$output_txt"

    echo '* SNPs called and snp.txt generated'
else
    touch "$output_txt"
    touch "$output_vcf"
fi
