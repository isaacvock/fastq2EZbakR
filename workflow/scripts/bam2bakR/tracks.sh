#!/bin/bash
set -euo pipefail

# Make mutation-stratified browser tracks (.tdf) from pre-sorted coordinate BAMs

cpus=$1
sample=$2
input_counts=$3
input_coord_bam=$4
input_scale=$5
mut_tracks=$6
genome_fasta=$7
WSL_b=$8
normalize=$9
pyscript=${10}
strand=${11}
output_sentinel=${12}

if [ "$strand" = "R" ]; then
    pos_star_strand="str2"
    min_star_strand="str1"
else
    pos_star_strand="str1"
    min_star_strand="str2"
fi

mkdir -p ./results/tracks
touch "$output_sentinel"

if [ "$normalize" = "True" ]; then
    normVal=$(awk -v sam="$sample" '$1 == sam {print $2}' "$input_scale")
else
    normVal='1'
    echo '* Tracks will not be normalized with scale factor. Setting normalization factor to 1'
fi

echo '* Creating files for each level of counting.'
python "$pyscript" \
    -i "$input_counts" \
    --mutType "$mut_tracks" \
    -s "./results/tracks/$sample"

# Per-sample chrom.sizes and temporary working directory prevent concurrent maketdf races
chrom_sizes="./results/tracks/${sample}.chrom.sizes"
igv_tmpdir=$(mktemp -d "./results/tracks/igv_tmp_${sample}_XXXXXX")
trap 'rm -rf "$igv_tmpdir" "$chrom_sizes"' EXIT

echo '* Making .chrom.sizes file'
samtools view -H "$input_coord_bam" \
    | awk -v OFS="\t" '$1 ~ /^@SQ/ {
        split($2, chr, ":")
        split($3, size, ":")
        print chr[2], size[2]
    }' > "$chrom_sizes"

muts=$(echo "$mut_tracks" | tr ',' ' ')

echo '* Making track headers'
for b in $muts; do
    if [ "$b" == "GA" ]; then
        colVal[0]='200,200,200'
        colVal[1]='120,188,230'
        colVal[2]='65,125,195'
        colVal[3]='36,110,182'
        colVal[4]='27,78,165'
        colVal[5]='18,50,120'
    else
        colVal[0]='200,200,200'
        colVal[1]='250,150,150'
        colVal[2]='250,0,0'
        colVal[3]='150,0,0'
        colVal[4]='100,0,0'
        colVal[5]='50,0,0'
    fi

    for count in $(seq 0 5); do
        echo "track type=bedGraph name=\" $sample $b $count plus \" description=\" $sample $b $count positive strand\" visibility=full autoScale=off windowingFunction=mean viewLimits=0:10 color=${colVal[$count]} altColor=10,10,10 priority=20" > "./results/tracks/${sample}.${b}.${count}.pos.bedGraph"
        echo "track type=bedGraph name=\" $sample $b $count minus \" description=\" $sample $b $count minus strand\" visibility=full negateValues=on autoScale=off windowingFunction=mean viewLimits=-10:0 color=10,10,10 altColor=${colVal[$count]} priority=20" > "./results/tracks/${sample}.${b}.${count}.min.bedGraph"
    done
done

append_star_signal() {
    local mut=$1
    local count=$2
    local sign=$3
    local star_strand=$4
    local out_strand=$5
    local in_bg="./results/tracks/${sample}_${mut}_${count}_Signal.Unique.${star_strand}.out.bg"
    local out_bg="./results/tracks/${sample}.${mut}.${count}.${out_strand}.bedGraph"

    if [ ! -s "$in_bg" ]; then
        echo "Warning: STAR signal file $in_bg is missing or empty; leaving $out_bg with no signal values." >&2
        return 0
    fi

    awk -v norm="$normVal" -v sign="$sign" 'BEGIN {OFS="\t"} $1 !~ /^track/ && NF >= 4 {print $1, $2, $3, sign * norm * $4}' "$in_bg" >> "$out_bg"
}

append_track_signals() {
    for b in $muts; do
        for count in $(seq 0 5); do
            append_star_signal "$b" "$count" 1 "$pos_star_strand" "pos"
            append_star_signal "$b" "$count" -1 "$min_star_strand" "min"
        done
    done
}

ensure_bedgraph_data() {
    local bg=$1

    if ! awk '$1 !~ /^track/ && NF >= 4 {found=1; exit} END {exit found ? 0 : 1}' "$bg"; then
        awk 'NF >= 2 {print $1 "\t0\t1\t0"; exit}' "$chrom_sizes" >> "$bg"
    fi
}

ensure_all_bedgraphs_have_data() {
    if [ ! -s "$chrom_sizes" ]; then
        echo "Error: chromosome sizes file $chrom_sizes is missing or empty." >&2
        exit 1
    fi

    for b in $muts; do
        for count in $(seq 0 5); do
            ensure_bedgraph_data "./results/tracks/${sample}.${b}.${count}.pos.bedGraph"
            ensure_bedgraph_data "./results/tracks/${sample}.${b}.${count}.min.bedGraph"
        done
    done
}

echo '* Filtering reads'
# shellcheck disable=SC2086
parallel -j 1 samtools view -@ "$cpus" \
                            -b \
                            -N "./results/tracks/${sample}_{1}_{2}_reads.txt" \
                            -o "./results/tracks/${sample}_{1}_{2}.bam" \
                            "$input_coord_bam" ::: $muts ::: $(seq 0 5)

if [ "$WSL_b" = "True" ]; then
    echo "Running STAR iteratively"
    for b in $muts; do
        for iter in $(seq 0 5); do
            STAR \
                --runMode inputAlignmentsFromBAM \
                --inputBAMfile "./results/tracks/${sample}_${b}_${iter}.bam" \
                --outWigType bedGraph \
                --outWigNorm None \
                --outWigStrand Stranded \
                --outFileNamePrefix "./results/tracks/${sample}_${b}_${iter}_"
        done
    done
else
    echo '* Making tracks with STAR'
    # shellcheck disable=SC2086
    parallel -j "$cpus" STAR \
                            --runMode inputAlignmentsFromBAM \
                            --inputBAMfile "./results/tracks/${sample}_{1}_{2}.bam" \
                            --outWigType bedGraph \
                            --outWigNorm None \
                            --outWigStrand Stranded \
                            --outFileNamePrefix "./results/tracks/${sample}_{1}_{2}_" ::: $muts ::: $(seq 0 5)
fi

echo '* Taking only unique components of tracks'
append_track_signals
rm -f ./results/tracks/"$sample"_*.bg

echo '* Make tdf files from tracks'
ensure_all_bedgraphs_have_data
rm -f ./results/tracks/"$sample".*.tdf

tracks_abs=$(cd ./results/tracks && pwd)
chrom_sizes_abs="${tracks_abs}/${sample}.chrom.sizes"

(
    cd "$igv_tmpdir"
    # shellcheck disable=SC2086
    parallel -j "$cpus" igvtools toTDF \
                            -f mean,max \
                            "${tracks_abs}/${sample}.{1}.{2}.{3}.bedGraph" \
                            "${tracks_abs}/${sample}.{1}.{2}.{3}.tdf" \
                            "$chrom_sizes_abs" ::: $muts ::: $(seq 0 5) ::: pos min
)

rm -f ./results/tracks/"$sample"_*.bam
rm -f ./results/tracks/"$sample"_*_reads.txt
rm -f ./results/tracks/"$sample".*.bedGraph
rm -f ./results/tracks/"$sample"_*.bai
rm -f ./results/tracks/"$sample"_*.out

echo "* TDF track files created for sample $sample."
