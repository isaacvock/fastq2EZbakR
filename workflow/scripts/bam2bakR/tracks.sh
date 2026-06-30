#!/bin/bash
set -euo pipefail
#
# To make filtered bam files and tracks with different numbers of mutations

## Set normalization value
#normVal='1'

## Dependencies
# Python
# Pysam
# Igvtools
# parallel
# Samtools
# Star

cpus=$1
sample=$2
input1=$3
input2=$4
input3=$5
mut_tracks=$6
genome_fasta=$7
WSL_b=$8
normalize=$9
pyscript=${10}
strand=${11}
output=${12}

if [ "$strand" = "R" ]; then
    pos_star_strand="str2"
    min_star_strand="str1"
else
    pos_star_strand="str1"
    min_star_strand="str2"
fi


    # Create ./results/tracks/
    touch "$output"

	if [ "$normalize" = "True" ]; then
		normVal=$(awk -v sam=$sample '$1 == sam {print $2}' ./results/normalization/scale)
	else
		normVal='1'
		echo '* Tracks will not be normalized with scale factor. Setting normalization factor to 1'
	fi



# Test if count files exist

    echo '* Creating files for each level of counting.'

    python $pyscript \
        -i $input1 \
        --mutType $mut_tracks \
        -s ./results/tracks/$sample


# Sort the starting .bam file by coordinate
echo '* Sorting bam files by coordiantes'

    samtools sort -@ "$cpus" -o ./results/tracks/"$sample"_sort.bam "$input2"
    samtools index -@ "$cpus" ./results/tracks/"$sample"_sort.bam


# Make .chrom.sizes file from .bam file header (alternative to .genome file for toTDF)

echo '* Making .chrom.sizes file'

    chrom_sizes=./results/tracks/$(echo ${genome_fasta##*/} | cut -f 1 -d '.')".chrom.sizes"
    if [ ! -f $chrom_sizes ]; then
            samtools view -H ./results/tracks/"$sample"_sort.bam \
                | awk -v OFS="\t" ' $1 ~ /^@SQ/ {split($2, chr, ":")
                                                 split($3, size, ":")
                                                 print chr[2], size[2]}' > "$chrom_sizes"
    fi


    muts=$(echo $mut_tracks | tr ',' ' ')

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

        # Make track headers
        for count in $(seq 0 5); do

            echo "track type=bedGraph name=\" $sample $b $count "plus \" description=\" $sample $b $count "positive strand\" visibility=full autoScale=off windowingFunction=mean viewLimits=0:10 color="${colVal[$count]} " altColor=10,10,10 priority=20" > ./results/tracks/"$sample"."$b"."$count".pos.bedGraph
            echo "track type=bedGraph name=\" $sample $b $count "minus \" description=\" $sample $b $count "minus strand\" visibility=full negateValues=on autoScale=off windowingFunction=mean viewLimits=-10:0 color=10,10,10 altColor="${colVal[$count]} " priority=20" > ./results/tracks/"$sample"."$b"."$count".min.bedGraph

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


    # Filter the reads
    echo '* Filtering reads'

    parallel -j 1 samtools view -@ "$cpus" \
                                     -b \
                                     -N ./results/tracks/"$sample"_{1}_{2}_reads.txt \
                                     -o ./results/tracks/"$sample"_{1}_{2}.bam \
                                     ./results/tracks/"$sample"_sort.bam ::: $muts \
                                                        ::: $(seq 0 5)


    if [ "$WSL_b" = "True" ]; then

        echo "Running STAR iteratively"

        for iter in $(seq 0 5); do

            STAR \
                --runMode inputAlignmentsFromBAM \
                --inputBAMfile ./results/tracks/"$sample"_"$muts"_"$iter".bam \
                --outWigType bedGraph \
                --outWigNorm None \
                --outWigStrand Stranded \
                --outFileNamePrefix ./results/tracks/"$sample"_"$muts"_"$iter"_

        done



        # Take only unique component of track
        append_track_signals

        rm -f ./results/tracks/"$sample"*.bg

    else
        # Make tracks
        echo '*making tracks with STAR'
        parallel -j "$cpus" STAR \
                                --runMode inputAlignmentsFromBAM \
                                --inputBAMfile ./results/tracks/"$sample"_{1}_{2}.bam \
                                --outWigType bedGraph \
                                --outWigNorm None \
                                --outWigStrand Stranded \
                                --outFileNamePrefix ./results/tracks/"$sample"_{1}_{2}_ ::: $muts \
                                                                         ::: $(seq 0 5)

        # Take only unique component of track
        echo '* Taking only unique components of tracks'
        append_track_signals

        rm -f ./results/tracks/"$sample"*.bg


    fi

    # parallel -j "$cpus" "bedtools genomecov \
    #                                 -ibam ${sample}_{1}_{2}.bam \
    #                                 -bg \
    #                                 -pc \
    #                                 -strand {3} \
    #                             | awk -v norm=${normVal} \
    #                                 '{print \$1, \$2, \$3, {3}1*norm*\$4}' \
    #                                 >> ${sample}.{1}.{2}.{4}.bedGraph" ::: $muts \
    #                                                                                  ::: $(seq 0 5) \
    #                                                                                  ::: + - \
    #                                                                                  :::+ pos min
    #



    # Make tdf files from the tracks
    echo '* Make tdf files from tracks'
    ensure_all_bedgraphs_have_data
    rm -f ./results/tracks/"$sample".*.tdf

    parallel -j "$cpus" igvtools toTDF \
                                -f mean,max \
                                ./results/tracks/"$sample".{1}.{2}.{3}.bedGraph \
                                ./results/tracks/"$sample".{1}.{2}.{3}.tdf \
                                "$chrom_sizes" ::: $muts \
                                             ::: $(seq 0 5) \
                                             ::: pos min

    rm -f igv.log


    ## Can comment out for debugging purposes
    rm -f ./results/tracks/"$sample"*.bam
    rm -f ./results/tracks/"$sample"*_reads.txt
    rm -f ./results/tracks/"$sample"*.bedGraph
    rm -f ./results/tracks/"$sample"*.bai
    rm -f ./results/tracks/"$sample"*.out

    # rm -f "$sample"*.chrom.sizes
    #rm -f igv*

    echo "* TDF track files created for sample $sample."
