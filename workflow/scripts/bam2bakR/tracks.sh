#!/bin/bash
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
    star_strands=(str2 str1)
else
    star_strands=(str1 str2)
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
        if [ $b == "GA" ]; then
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
        parallel -j "$cpus" "awk -v norm=${normVal} \
                                        '{print \$1, \$2, \$3, {3}1*norm*\$4}' \
                                        ./results/tracks/${sample}_{1}_{2}_Signal.Unique.{4}.out.bg \
                                        >> ./results/tracks/${sample}.{1}.{2}.{5}.bedGraph" ::: $muts \
                                                                           ::: $(seq 0 5) \
                                                                           ::: + - \
                                                                           :::+ "${star_strands[@]}" \
                                                                           :::+ pos min

        rm ./results/tracks/"$sample"*.bg

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
        parallel -j "$cpus" "awk -v norm=${normVal} \
                                        '{print \$1, \$2, \$3, {3}1*norm*\$4}' \
                                        ./results/tracks/${sample}_{1}_{2}_Signal.Unique.{4}.out.bg \
                                        >> ./results/tracks/${sample}.{1}.{2}.{5}.bedGraph" ::: $muts \
                                                                           ::: $(seq 0 5) \
                                                                           ::: + - \
                                                                           :::+ "${star_strands[@]}" \
                                                                           :::+ pos min

        rm ./results/tracks/"$sample"*.bg


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
    parallel -j "$cpus" igvtools toTDF \
                                -f mean,max \
                                ./results/tracks/"$sample".{1}.{2}.{3}.bedGraph \
                                ./results/tracks/"$sample".{1}.{2}.{3}.tdf \
                                "$chrom_sizes" ::: $muts \
                                             ::: $(seq 0 5) \
                                             ::: pos min

    rm -f igv.log


    ## Can comment out for debugging purposes
    rm ./results/tracks/"$sample"*.bam
    rm ./results/tracks/"$sample"*_reads.txt
    rm ./results/tracks/"$sample"*.bedGraph
    rm ./results/tracks/"$sample"*.bai
    rm ./results/tracks/"$sample"*.out

    # rm -f "$sample"*.chrom.sizes
    #rm -f igv*

    echo "* TDF track files created for sample $sample."
