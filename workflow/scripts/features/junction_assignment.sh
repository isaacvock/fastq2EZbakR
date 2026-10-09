#!/bin/bash

set -euo pipefail

# Source the paths and variables:
cpus=$1
sample=$2
input=$3
output=$4
output2=$5
pyscript=$6
awkscript=$7

# Create results/rsem_csv/
touch "$output2"

# Infer a decent fragment size
num_reads=$(samtools view -@ "$cpus" -c "$input")
if [[ "$num_reads" -eq 0 ]]; then
    printf 'qname,junction_start,junction_end\n' | pigz -p "$cpus" > "$output"
    exit 0
fi

fragment_size=$(echo "scale=0; $num_reads/$cpus" | bc)
fragment_size=$((fragment_size + 1))


# Spread the work load so all cpus are working at all times
    declare $(samtools view -@ "$cpus"  "$input" \
    			| wc -l \
    			| awk -v fragment_size="$fragment_size" \
                      -v cpus="$cpus" '{
                                        nFrag = int($1 / fragment_size)                         # Calculate to how many fragments must the .bam be split
                                        if ($1 % fragment_size != 0) {nFrag = nFrag + 1}        # If there are some left over reads then increase number of fractions

                                        nBatchRun = int(nFrag / cpus)                           # In how many batches will the fragments be processed
                                        if (nFrag % cpus != 0) { nBatchRun = nBatchRun + 1}     # If there are some leftover fragments that will not fill up all cpus, then incerease number of batches

                                        newFragmentNumber = nBatchRun * cpus

                                        newFragmentSize = int($1 / newFragmentNumber)
                                        if ($1 % newFragmentNumber != 0) {newFragmentSize = newFragmentSize + 2}    # if there are still some reads leftover, then increase fragment size by 1 read-pair

                                        print "newFragmentNumber="newFragmentNumber
                                        print "newFragmentSize="newFragmentSize
                                    }')

    echo "* The number of fragments for sample $sample is $newFragmentNumber"
    echo "* The fragment size is set to $newFragmentSize alignments per fragment"



    for f in $(seq 1 $newFragmentNumber); do
        samtools view -H "$input" > ./results/read_to_junctions/"$f"_"$sample".sam
    done &&


    samtools view "$input" \
    	| awk \
    		-v fragment_size="$newFragmentSize" \
    		-v sample="$sample" \
            -v directory="./results/read_to_junctions" \
    		-f "$awkscript" &&

    for f in $(seq 1 $newFragmentNumber); do
        samtools view -@ "$cpus" -o ./results/read_to_junctions/"$f"_"$sample"_frag.bam ./results/read_to_junctions/"$f"_"$sample".sam
        rm ./results/read_to_junctions/"$f"_"$sample".sam
    done &&

    echo "* Aligned .sam file fragmented for sample $sample"


# remove any remaining sam files:

    #rm -f *.sam

#!/bin/bash
# Main sript for mutation calling in samples

## Need to add mutation call script to scripts!

# Call mutations
    parallel -j $cpus "python $pyscript -b {1}" ::: ./results/read_to_junctions/*_"$sample"_frag.bam

    echo "** Junction information extracted for $sample"


    rm -f ./results/read_to_junctions/*_"$sample"_frag.bam

    echo "** fragmented bam files cleaned up for $sample"



# Combine output from fragments into single file
    # 1) _rsem.csv files
    awk 'FNR > 1 || NR == 1' ./results/read_to_junctions/*_"$sample"_frag_junctions.csv \
        | pigz -p $cpus > "$output"

    rm ./results/read_to_junctions/*_"$sample"_frag_junctions.csv



    echo "** Results fragments merged into final files"



	echo '* Cleaning up fragmented .bam files'
