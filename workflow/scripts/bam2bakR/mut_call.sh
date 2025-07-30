#!/bin/bash

# Source the paths and variables:
cpus=$1
sample=$2
input=$3
input_snp=$4
output=$5
output2=$6
minqual=$7
mut_tracks=$8
format=$9
strand=${10}
mutcnt=${11}
awkscript=${12}
mutpos=${13}

# Exit immediately if any command returns a non-zero status
set -e


# Create results/counts/
touch "$output2"

# Infer a decent fragment size
num_reads=$(samtools view -@ "$cpus" -c "$input")

# Get the available RAM in MB and convert it to GB
available_ram=$(free -m | awk '/^Mem:/{print $7}') # Memory in MB
available_ram=$(echo "$available_ram / 1024" | bc) # Convert to GB

# Calculate fragment_size
    # 20 million reads use roughly 1 TB of RAM in mutation counting
    # script when mutpos is True, so plan accordingly
fragment_size_ram_based=$(echo "scale=0; $available_ram / ($cpus * 0.00025)" | bc)


# Calculate fragment_size based on number_of_reads and CPUs
fragment_size_read_based=$(echo "scale=0; $num_reads / $cpus" | bc)

# Choose the minimum of the two calculated fragment sizes
if [ $fragment_size_ram_based -lt $fragment_size_read_based ]; then

    fragment_size=$fragment_size_ram_based

else

    fragment_size=$fragment_size_read_based
    
fi

# Make sure fragment size isn't 0, though if it were, you probably have bigger problems on your hand...
(( fragment_size++ ))


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
        samtools view -H "$input" > ./results/counts/"$f"_"$sample".sam
    done &&


    # If running mutation position counting, get bam file roughly coordinate sorted to 
    # try and get all reads covering a particular base close together, so as to reduce
    # disk space and RAM usage when running the mutation counting script.
samtools sort -@ "$cpus" "$input" | samtools sort -n -@ "$cpus" - | samtools view - \
    | awk \
        -v fragment_size="$newFragmentSize" \
        -v sample="$sample" \
        -f "$awkscript" &&

for f in $(seq 1 $newFragmentNumber); do
    samtools view -@ "$cpus" -o ./results/counts/"$f"_"$sample"_frag.bam ./results/counts/"$f"_"$sample".sam
    rm ./results/counts/"$f"_"$sample".sam
done &&

echo "* Aligned .sam file fragmented for sample $sample"

# remove any remaining sam files:

    #rm -f *.sam

#!/bin/bash
# Main sript for mutation calling in samples

## Need to add mutation call script to scripts!

shopt -s extglob # to deal with more specific regex below
# Call mutations
    parallel -j $cpus "python $mutcnt -b {1} \
                                              --mutType $mut_tracks \
                                              --minQual $minqual \
											  --SNPs "./results/snps/snp.txt" \
                                              --strandedness $strand \
                                              --mutPos \
                                              --reads $format" ::: ./results/counts/+([0-9])_"$sample"_frag.bam \


    echo "** Mutations called for sample $sample"


	rm -f ./results/counts/+([0-9])_"$sample"_frag.bam

	echo '* Cleaning up fragmented .bam files'


    # Pre-sort cU fragments
    parallel -j $cpus "tail -n +2 {1} \
                            | LC_COLLATE=C sort > {1.}_sort.csv" ::: ./results/counts/+([0-9])_"$sample"_frag_cU.csv
    rm ./results/counts/+([0-9])_"$sample"_frag_cU.csv                    

    # Combine pre-sorted fragments
    LC_COLLATE=C sort -m ./results/counts/+([0-9])_"$sample"_frag_cU_sort.csv > ./results/counts/"$sample"_cU_comb.csv 
    rm ./results/counts/+([0-9])_"$sample"_frag_cU_sort.csv

    # Get ammount of data
    cUsize=$(wc -l ./results/counts/"$sample"_cU_comb.csv | cut -d ' ' -f 1)

    # Split file to sorted fragments keeping the same position groups in the same file fragment 
    awk -v FS="," \
        -v cpus=$cpus \
        -v cUsize=$cUsize \
        -v sample=$sample \
        'BEGIN { 
                fragment_size = int(cUsize / cpus) + 1
                i = 1
        }
        NR == (fragment_size * i - 1 ) { x = $1","$2 } 

        NR < (fragment_size * i ) { print >> "./results/counts/"i"_"sample"_cU_comb.csv" }

        NR >= (fragment_size * i ) { if ($1","$2 == x) 
                                        { 
                                            print >> "./results/counts/"i"_"sample"_cU_comb.csv"
                                        } 
                                    else 
                                        {           
                                            i++
                                            print >> "./results/counts/"i"_"sample"_cU_comb.csv"
                                        }
            
        }' ./results/counts/"$sample"_cU_comb.csv

    rm ./results/counts/"$sample"_cU_comb.csv

    # Awk script for processing sorted cU files; assumes all same positions are grouped together in subsequent lines
        # Thanks to this we have to keep in memory only one genomic position at given time
    function awkProcessCU () {
        awk -v FS="," 'NR == 1 {
                    nuc = $1":"$2
            }
            $1":"$2 != nuc {
                    for (pos in trial) {
                        if (n[pos] >= 0)
                            print pos","trial[pos]","n[pos]
                    }
                    delete trial
                    delete n
                    nuc = $1":"$2
            }
            $1":"$2 == nuc {
                    trial[$1","$2","$3] += $4
                    n[$1","$2","$3] += $5
            }
            END {
                    for (pos in trial) {
                        if (n[pos] >= 0)
                            print pos","trial[pos]","n[pos]
                    }
            }' $1
    }

    export -f awkProcessCU

    # Add header and process sorted fragments in parallel
    cat <(echo "rname,gloc,tp,trials,n") \
        <(parallel -j $cpus awkProcessCU {1} ::: ./results/counts/+([0-9])_${sample}_cU_comb.csv) \
        | pigz -p $cpus > ./results/counts/"$sample"_cU.csv.gz

    rm ./results/counts/+([0-9])_${sample}_cU_comb.csv

    echo "Finished mutation position counts processing"