#!/bin/bash
set -euo pipefail

# Cute trick to deal with fact that there is an uncertain number of control sample
    # control_samples becomes array with all args
    # remove args that I know aren't the actual control_samples

cpus=$1
nsamps=$2
output_txt=$3
output_vcf=$4
mpileup_options="$5"
call_options="$6"
genome_fasta=$7

shift 8

control_samples=("$@")




if [[ "$nsamps" -gt 0 ]] 
then
    # Loop through control samples:
        declare -a bam_list=()

        for cs in ${control_samples[@]}; do
            name=$(echo "$cs" | cut -d '/' -f 3 | rev | cut -c8- | rev)

            sorted="./results/snps/"$name"_sort.bam"

            samtools sort -@ "$cpus" -o ./results/snps/"$name"_sort.bam "$cs"
            samtools index -@ "$cpus" ./results/snps/"$name"_sort.bam
            #bcftools mpileup --threads "$cpus" -f "$genome_fasta" ./results/snps/"$name"_sort.bam | bcftools call --threads "$cpus" -mv > ./results/snps/snp-"$name".vcf

            bam_list+=( "$sorted" )

        done

        #cat ./results/snps/*.vcf > $output_vcf
        #rm ./results/snps/snp-*

        echo "${bam_list[@]}" | tr ' ' '\n' > ./results/snps/bam.list

        # for cs in ${control_samples[@]}; do
        #     NAMES+=($(echo "$cs" | cut -d '/' -f 3 | rev | cut -c8- | rev))
        # done

    # Parallelize SNPs calling. Each chromosome in each .bam file is processed as separate job
        # Note: This approach does not give the same snp.txt result. In 2199712 SNPs there were 16 different.
        # {1} : path to fasta file
        # {2} : samtools view -H ${control_samples[0]}_sort.bam | awk ' $1 == "@SQ" {split($2,a,":"); print a[2]}' : Extracts chromosome names from .bam file header
        # {3} : ${control_samples[@]/%/_sort.bam}                                                    : Appends "_sort.bam" to the end of control names and prints
        # parallel -j "$cpus" "bcftools mpileup -f {1} \
        #                                            -r {2} {3} \
        #                         | bcftools call -mv" ::: $genome_fasta \
        #                                               ::: $(samtools view -H ./results/snps/${NAMES[0]}_sort.bam \
        #                                                             | awk ' $1 == "@SQ" {split($2,a,":"); print a[2]}') \
        #                                               ::: ./results/snps/*_sort.bam > $output_vcf

        bcftools mpileup --threads "$cpus" \
                         -f "$genome_fasta" \
                         -b ./results/snps/bam.list \
                         $mpileup_options \
                         -Ou \
        | bcftools call --threads "$cpus" $call_options -mv -Oz -o $output_vcf


        # Note: Easier and also fast option would be:  bcftools mpileup --threads $cpus -f $genome_fasta "$cs"_sort.bam | bcftools call --threads $cpus-mv > snp-"$cs".vcf


    # Clean this up for all possible mutations:
        # Note: Filtering now inludes cases where both alleles are mutated
        bcftools view -H "$output_vcf" | \
        awk '$1 !~ /^#/ && length($4) == 1 {if (length($5) == 1) {print $4":"$5":"$1":"$2}
                                            else if (length($5) == 3 && $5 ~ /,/) {split($5, mut, ",")
                                                                                   print $4":"mut[1]":"$1":"$2
                                                                                   print $4":"mut[2]":"$1":"$2}
                                            }' \
            | sort \
            | uniq > $output_txt

        echo '* SNPs called and snp.txt generated'
else
    touch "$output_txt"
    touch "$output_vcf"
fi


rm -f 0