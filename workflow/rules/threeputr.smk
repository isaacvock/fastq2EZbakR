"""
Rules to build a custom 3'UTR annotation from 3'-end data
"""


### For paired-end experiments, keep only the read that "matters"
# Read 1, first nucleotide represents the transcript end
rule get_informative_read:
    input:
        "results/sf_reads/{sample}.s.bam",
    output:
        "results/informative_read/{sample}_informative.bam",
    log:
        "logs/get_informative_read/{sample}.log",
    conda:
        "../envs/genomictools.yaml"
    threads: 8
    params:
        # Choose which mate is informative based on strandedness.
        #   - reverse-stranded: read 1 (flag 64) gives 3'-end info
        #   - forward-stranded: read 2 (flag 128) gives 3'-end info
        informative_flag=lambda wildcards: (
            "64" if config["strandedness"] == "reverse" else "128"
        ),
    shell:
        """
        samtools view -@ {threads} -hb -f {params.informative_flag} {input} -o {output} 1> {log} 2>&1
        """


### Convert to PAS location tracks
# Strandedness handling a bity wonky
rule bam_to_3pend_bg:
    input:
        fetch_informative_read,
    output:
        "results/bam2bg/{sample}_informative_{strand}.bg",
    log:
        "logs/bam2_3pend_bg/{sample}_{strand}.log",
    conda:
        "../envs/genomictools.yaml"
    params:
        strandedness=config.get("strandedness", "reverse"),
        coverage_cutoff=config.get("coverage_cutoff", 10),
    threads: 1
    shell:
        r"""
        T={params.coverage_cutoff}
        strand="{wildcards.strand}"
        library="{params.strandedness}"

        if [[ "$strand" == "plus" ]]; then
            if [[ "$library"  == "reverse" ]]; then
                strand_symbol="-"
            else
                strand_symbol="+"
            fi
        elif [[ "$strand" == "minus" ]]; then
            if [[ "$library"  == "reverse" ]]; then
                strand_symbol="+"
            else
                strand_symbol="-"
            fi
        fi

        genomeCoverageBed -ibam {input} -bg -strand "$strand_symbol" -5 \
        | awk -v T="$T" -v OFS='\t' '$4>=T' \
        | LC_COLLATE=C sort -k1,1 -k2,2n > {output} 2>{log}
        """


### Merge all bedgraphs into one
rule merge_3pend_bg:
    input:
        expand("results/bam2bg/{sample}_informative_{{strand}}.bg", sample=SAMP_NAMES),
    output:
        "results/merge_3pend_bg/merged_3pend_{strand}.bg",
    params:
        strandedness=config.get("strandedness", "reverse"),
    conda:
        "../envs/genomictools.yaml"
    log:
        "logs/merge_3pend_bg/{strand}.log",
    shell:
        r"""
        strand="{wildcards.strand}"
        library="{params.strandedness}"

        if [[ "$strand" == "plus" ]]; then
            if [[ "$library"  == "reverse" ]]; then
                strand_symbol="-"
            else
                strand_symbol="+"
            fi
        elif [[ "$strand" == "minus" ]]; then
            if [[ "$library"  == "reverse" ]]; then
                strand_symbol="+"
            else
                strand_symbol="-"
            fi
        fi

        bedtools unionbedg -i {input} \
        | awk 'BEGIN{{OFS="\t"}}{{s=0; for(i=4;i<=NF;i++) s+=$i; print $1,$2,$3,s}}' \
        | LC_COLLATE=C sort -k1,1 -k2,2n > {output} 2>{log}
        """


### Alternatively, identify regions of transcript end clusters, representing likely
### cleavage and polyadenylation sites
rule cluster_PAS_bedtools:
    input:
        "results/merge_3pend_bg/merged_3pend_{strand}.bg",
    output:
        "results/call_PAS/bedtools_clusters_{strand}.bg",
    log:
        "logs/call_PAS_bedtools/bedtools_cluster_{strand}.log",
    params:
        extra=config.get("bedtools_cluster_distance", 50),
    conda:
        "../envs/genomictools.yaml"
    threads: 1
    shell:
        """
        bedtools cluster -d {params.extra} -i {input} > {output} 2>{log}
        """


### Get coverage, start & end for all peak clusters
rule summarise_PAS_clusters:
    input:
        "results/call_PAS/bedtools_clusters_{strand}.bg",
    output:
        "results/summarise_PAS_clusters/summarise_clusters_{strand}.bg",
    log:
        "logs/summarise_PAS_clusters/{strand}.log",
    conda:
        "../envs/genomictools.yaml"
    threads: 1
    shell:
        """
        bedtools groupby -g 1,5 -c 2,3,4 -o min,max,sum -i {input} > {output} 2> {log}
        """


### Use peak calls to make UTR GTF
# Strand logic a bit confusing as + is - and - is + due to reverse
# strandedness of libraries
rule make_threepUTR_gtf:
    input:
        bed=expand(
            "results/summarise_PAS_clusters/summarise_clusters_{strand}.bg",
            strand=["minus", "plus"],
        ),
        gtf=config.get("annotation"),
        fasta=config.get("genome"),
    output:
        "annotations/threepUTR_annotation.gtf",
    log:
        "logs/make_threepUTR_gtf/make_threepUTR_gtf.log",
    conda:
        "../envs/Rbio.yaml"
    params:
        rscript=workflow.source_path("../scripts/threeputrs/make_3pgtf.R"),
        coverage=config.get("cluster_coverage", 20 * NUM_SAMPS),
        fxn=config.get("cluster_fxn", 0.0),
        extension=config.get("extension", 0),
        polyA=config.get("false_polyA_len", 7),
        CPA=config.get("require_CPA_site", False),
        only_annotated=config.get("only_annotated_threeputrs", False),
    threads: 1
    shell:
        """
        chmod +x {params.rscript}
        {params.rscript} \
            --bed_minus {input.bed[0]} \
            --bed_plus {input.bed[1]} \
            --gtf {input.gtf} \
            --fasta {input.fasta} \
            --output {output} \
            --extension {params.extension} \
            --min_coverage {params.coverage} \
            --false_polyA_len {params.polyA} \
            --require_CPA {params.CPA} \
            --only_annotated {params.only_annotated} \
            --min_fxn {params.fxn} 1> {log} 2>&1
        """
