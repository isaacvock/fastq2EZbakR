"""
Rules to build a 3'UTR annotation from a GTF or 3'-end data
"""


### Annotate terminal exons when de novo 3'-end calling is disabled
rule annotate_threepUTR_gtf:
    input:
        gtf=config["annotation"],
    output:
        "annotations/annotated_threepUTR_annotation.gtf",
    log:
        "logs/annotate_threepUTR_gtf/annotate_threepUTR_gtf.log",
    conda:
        "../envs/Rbio.yaml"
    threads: 1
    params:
        rscript=workflow.source_path("../scripts/threeputrs/annotate_3pgtf.R"),
    shell:
        """
        Rscript {params.rscript:q} --gtf {input.gtf:q} --output {output:q} 1> {log:q} 2>&1
        """


### For paired-end experiments, keep only the informative 3'-end mate
rule get_informative_read:
    input:
        lambda wc: (
            f"results/coord_bam/{wc.sample}.bam"
            if config.get("call_threeputrs", False)
            else f"results/sf_reads/{wc.sample}.s.bam"
        ),
    output:
        temp("results/informative_read/{sample}_informative.bam"),
    log:
        "logs/get_informative_read/{sample}.log",
    conda:
        "../envs/genomictools.yaml"
    threads: 8
    params:
        # Choose which mate is informative based on strandedness:
        #   - reverse-stranded: read 1 (flag 64) gives 3'-end info
        #   - forward-stranded: read 2 (flag 128) gives 3'-end info
        informative_flag=lambda wildcards: (
            "64" if config["strandedness"] == "reverse" else "128"
        ),
    shell:
        """
        samtools view -@ {threads} -hb -f {params.informative_flag} {input:q} -o {output:q} 1> {log:q} 2>&1
        """


### Convert coordinate-sorted BAM directly to PAS location tracks (zero redundant sort)
rule bam_to_3pend_bg:
    input:
        lambda wc: (
            f"results/informative_read/{wc.sample}_informative.bam"
            if config.get("PE", False)
            else f"results/coord_bam/{wc.sample}.bam"
        ),
    output:
        temp("results/bam2bg/{sample}_informative_{strand}.bg"),
    log:
        "logs/bam_to_3pend_bg/{sample}_{strand}.log",
    conda:
        "../envs/genomictools.yaml"
    threads: 1
    params:
        strandedness=config.get("strandedness", "reverse"),
        coverage_cutoff=config.get("coverage_cutoff", 10),
    shell:
        r"""
        T={params.coverage_cutoff}
        strand="{wildcards.strand}"
        library="{params.strandedness}"

        if [[ "$strand" == "plus" ]]; then
            if [[ "$library" == "reverse" ]]; then
                strand_symbol="-"
                pos_symbol="-5"
            else
                strand_symbol="+"
                pos_symbol="-3"
            fi
        elif [[ "$strand" == "minus" ]]; then
            if [[ "$library" == "reverse" ]]; then
                strand_symbol="+"
                pos_symbol="-5"
            else
                strand_symbol="-"
                pos_symbol="-3"
            fi
        fi

        genomeCoverageBed -ibam {input:q} -bg -strand "$strand_symbol" "$pos_symbol" 2> {log:q} \
        | awk -v T="$T" -v OFS='\t' '$4>=T' \
        | LC_COLLATE=C sort -k1,1 -k2,2n > {output:q}
        """


### Merge all bedgraphs into one (supports single-sample and multi-sample experiments)
rule merge_3pend_bg:
    input:
        expand("results/bam2bg/{sample}_informative_{{strand}}.bg", sample=SAMP_NAMES),
    output:
        temp("results/merge_3pend_bg/merged_3pend_{strand}.bg"),
    log:
        "logs/merge_3pend_bg/{strand}.log",
    conda:
        "../envs/genomictools.yaml"
    shell:
        r"""
        inputs=({input:q})
        if [[ "${{#inputs[@]}}" -eq 1 ]]; then
            LC_COLLATE=C sort -k1,1 -k2,2n "${{inputs[0]}}" > {output:q} 2> {log:q}
        else
            bedtools unionbedg -i {input:q} 2> {log:q} \
            | awk 'BEGIN{{OFS="\t"}}{{s=0; for(i=4;i<=NF;i++) s+=$i; print $1,$2,$3,s}}' \
            | LC_COLLATE=C sort -k1,1 -k2,2n > {output:q} 2>> {log:q}
        fi
        """


### Identify regions of transcript end clusters representing cleavage/polyadenylation sites
rule cluster_PAS_bedtools:
    input:
        "results/merge_3pend_bg/merged_3pend_{strand}.bg",
    output:
        temp("results/call_PAS/bedtools_clusters_{strand}.bg"),
    log:
        "logs/call_PAS_bedtools/bedtools_cluster_{strand}.log",
    conda:
        "../envs/genomictools.yaml"
    threads: 1
    params:
        extra=config.get("bedtools_cluster_distance", 100),
    shell:
        """
        bedtools cluster -d {params.extra} -i {input:q} > {output:q} 2> {log:q}
        """


### Get total endpoint support (depth x interval width), start & end for clusters
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
        awk 'BEGIN{{OFS="\t"}} {{$4 *= ($3 - $2); print}}' {input:q} \
        | bedtools groupby -g 1,5 -c 2,3,4 -o min,max,sum -i - > {output:q} 2> {log:q}
        """


### Use peak calls to make UTR GTF
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
    threads: 1
    params:
        rscript=workflow.source_path("../scripts/threeputrs/make_3pgtf.R"),
        coverage=config.get("cluster_coverage", 20 * NUM_SAMPS),
        fxn=config.get("cluster_fxn", 0.0),
        extension=config.get("extension", 0),
        polyA=config.get("false_polyA_len", 7),
        CPA=config.get("require_CPA_site", False),
        only_annotated=config.get("only_annotated_threeputrs", False),
    shell:
        """
        Rscript {params.rscript:q} \
            --bed_minus {input.bed[0]:q} \
            --bed_plus {input.bed[1]:q} \
            --gtf {input.gtf:q} \
            --fasta {input.fasta:q} \
            --output {output:q} \
            --extension {params.extension} \
            --min_coverage {params.coverage} \
            --false_polyA_len {params.polyA} \
            --require_CPA {params.CPA} \
            --only_annotated {params.only_annotated} \
            --min_fxn {params.fxn} 1> {log:q} 2>&1
        """
