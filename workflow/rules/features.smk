### THESE RULES PERTAIN TO THE ASSIGNMENT OF READS TO FEATURES


# Assign reads to genes
rule featurecounts_genes:
    input:
        samples="results/sf_reads/{sample}.s.bam",
        annotation=config["annotation"],
    output:
        multiext(
            "results/featurecounts_genes/{sample}",
            ".featureCounts",
            ".featureCounts.summary",
        ),
        temp("results/featurecounts_genes/{sample}.s.bam.featureCounts"),
    threads: 20
    params:
        strand=FC_STRAND,
        extra=config["fc_genes_extra"] + FC_GENES_PARAMS,
    log:
        "logs/featurecounts_genes/{sample}.log",
    wrapper:
        "v7.2.0/bio/subread/featurecounts"


# Assign reads to exons
rule featurecounts_exons:
    input:
        samples="results/sf_reads/{sample}.s.bam",
        annotation=config["annotation"],
    output:
        multiext(
            "results/featurecounts_exons/{sample}",
            ".featureCounts",
            ".featureCounts.summary",
            ".featureCounts.jcounts",
        ),
        temp("results/featurecounts_exons/{sample}.s.bam.featureCounts"),
    threads: 20
    params:
        strand=FC_STRAND,
        extra=config["fc_exons_extra"] + FC_EXONS_PARAMS,
    log:
        "logs/featurecounts_exons/{sample}.log",
    wrapper:
        "v7.2.0/bio/subread/featurecounts"


# Assign reads to transcripts
rule featurecounts_transcripts:
    input:
        samples="results/sf_reads/{sample}.s.bam",
        annotation=config["annotation"],
    output:
        multiext(
            "results/featurecounts_transcripts/{sample}",
            ".featureCounts",
            ".featureCounts.summary",
        ),
        temp("results/featurecounts_transcripts/{sample}.s.bam.featureCounts"),
    threads: 20
    params:
        strand=FC_STRAND,
        extra=config["fc_transcripts_extra"] + FC_TRANSCRIPTS_PARAMS,
    log:
        "logs/featurecounts_transcripts/{sample}.log",
    wrapper:
        "v7.2.0/bio/subread/featurecounts"


# Assign reads to exonic bins
rule featurecounts_exonbins:
    input:
        samples="results/sf_reads/{sample}.s.bam",
        annotation=config["flat_annotation"],
    output:
        multiext(
            "results/featurecounts_exonbins/{sample}",
            ".featureCounts",
            ".featureCounts.summary",
        ),
        temp("results/featurecounts_exonbins/{sample}.s.bam.featureCounts"),
    threads: 20
    params:
        strand=FC_STRAND,
        extra=config["fc_exonbins_extra"] + FC_EXONBINS_PARAMS,
    log:
        "logs/featurecounts_exonbins/{sample}.log",
    wrapper:
        "v7.2.0/bio/subread/featurecounts"


if TEC_MODE == "transcriptome":
    if config.get("run_rsem", True):

        # Get strand-compatible isoforms a read maps to from the RSEM transcriptome BAM
        rule read_to_transcripts:
            input:
                bam="results/rsem/{sample}.transcript.bam",
            output:
                table=temp("results/read_to_transcripts/{sample}.csv"),
            params:
                strandedness=config["strandedness"],
            log:
                "logs/read_to_transcripts/{sample}.log",
            conda:
                "../envs/full.yaml"
            threads: 4
            script:
                "../scripts/features/transcript_assignment.py"

    else:

        # Get strand-compatible isoforms a read maps to from the STAR transcriptome BAM
        rule read_to_transcripts:
            input:
                bam="results/align/{sample}-Aligned.toTranscriptome.out.bam",
            output:
                table=temp("results/read_to_transcripts/{sample}.csv"),
            params:
                strandedness=config["strandedness"],
            log:
                "logs/read_to_transcripts/{sample}.log",
            conda:
                "../envs/full.yaml"
            threads: 4
            script:
                "../scripts/features/transcript_assignment.py"

else:

    # Aligner-agnostic genome-BAM TEC assignment (HISAT-3N, HISAT2, bam2bakr, or tec_mode='genome')
    rule read_to_transcripts:
        input:
            bam="results/sf_reads/{sample}.s.bam",
            gtf=config["annotation"],
        output:
            table=temp("results/read_to_transcripts/{sample}.csv"),
        params:
            pythonscript=workflow.source_path(
                "../scripts/features/genome_tec_assignment.py"
            ),
            strandedness=config["strandedness"],
            end_overhang=config.get("tec_end_overhang", 2),
        log:
            "logs/read_to_transcripts/{sample}.log",
        conda:
            "../envs/full.yaml"
        threads: 8
        shell:
            """
            python {params.pythonscript:q} \
                -b {input.bam:q} \
                -g {input.gtf:q} \
                -o {output.table:q} \
                -s {params.strandedness:q} \
                -t {threads} \
                --end-overhang {params.end_overhang} 1> {log:q} 2>&1
            """


# Get deduplicated set of junctions each read/pair overlaps directly via pysam
rule read_to_junctions:
    input:
        bam="results/sf_reads/{sample}.s.bam",
    output:
        table=temp("results/read_to_junctions/{sample}.csv.gz"),
    params:
        pythonscript=workflow.source_path("../scripts/features/junction_assignment.py"),
    log:
        "logs/read_to_junctions/{sample}.log",
    threads: 8
    conda:
        "../envs/full.yaml"
    shell:
        """
        python {params.pythonscript:q} -b {input.bam:q} -o {output.table:q} -t {threads} 1> {log:q} 2>&1
        """


# Make junction annotation that featureCounts can assign reads with respect to
rule junction_annotation:
    input:
        config["annotation"],
    output:
        "junction_annotation/junctions.gtf",
    params:
        rscript=workflow.source_path("../scripts/features/junction_annotation.R"),
        extra=config["junction_annotation_params"],
    log:
        "logs/junction_annotation/junctions.log",
    threads: 1
    conda:
        "../envs/junctions.yaml"
    shell:
        r"""
        Rscript {params.rscript:q} -r {input:q} -o {output:q} {params.extra} 1> {log:q} 2>&1
        """


# Experimental: assign reads to exon-exon junctions through ad-hoc approach
rule featurecounts_eej:
    input:
        samples="results/sf_reads/{sample}.s.bam",
        annotation="junction_annotation/junctions.gtf",
    output:
        multiext(
            "results/featurecounts_eej/{sample}",
            ".featureCounts",
            ".featureCounts.summary",
        ),
        temp("results/featurecounts_eej/{sample}.s.bam.featureCounts"),
    threads: 20
    params:
        strand=FC_STRAND,
        extra=config["fc_eej_extra"] + FC_EEJ_PARAMS,
    log:
        "logs/featurecounts_eej/{sample}.log",
    wrapper:
        "v7.2.0/bio/subread/featurecounts"


# Experimental: assign reads to exon-intron junctions through ad-hoc approach
rule featurecounts_eij:
    input:
        samples="results/sf_reads/{sample}.s.bam",
        annotation="junction_annotation/junctions.gtf",
    output:
        multiext(
            "results/featurecounts_eij/{sample}",
            ".featureCounts",
            ".featureCounts.summary",
        ),
        temp("results/featurecounts_eij/{sample}.s.bam.featureCounts"),
    threads: 20
    params:
        strand=FC_STRAND,
        extra=config["fc_eij_extra"] + FC_EIJ_PARAMS,
    log:
        "logs/featurecounts_eij/{sample}.log",
    wrapper:
        "v7.2.0/bio/subread/featurecounts"


# Assign reads to 3'-UTRs using only the single 3'-end nucleotide of the informative mate
rule featurecounts_3utr:
    input:
        samples=fetch_informative_read,
        annotation=THREEPUTR_ANNOTATION,
    output:
        multiext(
            "results/featurecounts_3utr/{sample}",
            ".featureCounts",
            ".featureCounts.summary",
        ),
        temp(f"results/featurecounts_3utr/{{sample}}{FC_3UTR_CORE_EXT}"),
    conda:
        "../envs/genomictools.yaml"
    params:
        strand=FC_STRAND,
        extra=FC_3UTR_PARAMS,
    threads: 20
    log:
        "logs/featurecounts_3utr/{sample}.log",
    wrapper:
        "v7.2.0/bio/subread/featurecounts"
