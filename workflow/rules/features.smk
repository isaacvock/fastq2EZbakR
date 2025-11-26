### THESE RULES PERTAIN TO THE ASSIGNMENT OF READS TO FEATURES WITH FEATURECOUNTS


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
        strand=FC_STRAND,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
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
        strand=FC_STRAND,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
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
        strand=FC_STRAND,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
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
        strand=FC_STRAND,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
        extra=config["fc_exonbins_extra"] + FC_EXONBINS_PARAMS,
    log:
        "logs/featurecounts_exonbins/{sample}.log",
    wrapper:
        "v7.2.0/bio/subread/featurecounts"


if config.get("run_rsem", True):

    # Get the set of isoforms a read maps to from the transcriptome bam
    # TO-DO: No reason this can't be split up and multi-threaded
    rule read_to_transcripts:
        input:
            bam="results/rsem/{sample}.transcript.bam",
        output:
            table=temp("results/read_to_transcripts/{sample}.csv"),
        log:
            "logs/read_to_transcripts/{sample}.log",
        conda:
            "../envs/full.yaml"
        threads: 1
        script:
            "../scripts/features/transcript_assignment.py"

else:

    # Get the set of isoforms a read maps to from the transcriptome bam
    # TO-DO: No reason this can't be split up and multi-threaded
    rule read_to_transcripts:
        input:
            bam="results/align/{sample}-Aligned.toTranscriptome.out.bam",
        output:
            table=temp("results/read_to_transcripts/{sample}.csv"),
        log:
            "logs/read_to_transcripts/{sample}.log",
        conda:
            "../envs/full.yaml"
        threads: 1
        script:
            "../scripts/features/transcript_assignment.py"


# Get set of junctions a read overlaps
rule read_to_junctions:
    input:
        "results/sf_reads/{sample}.s.bam",
    output:
        temp("results/read_to_junctions/{sample}.csv.gz"),
        temp("results/read_to_junctions/{sample}_check.txt"),
    params:
        shellscript=workflow.source_path("../scripts/features/junction_assignment.sh"),
        pythonscript=workflow.source_path("../scripts/features/junction_assignment.py"),
        awkscript=workflow.source_path("../scripts/rsem_plus/fragment_sam_rsem.awk"),
    log:
        "logs/read_to_junctions/{sample}.log",
    threads: 32
    conda:
        "../envs/full.yaml"
    shell:
        """
        chmod +x {params.shellscript}
        chmod +x {params.pythonscript}
        chmod +x {params.awkscript}
        {params.shellscript} {threads} {wildcards.sample} {input} {output} {params.pythonscript} {params.awkscript}  1> {log} 2>&1
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
        chmod +x {params.rscript}
        {params.rscript} -r {input} -o {output} {params.extra} 1> {log} 2>&1
        """


# Experimental: assign reads to exon-exon junctions through ad-hoc not 100%-accurate approach
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
        strand=FC_STRAND,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
        extra=config["fc_eej_extra"] + FC_EEJ_PARAMS,
    log:
        "logs/featurecounts_eej/{sample}.log",
    wrapper:
        "v7.2.0/bio/subread/featurecounts"

# Experimental: assign reads to exon-intron junctions through ad-hoc not 100%-accurate approach
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
        strand=FC_STRAND,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
        extra=config["fc_eij_extra"] + FC_EIJ_PARAMS,
    log:
        "logs/featurecounts_eij/{sample}.log",
    wrapper:
        "v7.2.0/bio/subread/featurecounts"


# Assign reads to 3'-UTRs using only single nucleotide of read (3'-end position) for assignment
    # TO-DO: Might need to add a step for paired-end data where only informative read is extracted
    # and passed to this step. Not sure for paired-end data how featureCounts' read2pos works (will it
    # use both read's single base?)
rule featurecounts_3utr:
    input:
        samples="results/sf_reads/{sample}.s.bam",
        annotation=THREEPUTR_ANNOTATION,
    output:
        multiext(
            "results/featurecounts_3utr/{sample}",
            ".featureCounts",
            ".featureCounts.summary",
            ".featureCounts.jcounts",
        ),
        temp("results/featurecounts_3utr/{sample}.s.bam.featureCounts"),
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


