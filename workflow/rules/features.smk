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
    threads: 20
    params:
        strand=FC_STRAND,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
        extra=config["fc_genes_extra"] + FC_GENES_PARAMS,
    log:
        "logs/featurecounts_gene/{sample}.log",
    wrapper:
        "v3.0.2/bio/subread/featurecounts"


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
    threads: 20
    params:
        strand=FC_STRAND,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
        extra=config["fc_exons_extra"] + FC_EXONS_PARAMS,
    log:
        "logs/featurecounts_exons/{sample}.log",
    wrapper:
        "v3.0.2/bio/subread/featurecounts"


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
    threads: 20
    params:
        strand=FC_STRAND,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
        extra=config["fc_transcripts_extra"] + FC_TRANSCRIPTS_PARAMS,
    log:
        "logs/featurecounts_transcripts/{sample}.log",
    wrapper:
        "v3.0.2/bio/subread/featurecounts"


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
    threads: 20
    params:
        strand=FC_STRAND,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
        extra=config["fc_exonbins_extra"] + FC_EXONBINS_PARAMS,
    log:
        "logs/featurecounts_exonbins/{sample}.log",
    wrapper:
        "v3.0.2/bio/subread/featurecounts"


# Get the set of isoforms a read maps to from the transcriptome bam
# TO-DO: No reason this can't be split up and multi-threaded
rule read_to_transcripts:
    input:
        bam="results/align/{sample}-Aligned.toTranscriptome.out.bam",
    output:
        table="results/read_to_transcripts/{sample}.csv"
    log:
        "logs/read_to_transcripts/{sample}.log"
    conda:
        "../envs/full.yaml"
    threads: 1
    script:
        "../scripts/bam2bakR/transcript_assignment.py"
