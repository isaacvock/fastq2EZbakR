### THESE RULES PERTAIN TO THE ASSIGNMENT OF READS TO FEATURES WITH FEATURECOUNTS


# Assign reads to genes
rule featurecounts_genes:
    input:
        sam = "results/sf_reads/{sample}.s.bam",
        gtf = config["annotation"]
    output:
        multiext(
            "results/featurecounts_genes/{sample}",
            ".featureCounts",
            ".featureCounts.summary"
        ),
    threads: 20
    params:
        strand=FC_STRAND,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
        r_path= lambda w, output: os.path.splitext(output[0])[0],  # implicitly sets the --Rpath flag
        extra= config["fc_genes_extra"] + FC_GENES_PARAMS,
    log:
        "logs/{sample}.log",
    wrapper:
        "v3.0.2/bio/subread/featurecounts"


# Assign reads to exons
rule featurecounts_exons:
    input:
        sam = "results/sf_reads/{sample}.s.bam",
        gtf = config["annotation"]
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
        r_path=lambda w, output: os.path.splitext(output[0])[0],  # implicitly sets the --Rpath flag
        extra= config["fc_exons_extra"] + FC_EXONS_PARAMS,
    log:
        "logs/{sample}.log",
    wrapper:
        "v3.0.2/bio/subread/featurecounts"

# Assign reads to transcripts
rule featurecounts_transcripts:
    input:
        sam = "results/sf_reads/{sample}.s.bam",
        gtf = config["annotation"]
    output:
        multiext(
            "results/featurecounts_transcripts/{sample}",
            ".featureCounts",
            ".featureCounts.summary"
        ),
    threads: 20
    params:
        strand=FC_STRAND,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
        r_path=lambda w, output: os.path.splitext(output[0])[0],  # implicitly sets the --Rpath flag
        extra= config["fc_transcripts_extra"] + FC_TRANSCRIPTS_PARAMS,
    log:
        "logs/{sample}.log",
    wrapper:
        "v3.0.2/bio/subread/featurecounts"

# Assign reads to exonic bins
rule featurecounts_exonbins:
    input:
        sam = "results/sf_reads/{sample}.s.bam",
        gtf = config["flat_annotation"]
    output:
        multiext(
            "results/featurecounts_exonbins/{sample}",
            ".featureCounts",
            ".featureCounts.summary"
        ),
    threads: 20
    params:
        strand=FC_STRAND,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
        r_path=lambda w, output: os.path.splitext(output[0])[0],  # implicitly sets the --Rpath flag
        extra= config["fc_exonbins_extra"] + FC_EXONBINS_PARAMS,
    log:
        "logs/{sample}.log",
    wrapper:
        "v3.0.2/bio/subread/featurecounts"




