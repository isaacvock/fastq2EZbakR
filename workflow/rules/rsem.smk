# Create index for RSEM
rule RSEM_index:
    input:
        reference_genome=config["genome"],
        gtf=AandQ_ANNOTATION,
    output:
        seq="rsem_index/reference.seq",
        grp="rsem_index/reference.grp",
        ti="rsem_index/reference.ti",
        tfa="rsem_index/reference.transcripts.fa",
        idxfa="rsem_index/reference.idx.fa",
        n2g="rsem_index/reference.n2g.idx.fa",
    params:
        extra="--gtf {} {}".format(
            str(AandQ_ANNOTATION), str(config["rsem_index_params"])
        ),
    log:
        "logs/RSEM_index/prepare-reference.log",
    threads: 20
    wrapper:
        "v2.3.1/bio/rsem/prepare-reference"


# Run RSEM to quantify transcript abundances
rule RSEM:
    input:
        bam="results/align/{sample}-Aligned.toTranscriptome.out.bam",
        reference=multiext(
            "rsem_index/reference",
            ".grp",
            ".ti",
            ".transcripts.fa",
            ".seq",
            ".idx.fa",
            ".n2g.idx.fa",
        ),
    output:
        genes_results="results/rsem/{sample}.genes.results",
        isoforms_results="results/rsem/{sample}.isoforms.results",
        bam="results/rsem/{sample}.transcript.bam",
    params:
        # optional, specify if sequencing is paired-end
        paired_end=BAM2BAKR_FORMAT,
        # additional optional parameters to pass to rsem, for example,
        extra=config["rsem_quant_params"],
    log:
        "logs/RSEM/{sample}.log",
    conda:
        "../envs/rsem.yaml"
    threads: 20
    script:
        "../scripts/rsem/rsem-calc.py"
