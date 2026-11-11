# Flatten annotation with DEXseq script
rule flatten:
    input:
        gtf=config["annotation"],
    output:
        flatgtf="results/raw_flattened/flat_genome.gtf",
    log:
        "logs/get-genome.log",
    conda:
        "../envs/flatten.yaml"
    threads: 1
    script:
        "../scripts/flatten/dexseq_prepare_annotation.py"


# Add exon ID to flattened annotation
rule add_exon:
    input:
        "results/raw_flattened/flat_genome.gtf",
    output:
        config["flat_annotation"],
    log:
        "logs/add_exon.log",
    params:
        shellscript=workflow.source_path("../scripts/flatten/exon_ID.sh"),
    conda:
        "../envs/full.yaml"
    threads: 1
    shell:
        """
        bash {params.shellscript:q} {input:q} {output:q} 1> {log:q} 2>&1
        """
