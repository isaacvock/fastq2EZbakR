# Flatten annotation into DEXSeq-style exonic bins with unique exon_id attributes
rule flatten:
    input:
        gtf=config["annotation"],
    output:
        flatgtf=config["flat_annotation"],
    params:
        pythonscript=workflow.source_path("../scripts/flatten/flatten_annotation.py"),
    log:
        "logs/flatten/flatten.log",
    conda:
        "../envs/full.yaml"
    threads: 1
    shell:
        """
        python {params.pythonscript:q} -g {input.gtf:q} -o {output.flatgtf:q} 1> {log:q} 2>&1
        """
