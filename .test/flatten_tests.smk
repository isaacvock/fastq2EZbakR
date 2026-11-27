rule test_flatten_annotation:
    input:
        gtf="data/annotation/genome.gtf",
    output:
        "results/flatten_tests/flatten.ok",
    log:
        "logs/flatten_tests/flatten.log",
    conda:
        "../workflow/envs/full.yaml"
    threads: 1
    params:
        test_script=workflow.source_path("scripts/test_flatten_annotation.py"),
        flattener=workflow.source_path(
            "../workflow/scripts/flatten/flatten_annotation.py"
        ),
    shell:
        """
        python {params.test_script:q} {params.flattener:q} {input.gtf:q} {output:q} > {log:q} 2>&1
        """
