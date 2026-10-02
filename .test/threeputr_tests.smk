rule test_annotate_threepUTR_gtf:
    input:
        gtf="data/annotation/genome.gtf",
    output:
        "results/threeputr_tests/annotation.ok",
    log:
        "logs/threeputr_tests/annotation.log",
    conda:
        "../workflow/envs/Rbio.yaml"
    threads: 1
    params:
        rscript=workflow.source_path("scripts/test_annotate_3pgtf.R"),
        helper=workflow.source_path("../workflow/scripts/threeputrs/annotate_3pgtf.R"),
    shell:
        """
        Rscript {params.rscript:q} {params.helper:q} {input.gtf:q} {output:q} > {log:q} 2>&1
        """


rule test_make_threepUTR_gtf:
    output:
        "results/threeputr_tests/de_novo.ok",
    log:
        "logs/threeputr_tests/de_novo.log",
    conda:
        "../workflow/envs/Rbio.yaml"
    threads: 1
    params:
        rscript=workflow.source_path("scripts/test_make_3pgtf.R"),
        caller=workflow.source_path("../workflow/scripts/threeputrs/make_3pgtf.R"),
    shell:
        """
        Rscript {params.rscript:q} {params.caller:q} {output:q} > {log:q} 2>&1
        """


rule check_threeputr_output:
    input:
        gtf="annotations/annotated_threepUTR_annotation.gtf",
        cb="results/cB/cB.csv.gz",
        counts=expand(
            "results/featurecounts_3utr/{sample}.featureCounts",
            sample=["WT_1", "WT_2"],
        ),
        assignments=expand(
            "results/featurecounts_3utr/{sample}.s.bam.featureCounts",
            sample=["WT_1", "WT_2"],
        ),
    output:
        "results/threeputr_tests/integration.ok",
    log:
        "logs/threeputr_tests/integration.log",
    conda:
        "../workflow/envs/Rbio.yaml"
    threads: 1
    params:
        rscript=workflow.source_path("scripts/check_threeputr_output.R"),
    shell:
        """
        Rscript {params.rscript:q} {input.gtf:q} {input.cb:q} {output:q} \
            {input.counts:q} {input.assignments:q} > {log:q} 2>&1
        """
