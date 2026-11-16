rule test_snp_sample_selection:
    input:
        bams=expand(
            "data/bams/{name}.bam",
            name=["WT_nos4U", "WT_replicate_1", "WT_replicate_2"],
        ),
    output:
        "results/snp_tests/snp_selection.ok",
    log:
        "logs/snp_tests/snp_selection.log",
    conda:
        "../workflow/envs/full.yaml"
    threads: 2
    params:
        test_script=workflow.source_path("scripts/test_snp_sample_selection.py"),
        selector=workflow.source_path(
            "../workflow/scripts/bam2bakR/select_snp_samples.py"
        ),
        common_smk=workflow.source_path("../workflow/rules/common.smk"),
        bam_dir="data/bams",
    shell:
        """
        python {params.test_script:q} {params.selector:q} {params.common_smk:q} {params.bam_dir:q} {output:q} > {log:q} 2>&1
        """
