rule test_workflow_validator:
    output:
        "validation_checks/self_test.ok",
    log:
        "validation_checks/self_test.log",
    conda:
        "../workflow/envs/full.yaml"
    threads: 1
    params:
        validator=workflow.source_path("scripts/check_workflow_output.py"),
        configs_dir="configs",
    shell:
        """
        python {params.validator:q} \
            --self-test {params.configs_dir:q} \
            --output-sentinel {output:q} > {log:q} 2>&1
        """


rule validate_and_clean:
    output:
        "validation_checks/last_check.ok",
    log:
        "validation_checks/last_check.log",
    conda:
        "../workflow/envs/full.yaml"
    threads: 1
    params:
        validator=workflow.source_path("scripts/check_workflow_output.py"),
        target_config=config.get("target_config", "configs/transcript_config.yaml"),
        clean_flag="--clean" if config.get("clean", True) else "",
        clean_idx_flag="--clean-indices" if config.get("clean_indices", False) else "",
    shell:
        """
        python {params.validator:q} \
            --config {params.target_config:q} \
            --workdir . \
            {params.clean_flag} \
            {params.clean_idx_flag} \
            --output-sentinel {output:q} > {log:q} 2>&1
        """
