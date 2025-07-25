### Compile all QC into one report
rule multiqc:
    input:
        MULTIQC_INPUT,
    output:
        "results/multiqc/multiqc_report.html",
    params:
        extra=config.get("multiqc_extra", "")
    log:
        "logs/multiqc/multiqc.log",
    conda:
        "../envs/multiqc.yaml",
    script:
        "../scripts/multiQC.py"