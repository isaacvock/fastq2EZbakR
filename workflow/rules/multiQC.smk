### Compile all QC into one report
rule multiqc:
    input:
        expand(
            "results/fastqc/{sample}_r{read}_fastqc.{ext}",
            sample=SAMP_NAMES,
            read=READS,
            ext=["html", "zip"],
        ),
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