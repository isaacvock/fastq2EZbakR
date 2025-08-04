### Compile all QC into one report
rule multiqc:
    input:
        MULTIQC_INPUT,
    output:
        "results/multiqc/multiqc_report.html",
    params:
        extra=config.get("multiqc_extra", ""),
    log:
        "logs/multiqc/multiqc.log",
    wrapper:
        "v3.5.3/bio/multiqc"


# ### StackOverflow solution
# rule multiqc:
#     input:
#         MULTIQC_INPUT,
#     output:
#         "results/multiqc/multiqc_report.html",
#     params:
#         extra=config.get("multiqc_extra", ""),
#         output_dir=lambda wildcards, output: os.path.dirname(output[0]),
#         output_file_name=lambda wildcards, output: os.path.basename(output[0]),
#         input_directories=lambda wildcards, input: set(os.path.dirname(fp) for fp in input)
#     log:
#         "logs/multiqc/multiqc.log",
#     conda:
#         "../envs/multiqc.yaml"
#     shell:
#         "multiqc {params.extra} --force "
#         "-o {params.output_dir} -n {params.output_file_name} "
#         "{params.input_directories} &> {log}"