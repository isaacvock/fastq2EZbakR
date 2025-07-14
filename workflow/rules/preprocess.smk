## Trim adapters
if config["PE"]:

    # Trim with fastp
    rule fastp:
        input:
            sample=get_input_fastqs,
        output:
            trimmed=temp(
                [
                    "results/trimmed/{sample}.1.fastq",
                    "results/trimmed/{sample}.2.fastq",
                ]
            ),
            # Unpaired reads separately
            unpaired1=temp("results/trimmed/{sample}.u1.fastq"),
            unpaired2=temp("results/trimmed/{sample}.u2.fastq"),
            failed=temp("results/trimmed/{sample}.failed.fastq"),
            html="results/reports/{sample}.html",
            json="results/reports/{sample}.json",
        log:
            "logs/fastp/{sample}.log",
        params:
            adapters=config["fastp_adapters"],
            extra=config["fastp_parameters"],
        threads: 8
        wrapper:
            "v2.2.1/bio/fastp"

else:

    # Trim with fastp
    rule fastp:
        input:
            sample=get_input_fastqs,
        output:
            trimmed=temp("results/trimmed/{sample}.1.fastq"),
            failed=temp("results/trimmed/{sample}.1.failed.fastq"),
            html="results/reports/{sample}.1.html",
            json="results/reports/{sample}.1.json",
        log:
            "logs/fastp/{sample}.log",
        params:
            adapters=config["fastp_adapters"],
            extra=config["fastp_parameters"],
        threads: 8
        wrapper:
            "v2.2.1/bio/fastp"


if config["skip_trimming"] and is_gz:

    # Decompression with pigz cannot be parallelized, so force use of 1 thread
    rule unzip:
        input:
            fastqs=get_input_fastqs,
        output:
            unzipped_fqs=temp(
                expand("results/unzipped/{{sample}}.{read}.fastq", read=READS)
            ),
        log:
            "logs/unzip/{sample}.log",
        conda:
            "../envs/full.yaml"
        threads: 1
        script:
            "../scripts/preprocess/pigz.py"


# Run fastqc on trimmed fastqs
rule fastqc:
    input:
        get_fastqc_read,
    output:
        html="results/fastqc/{sample}_r{read}.html",
        zip="results/fastqc/{sample}_r{read}_fastqc.zip",
    log:
        "logs/fastqc/{sample}_r{read}.log",
    params:
        extra=config["fastqc_params"],
    resources:
        mem_mb=9000,
    threads: 4
    wrapper:
        "v2.2.1/bio/fastqc"
