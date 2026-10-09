### Little bit of Snakemake smoving here. If downloading
### fastq files, SAMP_NAMES is set equal to SRA accession IDs
### provided in config rather than the sample IDs specified under
### samples in the config. This allows for the {accession} wild
### card to be determined.
if TWO_FASTQS:

    rule download_fastq:
        output:
            "results/download_fastq/{accession}_1.fastq",
            "results/download_fastq/{accession}_2.fastq",
        params:
            extra=config["fasterq_dump_extras"],
        threads: 12
        wrapper:
            "v3.3.6/bio/sra-tools/fasterq-dump"

else:

    rule download_fastq:
        output:
            "results/download_fastq/{accession}.fastq",
        params:
            extra=config["fasterq_dump_extras"],
        threads: 12
        wrapper:
            "v3.3.6/bio/sra-tools/fasterq-dump"
