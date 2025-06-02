######################################################################################
##### ALIGNMENT WITH STAR
######################################################################################


# Make modified annotation if necessary
# For "regressing out" pre-mRNA signal
rule modify_annotation:
    input:
        gtf=config["annotation"],
    output:
        mod_gtf="results/modify_annotation/modified_annotation.gtf",
    params:
        rscript=workflow.source_path("../scripts/alignment/modify_annotation.R"),
    threads: 1
    conda:
        "../envs/modify.yaml"
    log:
        "logs/modify_annotation/modify_annotation.log",
    shell:
        """
        chmod +x {params.rscript}
        {params.rscript} -g {input.gtf} -o {output.mod_gtf} 1> {log} 2>&1
        """


if config["aligner"] == "star":
    if make_index:

        # Build STAR index
        rule index:
            input:
                fasta=config["genome"],
                gtf=AandQ_ANNOTATION,
            output:
                outdir=directory(config["indices"]),
                outfile=star_indices,
            params:
                extra=config["star_index_params"],
            log:
                "logs/index/star_index_genome.log",
            conda:
                "../envs/star.yaml"
            threads: 24
            script:
                "../scripts/alignment/star-index.py"

    # Align with STAR

    if config["PE"]:

        rule align:
            input:
                fq1=get_fastq_r1,
                fq2=get_fastq_r2,
                index=config["indices"],
            output:
                aln="results/align/{sample}.bam",
                sj="results/align/{sample}-SJ.out.tab",
                log="results/align/{sample}-Log.out",
                log_progress="results/align/{sample}-Log.progress.out",
                log_final="results/align/{sample}-Log.final.out",
                aln_tx="results/align/{sample}-Aligned.toTranscriptome.out.bam",
            log:
                "logs/align/{sample}_star.log",
            params:
                reads_per_gene=lambda wc: "GeneCounts" in config["star_align_params"],
                chim_junc=lambda wc: "--chimOutType Junctions"
                in config["star_align_params"],
                idx=lambda wc, input: input.index,
                tmpdir="data/tmp_star/{sample}",
                extra=STAR_EXTRA,
                out_reads_per_gene="results/align/{sample}-ReadsPerGene.out.tab",
                out_chim_junc="results/align/{sample}-Chimeric.out.junction",
            conda:
                "../envs/star.yaml"
            threads: 24
            script:
                "../scripts/alignment/star-align.py"

    else:

        rule align:
            input:
                fq1=get_fastq_r1,
                index=config["indices"],
            output:
                aln="results/align/{sample}.bam",
                sj="results/align/{sample}-SJ.out.tab",
                log="results/align/{sample}-Log.out",
                log_progress="results/align/{sample}-Log.progress.out",
                log_final="results/align/{sample}-Log.final.out",
                aln_tx="results/align/{sample}-Aligned.toTranscriptome.out.bam",
            log:
                "logs/align/{sample}_star.log",
            params:
                reads_per_gene=lambda wc: "GeneCounts" in config["star_align_params"],
                chim_junc=lambda wc: "--chimOutType Junctions"
                in config["star_align_params"],
                idx=lambda wc, input: input.index,
                tmpdir="data/tmp_star/{sample}",
                extra=STAR_EXTRA,
                out_reads_per_gene="results/align/{sample}-ReadsPerGene.out.tab",
                out_chim_junc="results/align/{sample}-Chimeric.out.junction",
            conda:
                "../envs/star.yaml"
            threads: 24
            script:
                "../scripts/alignment/star-align.py"


######################################################################################
##### ALIGNMENT WITH HISAT2
######################################################################################


if config["aligner"] == "hisat2":
    if make_index:
        ### Add annotated splice junctions to index
        if config["annotation"]:

            # Get exons from annotation using hisat2's custom python script
            rule get_exons:
                input:
                    annotation=config["annotation"],
                output:
                    "results/get_exons/exons.exon",
                log:
                    "logs/get_exons/exons.log",
                conda:
                    "../envs/hisat2.yaml"
                threads: 1
                shell:
                    "hisat2_extract_exons.py {input.annotation} 1> {output} 2> {log}"

            # Get splice sites from annotation using hisat2's custom python script
            rule get_ss:
                input:
                    annotation=config["annotation"],
                output:
                    "results/get_ss/splice_sites.ss",
                log:
                    "logs/get_ss/ss.log",
                conda:
                    "../envs/hisat2.yaml"
                threads: 1
                shell:
                    "hisat2_extract_splice_sites.py {input.annotation} 1> {output} 2> {log}"

            # Build hisat2 index
            rule index:
                input:
                    fasta=config["genome"],
                    annotation=config["annotation"],
                    ss="results/get_ss/splice_sites.ss",
                    exons="results/get_exons/exons.exon",
                output:
                    directory(config["indices"]),
                params:
                    prefix=HISAT2_BASE,
                    extra="{} {}".format(
                        "--ss results/get_ss/splice_sites.ss --exon results/get_exons/exons.exon",
                        config["hisat2_index_params"],
                    ),
                log:
                    "logs/index/hisat2_index.log",
                threads: 20
                wrapper:
                    "v2.6.0/bio/hisat2/index"

        else:

            # Build hisat2 index
            rule index:
                input:
                    fasta=config["genome"],
                    annotation=config["annotation"],
                    ss="results/get_ss/splice_sites.ss",
                    exons="results/get_exon/exons.exon",
                output:
                    directory(config["indices"]),
                params:
                    prefix=HISAT2_BASE,
                    extra=config["hisat2_index_params"],
                log:
                    "logs/index/hisat2_index.log",
                threads: 20
                wrapper:
                    "v2.6.0/bio/hisat2/index"

    # Align with hisat2
    rule align:
        input:
            reads=get_hisat2_reads,
            idx=config["indices"],
        output:
            "results/align/{sample}.bam",
        log:
            "logs/align/{sample}_hisat2.log",
        params:
            extra="{} {}".format(HISAT2_STRANDEDNESS, config["hisat2_align_params"]),
        threads: 20
        wrapper:
            "v2.6.0/bio/hisat2/align"
