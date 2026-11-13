### BAM FILTERING, MUTATION COUNTING, MERGING, AND TRACK GENERATION

if config.get("bam2bakr", False):
    if config.get("remove_tags", False):

        # Remove tags from bam files if requested
        rule remove_tags:
            input:
                input_bam=get_input_bams,
            output:
                output_bam=temp("results/remove_tags/{sample}_no_jI_jM.bam"),
            log:
                "logs/remove_tags/{sample}.log",
            conda:
                "../envs/full.yaml"
            script:
                "../scripts/bam2bakR/remove_tags.py"

        # Filter out multi-mappers and queryname-sort reads once
        rule sort_filter:
            input:
                "results/remove_tags/{sample}_no_jI_jM.bam",
            output:
                temp("results/sf_reads/{sample}.s.bam"),
            log:
                "logs/sort_filter/{sample}.log",
            params:
                shellscript=workflow.source_path("../scripts/bam2bakR/sort_filter.sh"),
            threads: 8
            conda:
                "../envs/full.yaml"
            shell:
                """
                bash {params.shellscript:q} {threads} {wildcards.sample:q} {input:q} {output:q} 1> {log:q} 2>&1
                """

    else:

        # Filter out multi-mappers and queryname-sort reads once
        rule sort_filter:
            input:
                get_input_bams,
            output:
                temp("results/sf_reads/{sample}.s.bam"),
            log:
                "logs/sort_filter/{sample}.log",
            params:
                shellscript=workflow.source_path("../scripts/bam2bakR/sort_filter.sh"),
            threads: 8
            conda:
                "../envs/full.yaml"
            shell:
                """
                bash {params.shellscript:q} {threads} {wildcards.sample:q} {input:q} {output:q} 1> {log:q} 2>&1
                """

else:

    # Filter out multi-mappers and queryname-sort reads once
    rule sort_filter:
        input:
            "results/align/{sample}.bam",
        output:
            temp("results/sf_reads/{sample}.s.bam"),
        log:
            "logs/sort_filter/{sample}.log",
        params:
            shellscript=workflow.source_path("../scripts/bam2bakR/sort_filter.sh"),
        threads: 8
        conda:
            "../envs/full.yaml"
        shell:
            """
            bash {params.shellscript:q} {threads} {wildcards.sample:q} {input:q} {output:q} 1> {log:q} 2>&1
            """


# Shared coordinate-sorted BAM produced on demand for call_snps, maketdf, and de novo 3'UTR calling
rule sort_coord_bam:
    input:
        "results/sf_reads/{sample}.s.bam",
    output:
        bam=temp("results/coord_bam/{sample}.bam"),
        bai=temp("results/coord_bam/{sample}.bam.bai"),
    log:
        "logs/sort_coord_bam/{sample}.log",
    threads: 8
    conda:
        "../envs/full.yaml"
    shell:
        """
        samtools sort -@ {threads} -o {output.bam:q} {input:q} 1> {log:q} 2>&1
        samtools index -@ {threads} {output.bam:q} 1>> {log:q} 2>&1
        """


# Calculate normalization scale factor to be applied to tracks
if NORMALIZE:

    rule normalize:
        input:
            NORMALIZATION_INPUT,
        output:
            "results/normalization/scale",
        log:
            "logs/normalize/normalize.log",
        threads: 1
        params:
            rscript=workflow.source_path("../scripts/bam2bakR/normalize.R"),
            spikename=config["spikename"],
            extra=NORMALIZATION_EXTRA,
        conda:
            "../envs/normalize.yaml"
        shell:
            r"""
            Rscript {params.rscript:q} {params.extra} --output {output:q} --spikename {params.spikename} 1> {log:q} 2>&1
            """

else:

    rule normalize:
        input:
            expand("results/sf_reads/{sample}.s.bam", sample=SAMP_NAMES),
        output:
            "results/normalization/scale",
        log:
            "logs/normalize/normalize.log",
        threads: 1
        conda:
            "../envs/full.yaml"
        shell:
            """
            touch {output:q}
            """


# Index genome fasta file for snp calling
rule genome_index:
    input:
        str(config["genome"]),
    output:
        get_index_name(),
    log:
        "logs/genome_index/genome-faidx.log",
    threads: 1
    conda:
        "../envs/index.yaml"
    script:
        "../scripts/bam2bakR/genome-faidx.py"


# Identify mutational anomalies / SNPs to mask when counting mutations
if SNP_BAMS:

    rule call_snps:
        input:
            fasta=str(config["genome"]),
            fai=get_index_name(),
            bams=SNP_BAMS,
            bais=SNP_BAIS,
        params:
            pythonscript=workflow.source_path("../scripts/bam2bakR/call_snps.py"),
            minqual=config["minqual"],
            min_depth=SNP_MIN_DEPTH,
            min_alt=config["snp_min_alt"],
            min_freq=SNP_MIN_FREQ,
            err_rate=SNP_ERR_RATE,
            p_cutoff=config["snp_p_cutoff"],
            mpileup_extra_arg=SNP_MPILEUP_EXTRA_ARG,
        output:
            txt="results/snps/snp.txt",
            tsv="results/snps/snp_counts.tsv",
        log:
            "logs/call_snps/ctl_samps.log",
        threads: 20
        conda:
            "../envs/full.yaml"
        shell:
            """
            python {params.pythonscript:q} \
                --bam {input.bams:q} \
                --fasta {input.fasta:q} \
                --fai {input.fai:q} \
                --output-txt {output.txt:q} \
                --output-tsv {output.tsv:q} \
                --threads {threads} \
                --min-qual {params.minqual} \
                --min-depth {params.min_depth} \
                --min-alt {params.min_alt} \
                --min-freq {params.min_freq} \
                --err-rate {params.err_rate} \
                --p-cutoff {params.p_cutoff} \
                {params.mpileup_extra_arg} 1> {log:q} 2>&1
            """

else:

    rule call_snps:
        params:
            pythonscript=workflow.source_path("../scripts/bam2bakR/call_snps.py"),
        output:
            txt="results/snps/snp.txt",
            tsv="results/snps/snp_counts.tsv",
        log:
            "logs/call_snps/ctl_samps.log",
        threads: 1
        conda:
            "../envs/full.yaml"
        shell:
            """
            python {params.pythonscript:q} \
                --output-txt {output.txt:q} \
                --output-tsv {output.tsv:q} 1> {log:q} 2>&1
            """


# Count mutations directly via parallel BGZF virtual-offset streaming
if config.get("mutpos", False):

    rule cnt_muts:
        input:
            bam="results/sf_reads/{sample}.s.bam",
            snps="results/snps/snp.txt",
        params:
            format=FORMAT,
            minqual=config["minqual"],
            mut_tracks=config["mut_tracks"],
            strand=STRAND,
            pythonscript=workflow.source_path("../scripts/bam2bakR/mut_call.py"),
        output:
            counts=temp("results/counts/{sample}_counts.csv.gz"),
            cU=temp("results/counts/{sample}_cU.csv.gz"),
        log:
            "logs/cnt_muts/{sample}.log",
        threads: 32
        conda:
            "../envs/full.yaml"
        shell:
            """
            python {params.pythonscript:q} \
                -b {input.bam:q} \
                -o {output.counts:q} \
                --cUOutput {output.cU:q} \
                -t {threads} \
                --mutType {params.mut_tracks:q} \
                --minQual {params.minqual} \
                --SNPs {input.snps:q} \
                --strandedness {params.strand:q} \
                --reads {params.format:q} \
                --mutPos 1> {log:q} 2>&1
            """

else:

    rule cnt_muts:
        input:
            bam="results/sf_reads/{sample}.s.bam",
            snps="results/snps/snp.txt",
        params:
            format=FORMAT,
            minqual=config["minqual"],
            mut_tracks=config["mut_tracks"],
            strand=STRAND,
            pythonscript=workflow.source_path("../scripts/bam2bakR/mut_call.py"),
        output:
            counts=temp("results/counts/{sample}_counts.csv.gz"),
        log:
            "logs/cnt_muts/{sample}.log",
        threads: 32
        conda:
            "../envs/full.yaml"
        shell:
            """
            python {params.pythonscript:q} \
                -b {input.bam:q} \
                -o {output.counts:q} \
                -t {threads} \
                --mutType {params.mut_tracks:q} \
                --minQual {params.minqual} \
                --SNPs {input.snps:q} \
                --strandedness {params.strand:q} \
                --reads {params.format:q} 1> {log:q} 2>&1
            """


# Merge mutation counts with feature assignments via streaming DuckDB views
if config.get("final_output", {}).get("arrow", False):

    rule merge_features_and_muts:
        input:
            get_merge_input,
        output:
            cBout=temp("results/merge_features_and_muts/{sample}_cB.csv"),
            cUPout=temp("results/merge_features_and_muts/{sample}_cUP.csv"),
            Arrowout="results/arrow_dataset/sample={sample}/part-0.parquet",
        params:
            genes_included=config.get("features").get("genes", True),
            exons_included=config.get("features").get("exons", True),
            exonbins_included=config.get("features").get("exonic_bins", False),
            transcripts_included=config.get("features").get("transcripts", False),
            bamfiletranscripts_included=config.get("features").get("tec", False),
            eej_included=config.get("features").get("eej", False),
            eij_included=config.get("features").get("eij", False),
            threeputr_included=config.get("features").get("threeputr", False),
            starjunc_included=config.get("features").get("junctions", False),
            rscript=workflow.source_path(
                "../scripts/bam2bakR/merge_features_and_muts.R"
            ),
            muttypes=config.get("mut_tracks", "TC"),
            annotation=config.get("annotation"),
            makecB=config.get("final_output").get("cB", True),
            makecUP=config.get("final_output").get("cUP", False),
            MaxMem=config.get("max_merge_mem", "8G"),
        log:
            "logs/merge_features_and_muts/{sample}.log",
        threads: 32
        conda:
            "../envs/merge.yaml"
        shell:
            """
            Rscript {params.rscript:q} \
                -g {params.genes_included} -e {params.exons_included} -b {params.exonbins_included} \
                -t {params.transcripts_included} --frombam {params.bamfiletranscripts_included} -s {wildcards.sample:q} \
                -j {params.eej_included} --threeputr {params.threeputr_included} --starjunc {params.starjunc_included} --eij {params.eij_included} \
                --annotation {params.annotation:q} -c {output.cBout:q} -m {params.muttypes:q} \
                --makecB {params.makecB} --makecUP {params.makecUP} --makeArrow TRUE \
                --cUPoutput {output.cUPout:q} --Arrowoutput {output.Arrowout:q} \
                --MaxMem {params.MaxMem:q} --threads {threads} 1> {log:q} 2>&1
            """

else:

    rule merge_features_and_muts:
        input:
            get_merge_input,
        output:
            cBout=temp("results/merge_features_and_muts/{sample}_cB.csv"),
            cUPout=temp("results/merge_features_and_muts/{sample}_cUP.csv"),
        params:
            genes_included=config.get("features").get("genes", True),
            exons_included=config.get("features").get("exons", True),
            exonbins_included=config.get("features").get("exonic_bins", False),
            transcripts_included=config.get("features").get("transcripts", False),
            bamfiletranscripts_included=config.get("features").get("tec", False),
            eej_included=config.get("features").get("eej", False),
            eij_included=config.get("features").get("eij", False),
            threeputr_included=config.get("features").get("threeputr", False),
            starjunc_included=config.get("features").get("junctions", False),
            rscript=workflow.source_path(
                "../scripts/bam2bakR/merge_features_and_muts.R"
            ),
            muttypes=config.get("mut_tracks", "TC"),
            annotation=config.get("annotation"),
            makecB=config.get("final_output").get("cB", True),
            makecUP=config.get("final_output").get("cUP", False),
            MaxMem=config.get("max_merge_mem", "8G"),
        log:
            "logs/merge_features_and_muts/{sample}.log",
        threads: 32
        conda:
            "../envs/merge.yaml"
        shell:
            """
            Rscript {params.rscript:q} \
                -g {params.genes_included} -e {params.exons_included} -b {params.exonbins_included} \
                -t {params.transcripts_included} --frombam {params.bamfiletranscripts_included} -s {wildcards.sample:q} \
                -j {params.eej_included} --threeputr {params.threeputr_included} --starjunc {params.starjunc_included} --eij {params.eij_included} \
                --annotation {params.annotation:q} -c {output.cBout:q} -m {params.muttypes:q} \
                --makecB {params.makecB} --makecUP {params.makecUP} --makeArrow FALSE \
                --cUPoutput {output.cUPout:q} \
                --MaxMem {params.MaxMem:q} --threads {threads} 1> {log:q} 2>&1
            """


# Combine per-sample site-specific mutation tables (cU)
rule makecU:
    input:
        cUins=expand(
            "results/counts/{sample}_cU.csv.gz",
            sample=SAMP_NAMES,
        ),
    output:
        mutpos="results/cB/mutpos.csv.gz",
        mutposfilter="results/cB/mutpos_filtered.csv.gz",
    params:
        shellscript=workflow.source_path("../scripts/bam2bakR/makecU.sh"),
        min_pos_coverage=config["min_pos_coverage"],
        max_pos_coverage=config["max_pos_coverage"],
    log:
        "logs/makecU/makecU.log",
    threads: MAKECB_THREADS
    conda:
        "../envs/full.yaml"
    shell:
        """
        bash {params.shellscript:q} {threads} \
            {params.min_pos_coverage} {output.mutpos:q} \
            {output.mutposfilter:q} {params.max_pos_coverage} \
            {input.cUins:q} 1> {log:q} 2>&1
        """


rule makecB:
    input:
        cBins=CBINPUT,
    output:
        cB="results/cB/cB.csv.gz",
    log:
        "logs/makecB/makecB.log",
    threads: 8
    conda:
        "../envs/full.yaml"
    shell:
        """
        {{ head -n 1 {input.cBins[0]:q}; tail -q -n +2 {input.cBins:q}; }} | pigz -p {threads} > {output.cB:q} 2> {log:q}
        """


rule makecUP:
    input:
        cUPins=CUPINPUT,
    output:
        cUP="results/cUP/cUP.csv.gz",
    log:
        "logs/makecUP/makecUP.log",
    threads: 8
    conda:
        "../envs/full.yaml"
    shell:
        """
        {{ head -n 1 {input.cUPins[0]:q}; tail -q -n +2 {input.cUPins:q}; }} | pigz -p {threads} > {output.cUP:q} 2> {log:q}
        """


# Make color-coded tracks
rule maketdf:
    input:
        counts="results/counts/{sample}_counts.csv.gz",
        bam="results/coord_bam/{sample}.bam",
        bai="results/coord_bam/{sample}.bam.bai",
        scale="results/normalization/scale",
    output:
        sentinel=temp("results/tracks/{sample}_success.txt"),
        tdfs=expand(
            "results/tracks/{{sample}}.{mut}.{id}.{strand}.tdf",
            mut=Mutation_Types,
            id=[0, 1, 2, 3, 4, 5],
            strand=["pos", "min"],
        ),
    params:
        shellscript=workflow.source_path("../scripts/bam2bakR/tracks.sh"),
        pythonscript=workflow.source_path("../scripts/bam2bakR/count_to_tracks.py"),
        mut_tracks=config["mut_tracks"],
        genome=config["genome"],
        WSL=config["WSL"],
        normalize=config["normalize"],
        strand=STRAND,
    log:
        "logs/maketdf/{sample}.log",
    threads: 10
    conda:
        "../envs/full.yaml"
    shell:
        """
        bash {params.shellscript:q} {threads} {wildcards.sample:q} {input.counts:q} {input.bam:q} {input.scale:q} {params.mut_tracks:q} {params.genome:q} {params.WSL} {params.normalize} {params.pythonscript:q} {params.strand:q} {output.sentinel:q} 1> {log:q} 2>&1
        """
