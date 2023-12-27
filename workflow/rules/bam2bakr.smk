### BIG PICTURE TO-DO
## 1) Make compatible with new config parameter names
## 2) Change names of rules and output directories to make more sense

### TO-DO
## 1) Clean up sort/filter function; maybe reduce number of output files
if config["bam2bakr"]:

    if config["remove_tags"]:
        
        # Remove tags from bam files that can break HTSeq
        rule remove_tags:
            input:
                input_bam=get_input_bams,
            output:
                output_bam="results/remove_tags/{sample}_no_jI_jM.bam",
            log:
                "logs/remove_tags/{sample}.log"
            conda:
                "../envs/full.yaml"
            script:
                "../scripts/bam2bakR/remove_tags.py"

        # Filter out multi-mappers and sort reads
        rule sort_filter:
            input:
                "results/remove_tags/{sample}_no_jI_jM.bam",
            output:
                "results/sf_reads/{sample}.s.sam",
                "results/sf_reads/{sample}_fixed_mate.bam",
                "results/sf_reads/{sample}.f.sam",
            log:
                "logs/sort_filter/{sample}.log"
            params: 
                shellscript=workflow.source_path("../scripts/bam2bakR/sort_filter.sh"),
                format=FORMAT
            threads: 8
            conda:
                "../envs/full.yaml"
            shell:
                """
                chmod +x {params.shellscript}
                {params.shellscript} {threads} {wildcards.sample} {input} {output} {params.format} 1> {log} 2>&1
                """

    else:

        # Filter out multi-mappers and sort reads
        rule sort_filter:
            input:
                get_input_bams
            output:
                "results/sf_reads/{sample}.s.sam",
                "results/sf_reads/{sample}_fixed_mate.bam",
                "results/sf_reads/{sample}.f.sam",
            log:
                "logs/sort_filter/{sample}.log"
            params: 
                shellscript=workflow.source_path("../scripts/bam2bakR/sort_filter.sh"),
                format=FORMAT
            threads: 8
            conda:
                "../envs/full.yaml"
            shell:
                """
                chmod +x {params.shellscript}
                {params.shellscript} {threads} {wildcards.sample} {input} {output} {params.format} 1> {log} 2>&1
                """


else:

    # Filter out multi-mappers and sort reads
    rule sort_filter:
        input:
            "results/align/{sample}.bam"
        output:
            "results/sf_reads/{sample}.s.bam",
            "results/sf_reads/{sample}_fixed_mate.bam",
            "results/sf_reads/{sample}.f.sam",
        log:
            "logs/sort_filter/{sample}.log"
        params: 
            shellscript=workflow.source_path("../scripts/bam2bakR/sort_filter.sh"),
            format=FORMAT
        threads: 8
        conda:
            "../envs/full.yaml"
        shell:
            """
            chmod +x {params.shellscript}
            {params.shellscript} {threads} {wildcards.sample} {input} {output} {params.format} 1> {log} 2>&1
            """



### TO-DO
## 1) Properly log standard out
# Calculate normalization scale factor to be applied to tracks        
if NORMALIZE:
    rule normalize:
        input:
            expand("results/sf_reads/{sample}.s.bam", sample = SAMP_NAMES)
        output:
            "results/normalization/scale"
        log:
            "logs/normalize/normalize.log"
        params:
            rscript=workflow.source_path("../scripts/bam2bakR/normalize.R"),
            spikename=config["spikename"]
        threads: 1
        conda:
            "../envs/full.yaml"
        shell:
            r"""
            chmod +x {params.rscript}
            {params.rscript} --dirs ./results/htseq/ --spikename {params.spikename} 1> {log} 2>&1
            mv scale {output}
            """
else:
    rule normalize:
        input:
            expand("results/sf_reads/{sample}.s.bam", sample = SAMP_NAMES)
        output:
            "results/normalization/scale"
        log:
            "logs/normalize/normalize.log"
        threads: 1
        conda:
            "../envs/full.yaml"
        shell:
            """
            touch {output}
            """

# Index genome fasta file for snp calling
rule genome_index:
    input:
        str(config["genome"])
    output:
        get_index_name()
    log:
        "logs/genome_index/genome-faidx.log",
    threads: 1
    conda:
        "../envs/index.yaml"
    script:
        "../scripts/bam2bakR/genome-faidx.py"

## TO-DO
# 1) Allow users to provide custom SNP file
# Identify SNPs to be accounted for when counting mutations
rule call_snps:
    input:
        str(config["genome"]),
        get_index_name(),
        expand("results/sf_reads/{ctl}.s.bam", ctl = CTL_NAMES)
    params:
        nctl = nctl,
        shellscript = workflow.source_path("../scripts/bam2bakR/call_snps.sh"),
    output:
        "results/snps/snp.txt",
        "results/snps/snp.vcf",
        temp("results/snps/mkdir.txt")
    log:
        "logs/call_snps/ctl_samps.log"
    threads: 20
    conda:
        "../envs/full.yaml"
    shell:
        """
        chmod +x {params.shellscript}
        {params.shellscript} {threads} {params.nctl} {output} {input} 1> {log} 2>&1
        """

# TO-DO:
# 1) Add mutation position optimizations and functionality
# Count mutations 
rule cnt_muts:
    input:
        "results/sf_reads/{sample}.s.bam",
        "results/snps/snp.txt"
    params:
        format = FORMAT,
        minqual = config["minqual"],
        mut_tracks = config["mut_tracks"],
        strand = STRAND,
        shellscript = workflow.source_path("../scripts/bam2bakR/mut_call.sh"),
        pythonscript = workflow.source_path("../scripts/bam2bakR/mut_call.py"),
        awkscript = workflow.source_path("../scripts/bam2bakR/fragment_sam.awk"),
        mutpos = config["mutpos"]
    output:
        "results/counts/{sample}_counts.csv.gz",
        temp("results/counts/{sample}_check.txt")
    log:
        "logs/cnt_muts/{sample}.log"
    threads: 32
    conda:
        "../envs/full.yaml"
    shell:
        """
        chmod +x {params.shellscript}
        chmod +x {params.pythonscript}
        chmod +x {params.awkscript}
        {params.shellscript} {threads} {wildcards.sample} {input} {output} {params.minqual} {params.mut_tracks} {params.format} {params.strand} {params.pythonscript} {params.awkscript} {params.mutpos} 1> {log} 2>&1
        """

# Merge mutation counts with feature assignment
rule merge_features_and_muts:
    input:
        get_merge_input
    output:
        "results/merge_features_and_muts/{sample}_counts.csv.gz"
    params:
        genes_included = config["features"]["genes"],
        exons_included = config["features"]["exons"],
        exonbins_included = config["features"]["exonbins"],
        transcripts_included = config["features"]["transcripts"],
        rscript = workflow.source_path("../scripts/bam2bakr/merge_features_and_muts.R")
    log:
        "logs/merge_features_and_muts/{sample}.log"
    threads: 1
    conda: 
        "../envs/full.yaml"
    shell:
        """
        chmod +x {params.rscript}

        {params.rscript} -g {params.gene_included} -e {params.exons_included} -b {params.exonbins_included} \
        -t {params.transcripts_included} -o ./results/merge_features_and_muts/{wildcards.sample}_counts.csv 1> {log} 2>&1

        pigz -p {threads} ./results/merge_features_and_muts/{wildcards.sample}_counts.csv
        """


# Make cB (and potentially cU) file
if config["mutpos"]:

    rule makecB:
        input:
            expand("results/merge_features_and_muts/{sample}_counts.csv.gz", sample=SAMP_NAMES)
        output:
            cB = "results/cB/cB.csv.gz",
            mutpos = "results/cB/mutpos.csv.gz",
            mutposfilter = "results/cB/mutpos_filtered.csv.gz",
        params:
            shellscript = workflow.source_path("../scripts/bam2bakR/master.sh"),
            keepcols = keepcols,
            mut_tracks = config["mut_tracks"],
            mut_pos = config["mutpos"],
            min_pos_coverage = config["min_pos_coverage"],
            max_pos_coverage = config["max_pos_coverage"],
            relative_counts_dir = COUNTS_DIR_RELATIVE,
        log:
            "logs/makecB/master.log"
        threads: 20
        conda:
            "../envs/full.yaml"
        shell:
            """
            chmod +x {params.shellscript}
            {params.shellscript} {threads} {output.cB} {params.keepcols} {params.mut_tracks} \
            {params.relative_counts_dir} {params.mut_pos} {params.min_pos_coverage} {output.mutpos} \
            {output.mutposfilter} {params.max_pos_coverage} 1> {log} 2>&1
            """

else:

    rule makecB:
        input:
            expand("results/merge_features_and_muts/{sample}_counts.csv.gz", sample=SAMP_NAMES)
        output:
            cB = "results/cB/cB.csv.gz"
        params:
            shellscript = workflow.source_path("../scripts/bam2bakR/master.sh"),
            keepcols = keepcols,
            mut_tracks = config["mut_tracks"],
            mut_pos = config["mutpos"],
            relative_counts_dir = COUNTS_DIR_RELATIVE,
        log:
            "logs/makecB/master.log"
        threads: 20
        conda:
            "../envs/full.yaml"
        shell:
            """
            chmod +x {params.shellscript}
            {params.shellscript} {threads} {output.cB} {params.keepcols} {params.mut_tracks} \
            {params.relative_counts_dir} {params.mut_pos} 1> {log} 2>&1
            """


# Make color-coded tracks
rule maketdf:
    input:
        "results/counts/{sample}_counts.csv.gz",
        "results/sf_reads/{sample}.s.bam",
	    "results/normalization/scale"
    output:
        temp("results/tracks/{sample}_success.txt"),
        expand("results/tracks/{{sample}}.{mut}.{id}.{strand}.tdf", mut=Mutation_Types, id=[0,1,2,3,4,5], strand = ['pos', 'min'])
    params:
        shellscript = workflow.source_path("../scripts/bam2bakR/tracks.sh"),
        pythonscript = workflow.source_path("../scripts/bam2bakR/count_to_tracks.py"),
        mut_tracks=config["mut_tracks"],
        genome=config["genome"],
        WSL=config["WSL"],
        normalize=config["normalize"]
    log:
        "logs/maketdf/{sample}.log"
    threads: 20
    conda:
        "../envs/full.yaml"
    shell:
        """
        chmod +x {params.shellscript}
        chmod +x {params.pythonscript}
        {params.shellscript} {threads} {wildcards.sample} {input} {params.mut_tracks} {params.genome} {params.WSL} {params.normalize} {params.pythonscript} {output} 1> {log} 2>&1
        """

