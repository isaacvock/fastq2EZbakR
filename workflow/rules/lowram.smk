### STRATEGY
# Step in pipeline whose RAM usage is a function of
# file size is the merge_features_and_muts rule. Also,
# cB files can similarly be created without using excessive
# RAM. To address these problems, I need to:
# 1) Sort the files to be merged by the read name column
# 2) Iterate through each file row by row and perform
# a left join with the mutation counts.
# 3) Sort files by columns that will be kept in cB
# 4) Iterate again this time counting the number of identical
# instances of the set of columns to keep, in order to
# summarize over these columns. The output of this can be passed
# as is to the current makecB script

### NO LONGER NECESSARY GIVEN DUCKDB FOR MERGING!!! (version 0.7.0)
# Keeping here for legacy reasons and backwards compatibility

### Sorting rules


# Sort by read name (qname)
rule sort_mutcounts_by_qname:
    input:
        "results/counts/{sample}_counts.csv.gz",
    output:
        sortout=temp("results/sort_mutcounts_by_qname/{sample}_counts.csv"),
        decomp=temp("results/sort_mutcounts_by_qname/{sample}_decompressed_counts.csv"),
    threads: 20
    conda:
        "../envs/full.yaml"
    log:
        "logs/sort_mutcounts_by_qname/{sample}.log",
    shell:
        """
        ### GOAL: Sort but preserve header
        gzip -d -c {input} > {output.decomp}

        head -n 1 {output.decomp} > {output.sortout}
        
        tail -n +2 {output.decomp} | sort --parallel={threads} -k1 -V >> {output.sortout}
        """


rule sort_fcgene_by_qname:
    input:
        "results/featurecounts_genes/{sample}.s.bam.featureCounts",
    output:
        temp("results/sort_fcgene_by_qname/{sample}.featureCounts"),
    threads: 20
    conda:
        "../envs/full.yaml"
    log:
        "logs/sort_fcgene_by_qname/{sample}.log",
    shell:
        """
        sort --parallel={threads} -k1 -V {input} > {output}
        """


rule sort_fcexon_by_qname:
    input:
        "results/featurecounts_exons/{sample}.s.bam.featureCounts",
    output:
        temp("results/sort_fcexon_by_qname/{sample}.featureCounts"),
    threads: 20
    conda:
        "../envs/full.yaml"
    log:
        "logs/sort_fcexon_by_qname/{sample}.log",
    shell:
        """
        sort --parallel={threads} -k1 -V {input} > {output}
        """


rule sort_junction_by_qname:
    input:
        "results/read_to_junctions/{sample}.csv.gz",
    output:
        sortout=temp("results/sort_junction_by_qname/{sample}_counts.csv"),
        decomp=temp("results/sort_junction_by_qname/{sample}_decompressed_counts.csv"),
    threads: 20
    conda:
        "../envs/full.yaml"
    log:
        "logs/sort_junction_by_qname/{sample}.log",
    shell:
        """
        ### GOAL: Sort but preserve header
        gzip -d -c {input} > {output.decomp}

        head -n 1 {output.decomp} > {output.sortout}
        
        tail -n +2 {output.decomp} | sort  --parallel={threads} -k1 -V >> {output.sortout}
        """


rule sort_fcee_by_qname:
    input:
        "results/featurecounts_eej/{sample}.s.bam.featureCounts",
    output:
        temp("results/sort_fcee_by_qname/{sample}.featureCounts"),
    threads: 20
    conda:
        "../envs/full.yaml"
    log:
        "logs/sort_fcee_by_qname/{sample}.log",
    shell:
        """
        sort --parallel={threads} -k1 -V {input} > {output}
        """


rule sort_fcei_by_qname:
    input:
        "results/featurecounts_eij/{sample}.s.bam.featureCounts",
    output:
        temp("results/sort_fcei_by_qname/{sample}.featureCounts"),
    threads: 20
    conda:
        "../envs/full.yaml"
    log:
        "logs/sort_fcei_by_qname/{sample}.log",
    shell:
        """
        sort --parallel={threads} -k1 -V {input} > {output}
        """


rule sort_fctranscript_by_qname:
    input:
        "results/featurecounts_transcripts/{sample}.s.bam.featureCounts",
    output:
        temp("results/sort_fctranscript_by_qname/{sample}.featureCounts"),
    threads: 20
    conda:
        "../envs/full.yaml"
    log:
        "logs/sort_fctranscript_by_qname/{sample}.log",
    shell:
        """
        sort --parallel={threads} -k1 -V {input} > {output}
        """


rule sort_fcexonbin_by_qname:
    input:
        "results/featurecounts_exonbins/{sample}.s.bam.featureCounts",
    output:
        temp("results/sort_fcexonbin_by_qname/{sample}.featureCounts"),
    threads: 20
    conda:
        "../envs/full.yaml"
    log:
        "logs/sort_fcexonbin_by_qname/{sample}.log",
    shell:
        """
        sort --parallel={threads} -k1 -V {input} > {output}
        """


rule sort_bamtranscript_by_qname:
    input:
        "results/read_to_transcripts/{sample}.csv",
    output:
        temp("results/sort_bamtranscript_by_qname/{sample}.csv"),
    threads: 20
    conda:
        "../envs/full.yaml"
    log:
        "logs/sort_bamtranscript_by_qname/{sample}.log",
    shell:
        """
        ### GOAL: Sort but preserve header
        head -n 1 {input} > {output}
        
        tail -n +2 {input} | sort --parallel={threads} -k1 -V >> {output}
        """


### MERGING AND SORTING FILES

if config["lowRAM"]:

    rule lowram_merge_features_and_counts:
        input:
            get_lowram_merge_input,
        output:
            "results/lowram_merge_features_and_counts/{sample}.csv",
        params:
            genes=config["features"]["genes"],
            exons=config["features"]["exons"],
            junctions=config["features"]["junctions"],
            exonic_bins=config["features"]["exonic_bins"],
            ee=config["features"]["eej"],
            ei=config["features"]["eij"],
            transcripts=config["features"]["transcripts"],
            bamfile_transcripts=config["strategies"]["Transcripts"],
        threads: 1
        conda:
            "../envs/full.yaml"
        log:
            "logs/lowram_merge_features_and_counts/{sample}.log",
        script:
            "../scripts/lowram/lowram_join.py"


rule sort_merged_files:
    input:
        "results/lowram_merge_features_and_counts/{sample}.csv",
    output:
        temp("results/sort_merged_files/{sample}.csv"),
    params:
        sortparams=SORTPARAMS,
    threads: 20
    conda:
        "../envs/full.yaml"
    log:
        "logs/sort_merged_files/{sample}.log",
    shell:
        """
        head -n 1 {input} > {output}
        
        tail -n +2 {input} | sort -t, --parallel={threads} {params.sortparams} >> {output}
        """


### SUMMARISE MERGED FILES


if config["final_output"]["arrow"]:

    rule lowram_summarise:
        input:
            "results/sort_merged_files/{sample}.csv",
        output:
            "results/arrow_dataset/sample={sample}/part-0.csv",
        params:
            cols_to_sum=COLS_TO_SUM,
        threads: 1
        conda:
            "../envs/full.yaml"
        log:
            "logs/lowram_summarise/{sample}.log",
        script:
            "../scripts/lowram/lowram_summarise.py"

else:

    rule lowram_summarise:
        input:
            "results/sort_merged_files/{sample}.csv",
        output:
            temp("results/lowram_summarise/{sample}.csv"),
        params:
            cols_to_sum=COLS_TO_SUM,
        threads: 1
        conda:
            "../envs/full.yaml"
        log:
            "logs/lowram_summarise/{sample}.log",
        script:
            "../scripts/lowram/lowram_summarise.py"
