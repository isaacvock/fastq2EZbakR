### BIG PICTURE TO-DO
## 1) Make compatible with new config parameter names
## 2) Change names of rules and output directories to make more sense

### TO-DO
## 1) Clean up sort/filter function; maybe reduce number of output files
if config.get("bam2bakr", False):
    if config.get("remove_tags", False):

        # Remove tags from bam files that can break HTSeq
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

        # Filter out multi-mappers and sort reads
        rule sort_filter:
            input:
                "results/remove_tags/{sample}_no_jI_jM.bam",
            output:
                "results/sf_reads/{sample}.s.bam",
            log:
                "logs/sort_filter/{sample}.log",
            params:
                shellscript=workflow.source_path("../scripts/bam2bakR/sort_filter.sh"),
            threads: 8
            conda:
                "../envs/full.yaml"
            shell:
                """
                chmod +x {params.shellscript}
                {params.shellscript} {threads} {wildcards.sample} {input} {output} 1> {log} 2>&1
                """

    else:

        # Filter out multi-mappers and sort reads
        rule sort_filter:
            input:
                get_input_bams,
            output:
                "results/sf_reads/{sample}.s.bam",
            log:
                "logs/sort_filter/{sample}.log",
            params:
                shellscript=workflow.source_path("../scripts/bam2bakR/sort_filter.sh"),
            threads: 8
            conda:
                "../envs/full.yaml"
            shell:
                """
                chmod +x {params.shellscript}
                {params.shellscript} {threads} {wildcards.sample} {input} {output} 1> {log} 2>&1
                """

else:

    # Filter out multi-mappers and sort reads
    rule sort_filter:
        input:
            "results/align/{sample}.bam",
        output:
            "results/sf_reads/{sample}.s.bam",
        log:
            "logs/sort_filter/{sample}.log",
        params:
            shellscript=workflow.source_path("../scripts/bam2bakR/sort_filter.sh"),
        threads: 8
        conda:
            "../envs/full.yaml"
        shell:
            """
            chmod +x {params.shellscript}
            {params.shellscript} {threads} {wildcards.sample} {input} {output} 1> {log} 2>&1
            """


### TO-DO
## 1) Properly log standard out
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
            chmod +x {params.rscript}
            {params.rscript} {params.extra} --output {output} --spikename {params.spikename} 1> {log} 2>&1
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
            touch {output}
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


## TO-DO
# 1) Allow users to provide custom SNP file
# Identify SNPs to be accounted for when counting mutations
rule call_snps:
    input:
        str(config["genome"]),
        get_index_name(),
        expand("results/sf_reads/{ctl}.s.bam", ctl=CTL_NAMES),
    params:
        nctl=nctl,
        shellscript=workflow.source_path("../scripts/bam2bakR/call_snps.sh"),
        call_params=config.get("bcftools_call_params", ""),
        mpileup_params=config.get("bcftools_mpileup_params", ""),
    output:
        "results/snps/snp.txt",
        "results/snps/snp.vcf",
    log:
        "logs/call_snps/ctl_samps.log",
    threads: 20
    conda:
        "../envs/full.yaml"
    shell:
        """
        chmod +x {params.shellscript}
        {params.shellscript} {threads} {params.nctl} {output} "{params.mpileup_params}" "{params.call_params}" {input} 1> {log} 2>&1
        """


# TO-DO:
# 1) Add mutation position optimizations and functionality
# Count mutations
rule cnt_muts:
    input:
        "results/sf_reads/{sample}.s.bam",
        "results/snps/snp.txt",
    params:
        format=FORMAT,
        minqual=config["minqual"],
        mut_tracks=config["mut_tracks"],
        strand=STRAND,
        shellscript=workflow.source_path("../scripts/bam2bakR/mut_call.sh"),
        pythonscript=workflow.source_path("../scripts/bam2bakR/mut_call.py"),
        awkscript=workflow.source_path("../scripts/bam2bakR/fragment_sam.awk"),
        mutpos=config["mutpos"],
    output:
        "results/counts/{sample}_counts.csv.gz",
        temp("results/counts/{sample}_check.txt"),
    log:
        "logs/cnt_muts/{sample}.log",
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


if not config.get("lowRAM", False):

    # Merge mutation counts with feature assignment
    # Bit of a cheap hack here where I didn't want to deal with dynamically
    # deciding the output, so I just create all files by default, with the ones
    # not requested by the user being temporary empty files.
    rule merge_features_and_muts:
        input:
            get_merge_input,
        output:
            output="results/merge_features_and_muts/{sample}_counts.csv.gz",
            cBout=temp("results/merge_features_and_muts/{sample}_cB.csv"),
            cUPout=temp("results/merge_features_and_muts/{sample}_cUP.csv"),
            Arrowout="results/arrow_dataset/sample={sample}/part-0.parquet",
            duckdb="results/merge_features_and_muts/duckdb/{sample}.duckdb",
        params:
            genes_included=config.get("features").get("genes", True),
            exons_included=config.get("features").get("exons", True),
            exonbins_included=config.get("features").get("exonic_bins", False),
            transcripts_included=config.get("features").get("transcripts", False),
            bamfiletranscripts_included=config.get("features").get("tec", False),
            eej_included=config.get("features").get("eej", False),
            eij_included=config.get("features").get("eij", False),
            threeputr_included=config.get("features").get("threeputr", False),
            sc_included=config.get("features").get("sc", False),
            starjunc_included=config.get("features").get("junctions", False),
            rscript=workflow.source_path(
                "../scripts/bam2bakR/merge_features_and_muts.R"
            ),
            muttypes=config.get("mut_tracks", "TC"),
            annotation=config.get("annotation"),
            makecB=config.get("final_output").get("cB", True),
            makecUP=config.get("final_output").get("cUP", False),
            makeArrow=config.get("final_output").get("arrow", True),
            MaxMem=config.get("max_merge_mem", "8G"),
        log:
            "logs/merge_features_and_muts/{sample}.log",
        threads: 32
        conda:
            "../envs/merge.yaml"
        shell:
            """
            chmod +x {params.rscript}

            {params.rscript} -g {params.genes_included} -e {params.exons_included} -b {params.exonbins_included} \
            -t {params.transcripts_included} --frombam {params.bamfiletranscripts_included} -o {output.output} -s {wildcards.sample} \
            -j {params.eej_included} --threeputr {params.threeputr_included} --sc {params.sc_included} --starjunc {params.starjunc_included} --eij {params.eij_included} \
            --annotation {params.annotation} -c {output.cBout} -m {params.muttypes} \
            --makecB {params.makecB} --makecUP {params.makecUP} --makeArrow {params.makeArrow} \
            --cUPout {output.cUPout} --Arrowout {output.Arrowout} --MaxMem {params.MaxMem} 1> {log} 2>&1
            """


# Make cB (and potentially cU) file
rule makecU:
    input:
        expand(
            "results/counts/{sample}_counts.csv.gz",
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
        chmod +x {params.shellscript}
        {params.shellscript} {threads} \
        {params.min_pos_coverage} {output.mutpos} \
        {output.mutposfilter} {params.max_pos_coverage} 1> {log} 2>&1
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
        ### GOAL: Concatenate but make sure that headers get removed before concatenation.

        head -n 1 {input.cBins[0]} > temp_header.txt
        
        # Prepare an empty, gzipped file for the output
        : > {output.cB}
        
        # Compress the header and add to the output file
        pigz -c temp_header.txt >> {output.cB}
        
        # Iterate over all files, decompress, skip headers, and append to the output file
        for file in {input.cBins}; do
            tail -n +2 ${{file}} | pigz -c >> {output.cB}
        done
        
        # Cleanup the temporary header file
        rm temp_header.txt
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
        ### GOAL: Concatenate but make sure that headers get removed before concatenation.

        head -n 1 {input.cUPins[0]} > temp_cUP_header.txt
        
        # Prepare an empty, gzipped file for the output
        : > {output.cUP}
        
        # Compress the header and add to the output file
        pigz -c temp_cUP_header.txt >> {output.cUP}
        
        # Iterate over all files, decompress, skip headers, and append to the output file
        for file in {input.cUPins}; do
            tail -n +2 ${{file}} | pigz -c >> {output.cUP}
        done
        
        # Cleanup the temporary header file
        rm temp_cUP_header.txt
        """


# Make color-coded tracks
rule maketdf:
    input:
        "results/counts/{sample}_counts.csv.gz",
        "results/sf_reads/{sample}.s.bam",
        "results/normalization/scale",
    output:
        temp("results/tracks/{sample}_success.txt"),
        expand(
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
    log:
        "logs/maketdf/{sample}.log",
    threads: 10
    conda:
        "../envs/full.yaml"
    shell:
        """
        chmod +x {params.shellscript}
        chmod +x {params.pythonscript}
        {params.shellscript} {threads} {wildcards.sample} {input} {params.mut_tracks} {params.genome} {params.WSL} {params.normalize} {params.pythonscript} {output} 1> {log} 2>&1
        """
