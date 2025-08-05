# fastq2EZbakR Tips and Tricks

This page is devoted to a number of small pieces of advice and details about running fastq2EZbakR that you may find useful. In particular, information relevant ot the following topics is provided:

1. [fastq2EZbakR parameter choices](#parameters)
1. [fastq2EZbakR runtime](#runtime)
1. [fastq2EZbakR RAM usage](#ram)
1. [Miscellaneous tips](#misc)

## Config parameter suggestions<a name="parameters"></a>

There are a lot of parameters you can adjust in the fastq2EbakR config file. To help you out, I have compiled a list of general scenarios and their optimal parameter choices.

### Alignment

For aligning NR-seq data with fastq2EZbakR (or more generally if you are doing the alignment yourself and passing the bam files to fastq2EZbakR), I would suggest the following:

1. Don't allow multi-mappers. As discussed [here](https://github.com/simonlabcode/bam2bakR/issues/26), fastq2EZbakR does nothing fancy with multi-mapping reads, and currently just keeps the primary alignment. As aligners like STAR choose the primary alignment randomly, this amounts to randomly choosing to assign a read to one of its multiple possible genomic origins. 
    - If using STAR (the suggested aligner at this point), this is achieved by adding `--outFilterMultimapNmax 1` to your `star_align_params` parameter.
2. Don't allow soft-clipping. As discussed in the section on [feature assignment](features.md), soft-clipped reads will get called "intronic" by fastq2EZbakR's exon assignment strategy, due to limitations of featureCounts. In general, it is best to trim the ends of your reads as necessary upstream of alignment, so as to eliminate the need for soft-clipping.
    - If using STAR, this is achieved by adding `--alignEndsType EndToEnd` to your `star_align_params` parameter.
3. Be lenient on the mismatch penalization. The defining feature of NR-seq reads is that they should have a number of T-to-C mismatches (and/or potentially G-to-A if using s^6^G) due to metabolic label incorporation and chemcial recoding of the label.
    - If using STAR, the following parameters have proven effective in benchmarking that I have performed using simulated and real data: `--outFilterMismatchNmax 20 --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4`. Again, these can be added to your `star_align_params parameter.

### Trimming

As it is suggested to disable soft-clipping when aligning reads for use in fastq2EZbakR, sufficient upstream trimming of reads is of paramount importance. Some suggestions are:

1. Be fairly stringent. fastq2EZbakR implements trimming with fastp, and I like enabling polyX trimming, quality trimming of ends of reads, and even some hard clipping (i.e., removal of a fixed number of bases after all other trimming is performed)
    - For most total RNA datasets, adding something like `--trim_poly_g --trim_poly_x --cut_tail --cut_front --trim_tail1 3 --trim_front1 3` to `fastp_parameters` for single-end data, and additionally `--trim_tail2 3 --trim_front2 3` for paried-end data, is suggested.
    - For 3'-end sequencing data (e.g., using the Lexogen SLAM-seq kit), it is suggested to be even more stringent, hard clipping the first 12 nucleotides with `trim_front1 12`.
    - Check your FastQC reports (which in fastq2EZbakR are generated post-trimming), to see if additional trimming may be necessary.
2. Specify your adapters. If you leave the `fastp_adapters` argument as `""` in your fastq2EZbakR config, fastp will attempt to auto-detect adapters. While this is reported to be decently robust for paired-end data, it is much sketchier for single-end data. In either case though, you are safer just specifying your adapters explicitly. 


### Feature assignment

A unique aspect of fastq2EZbakR is the set of features that reads can be assigned to (as discussed [here](features.md)). What feature assignment strategies should you turn on? TLDR: definitely "exons", probably "genes", and only any of the others if you are specifically interested in them. In more detail:

1. The first rule is to that it is best to keep things as simple as possible, i.e., turn on as few feature assignments strategy as necessary. The more features you include your cB, the less compressed it will be. Each row of a cB corresponds to a set of reads with identical information (sample of origin, mutation content, nucleotide content, feature assignments, etc.). The more features you assign reads to, the less reads will be in each group. Reads that get assigned to the same gene may not get assigned to the same set of isoforms, and genes assigned to the same set of isoforms may not be assigned to the same set of exon-exon junctions. Thus, more feature assignments will lead to larger cB files, which can make them harder to work with.
2. A basic, vanilla NR-seq analysis will likely only require the "exons" feature assignment (assigns reads to exonic regions of genes). If you have full-length NR-seq data (rather than 3'-end data), I like to also include the "genes" assignment as well. Looking at trends in intronic data (assigned to a gene but not an exonic region of that gene) can be useful, as for standard label times these should be almost completely labeled. For example, intronic reads can provide an alternative means by which to estimate the mutation rate in reads from new RNA (pnew in EZbakR) and to assess dropout. Thus, "exons" and "genes" are the two strategies you should almost definitely always turn on.
3. All of the other strategies are particularly useful for any sort of isoform-specific analysis. "tec" is specifially designed for a full-blown isoform-level analysis of full-length NR-seq data, but this type of analysis is challenging for a number of reasons (see discussion [here](https://www.biorxiv.org/content/10.1101/2025.03.12.642874v1) for example, namely the problem of getting an accurate annotation). The other strategies ("exonic_bins" and "junctions") are nice alternatives that don't provide full isoform-specificity but are robust to many of the problems that can plague an isoform-level analysis. 

### Output files

The main output of interest from fastq2EZbakR is a cB file. Recently though, I added the ability to also (or alternatively) output two slightly different versions of the same basic file. These can be selected via the "final_output" parameter in the config file. They are:

1. "arrow": This is identical to a cB file, but each sample (i.e., input fastq/bam file) gets a separate file organized into directories named in a specific manner. This output can be used as is by EZbakR via its Apache arrow back end. This option is ideal for datasets with lots of samples. In this case, EZbakR can load your data one sample at a time, and thus significanlty limit the amount of RAM required to analyze the data. This appraoch is thus far more scalable than the standard single cB file. I would go so far as to suggest that you should always include this as output in case the full cB file proves too unwieldy.
2. "cUP": cB stands for "counts Binomial". cUP stands for "counts U-content adjusted Poisson". In the original bakR paper, we introduced an alternative to binomial mixture modeling, termed U-content adjusted Poisson mixture modeling. In this strategy, instead of grouping reads by their mutation content and nucleotide content (e.g., TC and nT), reads are grouped by their mutation content, and the average nucleotide content of these reads is computed and passed to the model. This allows for significant compression of a standard cB file, as reads typically have only anywhere from 0 to 10 mutations of a given type but can have a much wider range of nucleotide contents. This can make a cB file less unwieldly to work with while also keeping the information in a single, easy to navigate file (vs. splitting up the data over many files as in the "arrow" output). Currently, EZbakR even by default summarises the cB down to a cUP to make certain analyses more efficient. Currently, there is technically no way to provide a cUP file as input to EZbakR, but that should change in the near future (as of 7/17/2025).

In summary, you should seriously consider including and using the "arrow" output. The "cUP" is a nice idea but you lose a bit of information that the "arrow" output includes (the precise distribution of nucleotide contents).


## Runtime<a name="runtime"></a>

Rough estimates for runtimes for a single fastq/bam file from human data, with 20 cores provided to Snakemake, are below:

1. 25 million read fastq file: 1-2 hours + time to index genome if necessary (usually ~0.5-1 hours)
1. 300 million read fastq file: 5-7 hours + time to index genome if necessary
1. 25 million read bam file: ~30 minutes
1. 300 million read bam file: 3-4 hours

In most cases, you will likely provide more than a single fastq/bam file as input. If you follow the steps in [Deployment](deploy.md) exactly, then that means most rules will run serially (i.e., not in parallel). Most steps of the pipeline run on a single sample worth of data at a time, so runtimes in this case will increase approximately linearly with the number of input files. **This is often not the best way to run fastq2EZbakR, or any Snakemake pipeline, though!**

More specifically, if you are running fastq2EZbakR in an environment where it is possible to request separate jobs to run in parallel (e.g., an institution's shared HPC clusters that utilize a scheduler like slurm or qsub), then you can significantly improve fastq2EZbakR runtimes by making full use of the computational resources at your disposal. The [Slurm/Yale HPC](slurm.md) section has some details relevant to this, but those instructions are mostly specific to the Yale HPC. I will try and write some more general documentation for this eventually, but the general idea is that you can provide Snakemake with information about your system's job scheduler. When you do this, Snakemake will be able to effectively run multiple input files in parallel through fastq2EZbakR. This has the potential to nearly flatten the number of input files vs. runtime function, dependent on how many jobs your particular system allows you to have running at the same time. Some documentation to check out relevant to this is:

1. If you are using version 8.0 of Snakemake or later, then there are a number of plugins to support optimized execution on a wide array of systems. Checkout the sidebar of this page for the list of plug-ins : https://snakemake.github.io/snakemake-plugin-catalog/
1. If you are using a version of Snakemake older than 8.0, then the same thing is accomplishable via built-in functionality documented here: https://snakemake.readthedocs.io/en/v7.32.3/executing/cluster.html#cluster-generic
1. If you are using version 7.29 of Snakemake or later, I highly suggest checking out the solution described on [John Blischak's relevant repo](https://github.com/jdblischak/smk-simple-slurm), which uses [profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) to implement a simple and elegant solution. Technically, this repo is specific to systems using the slurm scheduler, but the general architecture of the solution can be applied to other systems as well. **NOTE**: A lot changed in version 8.0 of Snakemake, meaning that this repo had to significantly alter the solution to accomodate newer versions of Snakemake. The pre- and post-8.0 solution is available as separate branches on this repo though (e.g., version 7 solution is [here](https://github.com/jdblischak/smk-simple-slurm/tree/v7)).


## RAM usage<a name="ram"></a>

The RAM requirements of fastq2EZbakR are mainly a function of two things:

1. Genome size 
1. Sequencing depth (i.e., number of total reads in your fastq/bam files). **NOTE**: This is largely not the case as of version 0.7.0 of fastq2EZbakR (released on 8/4/2025). See discussion below, and [here](https://github.com/isaacvock/fastq2EZbakR/releases/tag/v0.7.0).

RAM usage is mostly a function of the former. The main RAM bottlenecks are any rules that call STAR, where the size of data structures used by STAR to efficiently accomplish tasks like alignment are a funtion of genome size but **not** fastq file size. These are (`rule name`: description):

1. `align`: Alignment of reads to genome
1. `index`: Indexing a genome for the aligner
1. `maketdf`: Making colored sequencing tracks, which under the hood calls STAR for creating tracks (not alignment; turns out to also be fairly RAM-intensive).

For human data, RAM usage of these steps is typically in the following ballparks:

1. `align`: ~50 GB of RAM
1. `index`: ~120 GB of RAM. **NOTE**: these only need to be built once for a given genome + annotation combination and can be provided to fastq2EZbakR if you already have them.
1. `maketdf`: ~50 GB of RAM for T-to-C mutation content colored tracks

Before version 0.7.0 of fastq2EZbakR there was one non-STAR step that could also be fairly RAM-intensive: `merge_features_and_muts`. This rule joins the counts of mutations created by the `cnt_muts` rule with all of the feature assignment tables generated by their relevant rules. Previously, fastq2EZbakR would load the feature assignment and mutation count tables (roughly bam-sized files) fully into memory and join them with data.table. Starting in version 0.7.0, fastq2EZbakR now uses DuckDB to do this in a far RAM-friendlier manner. You can set the RAM-cap provided to DuckDB by setting the `MaxMem` parameter added to the config (default = 8 GB). Even setting this as low as it is by default leads to an improvement in speed relative to the older data.table-based strategy, thus it is proably fine not to tweak this parameter. Due to the changes in v0.7.0, the `lowRAM` config parameter and alternative merging strategy is now obsolote and now recommended.

On some occasions, the mutation counting step (`cnt_muts`) can also use a considerable amount of RAM. This is only the case if a large number of SNPs are called in the `call_snps` step of the pipeline. SNPs are only called if you provide a non-empty list to the `controls` parameter of the config.yaml, and by extension is only relevant if you have -label data for calling SNPs. Mutation counting is run in parallel via a somewhat naive strategy of splitting up the bam file and using GNU Parallel to run a custom Python script on the separate bam file chunks in parallel. This means that the entire SNP call text file is loaded into RAM in each GNU Parallel job. Thus, RAM usage in this step is roughly equal to the size of the snp.txt file created by fastq2EZbakR * the number of cores provided to this step. 

## Miscellaneous Tips<a name="misc"></a>

Below are a collection of miscellaneous Snakemake/fastq2EZbakR usage tips:

1. If you are using a scheduler to batch jobs in parallel, you will need to set the default amount of RAM requested for each batched jobs. Most steps of fastq2EZbakR require very little RAM, while some require a lot. See discussion of RAM-requirements above. You could thus request the amount of RAM required for the RAM-intensive jobs, but this would mean requesting way more RAM than is necessary for other jobs. A better strategy is to request the amount of RAM necessary for most jobs (which should be around 5-10 GB in most cases), and then specifically request more RAM for RAM-intensive jobs. This can be done by adding `--set-resources <rule name>:<resource name>=<amount> <rule name 2>:<resource name 2>=<amount 2> ...`. `<rule name>` refers to the name of the specific Snakemake "rule" that you are requesting a custom resource for. These are listed when you run Snakemake (or a Snakemake dry run with `snakemake -n`). `<resource name>` is the name of the resource that you want to set a custom amount for. These names are arbitrary and depend on what you call each particular resource in your custom profile or call to Snakemake. This is discussed more [here](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#resources). `<amount>` refers to the custom value you would like to set that resource to for that rule. For example, if you use a Snakemake profile like the one discussed in the Slurm/Yale deployment section of this website (repo [here](https://github.com/isaacvock/yale_profile)), then the amount of RAM requested is given the name "mem_mb". Thus, if I want to request 180 GB of RAM for the alignment index creation step, I would add: `--set-resources index:mem_mb=180000` to my call to Snakemake. You can also request custom job runtimes using this same strategy, which can also be useful in some cases.