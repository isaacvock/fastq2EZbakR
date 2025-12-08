## Release Notes

### Version 0.8.0

A new feature assignment strategy has been added: assignment to 3'-ends (`threeputr`). This strategy is only appropriate if one has 3'-end sequencing data (e.g., if you used the popular Quant-seq SLAM-seq kit).

This strategy can be run in two different ways:

1. User provides an annotation file (GTF) that includes features of type "3UTR", and assignment is done to these annotated 3'UTRs. The GTF must also contain a column called "utr_id" assigning an alphanumeric ID to each 3'-UTR.
2. A 3'-UTR annotation is generated de novo using the user's 3'-end sequencing data.

For de novo calling, 3'-ends are called with single nucleotide resolution from aligned reads, and these are clustered into called 3'-ends. Parameters in the config file control post-hoc filtering of these to remove 1) unannotated 3'-ends near polyA stretches that are likely instances of mispriming and 2) unannotated 3'-ends not near a cleavage consensus motif (AAUAAA). Alternatively, users can choose to filter out all called 3'-end clusters that don't overlap an annotated final exon. 

### Version 0.7.6

Small bug fix where an added space in the fastq file path name when downloading from the SRA would break SLURM profiles

### Version 0.7.5

Some quick bug fixes:

- Fixed list index errors plaguing the no trimming route
- Allow fastq files to have .fq or .fastq extensions (.fastq was previously only extension allowed). Files can of course be gzipped


### Version 0.7.4

Some small bug fixes:

1. If A-to-T mutations are tracked, DuckDB would throw an error due to the AT column of the counts file being interpreted as a SQL term.
2. Sometimes, DuckDB's auto-type-detection would fail for the mutation counts csv file, with the FR column (a string) being misparsed as a boolean. This is now fixed.


### Version 0.7.3

Some small bug fixes/improvements:

1. Using template-coordinate sorting to reduce disk space usage when setting `mutpos: True`. 
2. Fix bug in merge_features_and_muts when no features are set to True.
3. Small version change in environment to make bam2bakR fully ARM64 compatible (alignment with STAR will not be though in order to address Issue #38).

### Version 0.7.2

Now fully Mac-OS compatible (had to remove one linux-specific dependency that i added a while back for the mutpos branch). Might also be ARM compatible now, but if not users can install Rosetta and run within a Rosetta terminal (see #38).

[EDIT]: Should by ARM compatible unless user is running RSEM. RSEM has no ARM-compatible conda releases.

### Version 0.7.1

Cleaned up the mega-environment that the bam filtering, mutation counting, and final file creation rules run in. This means getting rid of lots of dependencies that no longer get used and pinning/updating versions of ones that do. Also, the normalization scale factor calculation also used this same environment, but this has now been spun off into its own environment as it is all R packages.

### Version 0.7.0

Changes:

1. A significant improvement to the computational performance (mainly RAM-usage, but also should be a bit faster). The merging of features and mutational data (rule merge_features_and_muts) was previously done through sequential loading of mutations/feature tables and joining them with data.table. Thus, bam file-sized tables were fully loaded into RAM, meaning that RAM-usage in this step scaled with sequencing depth, and could become untenable (> 100 GB) for deeply sequenced libraries. Switched to using DuckDB with a hard memory cutoff (set by MaxMem parameter added to config, 8GB by default) which not only solves the RAM problem, but also allows for some optimizations that means an overall speed-up of this step.
2. Had to tweak the MultiQC rule a bit to get it to run without the annoying conda issues. Easiest solution was to regress to an older wrapper version, used by a popular [RNA-seq snakemake pipeline](https://github.com/snakemake-workflows/rna-seq-star-deseq2)

MultiQC issues prevented me from also cleaning up the environments used throughout the pipeline (pinning versions and removing unneeded dependencies), but I hope to release a patch soon with these changes.

### Version 0.6.3

Added the option to use both exon and intron mapping reads for track normalization, via setting the use_exons_only parameter added to the config to False.

### Version 0.6.2

Added STAR alignment stats to MultiQC output

### Version 0.6.1

Added the missing bam2bakr parameter in the config.

### Version 0.6.0

Some bug fixes and new functionality:

1. SE fastp output is now more reasonable; there was a missing "\" before causing file names and directory names to be merged
2. Fixed a bug that could cause problems in cnt_muts if running multiple jobs in parallel, for example via HPC and a profile. Naive grepping of files could cause crossover between temp files from different samples, meaning that files being used by one instance could be deleted by another.
3. Added multiQC output. Currently only aggregating fastQC outputs for easier viewing, but hoping to expand this in the near future.


### Version 0.5.0

Changes:

- If RSEM is used for quantification, then assignment of reads to transcript equivalence class is done via the transcriptome alignment output by RSEM, as a more efficient solution to filtering out isoforms on the wrong strand. 
- Fixed a bug where junction feature assignment was not successfully auto-detected for adding the necessary SAM tags
- Slight update to config defaults

### Version 0.4.1

Fixes major bug that can randomly flip strandedness of paired-end libraries, introduced in v0.4.0.

### Version 0.4.0

1. SNP calling is now more streamlined (better use of bcftools native parallelization rather than hack involving GNU parallel)
2. Users can provide options to `bcftools mpileup` and `bcftools call` via new config parameters
3. Addressed #23 by implementing solution suggested by user (ld32)

### Version 0.3.0

New types of output can now be generated, specifically designed to facilitate analyses of larger datasets.

* In the config file, you will now see a parameter called "final_output", with three sub-parameters, "cB", "cUP", and "arrow". See [output.md] for more details as to how these differ and what combination can be set to `True`.
* A bug was also fixed where additional, optional fastp parameters were not getting passed to fastp.

### Version 0.2.0

Adds a handful of quality of life improvements meant to reduce the complexity of pipeline configuration. Namely:

* Whether or not the user has provided BAM or FASTQ files is auto-detected by looking for .bam at the end of all strings specified under `samples` in the config file.
* Old or overly experimental functionality has been removed
* TEC assignment flag has been moved from a unique parameter sub-heading ("strategies") to a new sub-heading of the "features" parameter

### Version 0.1.0

First official release of fastq2EZbakR. Also represents version used for EZbakR suite preprint analyses.