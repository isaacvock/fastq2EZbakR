## fastq2EZbakR Output

All output files will be placed in a directory named `results` that will be created the first time you run bam2bakR. The output of greatest interest, the gzipped cB.csv file, will be in `results/cB/`. The cB file and its columns is discussed in great detail on the [intro](index.md) page.

As of version 0.3.0 (released on 12/6/2024), there is a new parameter in the config file called `final_output` that has three sub-parameters (`cB`, `cUP`, and `arrow`).
Each takes a boolean argument (`True` or `False`) and detemrines which file types get created. With this addition, there are now two additional types of output 
you can optionally generate:

* A "counts U-content adjust Poisson" (cUP) file
    - This is identical to the cB file, except the nucleotide content columns have been replaced with an average over reads with identical mutational content. This offers significant compression relative to a cB file while still allowing accurate NR-seq analyses, as originally described in the [bakR paper](https://rnajournal.cshlp.org/content/29/7/958.full). 
    - This can be found in the `results/cUP/` directory if you set `cUP: True` in your config file.
* A partitioned dataset of parquet files
    - This splits the cB file into one file per sample (`results/arrow_dataset/sample=<sampleID>/part-0.parquet`), with each being identical in structure to a full cB file. The structure and file naming scheme is meant to make it convenient to work with using R's [arrow package](https://arrow.apache.org/docs/r/). For example, EZbakR can use this dataset to limit the amount of RAM required to analyze large, multi-sample datasets (see [EZbakR docs](https://isaacvock.github.io/EZbakR/articles/EstimateFractions.html#using-the-apache-arrow-backend) for details).
    - This can be found in the `results/arrow_dataset` directory if you set `arrow: True` in your config file.

You can set any combination of `final_output` options (`cB`, `cUP`, `arrow`) to `True` (at least one of the three needs to be set to `True`). If `mutpos: True`, site-specific mutation tables (`results/cB/mutpos.csv.gz` and `results/cB/mutpos_filtered.csv.gz`) will also be created.


**Universal fastq2EZbakR output**

Processed bam files:

* Sorted and filtered bam/sam files are in `results/sf_reads/`
    - These are passed to the mutation counting and feature assignment scripts.
    - Non-primary alignments and unaligned reads are filtered out. Reads are sorted by read name (and deleted automatically as temporary files unless `--notemp` is passed).

Mutation counting output:

* `<sampleID>_counts.csv.gz` files are in `results/counts/` (temporary unless `--notemp` is passed)
    - Each row of this table corresponds to a single read or read pair
    - The columns represent counts of every mutation type and every type of nucleotide
    - These can be useful for tracking down problems in the mutation counting. See [FAQs](faq.md) for details.
    - Mutation counting is accomplished with a parallel Python script that streams BGZF blocks directly from the filtered BAM file.

Merged feature assignment and mutation counting:

* Per-sample summarized tables (`results/merge_features_and_muts/<sampleID>_cB.csv` and `<sampleID>_cUP.csv`, temporary unless `--notemp` is passed) are streamed via DuckDB directly into `results/cB/cB.csv.gz`, `results/cUP/cUP.csv.gz`, and/or `results/arrow_dataset/sample=<sampleID>/part-0.parquet`.
    - If a read was not assigned to a particular feature type (i.e., exon or gene), then it will have the string "__no_feature" in the relevant column.


Colored sequencing tracks:

* If `make_tracks: True` (default), `.tdf` files that can be used to make the sequencing tracks colored by mutational content (described [here](tracks.md)) are in `results/tracks`
    - Each sample has 12 .tdf files per tracked mutation type, named like `<sampleID>.TC.<#>.<strand>.tdf`, where `<#>` represents a number from 0 to 5 (number of T-to-C mutations in the reads used to make that file) and `<strand>` is either `pos` (plus strand) or `min` (minus strand).

Mutational anomaly / single nucleotide polymorphism (SNP) calls:

* Mutational anomaly / SNP calls are written to `snp.txt` (`REF:ALT:CHR:POS` format) and a diagnostic `snp_counts.tsv` table (recording depth, alternate count, allele frequency, exact binomial p-value, and source) in the `results/snps/` directory.
    - If you did not provide `-s4U` `control_samples` and left `call_snps_without_ctl: False`, `snp.txt` will be empty. If `call_snps_without_ctl: True` with `snp_sample_strategy: "depth"` and a `snp_max_samples` cap, `selected_samples.txt` and `sample_depths.tsv` are also written here.
    - Flagged positions are excluded from both mutation counts and mutable reference nucleotide denominators (e.g., `nT`).

Normalization:

* Scale factors calculated using edgeR's TMM strategy are located in `results/normalization/scale`
    - This is a simple tab-delimited text file with two "columns", one corresponding to the sample ID, and the other corresponding to the scale factor
    - These will be used to scale the heights of the sequencing tracks in `results/tracks/`.

Feature assignment:

* Tables of read counts for each feature assigned via featureCounts (`genes`, `exons`, `transcripts`, `exonic_bins`, `3utr`, `eej`, and `eij`) will be located in a directory called `results/featurecounts_<feature>/`, where `<feature>` refers to the type of feature reads were assigned to.
    - If using feature assignment strategies that don't use featureCounts (`tec` and `junctions`), the intermediate per-read assignment tables are temporary (retained with `--notemp`), and their assignments are included in the final `cB` / `cUP` / `arrow` outputs.

**fastq2EZbakR output if providing fastq files (rather than bams) as input**

Alignment:

* BAM files will be located in the `results/align` directory.
* If using STAR this directory will include a number of additional files, including the transcriptome alignment bam file, an SJ.out.tab file, and various log files (or alignment summary files if using HISAT-3N).

FastQC:

* fastqc output will be stored in the `results/fastqc/` directory

MultiQC:

* MultiQC output will be stored in the `results/multiqc/` directory. It represents a compilation of FastQC output, as well as STAR or HISAT-3N alignment stats when applicable.

RSEM quantification

* If RSEM is used to quantify isoform abundances, its output will be located in the `results/rsem/` directory.

There are several temporary files that are deleted once the steps requiring them finish running.

**Temporary files**

There are several intermediate files produced which, by default, are deleted once the steps of the pipeline that use them as input have finished running. If you would like to retain this output, you can add `--notemp` to your call to `snakemake` to prevent deletion of these files.

