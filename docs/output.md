## fastq2EZbakR Output

All output files will be placed in a directory named `results` that will be created the first time you run bam2bakR. The output of greatest interest, the gzipped cB.csv file, will be in `results/cB/`. The cB file and its columns is discussed in great detail on the [intro](index.md) page.

**Universal fastq2EZbakR output**

Processed bam files:

* Sorted and filtered bam/sam files are in `results/sf_reads/`
  - These are passed to the mutation counting and feature assignment scripts.
  - Non-primary alignments and unaligned reads are filtered out. Reads are sorted by read name.

Mutation counting output:

* `<sampleID>_counts.csv.gz` files are in `results/counts/`
  - Each row of this table corresponds to a single read or read pair
  - The columns represent counts of every mutation type and every type of nucleotide
  - These can be useful for tracking down problems in the mutation counting. See [FAQs](faqs.md) for details.
  - Mutation counting is accomplished with a custom python script called by a custom shell script.

Merged feature assignment and mutation counting:

* Tables that have combined the exonic and gene feature assignment information with the mutation calling output are in `results/merge_feature_and_muts/<sampleID>_counts.csv.gz`. 
  - If a read was not assigned to a particular feature type (i.e., exon or gene), then it will have the string "__no_feature" in the relevant column.
* If `lowRAM` in the config.yaml file is set to `True`, then this merged output will be located in the `results/lowram_merge_features_and_counts` directory.


Colored sequencing tracks:

* .tdf files that can be used to make the sequencing tracks colored by mutational content (described [here](tracks.md)) are in `results/tracks`
  - Each sample has 12 .tdf files, named like `<sampleID>.TC.<#>.<strand>.tdf`, where `<#>` represents a number from 0 to 5 (number of T-to-C mutations in the reads used to make that file) and `<strand>` is either `pos` (plus strand) or `min` (minus strand).
  - Currently, if your library is reverse stranded (i.e., first read in a pair represents reverse complement of original RNA sequence), then the plus and minus strand tracks will be flipped. This does not change interpretation of the tracks, you just have to be aware of that when using an annotation to visually decide what reads are the 
  product of sense and antisense transcription.

Single nucleotide polymorphism (SNP) calls:

* SNP calls are in two formats (.txt and VCF)in the `results/snps/` directory.
  - If you did not have any -s<sup>4</sup>U control samples, then the .vcf file will not exist and the .txt file will be empty
  - These SNP calls are used to identify nucleotides which should be ignored for T-to-C mutation counting

Normalization:

* Scale factors calculated using edgeR's TMM strategy are located in `results/normalization/scale`
  - This is a simple tab-delimited text file with two "columns", one corresponding to the sample ID, and the other corresponding to the scale factor
  - These will be used to scale the heights of the sequencing tracks in `results/tracks/`.

Feature assignment:

* Tables of read counts for each feature assigned via featureCounts (genes, exons, exonic_bins, eej, and eij) will be located in a directory called `results/featurecounts_<feature>/`, where `<feature>` refers to the type of feature reads were assigned to.
  - If using feature assignment strategies that don't use featureCounts, all output produced by these are temporary. See below for details regaring what this means and how to keep these files if necessary. In either case, the read assignment information is also present in the merged tables discussed above.

**fastq2EZbakR output if providing fastq files (rather than bams) as input**

Alignment:

* BAM files will be located in the `results/align` directory. I
* If using STAR this directory will include a number of additional files, including the transcriptome alignment bam file, an SJ.out.tab file, and various log files

FastQC:

* fastqc output will be stored in the `results/fastqc/` directory

RSEM quantification

* If RSEM is used to quantify isoform abundances, its output will be located in the `results/rsem/` directory.

There are several temporary files that are deleted once the steps requiring them finish running.

**Temporary files**

There are several intermediate files produced which, by default, are deleted once the steps of the pipeline that use them as input have finished running. If you would like to retain this output, you can add `--notemp` to your call to `snakemake` to prevent deletion of these files.

