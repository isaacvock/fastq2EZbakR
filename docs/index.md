# Welcome to fastq2EZbakR!

fastq2EZbakR is a Snakemake implementation of the [TimeLapse pipeline](https://bitbucket.org/mattsimon9/timelapse_pipeline/src/master/) developed by the [Simon lab](https://simonlab.yale.edu/) at Yale. The contributors to the original pipeline are Matthew Simon, Jeremy Schofield, Martin Machyna, Lea Kiefer, and Joshua Zimmer. Despite its origins, a lot has changed since the initial creation of the pipeline, and fastq2EZbakR has a load of novel functionality. It is also an extension/rewrite of bam2bakR, doing pretty much everything it does and then some.

## Where to go

Step 1: Read the [Deployment](deploy.md) documentation to get up and running with fastq2EZbakR, or the [Slurm](slurm.md) documentation for deploying fastq2EZbakR on an HPC system with a slurm scheduler (e.g., Yale's HPC).

Step 2: Read the [Configuration](configuration.md) documentation to get details about all config parameters.

Stpe 3: Read about fastq2EZbakR's [generalized feature assignment](features.md), one of the major benefits of using fastq2EZbakR.

Step 4: Read the [Tips and Tricks](tips.md) documentation for helpful pointers that may make your life easier.

Step 5: Read about [Output](output.md) produced by fastq2EZbakR.

Step 6: Check out ancillary documentation about creating [tracks](tracks.md) and [FAQs](faq.md).

## What fastq2EZbakR does

The input to fastq2EZbakR is either fastq files or aligned bam files (the latter must have the not-always-standard MD tag). The main output of fastq2EZbakR is a so-called cB (counts binomial) file that will always include the following columns:

* sample - Sample name
* rname - Chromosome name
* sj - Logical: TRUE if read contains exon-exon spliced junction
* n - Number of reads which have the identical set of values described above

In addition, columns reporting mutation counts and nucleotide counts will be included. For a standard NR-seq dataset (s4U labeling), that means tracking T-to-C mutation counts (column name: TC) and the number of reference Ts covered by a read (column name: nT). Finally, reads will be assigned to a set of annotated features, and columns will be included based on which of these feature assignment strategies you have activated in your particular pipeline run. The possibilities include:

* GF: gene read was assigned to (any region of gene)
* XF: gene read was assigned to (only exonic regions of gene)
* exonic_bin: exonic bins as defined in DEXSeq paper
* transcripts: set of transcripts a read is compatible with (i.e, its transcript equivalence class)
* junction_start: 5' splice site of exon-exon junction (genomic coordinate)
* junction_end: 3' splice site of exon-exon junction (genomic coordinate)
* ei_junction_id: Numerical ID given to a given exon-intron junction
* ee_junction_id: Numerical ID given to a given exon-exon junction

See [Configuration](configuration.md) and [Feature Assignment](features.md) for details about feature assignment strategies and how to select which to use.