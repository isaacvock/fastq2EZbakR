# fastq2EZbakR

fastq2EZbakR is a Snakemake pipeline designed to process nucleotide recoding RNA-seq data (NR-seq, e.g., [TimeLapse-seq](https://www.nature.com/articles/nmeth.4582), [SLAM-seq](https://www.nature.com/articles/nmeth.4435), [TUC-seq](https://pubmed.ncbi.nlm.nih.gov/31768978/), etc.). fastq2EZbakR provides output readily compatible with [EZbakR](https://github.com/isaacvock/EZbakR), but is also designed to provide processed NR-seq data in a convenient form that you can work with as you see fit. fastq2EZbakR is the successor to [bam2bakR](https://github.com/simonlabcode/bam2bakR). See [our paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1013179) for a discussion of the motivation behind fastq2EZbakR and its companion R package [EZbakR](https://github.com/isaacvock/EZbakR).

**Documentation can be found here: [https://fastq2ezbakr.readthedocs.io/en/latest/](https://fastq2ezbakr.readthedocs.io/en/latest/)**


## Input

1) FASTQ or BAM files
    - If the latter, make sure the not-always-standard MD tag is included.
2) A genome sequence file (FASTA)
    - Used for alignment if FASTQ files are provided and SNP calling (if -label controls are provided) and sequencing track creation in all cases.
3) An annotation file (GTF format)
    - Used for alignment index building if FASTQ files are provided and feature assignment in all cases.

## Output

The main output is a counts Binomial file (a.k.a. a cB file; not a standard file type, it's just what we call the csv file created by fastq2EZbakR). The cB file is a convenient, compressed, and [tidy](https://vita.had.co.nz/papers/tidy-data.pdf) representation of one's NR-seq data. It contains information about:

* The mutation content of aligned reads (e.g., number of T-to-C mutations in a standard NR-seq experiment).
* The mutable nucleotide content of aligned reads (e.g., number of reference T's with which a read overlapped in a standard NR-seq experiment).
* The genomic features to which aligned reads were assigned. See [here](https://fastq2ezbakr.readthedocs.io/en/latest/features/) for details about fastq2EZbakR's uniquely flexible feature assignment strategies.

All other output produced by fastq2EZbakR is documented [here](https://fastq2ezbakr.readthedocs.io/en/latest/output/).

## Update (3/15/2025): Analyses of transcript isoforms

We recently (3/14/2025) put out [a preprint](https://www.biorxiv.org/content/10.1101/2025.03.12.642874v1) describing a method by which to analyze the kinetics of individual transcript isoforms using short read NR-seq data from total RNA. While this strategy is touched on a little bit in one of the EZbakR vignettes ([this one](https://isaacvock.github.io/EZbakR/articles/EstimateFractions.html#isoform-deconvolution)), I have also developed a full fastq-to-volcano plot walkthrough using real downsampled fastq files from that preprint so you can see how every step of the fastq2EZbakR and EZbakR pipeline needs to be configured/run for these analyses. The tutorial is [here](https://isaacvock.github.io/Isoform_Tutorial_Docs/), and the data used in that tutorial is [here](https://github.com/isaacvock/Isoform_Analysis_Tutorial). Over the next couple weeks I will be adding some extra details/analyses to this tutorial, but in its current form (as of 3/15/2025), all of the basics of performing isoform-level analyses are covered there. It also acts as a hand-on tutorial for all of the EZbakR-suite and can thus useful to checkout and try out even if you aren't interested in this particular analysis strategy.

## Why use fastq2EZbakR?

fastq2EZbakR provides a number of unique functionalites not found in other established NR-seq data processing tools. These include:

1. Flexible assignment of reads to [genomic features](https://fastq2ezbakr.readthedocs.io/en/latest/features/).
1. Quantification of any mutation type you are interested in. T-to-C mutation counting is the most common NR-seq application, but any combination of mutation types are fair game. 
1. A [tidy](https://vita.had.co.nz/papers/tidy-data.pdf), easy to work with representation of your mutational data in the form of the aforementioned cB file.
1. Optional site-specific mutation counting (as was used [here](https://acs.figshare.com/collections/Disulfide_Tethering_to_Map_Small_Molecule_Binding_Sites_Transcriptome-wide/7421963) for example). Has allowed fastq2EZbakR to support processing of non-NR-seq mutational probing RNA-seq datasets.
1. Optional automatic downloading and processing of published data available on the [Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra/docs/).

## How do you run fastq2EZbakR?

All of the steps necessary to deploy fastq2EZbakR are discussed in great detail in the official [documentation](https://fastq2ezbakr.readthedocs.io/en/latest/). Here, I will present a super succinct description of what needs to be done, with all necessary code included:

``` bash
### 
# PREREQUISITES: INSTALL MAMBA or CONDA AND GIT (only need to do once per system)
###

# CREATE ENVIRONMENT (only need to do once per system)
mamba create -c conda-forge -c bioconda --name deploy_snakemake snakemake snakedeploy

# CREATE AND NAVIGATE TO WORKING DIRECTORY (only need to do once per dataset)
mkdir path/to/working/directory
cd path/to/working/directory

# DEPLOY PIPELINE TO YOUR WORKING DIRECTORY (only need to do once per dataset)
conda activate deploy_snakemake
snakedeploy deploy-workflow https://github.com/isaacvock/fastq2EZbakR.git . --branch main

###
# EDIT CONFIG FILE (need to do once for each new dataset)
###

# RUN PIPELINE

# See [here](https://snakemake.readthedocs.io/en/stable/executing/cli.html) for details on all of the configurable parameters
snakemake --cores all --use-conda --rerun-triggers mtime --keep-going
```

## Alternative pipeline deployment strategy (to deal with 429 errors)

A recent [change to GitHub](https://github.blog/changelog/2025-05-08-updated-rate-limits-for-unauthenticated-requests/) has led to intermittent fastq2EZbakR crashing errors, when running fastq2EZbakR as detailed above. The errors will look something like:

```
Failed to open source file https://raw.githubusercontent.com/isaacvock/fastq2EZbakR/develop/workflow/scripts/bam2bakR/mut_call.sh
HTTPError: 429 Client Error: Too Many Requests for url: https://raw.githubusercontent.com/isaacvock/fastq2EZbakR/develop/workflow/scripts/bam2bakR/mut_call.sh
```

While simply rerunning the pipeline can often resolve these errors, a more consistent workaround is to clone the repo locally rather than use Snakedeploy. In this case, the deployment workflow looks like:

``` bash
### 
# PREREQUISITES: INSTALL MAMBA or CONDA AND GIT (only need to do once per system)
###

# CREATE ENVIRONMENT (only need to do once per system)
mamba create -c conda-forge -c bioconda --name snakemake snakemake

# (NEW STEP) CLONE REPO LOCALLY (only need to do once per system
git clone https://github.com/isaacvock/fastq2EZbakR.git

# CREATE AND NAVIGATE TO WORKING DIRECTORY (only need to do once per dataset)
mkdir path/to/working/directory
cd path/to/working/directory

# (MODIFIED STEP) "DEPLOY" PIPELINE TO YOUR WORKING DIRECTORY (only need to do once per dataset)
cp -r path/to/fastq2EZbakR/workflow/ ./
rp -r path/to/fastq2EZbakR/config/ ./

###
# EDIT CONFIG FILE (need to do once for each new dataset)
###

# RUN PIPELINE

# See [here](https://snakemake.readthedocs.io/en/stable/executing/cli.html) for details on all of the configurable parameters
snakemake --cores all --use-conda --rerun-triggers mtime --keep-going
```

The difference is that now instead of using Snakedeploy (which you no longer even need to install), you clone the repo locally and copy over the workflow and config directories to wherever you want to run the pipeline. The downside of this approach is that if I make changes to the pipeline that you want to adopt, you will have to manually "update" the pipeline by going to the fastq2EZbakR directory and running `git pull`, and then retransferring the workflow directory (and potentially the config directory too if the pipeline change involves new config parameters you would like to adjust).


