## Configuring fastq2EZbakR

In the `config/` directory you will find a file named `config.yaml`. If you open it in a text editor, you will see several parameters which you can alter. These parameters are split into three major sections, the first being those that are very important to check and alter as necessary, the second being parameters whose default values are worth assessing but that are designed to be reasonable in most settings. and the third that allow you to specify optional parameters in tools used by fastq2EZbakR (once again, with hopefully reasonable defaults set).

See [here](tips.md) for a summary of non-default parameters in the config that are worth tweaking for your particular use case.

### Parameters you need to set

At the top of the config file, you can specify whether or not you want to download FASTQ files from SRA:

``` yaml
download_fastqs: False
...
sra_accessions: []
```

If so, set `download_fastqs` to True, and specify the SRA accessions in `sra_accessions`. If `download_fastqs` is `False`, `sra_accessions` will be ignored. Downloading FASTQ files will take precedence over providing your own FASTQ files, meaning that the next parameter discussed (`samples`; path to files you are providing) will be ignored if `download_fastqs` is `True`.

If you are providing your own files, the location of these will be specified in the next required parameter:

``` yaml
samples:
    WT_1: data/fastq/WT_1
    WT_2: data/fastq/WT_2
    WT_ctl: data/fastq/WT_ctl
    KO_1: data/fastq/KO_1
    KO_2: data/fastq/KO_2
    KO_ctl: data/fastq/KO_ctl
```
`samples` is the set of \[sample ID\]:\[path\] pairs, where the paths are either to directories containing FASTQ files that you want to process, or to the individual BAM files. **NOTE: if providing FASTQ files, each directory must contain a single FASTQ file or pair of FASTQ files**. Delete the existing sample names and paths and add yours. The sample names in this example are `WT_1`, `WT_2`, `WT_ctl`, `KO_1`, `KO_2`, and `KO_ctl`. These are the sample names that will show up in the `sample` column of the output cB.csv file. The `:` is necessary to distinguish the sample name from what follows, the path to the relevant FASTQ-containing directory. The path can either be relative to the directory that you deployed to (i.e., `workdir` in this example), or an absolute path. In this example, the FASTQ files are located in a directory called `samples` that is inside of a directory called `data` located in `workdir`. 


The next parameter you have to set denotes the sample names of any -s<sup>4</sup>U control samples (i.e., samples that were not fed s<sup>4</sup>U or a similar metabolic label). This is only used to determine which samples should be analyzed for calling SNPs:

``` yaml
control_samples: ["WT_ctl", "KO_ctl"]
```

In this case, the samples named WT_ctl and KO_ctl are the -s<sup>4</sup>U control samples. If using downloaded FASTQs, this should be the relevant SRA accessions of the control samples.

Next is a boolean indicating whether or not your data is paired-end:

```yaml
PE: True
```

(This will likely get removed in later versions and inferred automatically from the number of FASTQs in one of your directories).

Below that is the path to the genome FASTA file that reads will be aligned to:

``` yaml
genome: data/genome/genome.fasta
```

This is followed by the path to the annotation GTF file:

``` yaml
annotation: data/annotation/genome.gtf
```

Make sure that the chromosome names denoted in the FASTA file are identical to what they are called in the GTF file. In addition, fastq2EZbakR assumes that your annotation has the following fields "gene_id" and "type", with "type" including entries "transcript" and "exon". This is pretty standard but is noted here for completeness.

You can then specify the aligner you would like to use:

```yaml
aligner: "star"
```
Currently, only STAR and HISAT2 are implemented, and I highly recommend using STAR. I used to advocate for HISAT-3N when aligning NR-seq data, but two independent papers ([here](https://pubmed.ncbi.nlm.nih.gov/38381903/) and [here](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03313-8)) showed that STAR is about as good (and in some cases even better) at aligning reads from an NR-seq experiment. In addition, HISAT-3N cannot be easily installed with conda, and is not as actively maintained as STAR. Finally, much of the cool functionality of fastq2EZbakR (e.g., assignment of reads to transcript equivalence classes and exon-exon splice junctions) are currently only possible with the output of STAR.


This is followed by the path to the alignment indices:

```yaml
indices: data/indices/star_index
```

You can either provide these yourself, or have fastq2EZbakR create them automatically. In the latter case, just be aware that index creation is a RAM and time intensive process. Indexing the human genome with STAR (and using the provided annotation to include splice junctions and exons in the indices) takes between 1 and 2 hours on a 12 core machine and requires around 100 GB of RAM.

Next you will specify the strandedness of your sequencing library (options: "reverse", "yes", and "no"):

```yaml
strandedness: "reverse"
```

The terminology used here is borrowed from HTSeq (though fastq2EZbakR now uses featureCounts in place of HTSeq). "reverse" means that the first read in a pair (or the only read if your library is single-end) represents the reverse complement of the sequenced RNA. "yes" means that the first read represents the original sequence of the sequenced RNA. "no" means that your library is unstranded, though it is not recommended that you use an unstranded library for NR-seq data (if you have to though, make sure to count both T-to-C and A-to-G mutations if doing standard s<sup>4</sup>U labeling).

Finally, what makes fastq2EZbakR special is all of the ways in which you can assign reads to features (discussed in more detail [here](features.md)). In the `features` section is where you specify which of these you want to use:

```yaml
features:
    genes: True
    exons: True
    tec: False
    exonic_bins: False
    junctions: False
    threeputr: False
    eej: False
    eij: False
```

These are:

* `genes`: Assignment of reads to the gene(s) to which they align. A read will be assigned to a gene if it overlaps with any part of the gene.
* `exons`: Assignment of reads to gene(s) to which they align. A read will ony be assigned to a gene if it strictly overlaps annotated exonic regions of that gene. **NOTE**: featureCounts has a slight suboptimality that makes it impossible to perform this assignment with 100% accuracy. This is because soft-clipped bases are counted as not overlapping any feature, so an arbitrary non-zero cutoff for the number of non-overlapping bases has to be set to avoid failing to assign all soft-clipped reads. The `Transcripts` strategy described above is a more rigorous way of determining if a read only aligned to annotated exonic regions.
* `tec`: Only compatible with using fastq2EZbakR for alignment, and using STAR as the aligner. Assigns reads to transcript equivalence classes (TECs), which is the set of transcript isoforms with which a read is compatible.
* `exonic_bins`: Assignment of reads to exonic bins, as defined in the original [DEXSeq paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3460195/). Best to pair with `Transcripts` or `exons` to filter out reads overlapping intronic regions, as this assignment strategy works like `genes` an will assign reads regardless of whether they overlap any non-exonic regions.
* `junctions`: Only compatible with STAR alignment. Uses the custom jI and jM tags to identify the set of exon-exon junctions a read overlaps. Thus, if you are providing your own BAM files, these tags need to have been included to use this feature assignment strategy.
* `threeputr`: Only compatible with 3'-end sequencing data. Uses either 1) annotated features of type "3UTR" in the provided annotation to assign reads to 3'-ends or 2) builds an annotation of 3'-ends from your data. See [here](features.md) for more details.
* `eej`: A hack that attempts to mimic `junctions` but in a way that does not require STAR alignment or custom tags. A custom annotation is created automatically that includes annotation of exon-exon junction reads, that if a read aligns to, and if the `sj` column always included in the output cB is TRUE, indicates that the read likely overlapped the respective exon-exon junction.
* `eei`: Similar to `eej`, hacky strategy to use featureCounts to assign reads to exon-intron junctions. Use this and `eej` with caution.


### Parameters you should probably double check

The next section of parameters have default settings that will work in a lot of cases, but that you should still double check to make sure they fit your particular use case:

* `mut_tracks`: Specifies the types of mutations you would like to track in the final cB file. Should be a comma separated string of strings of the form \[reference nucleotide\]\[mutated nucleotide\]. For example "TC" denotes that T-to-C mutations should be tracked in the cB file, and "TC,GA" denotes that T-to-C and G-to-A mutations should both be included in the cB file.
* `normalize`: Boolean; if True, then scale factors are calculated using edgeR that are applied to the .tdf sequencing tracks created by fastq2EZbakR.
* `spikename`: If `normalize` is True, and you have spike-ins, then you can specify a string common to all of the gene_ids in your provided annotation GTF. A custom R script will grep for this string when deciding what features to use for normalization purposes.
* `skip_trimming`: If True, then fastp trimming of FASTQs will be skipped.
* `fastp_adapters`: Arguments to specify adapter sequences for trimming by [fastp](https://github.com/OpenGene/fastp). fastp can automatically detect adapters in paired-end libraries, but its always best to specify these explicitly if you know them.
* `flat_annotation`: If `exonic_bins: True`, then this will be the path to and name of the DEXSeq flattened annotation created automatically by fastq2EZbakR.
* `minqual`: Minimum base quality for a mismatch to be called as a bona fide mutation.
* `WSL`: I have found that running fastq2EZbakR on the Windows subsystem for Linux encounters a weird bug where GNU parallel doesn't work for a single step of the pipeline (colored track creation). Thus, if on the WSL, this step needs to be run iteratively, which is what setting `WSL: True` will do.
* `lowRAM`: If True, this can significantly cut down on RAM usage of the mutation count and feature assignment merging step. This is done through sorting the tables to be merged and merging them with a custom Python script that iterates through the rows of all tables to be merged, never loading any of them fully into RAM. The downside of this strategy is additional runtime and temporary disk space usage. With `lowRAM` set to False, the RAM used by the merging step is a function of the sequencing depth (i.e., number of reads in your individual BAM files). This can cause problems for particularly deeply sequenced libraries. For example, a 300 million read library may require > 100 GB for this step.
* `mutpos`: If True, then an additional output will be created, called a cU file, that tracks the mutational content of all mutation types specified in `mut_tracks`. **NOTE**: this option requires a lot of disk space for large datasets.
* `run_rsem`: If True, and if FASTQs are provided and `aligner` is `"star"`, then RSEM will be used to estimate transcript isoform abundances.
* `final_output`: Specifies which output file should be created. Options are:
    - cB: Standard cB file create by default.
    - cUP: cB file but with the average nucleotide content over a set of reads with identical mutational content rather than tracking both the mutational and nucleotide content. See [output.md] for details. NOTE: Not currently compatible with `lowRAM: True`.
    - arrow: Dataset of cB-like files, one sample per file, which pairs well with the [arrow R frontend](https://arrow.apache.org/docs/r/) and [EZbakR](https://isaacvock.github.io/EZbakR/articles/EstimateFractions.html#using-the-apache-arrow-backend).

### Remaining parameters

The remaining parameters tune the behaviors of individual rules and are mostly optional parameters to the command line tools or custom scripts used. Read the comments associated with these parameters, and post an [Issue](https://github.com/isaacvock/fastq2EZbakR/issues) to the fastq2EZbakR Github if you have any questions about these parameters.