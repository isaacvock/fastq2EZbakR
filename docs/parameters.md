## Config parameter suggestions

There are a lot of parameters you can adjust in the fastq2EbakR config file. To help you out, I have compiled a list of general scenarios and their optimal parameter choices.

### Alignment

For aligning NR-seq data with fastq2EZbakR (or more generally if you are doing the alignment yourself and passing the bam files to fastq2EZbakR), I would suggest the following:

1. Don't allow multi-mappers. As discussed [here](https://github.com/simonlabcode/bam2bakR/issues/26), fastq2EZbakR does nothing fancy with multi-mapping reads, and currently just keeps the primary alignment. As aligners like STAR choose the primary alignment randomly, this amounts to randomly choosing to assign a read to one of its multiple possible genomic origins. 
    - If using STAR (the suggested aligner at this point), this is achieved by adding `--outFilterMultimapNmax 1` to your `star_align_params` parameter.
2. Don't allow soft-clipping. As discussed in the section on [feature assignment](features.md), soft-clipped reads will get called "intronic" by fastq2EZbakR's exon assignment strategy, due to limitations of featureCounts. In general, it is best to trim the ends of your reads as necessary upstream of alignment, so as to eliminate the need for soft-clipping.
    - If using STAR, this is achieved by adding `--alignEndsType EndToEnd` to your `star_align_params` parameter.
3. Be lenient on the mismatch penalization. The defining feature of NR-seq reads is that they should have a number of T-to-C mismatches (and/or potentially G-to-A if using s6G) due to metabolic label incorporation and chemcial recoding of the label.
    - If using STAR, the following parameters have proven effective in benchmarking that I have performed using simulated and real data: `--outFilterMismatchNmax 20 --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4`. Again, these can be added to your `star_align_params parameter.

### Trimming

As it is suggested to disable soft-clipping when aligning reads for use in fastq2EZbakR, sufficient upstream trimming of reads is of paramount importance. Some suggestions are:

1. Be fairly stringent. fastq2EZbakR implements trimming with fastp, and I like enabling polyX trimming, quality trimming of ends of reads, and even some hard clipping (i.e., removal of a fixed number of bases after all other trimming is performed)
    - For most total RNA datasets, adding something like `--trim_poly_g --trim_poly_x --cut_tail --cut_front --trim_tail1 3 --trim_front1 3` to `fastp_parameters` for single-end data, and additionally `--trim_tail2 3 --trim_front2 3` for paried-end data, is suggested.
    - For 3'-end sequencing data (e.g., using the Lexogen SLAM-seq kit), it is suggested to be even more stringent, hard clipping the first 12 nucleotides with `trim_front1 12`.
    - Check your FastQC reports (which in fastq2EZbakR are generated post-trimming), to see if additional trimming may be necessary.
2. Specify your adapters. If you leave the `fastp_adapters` argument as `""` in your fastq2EZbakR config, fastp will attempt to auto-detect adapters. While this is reported to be decently robust for paired-end data, it is much sketchier for single-end data. In either case though, you are safer just specifying your adapters explicitly. 


### Feature assignment


### Output files

