#!/usr/bin/env Rscript
# Load dependencies ------------------------------------------------------------

library(data.table)
library(dplyr)
library(rtracklayer)
library(optparse)
library(GenomicRanges)


### Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

option_list <- list(
    make_option(c("--bed_minus", type = "character"),
                help = "Path to - strand 3'-end coverage bedgraph file"),
    make_option(c("--bed_plus", type = "character"),
                help = "Path to + strand 3'-end coverage bedgraph file"),
    make_option(c("--gtf", type = "character"),
                help = "Path to input GTF file for gene annotation"),
    make_option(c("--output", type = "character"),
                help = "Path to output GTF file created by script"),
    make_option(c("--extension", type = "numeric"),
                help = "Length by which to extend end of genes when assessing peak overlap",
                default = 0),
    make_option(c("--min_coverage", type = "numeric"),
                help = "Minimum total coverage for a cluster to be called a UTR",
                default = 100),
    make_option(c("--min_fxn", type = "numeric"),
                help = "Minimum fraction usage of a UTR cluster to call it a UTR",
                default = 0)
)


opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

### Length by which to extend end of annotated genes when assigning 3'-end peaks
### to genes
extend_len <- round(opt$extension)
coverage_filter <- opt$min_coverage
fxn_filter <- opt$min_fxn

# Convert bed file to annotation -----------------------------------------------

### Load data

# Plus and minus are technically flipped due to RT
bed_minus <- fread(opt$bed_minus,
                  col.names = c("seqname", "cluster", "start", "end", "reads"))

bed_plus <- fread(opt$bed_plus,
                   col.names = c("seqname", "cluster", "start", "end", "reads"))

gtf <- rtracklayer::import(opt$gtf)

### My annotation is a bit wonky; some genes don't have an entry of type "gene"
### and some genes don't have an entry of type "transcript".
gene_gtf_df <- gtf %>%
  as_tibble() %>%
  dplyr::filter(type %in% c("gene", "transcript") &
                  !is.na(gene_id)) %>%
  dplyr::group_by(seqnames, strand, gene_id) %>%
  dplyr::summarise(
    start = min(start) - extend_len,
    end = max(end) + extend_len,
    type = "gene"
  ) %>%
  dplyr::filter(
    !grepl("-", gene_id)
  )


genes <- GenomicRanges::GRanges(
  seqnames = Rle(gene_gtf_df$seqnames),
  ranges = IRanges(
    start = gene_gtf_df$start,
    end = gene_gtf_df$end
  ),
  strand = gene_gtf_df$strand,
  gene_id = gene_gtf_df$gene_id
)


### Find genes that each UTR peak overlaps to annotate peaks

gr_plus <- GenomicRanges::GRanges(
  seqnames = Rle(bed_plus$seqname),
  ranges = IRanges(
    start = bed_plus$start,
    end = bed_plus$end
  ),
  strand = "+",
  coverage = bed_plus$reads
)

gr_minus <- GenomicRanges::GRanges(
  seqnames = Rle(bed_minus$seqname),
  ranges = IRanges(
    start = bed_minus$start,
    end = bed_minus$end
  ),
  strand = "-",
  coverage = bed_minus$reads
)

gr_3pseq <- c(gr_plus, gr_minus)


overlaps <- findOverlaps(
  gr_3pseq,
  genes,
  ignore.strand = FALSE
)

### Want to annotate clusters with genes but allow each cluster to map to
### multiple genes (as there are some overlapping gene edge cases)

df_3pseq <- as_tibble(gr_3pseq) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    queryHits = 1:dplyr::n()
  )


df_genes <- as_tibble(genes) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    subjectHits = 1:dplyr::n()
  ) %>%
  dplyr::select(
    gene_id, subjectHits
  )


clusters_annotated <- df_3pseq %>%
  dplyr::inner_join(
    overlaps %>%
      as_tibble(),
    by = "queryHits"
  ) %>%
  dplyr::inner_join(
    df_genes,
    by = "subjectHits"
  )


#### Build 3'-UTR annotation from filtered peaks

### Exploration suggests that high coverage but low fraction usage peaks
### could be interesting, so for now will stick to just a simple coverage filter

peaks_filtered <- clusters_annotated %>%
  dplyr::group_by(gene_id) %>%
  dplyr::mutate(
    fxn_usage = coverage / sum(coverage)
  ) %>%
  dplyr::filter(
    coverage >= coverage_filter &
    fxn_usage >= fxn_filter
  )


ThreePUTR_gr <- GenomicRanges::GRanges(
  seqnames = Rle(peaks_filtered$seqnames),
  ranges = IRanges(
    start = peaks_filtered$start,
    end = peaks_filtered$end
  ),
  strand = peaks_filtered$strand,
  gene_id = peaks_filtered$gene_id,
  type = "3UTR"
)

rtracklayer::export(
  ThreePUTR_gr,
  opt$output
)

