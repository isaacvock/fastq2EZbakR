#!/usr/bin/env Rscript
# Load dependencies ------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(rtracklayer)
  library(optparse)
  library(GenomicRanges)
  library(Biostrings)
  library(BSgenome)
  library(stringr)
})

### Helper functions -----------------------------------------------------------

### Get sequences
get_seq_windows <- function(genome, chr, start, end) {
  # chr, start, end are vectors of the same length
  if (is(genome, "BSgenome")) {
    # BSgenome method supports start/end
    Biostrings::getSeq(genome, names = chr, start = start, end = end)
  } else if (is(genome, "XStringSet") || is(genome, "DNAStringSet")) {
    # FASTA case: genome is an XStringSet/DNAStringSet
    if (is.null(names(genome))) {
      stop("Genome XStringSet must have names corresponding to chromosome names.")
    }
    
    idx <- match(chr, names(genome))
    if (anyNA(idx)) {
      stop("Some seqnames in peaks not found in genome FASTA: ",
           paste(unique(chr[is.na(idx)]), collapse = ", "))
    }
    
    full_seqs <- genome[idx]
    # subseq() is vectorised over XStringSet + start/end
    Biostrings::subseq(full_seqs, start = start, end = end)
  } else {
    stop("Unsupported 'genome' class: ", paste(class(genome), collapse = ", "))
  }
}



### CPA motifs (AAUAAA / AUUAAA) check near cleavage site
has_cpa_motif <- function(gr, genome, chr_lengths,
                          upstream = 50L,
                          downstream = 10L) {
  n <- length(gr)
  if (n == 0L) return(logical(0))
  
  res <- logical(n)
  seqnames_char <- as.character(GenomicRanges::seqnames(gr))
  strand_char   <- as.character(GenomicRanges::strand(gr))
  
  ## + strand: cleavage at end; search [end - upstream + 1, end + downstream]
  plus_idx <- which(strand_char == "+")
  if (length(plus_idx) > 0L) {
    site <- GenomicRanges::end(gr)[plus_idx]
    chr  <- seqnames_char[plus_idx]
    chr_len <- chr_lengths[chr]
    
    start_win <- pmax(1L, site - upstream + 1L)
    end_win   <- pmin(chr_len, site + downstream)
    
    seqs <- get_seq_windows(genome, chr, start_win, end_win)
    res[plus_idx] <-
      Biostrings::vcountPattern("AATAAA", seqs, fixed = TRUE) > 0L |
      Biostrings::vcountPattern("ATTAAA", seqs, fixed = TRUE) > 0L
  }
  
  ## - strand: cleavage at start; search around that, then reverse complement
  minus_idx <- which(strand_char == "-")
  if (length(minus_idx) > 0L) {
    site <- GenomicRanges::start(gr)[minus_idx]
    chr  <- seqnames_char[minus_idx]
    chr_len <- chr_lengths[chr]
    
    start_win <- pmax(1L, site - downstream)
    end_win   <- pmin(chr_len, site + upstream - 1L)
    
    seqs <- get_seq_windows(genome, chr, start_win, end_win)
    seqs_rc <- Biostrings::reverseComplement(seqs)
    res[minus_idx] <-
      Biostrings::vcountPattern("AATAAA", seqs_rc, fixed = TRUE) > 0L |
      Biostrings::vcountPattern("ATTAAA", seqs_rc, fixed = TRUE) > 0L
  }
  
  res
}


### Internal priming / false polyA filter
is_near_false_polyA <- function(gr, genome, chr_lengths,
                                tract_len = 7L,
                                flank_down = 25L) {
  n <- length(gr)
  if (n == 0L) return(logical(0))
  
  seqnames_char <- as.character(GenomicRanges::seqnames(gr))
  strand_char   <- as.character(GenomicRanges::strand(gr))
  
  # define cleavage site in genomic coords
  cleavage <- ifelse(
    strand_char == "+",
    GenomicRanges::end(gr),
    GenomicRanges::start(gr)
  )
  
  chr     <- seqnames_char
  chr_len <- chr_lengths[chr]
  
  if (anyNA(chr_len)) {
    stop("Some seqnames in peaks not found in genome FASTA: ",
         paste(unique(chr[is.na(chr_len)]), collapse = ", "))
  }

  # Strictly downstream in transcript orientation; the cleavage base is excluded.
  start_win <- pmax(1L, ifelse(strand_char == "+",
                              cleavage + 1L, cleavage - flank_down))
  end_win <- pmin(chr_len, ifelse(strand_char == "+",
                                 cleavage + flank_down, cleavage - 1L))
  valid <- which(start_win <= end_win)
  res <- logical(n)
  if (length(valid) == 0L) return(res)

  seqs <- get_seq_windows(genome, chr[valid], start_win[valid], end_win[valid])
  minus_idx <- which(strand_char[valid] == "-")
  if (length(minus_idx) > 0L) {
    seqs[minus_idx] <- Biostrings::reverseComplement(seqs[minus_idx])
  }

  polyA_pat <- Biostrings::DNAString(paste(rep("A", tract_len), collapse = ""))
  hitsA <- Biostrings::vcountPattern(polyA_pat, seqs, fixed = TRUE) > 0L

  # A-richness also catches interrupted tracts: >=7 As in the first 10 bases.
  # Short chromosome-edge windows do not qualify for this density criterion.
  rich <- logical(length(seqs))
  full_window <- which(width(seqs) >= 10L)
  if (length(full_window) > 0L) {
    rich[full_window] <- Biostrings::vcountPattern(
      "A", Biostrings::subseq(seqs[full_window], start = 1L, end = 10L),
      fixed = TRUE
    ) >= 7L
  }
  res[valid] <- hitsA | rich
  res
}



### Parse command line arguments -----------------------------------------------

option_list <- list(
  make_option(c("--bed_minus"),
              type = "character",
              help = "Path to - strand 3'-end coverage bedgraph file"),
  make_option(c("--bed_plus"),
              type = "character",
              help = "Path to + strand 3'-end coverage bedgraph file"),
  make_option(c("--gtf"),
              type = "character",
              help = "Path to input GTF file for gene annotation"),
  make_option(c("--fasta"),
              type = "character",
              help = "Genome FASTA required for enabled internal-priming flags or CPA checks"),
  make_option(c("--output"),
              type = "character",
              help = "Path to output GTF file created by script"),
  make_option(c("--extension"),
              type = "numeric",
              default = 0,
              help = "Length by which to extend end of genes when assessing peak overlap"),
  make_option(c("--min_coverage"),
              type = "numeric",
              default = 100,
              help = "Minimum total coverage for a cluster to be called a UTR"),
  make_option(c("--min_fxn"),
              type = "numeric",
              default = 0,
              help = "Minimum fraction usage of a UTR cluster to call it a UTR"),
  make_option(c("--false_polyA_len"),
              type = "numeric",
              default = 7,
              help = paste("Flag downstream A runs of this length within 25 bases or >=7 As",
                           "in the first 10 transcript-oriented bases; filter only unannotated",
                           "sites. <=0 disables both criteria.")),
  make_option(c("--require_CPA"),
              type = "logical",
              default = FALSE,
              help = "Whether to require AATAAA or ATTAAA near unannotated 3'-UTRs"),
  make_option(c("--only_annotated"),
              type = "logical",
              default = FALSE,
              help = "If TRUE, keep only 3'-ends overlapping annotated 3'-UTRs (last exons)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

### Params / convenience vars --------------------------------------------------

extend_len      <- round(opt$extension)
coverage_filter <- opt$min_coverage
fxn_filter      <- opt$min_fxn
false_polyA_len <- as.integer(opt$false_polyA_len)
only_annotated  <- isTRUE(opt$only_annotated)
require_CPA     <- isTRUE(opt$require_CPA)

# Convert bed file to annotation -----------------------------------------------

### Load data

# Plus and minus are technically flipped due to RT (as per your note),
# but here we trust that the input bedgraphs are already labeled appropriately.
bed_minus <- fread(
  opt$bed_minus,
  col.names = c("seqname", "cluster", "start", "end", "reads")
)

bed_plus <- fread(
  opt$bed_plus,
  col.names = c("seqname", "cluster", "start", "end", "reads")
)

gtf <- rtracklayer::import(opt$gtf)

### Build “gene body” GRanges (extended) for assigning peaks to genes ----------

gene_gtf_df <- gtf %>%
  as_tibble() %>%
  dplyr::filter(type %in% c("gene", "transcript") &
                  !is.na(gene_id)) %>%
  dplyr::group_by(seqnames, strand, gene_id) %>%
  dplyr::summarise(
    start = min(start) - extend_len,
    end   = max(end) + extend_len,
    type  = "gene",
    .groups = "drop"
  ) %>%
  dplyr::filter(
    !grepl("-", gene_id)
  )

genes <- GenomicRanges::GRanges(
  seqnames = Rle(gene_gtf_df$seqnames),
  ranges = IRanges(
    start = gene_gtf_df$start,
    end   = gene_gtf_df$end
  ),
  strand  = gene_gtf_df$strand,
  gene_id = gene_gtf_df$gene_id
)

### GRanges for 3'-end peaks ---------------------------------------------------

# BED uses zero-based starts; GRanges and GTF use one-based closed intervals.
gr_plus <- GenomicRanges::GRanges(
  seqnames = Rle(bed_plus$seqname),
  ranges = IRanges(
    start = bed_plus$start + 1L,
    end   = bed_plus$end
  ),
  strand   = "+",
  coverage = bed_plus$reads
)

gr_minus <- GenomicRanges::GRanges(
  seqnames = Rle(bed_minus$seqname),
  ranges = IRanges(
    start = bed_minus$start + 1L,
    end   = bed_minus$end
  ),
  strand   = "-",
  coverage = bed_minus$reads
)

gr_3pseq <- c(gr_plus, gr_minus)

### Assign peaks to genes ------------------------------------------------------

overlaps <- findOverlaps(
  gr_3pseq,
  genes,
  ignore.strand = FALSE
)

df_3pseq <- as_tibble(gr_3pseq) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    queryHits = dplyr::row_number()
  )

df_genes <- as_tibble(genes) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    subjectHits = dplyr::row_number()
  ) %>%
  dplyr::select(
    gene_id, subjectHits
  )

clusters_annotated <- df_3pseq %>%
  dplyr::inner_join(
    as_tibble(overlaps),
    by = "queryHits"
  ) %>%
  dplyr::inner_join(
    df_genes,
    by = "subjectHits"
  )

### Filter peaks by coverage and usage ----------------------------------------

peaks_filtered <- clusters_annotated %>%
  dplyr::group_by(gene_id) %>%
  dplyr::mutate(
    fxn_usage = coverage / sum(coverage)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(
    coverage >= coverage_filter &
      fxn_usage >= fxn_filter
  ) %>%
  dplyr::arrange(
    gene_id
  )

peak_genes <- peaks_filtered$gene_id

rle_out <- rle(peak_genes)
runs <- rle_out$lengths

suffixes <- sapply(
  runs,
  function (x) 1:x
) %>%
  unlist()

utr_ids <- paste0(peak_genes, "_utr_", suffixes)

ThreePUTR_gr <- GenomicRanges::GRanges(
  seqnames = Rle(peaks_filtered$seqnames),
  ranges = IRanges(
    start = peaks_filtered$start,
    end   = peaks_filtered$end
  ),
  strand  = peaks_filtered$strand,
  gene_id = peaks_filtered$gene_id,
  utr_id = utr_ids,
  type    = "3UTR"
)

### NEW: mark annotated vs de novo using last exon overlaps --------------------

last_exon_df <- gtf %>%
  dplyr::as_tibble() %>%
  dplyr::filter(
    type == "exon" & !is.na(transcript_id)
  ) %>%
  dplyr::group_by(
    transcript_id
  ) %>%
  dplyr::filter(
    (strand == "+" & end == max(end)) |
      (strand == "-" & start == min(start))
  )

last_exons <- GenomicRanges::GRanges(
  seqnames = last_exon_df$seqnames,
  ranges = IRanges(
    start = last_exon_df$start,
    end = last_exon_df$end
  ),
  strand = last_exon_df$strand,
  transcript_id = last_exon_df$transcript_id,
  gene_id = last_exon_df$gene_id,
  type = "last_exon"
)

annotated <- rep(FALSE, length(ThreePUTR_gr))

hits_le <- findOverlaps(ThreePUTR_gr, last_exons, ignore.strand = FALSE)

q_gene <- ThreePUTR_gr$gene_id[queryHits(hits_le)]
s_gene <- last_exons$gene_id[subjectHits(hits_le)]
same_gene <- !is.na(q_gene) & !is.na(s_gene) & q_gene == s_gene

if (any(same_gene)) {
  annotated[unique(queryHits(hits_le[same_gene]))] <- TRUE
}

ThreePUTR_gr$annotated <- annotated

### Assess sequence evidence, including priming flags in terminal exons --------

unannot_idx <- which(!ThreePUTR_gr$annotated)
check_priming <- false_polyA_len > 0L && length(ThreePUTR_gr) > 0L
check_CPA <- !only_annotated && require_CPA && length(unannot_idx) > 0L

# NA means not assessed, rather than evidence that a site passed a sequence check.
ThreePUTR_gr$has_CPA <- rep(NA, length(ThreePUTR_gr))
ThreePUTR_gr$false_polyA <- rep(NA, length(ThreePUTR_gr))

if (check_priming || check_CPA) {
  if (is.null(opt$fasta) || !file.exists(opt$fasta)) {
    stop("A FASTA file must be provided via --fasta for internal-priming flags ",
         "(including annotated sites) or enabled CPA checks.")
  }

  genome <- Biostrings::readDNAStringSet(opt$fasta)
  names(genome) <- sapply(str_split(names(genome), pattern = " "),
                          function(x) x[1])
  chr_lengths <- setNames(as.integer(width(genome)), names(genome))

  if (check_priming) {
    ThreePUTR_gr$false_polyA <- is_near_false_polyA(
      ThreePUTR_gr, genome, chr_lengths, tract_len = false_polyA_len
    )
  }
  if (check_CPA) {
    ThreePUTR_gr$has_CPA[unannot_idx] <- has_cpa_motif(
      ThreePUTR_gr[unannot_idx], genome, chr_lengths
    )
  }
}

### Filter candidates; annotated sites retain their sequence flags -------------

ThreePUTR_gr_filtered <- ThreePUTR_gr
if (only_annotated) {
  ThreePUTR_gr_filtered <- ThreePUTR_gr[ThreePUTR_gr$annotated]
} else {
  cpa_ok <- !require_CPA | ThreePUTR_gr$has_CPA
  priming_ok <- false_polyA_len <= 0L | !ThreePUTR_gr$false_polyA
  ThreePUTR_gr_filtered <- ThreePUTR_gr[
    ThreePUTR_gr$annotated | (cpa_ok & priming_ok)
  ]
}


### Export ---------------------------------------------------------------------

# rtracklayer's GTF exporter cannot serialize zero-row metadata reliably.
if (length(ThreePUTR_gr_filtered) == 0L) {
  writeLines(character(), opt$output)
} else {
  rtracklayer::export(
    ThreePUTR_gr_filtered,
    opt$output
  )
}
