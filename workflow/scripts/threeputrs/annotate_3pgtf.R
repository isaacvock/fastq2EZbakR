#!/usr/bin/env Rscript
# Load dependencies ------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(rtracklayer)
  library(optparse)
  library(GenomicRanges)
})


### Helper functions -----------------------------------------------------------

missing_value <- function(x) {
  is.na(x) | trimws(x) == ""
}

check_transcript_loci <- function(features) {
  inconsistent <- features %>%
    filter(!missing_value(gene_id), !missing_value(transcript_id)) %>%
    group_by(gene_id, transcript_id) %>%
    summarise(loci = n_distinct(seqnames, strand), .groups = "drop") %>%
    filter(loci > 1L)

  if (nrow(inconsistent) > 0L) {
    stop("Inconsistent chromosome or strand for transcript: ",
         inconsistent$transcript_id[1], call. = FALSE)
  }
}

fill_utr_ids <- function(features) {
  missing <- missing_value(features$utr_id)
  if (!any(missing)) return(features)

  if (any(missing & (missing_value(features$gene_id) |
                     !features$strand %in% c("+", "-")))) {
    stop("Cannot assign missing utr_id values: 3UTR records require a gene_id ",
         "and a '+' or '-' strand.", call. = FALSE)
  }

  # Multiple rows for a custom transcript represent one feature. Without a
  # transcript_id, only identical stranded intervals share a generated ID.
  grouped <- features %>%
    mutate(
      transcript_id = ifelse(missing_value(transcript_id), NA_character_, transcript_id),
      group_start = ifelse(is.na(transcript_id), start, NA_integer_),
      group_end = ifelse(is.na(transcript_id), end, NA_integer_)
    ) %>%
    group_by(gene_id, transcript_id, seqnames, strand, group_start, group_end) %>%
    mutate(group_id = cur_group_id()) %>%
    ungroup()

  groups <- grouped %>%
    group_by(group_id, gene_id, transcript_id, seqnames, strand, group_start, group_end) %>%
    summarise(
      needs_id = any(missing_value(utr_id)),
      existing_ids = list(unique(utr_id[!missing_value(utr_id)])),
      .groups = "drop"
    ) %>%
    arrange(gene_id, transcript_id, seqnames, strand, group_start, group_end)

  # Reserve all supplied IDs, including ones belonging to other transcripts.
  used <- new.env(hash = TRUE, parent = emptyenv())
  for (id in unique(features$utr_id[!missing])) assign(id, TRUE, envir = used)
  assigned <- rep(NA_character_, nrow(groups))
  previous_gene <- NULL
  suffix <- 0L

  for (i in which(groups$needs_id)) {
    ids <- groups$existing_ids[[i]]
    if (length(ids) > 1L) {
      stop("Cannot assign missing utr_id: multiple existing IDs for gene/transcript ",
           groups$gene_id[i], "/", groups$transcript_id[i], call. = FALSE)
    }
    if (length(ids) == 1L) {
      assigned[groups$group_id[i]] <- ids
      next
    }

    gene <- groups$gene_id[i]
    if (!identical(gene, previous_gene)) suffix <- 0L
    repeat {
      suffix <- suffix + 1L
      id <- paste0(gene, "_utr_", suffix)
      if (!exists(id, envir = used, inherits = FALSE)) break
    }
    assign(id, TRUE, envir = used)
    assigned[groups$group_id[i]] <- id
    previous_gene <- gene
  }

  features$utr_id[missing] <- assigned[grouped$group_id[missing]]
  features
}


### Augment the provided annotation -------------------------------------------

annotate_3pgtf <- function(input, output) {
  if (normalizePath(input, mustWork = TRUE) ==
      normalizePath(output, mustWork = FALSE)) {
    stop("Input and output GTF paths must differ.", call. = FALSE)
  }

  gtf <- rtracklayer::import(input)
  annotation <- as_tibble(gtf)
  for (column in c("gene_id", "transcript_id", "utr_id")) {
    if (!column %in% names(annotation)) annotation[[column]] <- NA_character_
  }
  annotation <- annotation %>%
    mutate(seqnames = as.character(seqnames), strand = as.character(strand))

  if (any(annotation$type == "3UTR")) {
    custom <- annotation %>% filter(type == "3UTR")
    check_transcript_loci(custom)
    custom <- fill_utr_ids(custom)
    annotation$utr_id[annotation$type == "3UTR"] <- custom$utr_id
    message("Reusing ", nrow(custom), " custom 3UTR records.")
  } else {
    exons <- annotation %>% filter(type == "exon")
    usable <- !missing_value(exons$gene_id) &
      !missing_value(exons$transcript_id) & exons$strand %in% c("+", "-")
    if (any(!usable)) {
      message("Skipping ", sum(!usable),
              " exon records with missing gene_id/transcript_id or unusable strand.")
    }
    exons <- exons[usable, ]
    if (nrow(exons) == 0L) {
      stop("No usable 3UTR features: provide custom 3UTR records or exons with ",
           "gene_id, transcript_id, and '+'/'-' strands.", call. = FALSE)
    }
    check_transcript_loci(exons)

    terminal <- exons %>%
      select(seqnames, start, end, strand, gene_id, transcript_id) %>%
      distinct() %>%
      mutate(
        terminal_boundary = ifelse(strand == "+", -end, start),
        other_boundary = ifelse(strand == "+", -start, end)
      ) %>%
      arrange(gene_id, transcript_id, terminal_boundary, other_boundary) %>%
      group_by(gene_id, transcript_id) %>%
      slice_head(n = 1L) %>%
      ungroup() %>%
      select(-terminal_boundary, -other_boundary) %>%
      mutate(type = "3UTR", source = "fastq2EZbakR", utr_id = NA_character_)

    terminal <- fill_utr_ids(terminal)
    annotation <- bind_rows(annotation, terminal)
    message("Added ", nrow(terminal), " terminal-exon 3UTR features.")
  }

  augmented <- makeGRangesFromDataFrame(annotation, keep.extra.columns = TRUE)
  rtracklayer::export(augmented, output, format = "gtf")
}


### Parse command line arguments -----------------------------------------------

if (sys.nframe() == 0L) {
  option_list <- list(
    make_option("--gtf", type = "character", help = "Path to input GTF file"),
    make_option("--output", type = "character", help = "Path to augmented GTF file")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  if (is.null(opt$gtf) || is.null(opt$output)) {
    stop("Both --gtf and --output are required.", call. = FALSE)
  }
  annotate_3pgtf(opt$gtf, opt$output)
}
