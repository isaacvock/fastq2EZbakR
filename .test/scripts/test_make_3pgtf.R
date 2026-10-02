#!/usr/bin/env Rscript
# Exercise the production CLI with small, mirrored genome/GTF/cluster fixtures.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
caller <- normalizePath(args[1])
suppressPackageStartupMessages({
  library(Biostrings)
  library(GenomicRanges)
  library(rtracklayer)
})

run_tests <- function() {
  test_dir <- tempfile("de-novo-threeputr-")
  dir.create(test_dir)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

  case <- function(name, downstream = "", annotated = TRUE, site = 100L,
                   size = 400L, upstream = "", motif = "", expected = FALSE) {
    list(name = name, downstream = downstream, annotated = annotated,
         site = site, size = size, upstream = upstream, motif = motif,
         expected = expected)
  }

  run <- function(cases, polyA = 7L, CPA = FALSE, only = FALSE,
                  fasta = TRUE, error = NULL) {
    gtf <- bed_plus <- bed_minus <- character()
    sequences <- DNAStringSet()
    expected_flags <- logical()
    expected_coords <- list()
    for (fixture in cases) {
      sequence <- paste(rep("C", fixture$size), collapse = "")
      if (nzchar(fixture$downstream) && fixture$site < fixture$size) {
        count <- min(nchar(fixture$downstream), fixture$size - fixture$site)
        substr(sequence, fixture$site + 1L, fixture$site + count) <-
          substr(fixture$downstream, 1L, count)
      }
      if (nzchar(fixture$upstream)) {
        substr(sequence, fixture$site - nchar(fixture$upstream),
               fixture$site - 1L) <- fixture$upstream
      }
      if (nzchar(fixture$motif)) {
        substr(sequence, fixture$site - 30L,
               fixture$site - 30L + nchar(fixture$motif) - 1L) <- fixture$motif
      }
      # The terminal exon covers the site unless a later exon is added.
      exons <- matrix(c(max(1L, fixture$site - 10L),
                        min(fixture$size, fixture$site + 10L)), ncol = 2L)
      if (exons[1, 1] > 5L) exons <- rbind(c(1L, 5L), exons)
      if (!fixture$annotated) {
        exons <- rbind(exons, c(fixture$site + 40L, fixture$site + 50L))
      }

      for (strand in c("+", "-")) {
        id <- paste0(fixture$name, if (strand == "+") "_plus" else "_minus")
        chr <- paste0("chr_", id)
        seq <- DNAString(sequence)
        site <- fixture$site
        exon_coords <- exons
        if (strand == "-") {
          seq <- reverseComplement(seq)
          site <- fixture$size - site + 1L
          exon_coords <- cbind(fixture$size - exons[, 2] + 1L,
                               fixture$size - exons[, 1] + 1L)
        }
        sequences <- c(sequences, setNames(DNAStringSet(seq), chr))
        row <- function(type, start, end) {
          attributes <- paste0('gene_id "', id, '";')
          if (type == "exon") {
            attributes <- paste(attributes, paste0('transcript_id "t_', id, '";'))
          }
          paste(chr, "fixture", type, start, end, ".", strand, ".",
                attributes, sep = "\t")
        }
        gtf <- c(gtf, row("gene", 1L, fixture$size))
        for (i in seq_len(nrow(exon_coords))) {
          gtf <- c(gtf, row("exon", exon_coords[i, 1], exon_coords[i, 2]))
        }
        bed <- paste(chr, 1L, site - 1L, site, 100L, sep = "\t")
        if (strand == "+") bed_plus <- c(bed_plus, bed) else bed_minus <- c(bed_minus, bed)
        expected_flags[id] <- fixture$expected
        expected_coords[[id]] <- site
      }
    }
    input_gtf <- file.path(test_dir, "input.gtf")
    plus <- file.path(test_dir, "plus.bg")
    minus <- file.path(test_dir, "minus.bg")
    genome <- file.path(test_dir, "genome.fa")
    output <- file.path(test_dir, "output.gtf")
    log <- file.path(test_dir, "caller.log")
    writeLines(gtf, input_gtf)
    writeLines(bed_plus, plus)
    writeLines(bed_minus, minus)
    writeXStringSet(sequences, genome)
    cli <- c(caller, "--bed_plus", plus, "--bed_minus", minus,
             "--gtf", input_gtf, "--output", output,
             "--min_coverage", "1", "--false_polyA_len", as.character(polyA),
             "--require_CPA", as.character(CPA), "--only_annotated", as.character(only))
    if (fasta) cli <- c(cli, "--fasta", genome)
    status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(cli),
                      stdout = log, stderr = log)
    if (!is.null(error)) {
      stopifnot(status != 0L, any(grepl(error, readLines(log))))
      return(invisible(NULL))
    }
    if (status != 0L) stop(paste(readLines(log), collapse = "\n"))
    stopifnot(file.exists(output))
    if (file.info(output)$size == 0L) return(GRanges())
    result <- rtracklayer::import(output)
    stopifnot(all(start(result) == end(result)),
              all(start(result) == unlist(expected_coords[result$gene_id])))
    if (polyA > 0L) {
      stopifnot(!anyNA(result$false_polyA),
                identical(as.logical(result$false_polyA),
                          unname(expected_flags[result$gene_id])))
    } else {
      stopifnot(is.null(result$false_polyA) || all(is.na(result$false_polyA)))
    }
    result
  }

  annotated <- list(
    case("clean"),
    case("run", "CCCAAAAAAACCCCC", expected = TRUE),
    case("density", "AAACAACACA", expected = TRUE),
    case("below", "AAACAACACC"),
    case("t_run", "TTTTTTTTTT"),
    case("upstream", upstream = "AAAAAAA"),
    case("ambiguous", "AAANAAANNC"),
    case("ambiguous_rich", "AAACAACANA", expected = TRUE),
    case("late_density", "CCCCCCCCCCAAACAACACA"),
    case("short_density", "AAACAAAA", site = 92L, size = 100L),
    case("short_run", "AAAAAAA", site = 93L, size = 100L, expected = TRUE),
    case("empty_window", site = 400L),
    case("first_base", site = 1L)
  )
  result <- run(annotated)
  stopifnot(length(result) == 2L * length(annotated), all(as.logical(result$annotated)))
  # Custom tract length affects A runs but not the fixed density criterion.
  custom <- run(list(case("custom", "CCCCCCCCCCAAAAA", expected = TRUE)), polyA = 5L)
  stopifnot(length(custom) == 2L)
  custom <- run(list(case("density", "AAACAACACA", expected = TRUE)), polyA = 20L)
  stopifnot(length(custom) == 2L)

  novel <- list(
    case("novel_run", "AAAAAAA", annotated = FALSE, expected = TRUE),
    case("novel_density", "AAACAACACA", annotated = FALSE, expected = TRUE),
    case("novel_clean", annotated = FALSE)
  )
  result <- run(c(annotated, novel), only = TRUE)
  stopifnot(length(result) == 2L * length(annotated), all(as.logical(result$annotated)))
  result <- run(c(list(case("baseline")), novel))
  stopifnot(setequal(result$gene_id, c("baseline_plus", "baseline_minus",
                                     "novel_clean_plus", "novel_clean_minus")))
  clean_novel <- result[grepl("novel_clean", result$gene_id)]
  stopifnot(!any(as.logical(clean_novel$annotated)))

  result <- run(c(annotated, novel), polyA = 0L, fasta = FALSE)
  stopifnot(length(result) == 2L * (length(annotated) + length(novel)))
  result <- run(list(case("baseline")), polyA = -1L, only = TRUE, fasta = FALSE)
  stopifnot(length(result) == 2L)
  run(list(case("baseline")), only = TRUE, fasta = FALSE, error = "FASTA")

  cpa <- list(case("annotated_no_motif"),
              case("canonical", annotated = FALSE, motif = "AATAAA"),
              case("variant", annotated = FALSE, motif = "ATTAAA"),
              case("no_motif", annotated = FALSE))
  result <- run(cpa, polyA = 0L, CPA = TRUE)
  stopifnot(length(result) == 6L, !any(grepl("^no_motif", result$gene_id)),
            all(as.logical(result$has_CPA[!as.logical(result$annotated)])))
  run(cpa, polyA = 0L, CPA = TRUE, fasta = FALSE, error = "FASTA")
  result <- run(list(case("novel", annotated = FALSE)), polyA = 0L,
                only = TRUE, fasta = FALSE)
  stopifnot(length(result) == 0L)
  result <- run(list(case("novel", "AAAAAAA", annotated = FALSE, expected = TRUE)))
  stopifnot(length(result) == 0L)

  writeLines("De novo coordinates, terminal exons, priming flags, CPA motifs and export paths pass.",
             args[2])
}

run_tests()
