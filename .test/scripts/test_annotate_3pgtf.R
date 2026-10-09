#!/usr/bin/env Rscript
# Small regression fixtures; run in the workflow's Rbio environment.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
source(args[1])

run_tests <- function() {
  test_dir <- tempfile("threeputr-")
  dir.create(test_dir)
  on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)
  input <- file.path(test_dir, "input.gtf")
  output <- file.path(test_dir, "output.gtf")

  row <- function(start, end, strand, gene = "g", tx = "t", type = "exon",
                  extra = "", chr = "chr1") {
    attributes <- c(
      if (!is.na(gene)) paste0('gene_id "', gene, '";'),
      if (!is.na(tx)) paste0('transcript_id "', tx, '";'),
      extra
    )
    paste(chr, "fixture", type, start, end, ".", strand, ".",
          paste(attributes, collapse = " "), sep = "\t")
  }
  run <- function(rows) {
    writeLines(rows, input)
    before <- tools::md5sum(input)
    annotate_3pgtf(input, output)
    stopifnot(identical(before, tools::md5sum(input)))
    as_tibble(rtracklayer::import(output))
  }
  expect_error <- function(rows, pattern) {
    writeLines(rows, input)
    error <- tryCatch({
      annotate_3pgtf(input, output)
      NULL
    }, error = identity)
    stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
  }
  utrs <- function(x) {
    x %>% filter(type == "3UTR") %>%
      select(gene_id, transcript_id, start, end, strand, utr_id) %>%
      arrange(gene_id, transcript_id)
  }

  rows <- c(
    row(10, 60, "+", "plus", NA, "gene", 'gene_name "Alpha beta"; note "keep me";'),
    row(10, 20, "+", "plus", "a"),
    row(40, 60, "+", "plus", "a"),
    row(40, 60, "+", "plus", "a"),
    row(40, 60, "+", "plus", "b"),
    row(100, 120, "-", "minus", "m"),
    row(150, 200, "-", "minus", "m"),
    row(300, 310, "+", "single", "s"),
    row(400, 450, "+", "tie_plus", "p"),
    row(410, 450, "+", "tie_plus", "p"),
    row(500, 550, "-", "tie_minus", "n"),
    row(500, 540, "-", "tie_minus", "n"),
    row(600, 610, "+", NA, "missing_gene"),
    row(600, 610, "+", "missing_tx", NA),
    row(600, 610, ".", "unstranded", "u")
  )
  messages <- capture.output(result <- run(rows), type = "message")
  stopifnot(any(grepl("Skipping 3 exon records", messages)))
  original <- as_tibble(rtracklayer::import(input)) %>%
    mutate(across(where(is.factor), as.character))
  preserved <- result[seq_len(nrow(original)), names(original)] %>%
    mutate(across(where(is.factor), as.character))
  difference <- all.equal(original, preserved)
  if (!isTRUE(difference)) stop(paste(difference, collapse = "\n"))
  features <- utrs(result)
  stopifnot(
    nrow(features) == 6L,
    !anyDuplicated(features$utr_id),
    features$start[features$gene_id == "minus"] == 100L,
    features$end[features$gene_id == "minus"] == 120L,
    all(features$start[features$gene_id == "plus"] == 40L),
    identical(features$utr_id[features$gene_id == "plus"], c("plus_utr_1", "plus_utr_2")),
    features$start[features$gene_id == "single"] == 300L,
    features$start[features$gene_id == "tie_plus"] == 410L,
    features$end[features$gene_id == "tie_minus"] == 540L
  )
  set.seed(1)
  stopifnot(identical(features, utrs(run(sample(rows)))))

  # The augmented file can be supplied again without adding duplicate features.
  second_output <- file.path(test_dir, "second.gtf")
  annotate_3pgtf(output, second_output)
  stopifnot(identical(as_tibble(rtracklayer::import(output)),
                      as_tibble(rtracklayer::import(second_output))))

  custom <- run(c(
    row(1, 10, "+", tx = "a", type = "3UTR", extra = 'utr_id "g_utr_1"; note "custom";'),
    row(20, 30, "+", tx = "a", type = "3UTR"),
    row(40, 50, "+", tx = "b", type = "3UTR"),
    row(60, 70, "-", tx = NA, type = "3UTR"),
    row(60, 70, "-", tx = NA, type = "3UTR"),
    row(80, 90, "+", tx = "unused")
  ))
  stopifnot(
    nrow(custom) == 6L,
    identical(custom$utr_id[1:5], c("g_utr_1", "g_utr_1", "g_utr_2", "g_utr_3", "g_utr_3")),
    custom$note[1] == "custom",
    is.na(custom$utr_id[6])
  )
  # Fully specified custom IDs require no exon annotation or identifier changes.
  custom <- run(row(1, 10, "+", gene = NA, tx = NA, type = "3UTR",
                    extra = 'utr_id "provided";'))
  stopifnot(nrow(custom) == 1L, custom$utr_id == "provided")

  expect_error(row(1, 10, "+", type = "gene"), "No usable 3UTR")
  expect_error(row(1, 10, "+", tx = NA), "No usable 3UTR")
  expect_error(c(row(1, 10, "+"), row(20, 30, "-")), "Inconsistent chromosome or strand")
  expect_error(c(row(1, 10, "+"), row(20, 30, "+", chr = "chr2")),
               "Inconsistent chromosome or strand")
  expect_error(row(1, 10, "+", gene = NA, type = "3UTR"), "require a gene_id")
  expect_error(row(1, 10, ".", type = "3UTR"), "require a gene_id")
  expect_error(c(
    row(1, 10, "+", type = "3UTR", extra = 'utr_id "first";'),
    row(20, 30, "+", type = "3UTR", extra = 'utr_id "second";'),
    row(40, 50, "+", type = "3UTR")
  ), "multiple existing IDs")
  error <- tryCatch(annotate_3pgtf(input, input), error = identity)
  stopifnot(inherits(error, "error"), grepl("paths must differ", conditionMessage(error)))

  # Exercise the actual CLI against the repository's exon-only annotation.
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    c(shQuote(args[1]), "--gtf", shQuote(args[2]), "--output", shQuote(output)))
  stopifnot(status == 0L)
  reference <- as_tibble(rtracklayer::import(args[2]))
  expected <- reference %>% filter(type == "exon") %>% distinct(gene_id, transcript_id)
  generated <- as_tibble(rtracklayer::import(output)) %>% filter(type == "3UTR")
  stopifnot(nrow(generated) == nrow(expected), !anyDuplicated(generated$utr_id))
}

run_tests()
writeLines("Annotation regression tests passed.", args[3])
