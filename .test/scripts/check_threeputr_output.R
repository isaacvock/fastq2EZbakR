#!/usr/bin/env Rscript
# Validate the annotation-only workflow through featureCounts and the final cB.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 7L)
stopifnot(!any(file.exists(c(
  "results/informative_read", "results/bam2bg", "results/merge_3pend_bg",
  "results/call_PAS", "results/summarise_PAS_clusters", "annotations/threepUTR_annotation.gtf"
))))
gtf <- rtracklayer::import(args[1])
ids <- gtf$utr_id[gtf$type == "3UTR"]
stopifnot(length(ids) > 0L, !anyNA(ids), !anyDuplicated(ids))

counted_ids <- character()
assigned_ids <- character()
for (path in args[4:5]) {
  counts <- read.delim(path, comment.char = "#", check.names = FALSE)
  stopifnot(setequal(counts$Geneid, ids))
  counted_ids <- union(counted_ids, counts$Geneid[counts[[ncol(counts)]] > 0])
}
for (path in args[6:7]) {
  assignments <- read.delim(path, header = FALSE, quote = "", comment.char = "")
  assigned <- assignments[[4]][assignments[[3]] > 0]
  assigned_ids <- union(assigned_ids, unlist(strsplit(assigned, ",", fixed = TRUE)))
}
stopifnot(length(counted_ids) > 0L, length(assigned_ids) > 0L,
          all(assigned_ids %in% ids), setequal(counted_ids, assigned_ids))

cb <- read.csv(gzfile(args[2]), check.names = FALSE)
stopifnot("threepUTR" %in% names(cb))
features <- cb$threepUTR[!is.na(cb$threepUTR) & cb$threepUTR != "__no_feature"]
cb_ids <- unique(unlist(strsplit(features, "+", fixed = TRUE)))
stopifnot(length(cb_ids) > 0L, all(cb_ids %in% assigned_ids))
writeLines("Generated 3UTR IDs reach featureCounts assignments and the final cB.", args[3])
