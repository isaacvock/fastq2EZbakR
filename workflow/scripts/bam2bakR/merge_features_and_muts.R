#!/usr/bin/env Rscript
### PURPOSE OF THIS SCRIPT
## Merge feature assignment and mutation counts tables via streaming DuckDB views

# Load dependencies ------------------------------------------------------------

library(optparse)
library(duckdb)
library(DBI)
library(glue)
library(dplyr)
library(utils)


# Process CLI options ----------------------------------------------------------

option_list <- list(
  make_option(c("-g", "--genes"), type = "logical",
              default = FALSE,
              help = "Whether reads were assigned to genes"),
  make_option(c("-e", "--exons"), type = "logical",
              default = FALSE,
              help = "Whether reads were assigned to exons"),
  make_option(c("-b", "--exonbins"), type = "logical",
              default = FALSE,
              help = "Whether reads were assigned to exonbins"),
  make_option(c("-t", "--transcripts"), type = "logical",
              default = FALSE,
              help = "Whether reads were assigned to transcripts via featureCounts"),
  make_option(c("--frombam"), type = "logical",
              default = FALSE,
              help = "Whether reads were assigned to transcripts from the transcriptome BAM"),
  make_option(c("--starjunc"), type = "logical",
              default = FALSE,
              help = "Whether reads were assigned to junctions via STAR jI/jM tags"),
  make_option(c("--threeputr"), type = "logical",
              default = FALSE,
              help = "Whether reads were assigned to 3'-UTRs"),
  make_option(c("-j", "--eej"), type = "logical",
              default = FALSE,
              help = "Whether reads were assigned to exon-exon junctions"),
  make_option(c("--eij"), type = "logical",
              default = FALSE,
              help = "Whether reads were assigned to exon-intron junctions"),
  make_option(c("-o", "--output"), type = "character",
              default = NULL,
              help = "Deprecated; retained for CLI compatibility"),
  make_option(c("-c", "--cBoutput"), type = "character",
              default = NULL,
              help = "Path to per-sample cB CSV output"),
  make_option(c("-m", "--muttypes"), type = "character",
              help = "Comma-separated mutation types to retain in cB/cUP"),
  make_option(c("-s", "--sample"), type = "character",
              help = "Sample name"),
  make_option(c("--annotation"), type = "character",
              default = NULL,
              help = "Path to annotation GTF file"),
  make_option(c("--makecB"), type = "logical",
              default = FALSE,
              help = "Whether to create a summarized cB CSV output"),
  make_option(c("--makecUP"), type = "logical",
              default = FALSE,
              help = "Whether to create a summarized cUP CSV output"),
  make_option(c("--makeArrow"), type = "logical",
              default = FALSE,
              help = "Whether to create a Hive-partitioned Parquet file output"),
  make_option(c("--cUPoutput"), type = "character",
              default = NULL,
              help = "Path to per-sample cUP CSV output"),
  make_option(c("--Arrowoutput"), type = "character",
              default = NULL,
              help = "Path to Arrow Parquet output"),
  make_option(c("--MaxMem"), type = "character",
              default = "8GB",
              help = "Max RAM usage to enforce in DuckDB"),
  make_option(c("--threads"), type = "integer",
              default = 1L,
              help = "Number of DuckDB execution threads")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

message("opt contents:\n",
        paste(capture.output(str(opt)), collapse = "\n"))

sql_path <- function(p) {
  gsub("'", "''", gsub("\\\\", "/", p))
}


# Ephemeral DuckDB setup (spills safely to disk if MaxMem is exceeded, cleaned on exit) ----

tmp_db_dir <- tempfile(pattern = paste0("duckdb_", opt$sample, "_"))
dir.create(tmp_db_dir, recursive = TRUE, showWarnings = FALSE)
db_path <- file.path(tmp_db_dir, paste0(opt$sample, ".duckdb"))

con <- dbConnect(duckdb(), dbdir = db_path, read_only = FALSE)
on.exit({
  try(dbDisconnect(con, shutdown = TRUE), silent = TRUE)
  unlink(tmp_db_dir, recursive = TRUE, force = TRUE)
}, add = TRUE)

dbExecute(con, glue("SET memory_limit='{opt$MaxMem}';"))
dbExecute(con, glue("SET threads={max(1L, as.integer(opt$threads))};"))
dbExecute(con, glue("SET temp_directory='{sql_path(tmp_db_dir)}';"))


register_feat <- function(view, path, feat_col) {
  sql <- glue("
    CREATE OR REPLACE TEMP VIEW {`view`} AS
      SELECT qname,
             REPLACE({feat_col}, ',', '+') AS {feat_col}
      FROM read_csv_auto('{sql_path(path)}',
                         header = FALSE,
                         delim  = '\t',
                         columns = {{'qname':'VARCHAR', 'status':'VARCHAR', 'nhits':'INT', '{feat_col}':'VARCHAR'}})
      WHERE nhits > 0;
  ")
  dbExecute(con, sql)
}


#### Register mutation counts view ####

counts_path <- glue("results/counts/{opt$sample}_counts.csv.gz")
dbExecute(con, glue("
  CREATE OR REPLACE TEMP VIEW muts AS
    SELECT * FROM read_csv_auto(
      '{sql_path(counts_path)}',
      header = TRUE,
      types = {{'qname': 'VARCHAR', 'rname': 'VARCHAR', 'FR': 'VARCHAR', 'sj': 'VARCHAR'}}
    );
"))


#### Register feature views ####

feature_vect <- character()

if (isTRUE(opt$genes)) {
  cat("Registering genes view...\n")
  register_feat("genes", glue("./results/featurecounts_genes/{opt$sample}.s.bam.featureCounts"), "GF")
  feature_vect <- c(feature_vect, "GF")
}

if (isTRUE(opt$exons)) {
  cat("Registering exons view...\n")
  register_feat("exons", glue("./results/featurecounts_exons/{opt$sample}.s.bam.featureCounts"), "XF")
  feature_vect <- c(feature_vect, "XF")
}

if (isTRUE(opt$transcripts)) {
  cat("Registering featureCounts transcripts view...\n")
  register_feat("fc_transcripts", glue("./results/featurecounts_transcripts/{opt$sample}.s.bam.featureCounts"), "transcripts")
  feature_vect <- c(feature_vect, "transcripts")
}

if (isTRUE(opt$exonbins)) {
  cat("Registering exonbins view...\n")
  register_feat("exonbins", glue("./results/featurecounts_exonbins/{opt$sample}.s.bam.featureCounts"), "exon_bin")
  feature_vect <- c(feature_vect, "exon_bin")
}

if (isTRUE(opt$frombam)) {
  cat("Registering TEC view...\n")
  tec_path <- glue("./results/read_to_transcripts/{opt$sample}.csv")
  dbExecute(con, glue("
    CREATE OR REPLACE TEMP VIEW transcripts AS
      SELECT qname,
             bamfile_transcripts
      FROM read_csv_auto(
        '{sql_path(tec_path)}',
        header = TRUE,
        types = {{'qname': 'VARCHAR', 'bamfile_transcripts': 'VARCHAR'}}
      );
  "))
  feature_vect <- c(feature_vect, "TEC")
}

if (isTRUE(opt$starjunc)) {
  cat("Registering junctions view...\n")
  junc_path <- glue("./results/read_to_junctions/{opt$sample}.csv.gz")
  dbExecute(con, glue("
    CREATE OR REPLACE TEMP VIEW starjunc AS
      SELECT qname,
             junction_start,
             junction_end
      FROM read_csv_auto(
        '{sql_path(junc_path)}',
        header = TRUE,
        types = {{'qname': 'VARCHAR', 'junction_start': 'VARCHAR', 'junction_end': 'VARCHAR'}}
      );
  "))
  feature_vect <- c(feature_vect, "junction_start", "junction_end")
}

if (isTRUE(opt$threeputr)) {
  cat("Registering 3'-UTR view...\n")
  pe_3utr <- glue("./results/featurecounts_3utr/{opt$sample}_informative.bam.featureCounts")
  se_3utr <- glue("./results/featurecounts_3utr/{opt$sample}.s.bam.featureCounts")
  utr_path <- if (file.exists(pe_3utr)) pe_3utr else se_3utr
  register_feat("threeputr", utr_path, "threepUTR")
  feature_vect <- c(feature_vect, "threepUTR")
}

if (isTRUE(opt$eej)) {
  cat("Registering eej view...\n")
  register_feat("eej", glue("./results/featurecounts_eej/{opt$sample}.s.bam.featureCounts"), "ee_junction_id")
  feature_vect <- c(feature_vect, "ee_junction_id")
}

if (isTRUE(opt$eij)) {
  cat("Registering eij view...\n")
  register_feat("eij", glue("./results/featurecounts_eij/{opt$sample}.s.bam.featureCounts"), "ei_junction_id")
  feature_vect <- c(feature_vect, "ei_junction_id")
}


#### Construct non-materialized LEFT JOIN view ####

feat_catalogue <- list(
  genes = list(
    flag   = "genes",
    select = "COALESCE(g.GF, '__no_feature') AS GF",
    join   = "LEFT JOIN genes g USING (qname)"
  ),
  exons = list(
    flag   = "exons",
    select = "COALESCE(e.XF, '__no_feature') AS XF",
    join   = "LEFT JOIN exons e USING (qname)"
  ),
  transcripts = list(
    flag   = "transcripts",
    select = "COALESCE(tr.transcripts, '__no_feature') AS transcripts",
    join   = "LEFT JOIN fc_transcripts tr USING (qname)"
  ),
  exonbins = list(
    flag   = "exonbins",
    select = "COALESCE(eb.exon_bin, '__no_feature') AS exon_bin",
    join   = "LEFT JOIN exonbins eb USING (qname)"
  ),
  frombam = list(
    flag   = "frombam",
    select = "COALESCE(t.bamfile_transcripts, '__no_feature') AS TEC",
    join   = "LEFT JOIN transcripts t USING (qname)"
  ),
  starjunc = list(
    flag   = "starjunc",
    select = "COALESCE(sj.junction_start, '__no_feature') AS junction_start,\n         COALESCE(sj.junction_end, '__no_feature') AS junction_end",
    join   = "LEFT JOIN starjunc sj USING (qname)"
  ),
  threeputr = list(
    flag   = "threeputr",
    select = "COALESCE(tpu.threepUTR, '__no_feature') AS threepUTR",
    join   = "LEFT JOIN threeputr tpu USING (qname)"
  ),
  eej = list(
    flag   = "eej",
    select = "COALESCE(eej.ee_junction_id, '__no_feature') AS ee_junction_id",
    join   = "LEFT JOIN eej eej USING (qname)"
  ),
  eij = list(
    flag   = "eij",
    select = "COALESCE(eij.ei_junction_id, '__no_feature') AS ei_junction_id",
    join   = "LEFT JOIN eij eij USING (qname)"
  )
)

select_fragments <- c("m.*")
join_fragments   <- character()

for (feat in feat_catalogue) {
  if (isTRUE(opt[[feat$flag]])) {
    select_fragments <- c(select_fragments, feat$select)
    join_fragments   <- c(join_fragments,   feat$join)
  }
}

if (length(join_fragments) == 0L) {
  sql <- glue("
    CREATE OR REPLACE TEMP VIEW merged AS
    SELECT {paste(select_fragments, collapse = ',\n         ')}
    FROM   muts m;
  ")
} else {
  sql <- glue("
    CREATE OR REPLACE TEMP VIEW merged AS
    SELECT {paste(select_fragments, collapse = ',\n         ')}
    FROM   muts m
    {paste(join_fragments, collapse = '\n    ')};
  ")
}

cat(sql, "\n")
dbExecute(con, sql)


#### Summarize cB (shared by CSV cB and Hive Parquet outputs) ####

mut_cols     <- trimws(strsplit(opt$muttypes, ",")[[1]])
base_cols    <- unique(paste0("n", substr(mut_cols, 1, 1)))
feature_cols <- feature_vect

quoted_mut_cols  <- DBI::dbQuoteIdentifier(con, mut_cols)
quoted_base_cols <- DBI::dbQuoteIdentifier(con, base_cols)
quoted_feat_cols <- if (length(feature_cols) > 0L) DBI::dbQuoteIdentifier(con, feature_cols) else character()
sel_cols         <- c(quoted_mut_cols, quoted_base_cols)
sample_lit       <- DBI::dbQuoteString(con, opt$sample)

if (isTRUE(opt$makecB) || isTRUE(opt$makeArrow)) {
  if (length(quoted_feat_cols) == 0L) {
    dbExecute(con, glue("
      CREATE OR REPLACE TEMP TABLE cB AS
        SELECT {sample_lit} AS sample,
               rname, sj,
               {paste(sel_cols, collapse = ', ')},
               COUNT(*) AS n
        FROM   merged
        GROUP  BY ALL;
    "))
  } else {
    dbExecute(con, glue("
      CREATE OR REPLACE TEMP TABLE cB AS
        SELECT {sample_lit} AS sample,
               rname, sj, {paste(quoted_feat_cols, collapse = ', ')},
               {paste(sel_cols, collapse = ', ')},
               COUNT(*) AS n
        FROM   merged
        GROUP  BY ALL;
    "))
  }
}


#### Write per-sample cB ####

if (!is.null(opt$cBoutput) && nzchar(opt$cBoutput)) {
  dir.create(dirname(opt$cBoutput), recursive = TRUE, showWarnings = FALSE)
  if (isTRUE(opt$makecB)) {
    dbExecute(con, glue("
      COPY cB TO '{sql_path(opt$cBoutput)}' (FORMAT CSV, HEADER);
    "))
  } else {
    write.csv(tibble(), file = opt$cBoutput)
  }
}


#### Write Hive-partitioned Arrow Parquet dataset ####
# Retains `sample` column inside part-0.parquet and Snappy compression
# for 100% backwards compatibility with all versions of EZbakR.

if (isTRUE(opt$makeArrow) && !is.null(opt$Arrowoutput) && nzchar(opt$Arrowoutput)) {
  dir.create(dirname(opt$Arrowoutput), recursive = TRUE, showWarnings = FALSE)
  dbExecute(con, glue("
    COPY cB TO '{sql_path(opt$Arrowoutput)}' (FORMAT PARQUET, COMPRESSION snappy);
  "))
}


#### Write per-sample cUP ####

if (!is.null(opt$cUPoutput) && nzchar(opt$cUPoutput)) {
  dir.create(dirname(opt$cUPoutput), recursive = TRUE, showWarnings = FALSE)
  if (isTRUE(opt$makecUP)) {
    group_cols      <- c(DBI::dbQuoteIdentifier(con, "rname"), quoted_feat_cols, quoted_mut_cols)
    group_by_clause <- paste(group_cols, collapse = ", ")
    avg_fragments   <- paste(
      sprintf("AVG(%s) AS %s", quoted_base_cols, quoted_base_cols),
      collapse = ",\n           "
    )

    select_clause <- glue("
      {sample_lit} AS sample,
      {group_by_clause},
      {avg_fragments},
      COUNT(*) AS n
    ")

    dbExecute(con, glue("
      CREATE OR REPLACE TEMP TABLE cUP AS
      SELECT {select_clause}
      FROM   merged
      GROUP  BY {group_by_clause};
    "))

    dbExecute(con, glue("
      COPY cUP TO '{sql_path(opt$cUPoutput)}' (FORMAT CSV, HEADER);
    "))
  } else {
    write.csv(tibble(), file = opt$cUPoutput)
  }
}
