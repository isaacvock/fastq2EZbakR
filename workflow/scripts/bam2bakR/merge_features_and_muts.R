#!/usr/bin/env Rscript
### PURPOSE OF THIS SCRIPT
## Merge feature assignment and mutation counts tables

# Load dependencies ------------------------------------------------------------

library(optparse)
library(duckdb)
library(DBI)
library(glue)
library(dplyr)
library(utils)


# Process CLI options ----------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

option_list <- list(
  make_option(c("-g", "--genes", type = "logical"),
              default = "FALSE",
              help = "Whether reads were assigned to genes"),
  make_option(c("-e", "--exons", type = "logical"),
              default = "FALSE",
              help = "Whether reads were assigned to exons"),
  make_option(c("-b", "--exonbins", type = "logical"),
              default = "FALSE",
              help = "Whether reads were assigned to exonbins"),
  make_option(c("-t", "--transcripts", type = "logical"),
              default = "FALSE",
              help = "Whether reads were assigned to transcripts"),
  make_option(c("--frombam", type = "logical"),
              default = "FALSE",
              help = "Whether reads were assigned to transcripts from the 
              transcriptome aligned bam file directly. This is more accurate
              than featureCounts based transcript isoform assignment as
              the latter does not account for the splice junctions a read
              is mapped across."),
  make_option(c("--starjunc", type = "logical"),
              default = "FALSE",
              help = "Whether reads were assigned to junctions via relevant STAR tags."),              
  make_option(c("-j", "--eej", type = "logical"),
              default = "FALSE",
              help = "Whether reads were assigned to exon-exon junctions"),
  make_option(c("--eij", type = "logical"),
              default = "FALSE",
              help = "Whether reads were assigned to exon-intron junctions"),
  make_option(c("-o", "--output", type = "character"),
              help = "Path to full mutation counts/feature assignment output."),
  make_option(c("-c", "--cBoutput", type = "character"),
              help = "Path to cB output; same as full output with some columns averaged out."),
  make_option(c("-m", "--muttypes", type = "character"),
              help = "String of comma separated mutation types to keep in cBs."),
  make_option(c("-s", "--sample", type = "character"),
              help = "Sample name"),
  make_option(c("--annotation", type = "character"),
              help = "Path to annotation GTF file."),
  make_option(c("--makecB", "logical"),
              default = "FALSE",
              help = "Whether to create a summarized cB file as output."),
  make_option(c("--makecUP", "logical"),
              default = "FALSE",
              help = "Whether to create a summarized cUP file as output."),
  make_option(c("--makeArrow", "logical"),
              default = "FALSE",
              help = "Whether to create a summarized arrow parquet file as output."),
  make_option(c("--cUPoutput", type = "character"),
              help = "Path to cUP output; same as full output with some columns averaged out."),
  make_option(c("--Arrowoutput", type = "character"),
              help = "Path to arrow parquet file output; same as full output with some columns averaged out."),
  make_option(c("--MaxMem", type = "character"),
              default = "8GB",
              help = "Max RAM usage to (mostly) enforce in DuckDB")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser) # Load options from command line.
# 
# opt <- list(
#   genes = TRUE,
#   exons = TRUE,
#   exonbins = TRUE,
#   frombam = TRUE,
#   starjunc = TRUE,
#   eej = FALSE,
#   eij = FALSE,
#   output = "C:/Users/isaac/Documents/Simon_Lab/Sandbox/fastq2EZbakR/Data/new_output_test.csv.gz",
#   cBoutput = "C:/Users/isaac/Documents/Simon_Lab/Sandbox/fastq2EZbakR/Data/new_cB_test.csv.gz",
#   muttypes = "TC",
#   sample = "DMSO_8hr_1",
#   makecB = TRUE,
#   makecUP = TRUE,
#   makeArrow = TRUE,
#   cUPout = "C:/Users/isaac/Documents/Simon_Lab/Sandbox/fastq2EZbakR/Data/new_cUP_test.csv.gz",
#   ArrowOutput = "C:/Users/isaac/Documents/Simon_Lab/Sandbox/fastq2EZbakR/Data/sample=DMSO_8hr_1/new_cB_test.parquet"
# )

message("opt contents:\n",
        paste(capture.output(str(opt)), collapse = "\n"))

# DuckDB strategy --------------------------------------------------------------


dir.create("./results/merge_features_and_muts/duckdb", showWarnings = FALSE)
con <- dbConnect(duckdb(), dbdir = glue("./results/merge_features_and_muts/duckdb/{opt$sample}.duckdb"), read_only = FALSE)
dbExecute(con, glue("SET memory_limit='{opt$MaxMem}';"))


register_feat <- function(view, path, feat_col)
{
  sql <- glue("
    CREATE OR REPLACE VIEW {`view`} AS
      SELECT qname,
             REPLACE({feat_col}, ',', '+') AS {feat_col}
      FROM read_csv_auto('{path}',
                         header = FALSE,
                         delim  = '\t', 
                         columns = {{'qname':'VARCHAR', 'status':'VARCHAR', 'nhits':'INT', '{feat_col}':'VARCHAR'}})
      WHERE nhits > 0;
  ")
  dbExecute(con, sql)
}


#### "Load" mutation counts ####

dbExecute(con, glue("
  CREATE OR REPLACE VIEW muts AS
    SELECT * FROM read_csv_auto(
      'results/counts/{opt$sample}_counts.csv.gz', header = TRUE, types = {{'FR': 'VARCHAR'}});
"))


#### "Load" feature tables ####

feature_vect <- c()


if(opt$genes){
  
  cat("Making genes table...")
  
  register_feat("genes", glue("./results/featurecounts_genes/{opt$sample}.s.bam.featureCounts"), "GF")
  
  feature_vect <- c(feature_vect, "GF")
  
}


if(opt$frombam){
  
  cat("Making TEC table...")
  
  
  dbExecute(con, glue("
  CREATE OR REPLACE VIEW transcripts AS
    SELECT qname,
           bamfile_transcripts
    FROM read_csv_auto(
      './results/read_to_transcripts/{opt$sample}.csv', header = TRUE);
  "))
  
  feature_vect <- c(feature_vect, "TEC")
}


if(opt$exons){
  
  cat("Making exons table...")
  
  
  register_feat("exons", glue("./results/featurecounts_exons/{opt$sample}.s.bam.featureCounts"), "XF")
  
  feature_vect <- c(feature_vect, "XF")
  
}


if(opt$exonbins){
  
  cat("Making exonbins table...")
  
  
  register_feat("exonbins", glue("results/featurecounts_exonbins/{opt$sample}.s.bam.featureCounts"), "exon_bin")
  
  feature_vect <- c(feature_vect, "exon_bin")
  
}

if(opt$eej){
  
  cat("Making eej table...")
  
  
  register_feat("eej", glue("results/featurecounts_eej/{opt$sample}.s.bam.featureCounts"), "ee_junction_id")
  
  feature_vect <- c(feature_vect, "ee_junction_id")
  
}

if(opt$eij){
  
  cat("Making eij table...")
  
  
  register_feat("eij", glue("results/featurecounts_eij/{opt$sample}.s.bam.featureCounts"), "ei_junction_id")
  
  feature_vect <- c(feature_vect, "ei_junction_id")
  
}


if(opt$starjunc){
  
  cat("Making junctions table...")
  
  
  dbExecute(con, glue("
  CREATE OR REPLACE VIEW starjunc AS
    SELECT qname,
           junction_start,
           junction_end
    FROM read_csv_auto(
      'results/read_to_junctions/{opt$sample}.csv.gz', header = TRUE);
  "))
  
  feature_vect <- c(feature_vect, "junction_start", "junction_end")
  
}

#### Construct query ####

feat_catalogue <- list(
  genes = list(
    flag   = "genes",
    select = "COALESCE(g.GF,        '__no_feature') AS GF",
    join   = "LEFT JOIN genes        g  USING (qname)"
  ),
  exons = list(
    flag   = "exons",
    select = "COALESCE(e.XF,        '__no_feature') AS XF",
    join   = "LEFT JOIN exons        e  USING (qname)"
  ),
  exonbins = list(
    flag   = "exonbins",
    select = "COALESCE(eb.exon_bin, '__no_feature') AS exon_bin",
    join   = "LEFT JOIN exonbins     eb USING (qname)"
  ),
  transcripts = list(
    flag   = "frombam",                                # because opt$frombam
    select = "COALESCE(t.bamfile_transcripts, '__no_feature') AS TEC",
    join   = "LEFT JOIN transcripts   t  USING (qname)"
  ),
  starjunc = list(
    flag   = "starjunc",
    select = "COALESCE(sj.junction_start, '__no_feature') AS junction_start,
              COALESCE(sj.junction_end,   '__no_feature') AS junction_end",
    join   = "LEFT JOIN starjunc      sj USING (qname)"
  ),
  eej = list(
    flag   = "eej",                                
    select = "COALESCE(eej.ee_junction_id, '__no_feature') AS ee_junction_id",
    join   = "LEFT JOIN eej   eej  USING (qname)"
  ),
  eij = list(
    flag   = "eij",                                
    select = "COALESCE(eij.ei_junction_id, '__no_feature') AS ei_junction_id",
    join   = "LEFT JOIN eij   eij  USING (qname)"
  )
)


select_fragments <- c("m.*")   # always start with all mutation columns
join_fragments   <- character()

for (feat in feat_catalogue) {
  if (opt[[feat$flag]]) {
    select_fragments <- c(select_fragments, feat$select)
    join_fragments   <- c(join_fragments,   feat$join)
  }
}



#### LEFT JOIN ####

if(length(join_fragments) == 0){

  sql <- glue("
    CREATE OR REPLACE TABLE merged AS
    SELECT {paste(select_fragments, collapse = ',\n         ')}
    FROM   muts m
  ;")

}else{

  sql <- glue("
    CREATE OR REPLACE TABLE merged AS
    SELECT {paste(select_fragments, collapse = ',\n         ')}
    FROM   muts m
    {paste(join_fragments, collapse = '\n  ')}
  ;")
}



# Print for debugging
cat(sql)

dbExecute(con, sql)


# Write to csv
dbExecute(con, glue("
  COPY merged TO '{opt$output}' (FORMAT CSV, HEADER 1);
"))


#### cB table creation ####
mut_cols  <- strsplit(opt$muttypes, ",")[[1]]
base_cols <- unique(paste0('n', substr(mut_cols, 1, 1)))
feature_cols <- feature_vect

sel_cols <- DBI::dbQuoteIdentifier(con, c(mut_cols, base_cols))

cat(opt$muttypes)
cat(mut_cols)
cat(base_cols)
cat(opt$makecB)
cat(opt$makeArrow)

# Create summarized cB
if(length(join_fragments) == 0){

  dbExecute(con, glue("
  CREATE OR REPLACE TABLE cB AS
    SELECT '{opt$sample}'      AS sample,
            rname, sj,
            {paste(sel_cols, collapse = ',')},
            COUNT(*)            AS n
    FROM   merged
    GROUP  BY ALL;
  "))

}else{

  dbExecute(con, glue("
  CREATE OR REPLACE TABLE cB AS
    SELECT '{opt$sample}'      AS sample,
            rname, sj, {paste(feature_cols, collapse = ',')},
            {paste(c(mut_cols, base_cols), collapse = ',')}
            COUNT(*)            AS n
    FROM   merged
    GROUP  BY ALL;
  "))

}


#### Write to cB ####

if(opt$makecB){
  
  dir.create(dirname(opt$cBoutput), showWarnings = FALSE)
    
  # Write to csv
  dbExecute(con, glue("
    COPY cB TO '{opt$cBoutput}' (FORMAT CSV, HEADER);
  "))
  
    
}else{
  
  write.csv(tibble(),
            file = opt$cBoutput)
  
}


#### Write to arrow dataset ####

if(opt$makeArrow){
  
  dir.create(dirname(opt$Arrowoutput), showWarnings = FALSE)
  
  # 2b → Parquet
  dbExecute(con, glue("
    COPY cB TO '{opt$Arrowoutput}' (FORMAT PARQUET, COMPRESSION snappy);
  "))
  
  
}else{
  
  write.csv(tibble(),
                opt$Arrowoutput)
  
}


#### Write to cUP ####

if(opt$makecUP){
  
  group_cols        <- c("rname", feature_cols, mut_cols)   # columns to GROUP BY
  group_by_clause   <- paste(group_cols, collapse = ", ")
  
  avg_fragments     <- paste(sprintf("AVG(%s) AS %s", base_cols, base_cols),
                             collapse = ",\n           ")
  
  select_clause <- glue("
    '{opt$sample}' AS sample,
    {paste(group_cols, collapse = ', ')},
    {avg_fragments},
    COUNT(*) AS n
  ")
  
  # Create summarized cUP
  dbExecute(con, glue("
    CREATE OR REPLACE TABLE cUP AS
    SELECT {select_clause}
    FROM   merged
    GROUP  BY {group_by_clause};
  "))

  dbExecute(con, glue("
    COPY cUP TO '{opt$cUPoutput}' (FORMAT CSV, HEADER);
  "))
  
}else{
  
  write.csv(
    tibble(),
    file = opt$cUPoutput
  )
  
}


dbDisconnect(con, shutdown = TRUE)

