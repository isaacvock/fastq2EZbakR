#!/usr/bin/env Rscript
# R script for normalizing TimeLapse data
# July 17, 2017
# Oct 30, 2018 updated to prevent scale file with NaN/Inf entries.
# May 12, 2021 rewrite by Martin to allow stand-alone run
#              processing of multiple file locations at once (in case normalization between several TimeLapse experiments is needed)
#              option to calculate normalization factors for spikeins only
# Simon Lab
# Output a scale file for making tracks.

args = commandArgs(trailingOnly = TRUE)


# Read parameters
    library(optparse)

    option_list <- list(
        make_option(c("-d", "--dirs", type="character"),
                     default = ".",
                     help = "comma separated list of directories with mature.*.txt files"),
        make_option(c("-s", "--spikename", type="character"),
                     default = "",
                     help = 'character string that is common identifier of spikeins gene_name'),
        make_option(c("-o", "--output", type="character"),
                    default = "./results/normalization/scale",
                    help = 'File to output scale factors to'),
        make_option(c("-e", "--echocode", type="logical"),
                     default = "FALSE",
                     help = 'print R code to stdout'),
        make_option(c("--use_genes"),
                     action = "store_false",
                     default = TRUE,
                     dest= "exon",
                     help = 'If included, reads from exons and introns will be used for normalization'))

    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser) # Load options from command line.

# Set R code printing for debug mode
options(echo = as.logical(opt$echocode))


# Load libraries:
library(dplyr)
library(readr)
library(tidyr)
library(edgeR)
library(MASS)
library(data.table)


# Set date

date <- Sys.Date()
date

# Define the data frame

master <- tibble()



# Screw it, going to hardcode directory cause I can
if(opt$exon){
    dirs <- paste0(getwd(), "/results/featurecounts_exons/")
}else{
    dirs <- paste0(getwd(), "/results/featurecounts_genes/")
}

# Currently will not work alone due to fact that CORE files also have .featureCounts
samplefiles <- list.files(path = dirs,
                          pattern = "\\.featureCounts$",
                          recursive = FALSE)

print(dirs)
print(samplefiles)

# Remove the CORE files
samplefiles <- samplefiles[!grepl("\\.bam\\.featureCounts$", samplefiles)]

print(samplefiles)

# Create full directory path
samplenames <- paste0(dirs, samplefiles)

# Actual sample name wildcards
snames <- gsub(".featureCounts", "", samplefiles)


# Loop through the samples to load them:

for (i in seq_along(samplenames)){

    df <- fread(samplenames[i],
                sep = "\t")
    
    colnames(df) <- c("Geneid", "Chr", "Start", "End",
                      "Strand", "Length", "count")

    df$sample <- paste0(snames[i])

    master <- rbind(master, df)

}

rm(df)


m <- master %>%
    dplyr::select(-Chr, -Start, -End, -Strand, -Length) %>%
    pivot_wider(names_from = sample, values_from = count)

# Filter only for spikeins. gene_name (or other identifier must contain universal character string opt$spikename)
if (opt$spikename != "") {
    m <- m %>%
            filter(grepl(opt$spikename, Geneid))

    print(paste0("* Calculating normalization factors from spikeins with name *", opt$spikename))
}

# Turn into a count matrix accceptable for edgeR
rnames <- m[,1]
m <- as.matrix(m[,-1])
storage.mode(m) <- "numeric"
rownames(m) <- unlist(rnames)


### Normalize by edgeR:

# Add dummy IDs (don't need them yet)
ers <- rep("X", time = dim(m)[2])
ers <- factor(ers)

# Filter to only data-containing genes:
keep <- rowSums(m > 0) > 1
erde <- DGEList(m[keep,], group = ers)
erde <- calcNormFactors(erde, method = "upperquartile")
scale <- erde$samples$norm.factors*erde$samples$lib.size

if (sum(is.finite(scale)) == length(scale)){

  x <- mean(scale)/scale
  
  sdf <- tibble(sample = snames,
                scale = x)
  

  write.table(sdf, file = opt$output,
                   quote = FALSE,
                   row.names = FALSE,
                   col.names = FALSE)

} else {
  
  warning('!!no acceptable scale values!! Reverting to scale factor of 1')
  
  x <- rep(1, times = length(snames))
  
  sdf <- tibble(sample = snames,
                scale = x)
  

  write.table(sdf, file = opt$output,
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)
}

