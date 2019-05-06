#!/usr/bin/env Rscript

##########################################################################
##        R Script for 1) counting number of MT snps on each platform   ##
##        and 2) Writing out a list of platforms with MT SNPs           ##
##########################################################################

# Load Tidyverse
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(readr))


## Read in Arguments
# 1: Current working directory
args <- commandArgs(trailingOnly = T)
indir <- args[1]
if (length(args) == 1) {
  Nsnps_Mt_platforms <- sprintf("%s/Nsnps_Mt_platforms.txt", indir)
  Mt_platforms <- sprintf("%s/Mt_platforms.txt", indir)
} else {
  Nsnps_Mt_platforms <- args[2]
  Mt_platforms <- args[3]
}

# Read in file paths of individule platform strand files
fpath <- sprintf("%s/data/platforms", rwd)
b37 <- list.dirs(path = fpath, full.names = T, recursive = F)

## for each strand file
#   read in file and rename columns
#   if MT is present in CHROM column
#     count number of MT SNPs
#   else create empty tibble
#   write out
b37.ls <- function(x) {
  message(sprintf("IN PROGRESS: %s", x))
  fname <- sprintf("%s/platform.strand", x)
  df <- read_tsv(fname, col_names = F, col_types = "ccnncc") %>%
    rename(SNP = X1, chr = X2, pos = X3, match = X4, strand = X5)
  plat <- stringr::str_split_fixed(x, "platforms/", n = 2)[, 2]
  if ("MT" %in% df$chr) {
    out <- count(df, chr) %>%
      filter(chr == "MT") %>%
      mutate(platform = plat)
  } else {
    out <- tibble(chr = "MT", n = 0, platform = plat)
  }
  message("\t DONE \n")
  out
}

# Gather counts into tibble
platforms <- b37 %>%
  lapply(b37.ls) %>%
  bind_rows

# Write out tibble of number of MT snps on each platform
platforms %>%
  arrange(-n) %>%
  write_tsv(Nsnps_Mt_platforms)
message("Written:\t%s", Nsnps_Mt_platforms)

# Write out list of platforms containing > 1 MT SNPs
platforms %>%
  filter(n > 1) %>%
  pull(platform) %>%
  write(Mt_platforms, sep = "\n")
message("Written:\t%s", Mt_platforms)
