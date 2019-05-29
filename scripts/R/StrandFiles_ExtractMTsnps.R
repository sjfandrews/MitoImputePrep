##########################################################################
##        R Script for extracting a list of MT SNPs on each platform    ##
##########################################################################

## Load Tidyverse
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))

## Arguments
# 1: Stand file
# 2: Output file
args <- commandArgs(trailingOnly = T)
strand.file <- args[1]
out.file <- args[2]

## Extract MT snps from strand file
#     Read in strand file (only chrom and pos)
#     Filter for MT in CHROM
#     Arrange tibble by POS
#     write out file
message(sprintf("IN PROGRESS: %s", strand.file))
read_tsv(strand.file, col_names = c("CHROM", "POS"), col_types = "-cn---") %>%
  filter(CHROM == "MT") %>%
  arrange(POS) %>%
  write_tsv(out.file, col_names = F)
