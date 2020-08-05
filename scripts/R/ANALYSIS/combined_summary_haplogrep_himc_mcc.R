library(tidyverse)

args = commandArgs(trailingOnly = TRUE) # Set arguments from the command line
khap_run = args[1]
maf_run  = args[2]

if (is.na(khap_run) || is.null(khap_run)) {
  khap_run = T
}

if (is.na(maf_run) || is.null(maf_run)) {
  maf_run = T
}

if (khap_run == toupper("TRUE") || khap_run == toupper("T")) {
  khap_run = T
}

if (maf_run == toupper("TRUE") || maf_run == toupper("T")) {
  maf_run = T
}

message("")
message("")
if (khap_run == T) {
  message("COMBINING KHAP")
}
if (maf_run == T) {
  message("COMBINING MAF")
}
message("")
message("")

#stop("TEST!")

wd = "/g/data1a/te53/MitoImpute/data/STRANDS/"

out_file1 = "/g/data1a/te53/MitoImpute/analyses/combined_summaries/MAF_combined.csv"
out_file2 = "/g/data1a/te53/MitoImpute/analyses/combined_summaries/kHAP_combined.csv"

strands = read_tsv("~/GitCode/MitoImputePrep/scripts/INFORMATION_LISTS/b37_platforms.txt", col_names = F)

mcmc.dir = "MCMC_Experiments"
mcmc.var = c("MCMC1", "MCMC5", "MCMC10", "MCMC20", "MCMC30")
khap.dir = "kHAP_Experiments"
khap.var = c("kHAP100", "kHAP250", "kHAP500", "kHAP1000", "kHAP2500", "kHAP5000", "kHAP10000", "kHAP20000", "kHAP30000")
maf.dir = "MAF_Experiments"
maf.panel = c("ReferencePanel_v2", "ReferencePanel_v4", "ReferencePanel_v3")
maf.var = c("MAF1%", "MAF0.5%", "MAF0.1%")

# MAF EXPERIMENTS

if (maf_run == T) {
  for (i in 1:length(maf.panel)) {
    # SET REFERENCE PANEL
    ref_pan = maf.panel[i]
    maf     = maf.var[i]
    
    for (j in 1:nrow(strands)) {
      # SET STRAND
      strand = strands$X1[j]
      
      # ECHO INFO
      cat("WORKING ON ", ref_pan, " | ", strand, "\r")
      
      # SET SUMMARY FILE
      #summary_file = paste0(wd, strand, "/", ref_pan, "/MCMC1/chrMT_1kg_", strand, "_imputed_MCMC1_SUMMARY.csv")
      summary_file = paste0(wd, strand, "/", ref_pan, "/MCMC1/chrMT_1kg_", strand, "_imputed_MCMC1_SUMMARY.tsv")
      
      if (i == 1 && j == 1) {
        # READ IN SUMMARY FILE AS COMBINED SUMMARY IF FIRST IN LIST
        combined_summary = read_csv(summary_file, col_names = T)
      } else {
        # OTHERWISE CHECK TO SEE IF FILE EXISTS
        if (file.exists(summary_file)) {
          # IF IT FOES, READ IT IN, THEN COMBINED
          tmp_summary = read_csv(summary_file, col_names = T)
          combined_summary = bind_rows(combined_summary, tmp_summary)
        } else {
          # IF NOT, CREATE AN EMPTY DATA FRAME WITH THE ARRAY NAME AND IMPUTED AS FALSE
          tmp_empty_summary = as_tibble(data.frame(matrix(NA, nrow = 1, ncol = ncol(combined_summary))))
          names(tmp_empty_summary) = names(combined_summary)
          tmp_empty_summary$array[1] = strand
          tmp_empty_summary$imputed[1] = FALSE
          tmp_empty_summary$refpan_maf = maf
          tmp_empty_summary$mcmc = "MCMC1"
          tmp_empty_summary$k_hap = "kHAP500"
          combined_summary = bind_rows(combined_summary, tmp_empty_summary)
        }
        
      }
      
    }
    
  }
}

write_csv(x = combined_summary, path = out_file1)

if (file.exists(out_file1)) {
  message(paste0("OUTPUT FILE WRITTEN TO:  ", out_file1))
}

print(combined_summary)

# KHAP EXPERIMENTS

if (khap_run == T) {
  for (i in 1:length(khap.var)) {
    # SET REFERENCE PANEL
    khap =  khap.var[i]
    
    for (j in 1:nrow(strands)) {
      # SET STRAND
      strand = strands$X1[j]
      
      # ECHO INFO
      cat("WORKING ON ", khap, " | ", strand, "\r")
      
      # SET SUMMARY FILE
      #summary_file = paste0(wd, strand, "/kHAP_Experiments/", khap, "/chrMT_1kg_", strand, "_imputed_", khap, "_SUMMARY.csv")
      summary_file = paste0(wd, strand, "/kHAP_Experiments/", khap, "/chrMT_1kg_", strand, "_imputed_", khap, "_SUMMARY.tsv")
      
      if (i == 1 && j == 1) {
        # READ IN SUMMARY FILE AS COMBINED SUMMARY IF FIRST IN LIST
        combined_summary = read_csv(summary_file, col_names = T)
      } else {
        # OTHERWISE CHECK TO SEE IF FILE EXISTS
        if (file.exists(summary_file)) {
          # IF IT FOES, READ IT IN, THEN COMBINED
          tmp_summary = read_csv(summary_file, col_names = T)
          combined_summary = bind_rows(combined_summary, tmp_summary)
        } else {
          # IF NOT, CREATE AN EMPTY DATA FRAME WITH THE ARRAY NAME AND IMPUTED AS FALSE
          tmp_empty_summary = as_tibble(data.frame(matrix(NA, nrow = 1, ncol = ncol(combined_summary))))
          names(tmp_empty_summary) = names(combined_summary)
          tmp_empty_summary$array[1] = strand
          tmp_empty_summary$imputed[1] = FALSE
          tmp_empty_summary$refpan_maf = "MAF1%"
          tmp_empty_summary$mcmc = "MCMC1"
          tmp_empty_summary$k_hap = khap
          combined_summary = bind_rows(combined_summary, tmp_empty_summary)
        }
        
      }
      
    }
    
  }
}


write_csv(x = combined_summary, path = out_file2)

if (file.exists(out_file1)) {
  message(paste0("OUTPUT FILE WRITTEN TO:  ", out_file2))
}