library(tidyverse)

start.time = proc.time() # Start the timer!

### SPECIFY ARGUMENTS!
message("USAGE HINTS")
message("ARGUMENT 1:  IMPUTE2 INFO FILE")
message("ARGUMENT 1:  IMPUTE2 CUTOFF")

#impute_file = "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/BDCHP-1X10-HUMANHAP240S_11216501_A-b37/ReferencePanel_v1-unique_0.01/kHAP100_uniqueseqs/chrMT_1kg_BDCHP-1X10-HUMANHAP240S_11216501_A-b37_imputed_kHAP100"
#cutoff = 0.3

args = commandArgs(trailingOnly = TRUE) # Set arguments from the command line

impute_file = args[1] # WGS VCF file
cutoff      = args[2]

if (is.null(cutoff) || is.na(cutoff)) {
  cutoff = 0.3
  message(paste0("CUTOFF NOT SPECIFIED ... DEFAULTING TO: ", cutoff))
} else {
  cutoff = as.numeric(cutoff)
}

info_file = paste0(impute_file, "_info")

message("")
message(paste0("IMPUTE2 FILE:        ", impute_file))
message(paste0("IMPUTE2 INFO FILE:   ", info_file))
message(paste0("IMPUTE2 CUTOFF:      >= ", cutoff))
message("")

out_impute_retained_file = paste0(impute_file, "_cutoffRetained")
out_impute_discard_file  = paste0(impute_file, "_cutoffDiscarded")

out_info_retained_file = paste0(info_file, "_cutoffRetained")
out_info_discard_file  = paste0(info_file, "_cutoffDiscarded")

impute_df = read_delim(impute_file, delim = " ", col_names = F)
info_df   = read_delim(info_file, delim = " ", col_names = T)


message("")
info_df_retained = info_df %>%
  filter(info >= cutoff) %>%
  write_delim(path = out_info_retained_file, delim = " ")
message(paste0("RETAINED INFO FILE:       ", out_info_retained_file))

#info_df_discard  = info_df %>%
#  filter(info < cutoff) %>%
#  write_delim(path = out_info_discard_file, delim = " ")
#message(paste0("DISCARD INFO FILE:      >= ", out_info_retained_file))

impute_df %>%
  filter(impute_df$X3 %in% info_df_retained$position) %>%
  write_delim(path = out_impute_retained_file, delim = " ", col_names = F)
message(paste0("RETAINED IMPUTE2 FILE:      ", out_impute_retained_file))

#impute_df %>%
  #filter(impute_df$X3 %in% info_df_discard$position) %>%
  #write_delim(path = out_info_discard_file, delim = " ")
message("")



# END!
