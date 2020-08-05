library(tidyverse)

printTime = function(end.time) {
  if (end.time[3] < 60) {
    out.time = end.time[3]
    return(paste0('Process completed in ',round(out.time,2),' seconds'))
  }
  if (end.time[3] >= 60 & end.time[3] < 3600) {
    out.time = end.time[3] / 60
    return(paste0('Process completed in ',round(out.time,2),' minutes'))
  }
  if (end.time[3] >= 3600 & end.time[3] < 86400) {
    out.time = (end.time[3] / 60) / 60
    return(paste0('Process completed in ',round(out.time,2),' hours'))
  }
  if (end.time[3] >= 86400) {
    out.time = ((end.time[3] / 60) / 60) / 24
    return(paste0('Process completed in ',round(out.time,2),' days'))
  }
}

start.time = proc.time() # Start the timer!

### SPECIFY ARGUMENTS!
message("USAGE HINTS")
message("ARGUMENT 1:  HiMC FILE")
message("ARGUMENT 2:  HAPLOGREP FILE")
message("ARGUMENT 3:  INFO FILE")
message("ARGUMENT 4:  INFO (CUTOFF) FILE")
message("ARGUMENT 5:  ARRAY SNPS FILE")
message("ARGUMENT 6:  MCC FILE")
message("ARGUMENT 7:  MCC (CUTOFF) FILE")


args = commandArgs(trailingOnly = TRUE) # Set arguments from the command line
strand           = args[1]
mcmc             = args[2]
maf              = args[3]
k_hap            = args[4]
himc_file        = args[5]
haplogrep_file   = args[6]
info_file        = args[7]
info_cutoff_file = args[8]
array_snps_file  = args[9]
mcc_file         = args[10]
mcc_cutoff_file  = args[11]


#himc_file        = "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/Human610-Quadv1_B-b37/ReferencePanel_v1-unique_0.01/kHAP100_uniqueseqs/chrMT_1kg_Human610-Quadv1_B-b37_imputed_kHAP100_HiMC_haplogroups.csv"
#haplogrep_file   = "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/Human610-Quadv1_B-b37/ReferencePanel_v1-unique_0.01/kHAP100_uniqueseqs/chrMT_1kg_Human610-Quadv1_B-b37_imputed_kHAP100_HaploGrep_haplogroups.csv"
#info_file        = "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/Human610-Quadv1_B-b37/ReferencePanel_v1-unique_0.01/kHAP100_uniqueseqs/chrMT_1kg_Human610-Quadv1_B-b37_imputed_kHAP100_info"
#info_cutoff_file = "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/Human610-Quadv1_B-b37/ReferencePanel_v1-unique_0.01/kHAP100_uniqueseqs/chrMT_1kg_Human610-Quadv1_B-b37_imputed_kHAP100_info_cutoffRetained"
#array_snps_file  = "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/Human610-Quadv1_B-b37/Human610-Quadv1_B-b37_MT_snps.txt"
#mcc_file         = "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/Human610-Quadv1_B-b37/ReferencePanel_v1-unique_0.01/kHAP100_uniqueseqs/chrMT_1kg_Human610-Quadv1_B-b37_imputed_kHAP100_imputed_MCC.csv"
#mcc_cutoff_file  = "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/Human610-Quadv1_B-b37/ReferencePanel_v1-unique_0.01/kHAP100_uniqueseqs/chrMT_1kg_Human610-Quadv1_B-b37_imputed_kHAP100_cutoffRetained_imputed_MCC.csv"

out_file = sub(pattern = "_HaploGrep_haplogroups.csv", replacement = "_SUMMARY.csv", x = haplogrep_file)

#strand  = "Human610-Quadv1_B-b37"
#mcmc    = "MCMC1"
#maf     = "MAF1%"
#k_hap   = "kHAP100"
imputed = TRUE
cutoff  = 0.3

message("")
message("INPUTS PARAMETERS")
message(paste0("INPUT ARRAY:          ", strand))
message(paste0("MCMC PARAMETER:       ", mcmc))
message(paste0("MAF PARAMETER:        ", maf))
message(paste0("KHAP PARAMETER:       ", k_hap))
message(paste0("INFO CUTOFF:           ≥ ", cutoff))
message(paste0("HiMC FILE:            ", himc_file))
message(paste0("HAPLOGREP FILE:       ", haplogrep_file))
message(paste0("INFO FILE:            ", info_file))
message(paste0("INFO (CUTOFF) FILE:   ", info_cutoff_file))
message(paste0("MCC FILE:             ", mcc_file))
message(paste0("MCC (CUTOFF) FILE:    ", mcc_cutoff_file))
message("")
message(paste0("OUTPUT FILE:          ", out_file))
message("")

himc_df        = read_csv(himc_file, col_names = T)
#haplogrep_df   = read_csv(haplogrep_file, col_names = T)
haplogrep_df   = read_tsv(haplogrep_file, col_names = T)
info_df        = read_delim(info_file, delim = " ", col_names = T)
info_cutoff_df = read_delim(info_cutoff_file, delim = " ", col_names = T)
array_snps     = read_tsv(array_snps_file, col_names = F)
mcc_df         = read_csv(mcc_file)
mcc_cutoff_df  = read_csv(mcc_cutoff_file)

summary_df = as_tibble(data.frame(array = strand, mcmc = mcmc, refpan_maf = maf, k_hap = k_hap, imputed = imputed, info_cutoff = cutoff)) %>%
  mutate(n_snps_array                                      = nrow(array_snps),
         n_snps_imputed                                    = nrow(info_df),
         n_snps_cutoff_imputed                             = nrow(info_cutoff_df),
         n_type_0                                          = nrow(info_df %>% filter(type == 0)),
         n_type_1                                          = nrow(info_df %>% filter(type == 1)),
         n_type_2                                          = nrow(info_df %>% filter(type == 2)),
         n_type_3                                          = nrow(info_df %>% filter(type == 3)),
         n_type_0_cutoff                                   = nrow(info_cutoff_df %>% filter(type == 0)),
         n_type_1_cutoff                                   = nrow(info_cutoff_df %>% filter(type == 1)),
         n_type_2_cutoff                                   = nrow(info_cutoff_df %>% filter(type == 2)),
         n_type_3_cutoff                                   = nrow(info_cutoff_df %>% filter(type == 3)),
         # INFO
         q1_info                                           = quantile(info_df$info, na.rm = T, probs = 0.25),
         mean_info                                         = mean(info_df$info, na.rm = T),
         sd_info                                           = sd(info_df$info, na.rm = T),
         median_info                                       = median(info_df$info, na.rm = T),
         q3_info                                           = quantile(info_df$info, na.rm = T, probs = 0.75),
         # INFO CUTOFF
         q1_info_cutoff                                    = quantile(info_cutoff_df$info, na.rm = T, probs = 0.25),
         mean_info_cutoff                                  = mean(info_cutoff_df$info, na.rm = T),
         sd_info_cutoff                                    = sd(info_cutoff_df$info, na.rm = T),
         median_info_cutoff                                = median(info_cutoff_df$info, na.rm = T),
         q3_info_cutoff                                    = quantile(info_cutoff_df$info, na.rm = T, probs = 0.75),
         # MAF
         q1_maf                                            = quantile(info_df$exp_freq_a1, na.rm = T, probs = 0.25),
         mean_maf                                          = mean(info_df$exp_freq_a1, na.rm = T),
         sd_maf                                            = sd(info_df$exp_freq_a1, na.rm = T),
         median_maf                                        = median(info_df$exp_freq_a1, na.rm = T),
         q3_maf                                            = quantile(info_df$exp_freq_a1, na.rm = T, probs = 0.75),
         # MAF CUTOFF
         q1_maf_cutoff                                     = quantile(info_cutoff_df$exp_freq_a1, na.rm = T, probs = 0.25),
         mean_maf_cutoff                                   = mean(info_cutoff_df$exp_freq_a1, na.rm = T),
         sd_maf_cutoff                                     = sd(info_cutoff_df$exp_freq_a1, na.rm = T),
         median_maf_cutoff                                 = median(info_cutoff_df$exp_freq_a1, na.rm = T),
         q3_maf_cutoff                                     = quantile(info_cutoff_df$exp_freq_a1, na.rm = T, probs = 0.75),
         # MCC
         q1_mcc                                            = quantile(mcc_df$mcc, na.rm = T, probs = 0.25),
         mean_mcc                                          = mean(mcc_df$mcc, na.rm = T),
         sd_mcc                                            = sd(mcc_df$mcc, na.rm = T),
         median_mcc                                        = median(mcc_df$mcc, na.rm = T),
         q3_mcc                                            = quantile(mcc_df$mcc, na.rm = T, probs = 0.75),
         # MCC CUTOFF
         q1_mcc_cutoff                                     = quantile(mcc_cutoff_df$mcc, na.rm = T, probs = 0.25),
         mean_mcc_cutoff                                   = mean(mcc_cutoff_df$mcc, na.rm = T),
         sd_mcc_cutoff                                     = sd(mcc_cutoff_df$mcc, na.rm = T),
         median_mcc_cutoff                                 = median(mcc_cutoff_df$mcc, na.rm = T),
         q3_mcc_cutoff                                     = quantile(mcc_cutoff_df$mcc, na.rm = T, probs = 0.75),
         # CONCORDANCE
         q1_concordance                                    = quantile(mcc_df$concodance, na.rm = T, probs = 0.25),
         mean_concordance                                  = mean(mcc_df$concodance, na.rm = T),
         sd_concordance                                    = sd(mcc_df$concodance, na.rm = T),
         median_concordance                                = median(mcc_df$concodance, na.rm = T),
         q3_concordance                                    = quantile(mcc_df$concodance, na.rm = T, probs = 0.75),
         # CONCORDANCE CUTOFF
         q1_concordance_cutoff                             = quantile(mcc_cutoff_df$concodance, na.rm = T, probs = 0.25),
         mean_concordance_cutoff                           = mean(mcc_cutoff_df$concodance, na.rm = T),
         sd_concordance_cutoff                             = sd(mcc_cutoff_df$concodance, na.rm = T),
         median_concordance_cutoff                         = median(mcc_cutoff_df$concodance, na.rm = T),
         q3_concordance_cutoff                             = quantile(mcc_cutoff_df$concodance, na.rm = T, probs = 0.75),
         # CERTAINTY
         q1_certainty                                      = quantile(mcc_df$certainty, na.rm = T, probs = 0.25),
         mean_certainty                                    = mean(mcc_df$certainty, na.rm = T),
         sd_certainty                                      = sd(mcc_df$certainty, na.rm = T),
         median_certainty                                  = median(mcc_df$certainty, na.rm = T),
         q3_certainty                                      = quantile(mcc_df$certainty, na.rm = T, probs = 0.75),
         # CERTAINTY CUTOFF
         q1_certainty_cutoff                               = quantile(mcc_cutoff_df$certainty, na.rm = T, probs = 0.25),
         mean_certainty_cutoff                             = mean(mcc_cutoff_df$certainty, na.rm = T),
         sd_certainty_cutoff                               = sd(mcc_cutoff_df$certainty, na.rm = T),
         median_certainty_cutoff                           = median(mcc_cutoff_df$certainty, na.rm = T),
         q3_certainty_cutoff                               = quantile(mcc_cutoff_df$certainty, na.rm = T, probs = 0.75),
         # HiMC CONCORDANCE TYPED
         q1_himc_concordance_typed                         = quantile(himc_df$HiMC_typed_match, na.rm = T, probs = 0.25),
         mean_himc_concordance_typed                       = mean(himc_df$HiMC_typed_match, na.rm = T),
         sd_himc_concordance_typed                         = sd(himc_df$HiMC_typed_match, na.rm = T),
         median_himc_concordance_typed                     = median(as.numeric(himc_df$HiMC_typed_match), na.rm = T),
         q3_himc_concordance_typed                       = quantile(himc_df$HiMC_typed_match, na.rm = T, probs = 0.75),
         # HiMC CONCORDANCE TYPED MACRO
         q1_himc_concordance_typed_macro                   = quantile(himc_df$HiMC_typed_macro_match, na.rm = T, probs = 0.25),
         mean_himc_concordance_typed_macro                 = mean(himc_df$HiMC_typed_macro_match, na.rm = T),
         sd_himc_concordance_typed_macro                   = sd(himc_df$HiMC_typed_macro_match, na.rm = T),
         median_himc_concordance_typed_macro               = median(as.numeric(himc_df$HiMC_typed_macro_match), na.rm = T),
         q3_himc_concordance_typed_macro                   = quantile(himc_df$HiMC_typed_macro_match, na.rm = T, probs = 0.75),
         # HiMC CONCORDANCE IMPUTED
         q1_himc_concordance_imputed                       = quantile(himc_df$HiMC_imputed_match, na.rm = T, probs = 0.25),
         mean_himc_concordance_imputed                     = mean(himc_df$HiMC_imputed_match, na.rm = T),
         sd_himc_concordance_imputed                       = sd(himc_df$HiMC_imputed_match, na.rm = T),
         median_himc_concordance_imputed                   = median(as.numeric(himc_df$HiMC_imputed_match), na.rm = T),
         q3_himc_concordance_imputed                       = quantile(himc_df$HiMC_imputed_match, na.rm = T, probs = 0.75),
         # HiMC CONCORDANCE IMPUTED CUTOFF
         q1_himc_concordance_imputed_cutoff                = quantile(himc_df$HiMC_imputed_cutoff_match, na.rm = T, probs = 0.25),
         mean_himc_concordance_imputed_cutoff              = mean(himc_df$HiMC_imputed_cutoff_match, na.rm = T),
         sd_himc_concordance_imputed_cutoff                = sd(himc_df$HiMC_imputed_cutoff_match, na.rm = T),
         median_himc_concordance_imputed_cutoff            = median(as.numeric(himc_df$HiMC_imputed_cutoff_match), na.rm = T),
         q3_himc_concordance_imputed_cutoff                = quantile(himc_df$HiMC_imputed_cutoff_match, na.rm = T, probs = 0.75),
         # HiMC CONCORDANCE IMPUTED MACRO
         q1_himc_concordance_imputed_macro                 = quantile(himc_df$HiMC_imputed_macro_match, na.rm = T, probs = 0.25),
         mean_himc_concordance_imputed_macro               = mean(himc_df$HiMC_imputed_macro_match, na.rm = T),
         sd_himc_concordance_imputed_macro                 = sd(himc_df$HiMC_imputed_macro_match, na.rm = T),
         median_himc_concordance_imputed_macro             = median(as.numeric(himc_df$HiMC_imputed_macro_match), na.rm = T),
         q3_himc_concordance_imputed_macro                 = quantile(himc_df$HiMC_imputed_macro_match, na.rm = T, probs = 0.75),
         # HiMC CONCORDANCE IMPUTED MACRO CUTOFF
         q1_himc_concordance_imputed_macro_cutoff          = quantile(himc_df$HiMC_imputed_cutoff_macro_match, na.rm = T, probs = 0.25),
         mean_himc_concordance_imputed_macro_cutoff        = mean(himc_df$HiMC_imputed_cutoff_macro_match, na.rm = T),
         sd_himc_concordance_imputed_macro_cutoff          = sd(himc_df$HiMC_imputed_cutoff_macro_match, na.rm = T),
         median_himc_concordance_imputed_macro_cutoff      = median(as.numeric(himc_df$HiMC_imputed_cutoff_macro_match), na.rm = T),
         q3_himc_concordance_imputed_macro_cutoff          = quantile(himc_df$HiMC_imputed_cutoff_macro_match, na.rm = T, probs = 0.75),
         # HAPLOGREP CONCORDANCE TYPED
         q1_haplogrep_concordance_typed                    = quantile(haplogrep_df$HaploGrep_typed_match, na.rm = T, probs = 0.25),
         mean_haplogrep_concordance_typed                  = mean(haplogrep_df$HaploGrep_typed_match, na.rm = T),
         sd_haplogrep_concordance_typed                    = sd(haplogrep_df$HaploGrep_typed_match, na.rm = T),
         median_haplogrep_concordance_typed                = median(as.numeric(haplogrep_df$HaploGrep_typed_match), na.rm = T),
         q3_haplogrep_concordance_typed                    = quantile(haplogrep_df$HaploGrep_typed_match, na.rm = T, probs = 0.75),
         # HAPLOGREP CONCORDANCE TYPED MACRO
         q1_haplogrep_concordance_typed_macro              = quantile(haplogrep_df$HaploGrep_typed_macro_match, na.rm = T, probs = 0.25),
         mean_haplogrep_concordance_typed_macro            = mean(haplogrep_df$HaploGrep_typed_macro_match, na.rm = T),
         sd_haplogrep_concordance_typed_macro              = sd(haplogrep_df$HaploGrep_typed_macro_match, na.rm = T),
         median_haplogrep_concordance_typed_macro          = median(as.numeric(haplogrep_df$HaploGrep_typed_macro_match), na.rm = T),
         q3_haplogrep_concordance_typed_macro              = quantile(haplogrep_df$HaploGrep_typed_macro_match, na.rm = T, probs = 0.75),
         # HAPLOGREP CONCORDANCE IMPUTED
         q1_haplogrep_concordance_imputed                  = quantile(haplogrep_df$HaploGrep_imputed_match, na.rm = T, probs = 0.25),
         mean_haplogrep_concordance_imputed                = mean(haplogrep_df$HaploGrep_imputed_match, na.rm = T),
         sd_haplogrep_concordance_imputed                  = sd(haplogrep_df$HaploGrep_imputed_match, na.rm = T),
         median_haplogrep_concordance_imputed              = median(as.numeric(haplogrep_df$HaploGrep_imputed_match), na.rm = T),
         q3_haplogrep_concordance_imputed                  = quantile(haplogrep_df$HaploGrep_imputed_match, na.rm = T, probs = 0.75),
         # HAPLOGREP CONCORDANCE IMPUTED CUTOFF
         q1_haplogrep_concordance_imputed_cutoff           = quantile(haplogrep_df$HaploGrep_imputed_cutoff_match, na.rm = T, probs = 0.25),
         mean_haplogrep_concordance_imputed_cutoff         = mean(haplogrep_df$HaploGrep_imputed_cutoff_match, na.rm = T),
         sd_haplogrep_concordance_imputed_cutoff           = sd(haplogrep_df$HaploGrep_imputed_cutoff_match, na.rm = T),
         median_haplogrep_concordance_imputed_cutoff       = median(as.numeric(haplogrep_df$HaploGrep_imputed_cutoff_match), na.rm = T),
         q3_haplogrep_concordance_imputed_cutoff           = quantile(haplogrep_df$HaploGrep_imputed_cutoff_match, na.rm = T, probs = 0.75),
         # HAPLOGREP CONCORDANCE IMPUTED MACRO
         q1_haplogrep_concordance_imputed_macro            = quantile(haplogrep_df$HaploGrep_imputed_macro_match, na.rm = T, probs = 0.25),
         mean_haplogrep_concordance_imputed_macro          = mean(haplogrep_df$HaploGrep_imputed_macro_match, na.rm = T),
         sd_haplogrep_concordance_imputed_macro            = sd(haplogrep_df$HaploGrep_imputed_macro_match, na.rm = T),
         median_haplogrep_concordance_imputed_macro        = median(as.numeric(haplogrep_df$HaploGrep_imputed_macro_match), na.rm = T),
         q3_haplogrep_concordance_imputed_macro            = quantile(haplogrep_df$HaploGrep_imputed_macro_match, na.rm = T, probs = 0.75),
         # HAPLOGREP CONCORDANCE IMPUTED MACRO CUTOFF
         q1_haplogrep_concordance_imputed_macro_cutoff     = quantile(haplogrep_df$HaploGrep_imputed_cutoff_macro_match, na.rm = T, probs = 0.25),
         mean_haplogrep_concordance_imputed_macro_cutoff   = mean(haplogrep_df$HaploGrep_imputed_cutoff_macro_match, na.rm = T),
         sd_haplogrep_concordance_imputed_macro_cutoff     = sd(haplogrep_df$HaploGrep_imputed_cutoff_macro_match, na.rm = T),
         median_haplogrep_concordance_imputed_macro_cutoff = median(as.numeric(haplogrep_df$HaploGrep_imputed_cutoff_macro_match), na.rm = T),
         q3_haplogrep_concordance_imputed_macro_cutoff     = quantile(haplogrep_df$HaploGrep_imputed_cutoff_macro_match, na.rm = T, probs = 0.75),
         # HAPLOGREP QUALITY TRUTH
         q1_haplogrep_quality_truth                        = quantile(haplogrep_df$HaploGrep_Quality_truth, na.rm = T, probs = 0.25),
         mean_haplogrep_quality_truth                      = mean(haplogrep_df$HaploGrep_Quality_truth, na.rm = T),
         sd_haplogrep_quality_truth                        = sd(haplogrep_df$HaploGrep_Quality_truth, na.rm = T),
         median_haplogrep_quality_truth                    = median(as.numeric(haplogrep_df$HaploGrep_Quality_truth), na.rm = T),
         q3_haplogrep_quality_truth                        = quantile(haplogrep_df$HaploGrep_Quality_truth, na.rm = T, probs = 0.75),
         # HAPLOGREP QUALITY TYPED
         q1_haplogrep_quality_typed                        = quantile(haplogrep_df$HaploGrep_Quality_typed, na.rm = T, probs = 0.25),
         mean_haplogrep_quality_typed                      = mean(haplogrep_df$HaploGrep_Quality_typed, na.rm = T),
         sd_haplogrep_quality_typed                        = sd(haplogrep_df$HaploGrep_Quality_typed, na.rm = T),
         median_haplogrep_quality_typed                    = median(as.numeric(haplogrep_df$HaploGrep_Quality_typed), na.rm = T),
         q3_haplogrep_quality_typed                        = quantile(haplogrep_df$HaploGrep_Quality_typed, na.rm = T, probs = 0.75),
         # HAPLOGREP QUALITY IMPUTED
         q1_haplogrep_quality_imputed                      = quantile(haplogrep_df$HaploGrep_Quality_imputed, na.rm = T, probs = 0.25),
         mean_haplogrep_quality_imputed                    = mean(haplogrep_df$HaploGrep_Quality_imputed, na.rm = T),
         sd_haplogrep_quality_imputed                      = sd(haplogrep_df$HaploGrep_Quality_imputed, na.rm = T),
         median_haplogrep_quality_imputed                  = median(as.numeric(haplogrep_df$HaploGrep_Quality_imputed), na.rm = T),
         q3_haplogrep_quality_imputed                      = quantile(haplogrep_df$HaploGrep_Quality_imputed, na.rm = T, probs = 0.75),
         # HAPLOGREP QUALITY IMPUTED CUTOFF
         q1_haplogrep_quality_imputed_cutoff               = quantile(haplogrep_df$HaploGrep_Quality_imputed_cutoff, na.rm = T, probs = 0.25),
         mean_haplogrep_quality_imputed_cutoff             = mean(haplogrep_df$HaploGrep_Quality_imputed_cutoff, na.rm = T),
         sd_haplogrep_quality_imputed_cutoff               = sd(haplogrep_df$HaploGrep_Quality_imputed_cutoff, na.rm = T),
         median_haplogrep_quality_imputed_cutoff           = median(as.numeric(haplogrep_df$HaploGrep_Quality_imputed_cutoff), na.rm = T),
         q3_haplogrep_quality_imputed_cutoff               = quantile(haplogrep_df$HaploGrep_Quality_imputed_cutoff, na.rm = T, probs = 0.75),
         # HAPLOGREP DISTANCE DL TYPED
         q1_haplogrep_distance_dl_typed                    = quantile(haplogrep_df$HaploGrep_typed_dl, na.rm = T, probs = 0.25),
         mean_haplogrep_distance_dl_typed                  = mean(haplogrep_df$HaploGrep_typed_dl, na.rm = T),
         sd_haplogrep_distance_dl_typed                    = sd(haplogrep_df$HaploGrep_typed_dl, na.rm = T),
         median_haplogrep_distance_dl_typed                = median(as.numeric(haplogrep_df$HaploGrep_typed_dl), na.rm = T),
         q3_haplogrep_distance_dl_typed                    = quantile(haplogrep_df$HaploGrep_typed_dl, na.rm = T, probs = 0.75),
         # HAPLOGREP DISTANCE DL IMPUTED
         q1_haplogrep_distance_dl_imputed                  = quantile(haplogrep_df$HaploGrep_imputed_dl, na.rm = T, probs = 0.25),
         mean_haplogrep_distance_dl_imputed                = mean(haplogrep_df$HaploGrep_imputed_dl, na.rm = T),
         sd_haplogrep_distance_dl_imputed                  = sd(haplogrep_df$HaploGrep_imputed_dl, na.rm = T),
         median_haplogrep_distance_dl_imputed              = median(as.numeric(haplogrep_df$HaploGrep_imputed_dl), na.rm = T),
         q3_haplogrep_distance_dl_imputed                  = quantile(haplogrep_df$HaploGrep_imputed_dl, na.rm = T, probs = 0.75),
         # HAPLOGREP DISTANCE DL IMPUTED CUTOFF
         q1_haplogrep_distance_dl_imputed_cutoff           = quantile(haplogrep_df$HaploGrep_imputed_cutoff_dl, na.rm = T, probs = 0.25),
         mean_haplogrep_distance_dl_imputed_cutoff         = mean(haplogrep_df$HaploGrep_imputed_cutoff_dl, na.rm = T),
         sd_haplogrep_distance_dl_imputed_cutoff           = sd(haplogrep_df$HaploGrep_imputed_cutoff_dl, na.rm = T),
         median_haplogrep_distance_dl_imputed_cutoff       = median(as.numeric(haplogrep_df$HaploGrep_imputed_cutoff_dl), na.rm = T),
         q3_haplogrep_distance_dl_imputed_cutoff           = quantile(haplogrep_df$HaploGrep_imputed_cutoff_dl, na.rm = T, probs = 0.75),
         # HAPLOGREP DISTANCE LV TYPED
         q1_haplogrep_distance_lv_typed                    = quantile(haplogrep_df$HaploGrep_typed_lv, na.rm = T, probs = 0.25),
         mean_haplogrep_distance_lv_typed                  = mean(haplogrep_df$HaploGrep_typed_lv, na.rm = T),
         sd_haplogrep_distance_lv_typed                    = sd(haplogrep_df$HaploGrep_typed_lv, na.rm = T),
         median_haplogrep_distance_lv_typed                = median(as.numeric(haplogrep_df$HaploGrep_typed_lv), na.rm = T),
         q3_haplogrep_distance_lv_typed                    = quantile(haplogrep_df$HaploGrep_typed_lv, na.rm = T, probs = 0.75),
         # HAPLOGREP DISTANCE LV IMPUTED
         q1_haplogrep_distance_lv_imputed                  = quantile(haplogrep_df$HaploGrep_imputed_lv, na.rm = T, probs = 0.25),
         mean_haplogrep_distance_lv_imputed                = mean(haplogrep_df$HaploGrep_imputed_lv, na.rm = T),
         sd_haplogrep_distance_lv_imputed                  = sd(haplogrep_df$HaploGrep_imputed_lv, na.rm = T),
         median_haplogrep_distance_lv_imputed              = median(as.numeric(haplogrep_df$HaploGrep_imputed_lv), na.rm = T),
         q3_haplogrep_distance_lv_imputed                  = quantile(haplogrep_df$HaploGrep_imputed_lv, na.rm = T, probs = 0.75),
         # HAPLOGREP DISTANCE LV IMPUTED CUTOFF
         q1_haplogrep_distance_lv_imputed_cutoff           = quantile(haplogrep_df$HaploGrep_imputed_cutoff_lv, na.rm = T, probs = 0.25),
         mean_haplogrep_distance_lv_imputed_cutoff         = mean(haplogrep_df$HaploGrep_imputed_cutoff_lv, na.rm = T),
         sd_haplogrep_distance_lv_imputed_cutoff           = sd(haplogrep_df$HaploGrep_imputed_cutoff_lv, na.rm = T),
         median_haplogrep_distance_lv_imputed_cutoff       = median(as.numeric(haplogrep_df$HaploGrep_imputed_cutoff_lv), na.rm = T),
         q3_haplogrep_distance_lv_imputed_cutoff           = quantile(haplogrep_df$HaploGrep_imputed_cutoff_lv, na.rm = T, probs = 0.75),
         # HAPLOGREP DISTANCE JC TYPED
         q1_haplogrep_distance_jc_typed                    = quantile(haplogrep_df$HaploGrep_typed_jc, na.rm = T, probs = 0.25),
         mean_haplogrep_distance_jc_typed                  = mean(haplogrep_df$HaploGrep_typed_jc, na.rm = T),
         sd_haplogrep_distance_jc_typed                    = sd(haplogrep_df$HaploGrep_typed_jc, na.rm = T),
         median_haplogrep_distance_jc_typed                = median(as.numeric(haplogrep_df$HaploGrep_typed_jc), na.rm = T),
         q3_haplogrep_distance_jc_typed                    = quantile(haplogrep_df$HaploGrep_typed_jc, na.rm = T, probs = 0.75),
         # HAPLOGREP DISTANCE JC IMPUTED
         q1_haplogrep_distance_jc_imputed                  = quantile(haplogrep_df$HaploGrep_imputed_jc, na.rm = T, probs = 0.25),
         mean_haplogrep_distance_jc_imputed                = mean(haplogrep_df$HaploGrep_imputed_jc, na.rm = T),
         sd_haplogrep_distance_jc_imputed                  = sd(haplogrep_df$HaploGrep_imputed_jc, na.rm = T),
         median_haplogrep_distance_jc_imputed              = median(as.numeric(haplogrep_df$HaploGrep_imputed_jc), na.rm = T),
         q3_haplogrep_distance_jc_imputed                  = quantile(haplogrep_df$HaploGrep_imputed_jc, na.rm = T, probs = 0.75),
         # HAPLOGREP DISTANCE JC IMPUTED
         q1_haplogrep_distance_jc_imputed_cutoff           = quantile(haplogrep_df$HaploGrep_imputed_cutoff_jc, na.rm = T, probs = 0.25),
         mean_haplogrep_distance_jc_imputed_cutoff         = mean(haplogrep_df$HaploGrep_imputed_cutoff_jc, na.rm = T),
         sd_haplogrep_distance_jc_imputed_cutoff           = sd(haplogrep_df$HaploGrep_imputed_cutoff_jc, na.rm = T),
         median_haplogrep_distance_jc_imputed_cutoff       = median(as.numeric(haplogrep_df$HaploGrep_imputed_cutoff_jc), na.rm = T),
         q3_haplogrep_distance_jc_imputed_cutoff           = quantile(haplogrep_df$HaploGrep_imputed_cutoff_jc, na.rm = T, probs = 0.75)
         ) %>%
  write_csv(path = out_file)

summary_df

message("")
message(printTime(proc.time() - start.time))
message("")

# END!

