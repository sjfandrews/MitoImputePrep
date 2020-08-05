library(tidyverse)
#ibrary(HiMC); data(nodes)
library(stringdist)

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
message("ARGUMENT 1:  HAPLOGREP FILE FOR FOR FULL/WGS DATA")
message("ARGUMENT 2:  HAPLOGREP FILE FOR FOR GENOTYPED DATA")
message("ARGUMENT 3:  HAPLOGREP FILE FOR FOR IMPUTED DATA")
message("ARGUMENT 4:  HAPLOGREP FILE FOR FOR IMPUTED (with info score cutoff) DATA")
message("ARGUMENT 5:  OUTPUT FILE")

#impute_file = "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/BDCHP-1X10-HUMANHAP240S_11216501_A-b37/ReferencePanel_v1-unique_0.01/kHAP100_uniqueseqs/chrMT_1kg_BDCHP-1X10-HUMANHAP240S_11216501_A-b37_imputed_kHAP100"
#cutoff = 0.3

args = commandArgs(trailingOnly = TRUE) # Set arguments from the command line

#full_1kGP_file           = "~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/1000genomes_mtDNA_haplogrep.txt"
#typed_1kGP_file          = "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/Human610-Quadv1_B-b37/ReferencePanel_v1-unique_0.01/chrMT_1kg_Human610-Quadv1_B-b37_diploid_haplogrep.txt"
#imputed_1kGP_file        = "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/Human610-Quadv1_B-b37/ReferencePanel_v1-unique_0.01/kHAP100_uniqueseqs/chrMT_1kg_Human610-Quadv1_B-b37_imputed_kHAP100_haplogrep.txt"
#imputed_1kGP_cutoff_file = "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/Human610-Quadv1_B-b37/ReferencePanel_v1-unique_0.01/kHAP100_uniqueseqs/chrMT_1kg_Human610-Quadv1_B-b37_imputed_kHAP100_cutoffRetained_haplogrep.txt"

full_1kGP_file           = args[1] # HAPLOGREP FILE FOR FOR FULL/WGS DATA
typed_1kGP_file          = args[2] # HAPLOGREP FILE FOR FOR GENOTYPED DATA
imputed_1kGP_file        = args[3] # HAPLOGREP FILE FOR FOR IMPUTED DATA
imputed_1kGP_cutoff_file = args[4] # HAPLOGREP FILE FOR FOR IMPUTED (with info score cutoff) DATA
out_file                 = args[5] # OUTPUT FILE

message("")
message("INPUTS PARAMETERS")
message(paste0("TRUTH HAPLOGREP FILE:               ", full_1kGP_file))
message(paste0("GENOTYPED HAPLOGREP FILE:           ", typed_1kGP_file))
message(paste0("IMPUTED HAPLOGREP FILE:             ", imputed_1kGP_file))
message(paste0("IMPUTED (CUTOFF) HAPLOGREP FILE:    ", imputed_1kGP_cutoff_file))

full_1kGP           = read_delim(full_1kGP_file, delim = "\t", col_names = T)
typed_1kGP          = read_delim(typed_1kGP_file, delim = "\t", col_names = T)
imputed_1kGP        = read_delim(imputed_1kGP_file, delim = "\t", col_names = T)
imputed_1kGP_cutoff = read_delim(imputed_1kGP_cutoff_file, delim = "\t", col_names = T)

if (is.na(out_file) || is.null(out_file)) {
  message("OUTPUT FILE NOT DETECTED")
  #out_file = sub(pattern = "_cutoffRetained_haplogrep.txt", replacement = "_HaploGrep_haplogroups.csv", x = imputed_1kGP_cutoff_file)
  out_file = sub(pattern = "_cutoffRetained_haplogrep.txt", replacement = "_HaploGrep_haplogroups.tsv", x = imputed_1kGP_cutoff_file)
  message(paste0("DEFAULTING TO:  ", out_file))
}

message("")
message(paste0("OUTPUT FILE:          ", out_file))
message("")

full_1kGP = full_1kGP %>%
  #mutate(SampleID = unlist(str_split(string = SampleID, "_"))[1]) %>%
  select(SampleID, Range, Haplogroup, Rank, Quality)

hg_exceptions = c("^L", "HV", "JT")

full_1kGP = full_1kGP %>%
  mutate(Macrohaplogroup = if_else(str_detect(Haplogroup, hg_exceptions),
                                   substr(Haplogroup, start = 1, stop = 2),
                                   substr(Haplogroup, start = 1, stop = 1))) %>%
  select(-c("Range"))
typed_1kGP = typed_1kGP %>%
  mutate(Macrohaplogroup = if_else(str_detect(Haplogroup, hg_exceptions),
                                   substr(Haplogroup, start = 1, stop = 2),
                                   substr(Haplogroup, start = 1, stop = 1))) %>%
  select(-c("Range"))
imputed_1kGP = imputed_1kGP %>%
  mutate(Macrohaplogroup = if_else(str_detect(Haplogroup, hg_exceptions),
                                   substr(Haplogroup, start = 1, stop = 2),
                                   substr(Haplogroup, start = 1, stop = 1))) %>%
  select(-c("Range"))
imputed_1kGP_cutoff = imputed_1kGP_cutoff %>%
  mutate(Macrohaplogroup = if_else(str_detect(Haplogroup, hg_exceptions),
                                   substr(Haplogroup, start = 1, stop = 2),
                                   substr(Haplogroup, start = 1, stop = 1))) %>%
  select(-c("Range"))

joined_table = full_1kGP %>%
  full_join(typed_1kGP, by = "SampleID", suffix = c("_truth", "_typed")) %>%
  full_join(imputed_1kGP, by = "SampleID") %>%
  full_join(imputed_1kGP_cutoff, by = "SampleID", suffix = c("_imputed", "_imputed_cutoff")) %>%
  select(-contains("full_path")) %>%
  rename_at(vars(starts_with("Haplogroup")),
            funs(str_replace(., "Haplogroup", "HaploGrep_Haplogroup"))) %>%
  rename_at(vars(starts_with("Macrohaplogroup")),
            funs(str_replace(., "Macrohaplogroup", "HaploGrep_Macrohaplogroup"))) %>%
  rename_at(vars(starts_with("Rank")),
            funs(str_replace(., "Rank", "HaploGrep_Rank"))) %>%
  rename_at(vars(starts_with("Quality")),
            funs(str_replace(., "Quality", "HaploGrep_Quality"))) %>%
  mutate(HaploGrep_typed_match                = HaploGrep_Haplogroup_truth == HaploGrep_Haplogroup_typed,
         HaploGrep_typed_macro_match          = HaploGrep_Macrohaplogroup_truth == HaploGrep_Macrohaplogroup_typed,
         HaploGrep_imputed_match              = HaploGrep_Haplogroup_truth == HaploGrep_Haplogroup_imputed,
         HaploGrep_imputed_macro_match        = HaploGrep_Macrohaplogroup_truth == HaploGrep_Macrohaplogroup_imputed,
         HaploGrep_imputed_cutoff_match       = HaploGrep_Haplogroup_truth == HaploGrep_Haplogroup_imputed_cutoff,
         HaploGrep_imputed_cutoff_macro_match = HaploGrep_Macrohaplogroup_truth == HaploGrep_Macrohaplogroup_imputed_cutoff) %>%
  mutate(HaploGrep_typed_dl = stringdist(HaploGrep_Haplogroup_truth, HaploGrep_Haplogroup_typed, method = "dl"),
         HaploGrep_typed_lv = stringdist(HaploGrep_Haplogroup_truth, HaploGrep_Haplogroup_typed, method = "lv"),
         HaploGrep_typed_jc = stringdist(HaploGrep_Haplogroup_truth, HaploGrep_Haplogroup_typed, method = "jaccard"),
         HaploGrep_imputed_dl = stringdist(HaploGrep_Haplogroup_truth, HaploGrep_Haplogroup_imputed, method = "dl"),
         HaploGrep_imputed_lv = stringdist(HaploGrep_Haplogroup_truth, HaploGrep_Haplogroup_imputed, method = "lv"),
         HaploGrep_imputed_jc = stringdist(HaploGrep_Haplogroup_truth, HaploGrep_Haplogroup_imputed, method = "jaccard"),
         HaploGrep_imputed_cutoff_dl = stringdist(HaploGrep_Haplogroup_truth, HaploGrep_Haplogroup_imputed_cutoff, method = "dl"),
         HaploGrep_imputed_cutoff_lv = stringdist(HaploGrep_Haplogroup_truth, HaploGrep_Haplogroup_imputed_cutoff, method = "lv"),
         HaploGrep_imputed_cutoff_jc = stringdist(HaploGrep_Haplogroup_truth, HaploGrep_Haplogroup_imputed_cutoff, method = "jaccard")) %>%
  #write_csv(path = out_file)
  write_tsv(path = out_file)
  
if (file.exists(out_file)) {
  message(paste0("FILE SUCCESSFULLY WRITTEN TO: ", out_file))
} else {
  message("FILE APPEARS TO HAVE FAILED TO HAVE WRITTEN")
}


joined_table 

message("")
message(printTime(proc.time() - start.time))
message("")
