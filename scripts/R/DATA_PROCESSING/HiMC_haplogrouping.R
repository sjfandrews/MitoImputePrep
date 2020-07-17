library(tidyverse)
library(HiMC); data(nodes)
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
message("ARGUMENT 1:  MAP/PED PREFIX FOR FULL/WGS DATA")
message("ARGUMENT 2:  MAP/PED PREFIX FOR GENOTYPED DATA")
message("ARGUMENT 3:  MAP/PED PREFIX FOR IMPUTED DATA")
message("ARGUMENT 4:  MAP/PED PREFIX FOR IMPUTED (with info score cutoff) DATA")
message("ARGUMENT 5:  OUTPUT FILE")

#impute_file = "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/BDCHP-1X10-HUMANHAP240S_11216501_A-b37/ReferencePanel_v1-unique_0.01/kHAP100_uniqueseqs/chrMT_1kg_BDCHP-1X10-HUMANHAP240S_11216501_A-b37_imputed_kHAP100"
#cutoff = 0.3

args = commandArgs(trailingOnly = TRUE) # Set arguments from the command line

full_1kGP_file_prefix           = args[1] # MAP/PED PREFIX FOR FULL/WGS DATA
typed_1kGP_file_prefix          = args[2] # MAP/PED PREFIX FOR GENOTYPED DATA
imputed_1kGP_file_prefix        = args[3] # MAP/PED PREFIX FOR IMPUTED DATA
imputed_1kGP_cutoff_file_prefix = args[4] # MAP/PED PREFIX FOR IMPUTED (with info score cutoff) DATA
out_file                        = args[5] # OUTPUT FILE

full_1kGP_map_file = paste0(full_1kGP_file_prefix, ".map")
full_1kGP_ped_file = paste0(full_1kGP_file_prefix, ".ped")

typed_1kGP_map_file = paste0(typed_1kGP_file_prefix, ".map")
typed_1kGP_ped_file = paste0(typed_1kGP_file_prefix, ".ped")

imputed_1kGP_map_file = paste0(imputed_1kGP_file_prefix, ".map")
imputed_1kGP_ped_file = paste0(imputed_1kGP_file_prefix, ".ped")

imputed_1kGP_map_cutoff_file = paste0(imputed_1kGP_cutoff_file_prefix, ".map")
imputed_1kGP_ped_cutoff_file = paste0(imputed_1kGP_cutoff_file_prefix, ".ped")

if (is.na(out_file) || is.null(out_file)) {
  message("OUTPUT FILE NOT DETECTED")
  out_file = sub(pattern = "_cutoffRetained.ped", replacement = "_HiMC_haplogroups.csv", x = imputed_1kGP_ped_cutoff_file)
  message(paste0("DEFAULTING TO:  ", out_file))
}

full_1kGP = generate_snp_data(full_1kGP_map_file,
                              full_1kGP_ped_file)
truth_table = as_tibble(HiMC::getClassifications(full_1kGP))
#truth.table = as_tibble(truth.table)

typed_1kGP = generate_snp_data(typed_1kGP_map_file,
                               typed_1kGP_ped_file)
typed_table = as_tibble(HiMC::getClassifications(typed_1kGP))

imputed_1kGP = generate_snp_data(imputed_1kGP_map_file,
                                 imputed_1kGP_ped_file)
imputed_table = as_tibble(HiMC::getClassifications(imputed_1kGP))

imputed_cutoff_1kGP = generate_snp_data(imputed_1kGP_map_cutoff_file,
                                        imputed_1kGP_ped_cutoff_file)
imputed_cutoff_table = as_tibble(HiMC::getClassifications(imputed_cutoff_1kGP))

hg_exceptions = c("^L", "HV", "JT")

truth_table = truth_table %>%
  mutate(haplogroup = na_if(haplogroup, "unclassified"),
         macrohaplogroup = if_else(str_detect(haplogroup, hg_exceptions),
                                   substr(haplogroup, start = 1, stop = 2),
                                   substr(haplogroup, start = 1, stop = 1)))
typed_table = typed_table %>%
  mutate(haplogroup = na_if(haplogroup, "unclassified"),
         macrohaplogroup = if_else(str_detect(haplogroup, hg_exceptions),
                                   substr(haplogroup, start = 1, stop = 2),
                                   substr(haplogroup, start = 1, stop = 1)))
imputed_table = imputed_table %>%
  mutate(haplogroup = na_if(haplogroup, "unclassified"),
         macrohaplogroup = if_else(str_detect(haplogroup, hg_exceptions),
                                   substr(haplogroup, start = 1, stop = 2),
                                   substr(haplogroup, start = 1, stop = 1)))
imputed_cutoff_table = imputed_cutoff_table %>%
  mutate(haplogroup = na_if(haplogroup, "unclassified"),
         macrohaplogroup = if_else(str_detect(haplogroup, hg_exceptions),
                                   substr(haplogroup, start = 1, stop = 2),
                                   substr(haplogroup, start = 1, stop = 1)))

joined_table = truth_table %>%
  full_join(typed_table, by = "Individual", suffix = c("_truth", "_typed")) %>%
  full_join(imputed_table, by = "Individual") %>%
  full_join(imputed_cutoff_table, by = "Individual", suffix = c("_imputed", "_imputed_cutoff")) %>%
  select(-contains("full_path")) %>%
  rename_at(vars(starts_with("haplogroup")),
            funs(str_replace(., "haplogroup", "HiMC_haplogroup"))) %>%
  rename_at(vars(starts_with("macrohaplogroup")),
            funs(str_replace(., "macrohaplogroup", "HiMC_macrohaplogroup"))) %>%
  mutate(HiMC_typed_match                = HiMC_haplogroup_truth == HiMC_haplogroup_typed,
         HiMC_typed_macro_match          = HiMC_macrohaplogroup_truth == HiMC_macrohaplogroup_typed,
         HiMC_imputed_match              = HiMC_haplogroup_truth == HiMC_haplogroup_imputed,
         HiMC_imputed_macro_match        = HiMC_macrohaplogroup_truth == HiMC_macrohaplogroup_imputed,
         HiMC_imputed_cutoff_match       = HiMC_haplogroup_truth == HiMC_haplogroup_imputed_cutoff,
         HiMC_imputed_cutoff_macro_match = HiMC_macrohaplogroup_truth == HiMC_macrohaplogroup_imputed_cutoff) %>%
  mutate(HiMC_typed_dl = stringdist(HiMC_haplogroup_truth, HiMC_haplogroup_typed, method = "dl"),
         HiMC_typed_lv = stringdist(HiMC_haplogroup_truth, HiMC_haplogroup_typed, method = "lv"),
         HiMC_typed_jc = stringdist(HiMC_haplogroup_truth, HiMC_haplogroup_typed, method = "jaccard"),
         HiMC_imputed_dl = stringdist(HiMC_haplogroup_truth, HiMC_haplogroup_imputed, method = "dl"),
         HiMC_imputed_lv = stringdist(HiMC_haplogroup_truth, HiMC_haplogroup_imputed, method = "lv"),
         HiMC_imputed_jc = stringdist(HiMC_haplogroup_truth, HiMC_haplogroup_imputed, method = "jaccard"),
         HiMC_imputed_cutoff_dl = stringdist(HiMC_haplogroup_truth, HiMC_haplogroup_imputed_cutoff, method = "dl"),
         HiMC_imputed_cutoff_lv = stringdist(HiMC_haplogroup_truth, HiMC_haplogroup_imputed_cutoff, method = "lv"),
         HiMC_imputed_cutoff_jc = stringdist(HiMC_haplogroup_truth, HiMC_haplogroup_imputed_cutoff, method = "jaccard")) %>%
  write_csv(path = out_file)

if (file.exists(out_file)) {
  message(paste0("FILE SUCCESSFULLY WRITTEN TO: ", out_file))
} else {
  message("FILE APPEARS TO HAVE FAILED TO HAVE WRITTEN")
}


joined_table 

message("")
message(printTime(proc.time() - start.time))
message("")
