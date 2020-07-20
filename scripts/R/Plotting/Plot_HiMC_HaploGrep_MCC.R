library(tidyverse)

combined_summary_file = "/Volumes/TimMcInerney/MitoImpute/analyses/combined_summaries/MAF_combined.csv"

mcmc.dir = "MCMC_Experiments"
mcmc.var = c("MCMC1", "MCMC5", "MCMC10", "MCMC20", "MCMC30")
khap.dir = "kHAP_Experiments"
khap.var = c("kHAP100", "kHAP250", "kHAP500", "kHAP1000", "kHAP2500", "kHAP5000", "kHAP10000", "kHAP20000", "kHAP30000")
maf.dir = "MAF_Experiments"
maf.panel = c("ReferencePanel_v2", "ReferencePanel_v4", "ReferencePanel_v3")
maf.var = c("MAF1%", "MAF0.5%", "MAF0.1%")

combined_summary = read_csv(combined_summary_file)

combined_summary %>% 
  mutate(array = as.factor(array),
         refpan_maf = factor(refpan_maf,
                             levels = maf.var,
                             labels = maf.var,
                             ordered = T))

combined_summary_maf = combined_summary %>%
  select(c(1:17), contains("mean")) %>%
  mutate(himc_diff              = mean_himc_concordance_imputed              - mean_himc_concordance_typed,
         himc_cutoff_diff       = mean_himc_concordance_imputed_cutoff       - mean_himc_concordance_typed,
         himc_macro_diff        = mean_himc_concordance_imputed_macro        - mean_himc_concordance_typed_macro,
         himc_macro_cutoff_diff = mean_himc_concordance_imputed_macro_cutoff - mean_himc_concordance_typed_macro,
         # HAPLOGREP CONCORDANCE
         haplogrep_diff              = mean_haplogrep_concordance_imputed              - mean_haplogrep_concordance_typed,
         haplogrep_cutoff_diff       = mean_haplogrep_concordance_imputed_cutoff       - mean_haplogrep_concordance_typed,
         haplogrep_macro_diff        = mean_haplogrep_concordance_imputed_macro        - mean_haplogrep_concordance_typed_macro,
         haplogrep_macro_cutoff_diff = mean_haplogrep_concordance_imputed_macro_cutoff - mean_haplogrep_concordance_typed_macro,
         # HAPLOGREP QUALITY - IMPUTED | TYPED
         haplogrep_quality_diff              = mean_haplogrep_quality_imputed              - mean_haplogrep_quality_typed,
         haplogrep_quality_cutoff_diff       = mean_haplogrep_quality_imputed_cutoff       - mean_haplogrep_quality_typed,
         # HAPLOGREP QUALITY - TRUTH | TYPED + IMPUTED
         haplogrep_quality_diff_truth_typed                = mean_haplogrep_quality_truth       - mean_haplogrep_quality_typed,
         haplogrep_quality_diff_truth_imputed              = mean_haplogrep_quality_truth       - mean_haplogrep_quality_imputed,
         haplogrep_quality_diff_truth_imputed_cutoff       = mean_haplogrep_quality_truth       - mean_haplogrep_quality_imputed_cutoff,
         # HAPLOGREP DISTANCE DL
         haplogrep_distance_dl_diff              = mean_haplogrep_distance_dl_imputed              - mean_haplogrep_distance_dl_typed,
         haplogrep_distance_dl_cutoff_diff       = mean_haplogrep_distance_dl_imputed_cutoff       - mean_haplogrep_distance_dl_typed,
         # HAPLOGREP DISTANCE LV
         haplogrep_distance_lv_diff              = mean_haplogrep_distance_lv_imputed              - mean_haplogrep_distance_lv_typed,
         haplogrep_distance_lv_cutoff_diff       = mean_haplogrep_distance_lv_imputed_cutoff       - mean_haplogrep_distance_lv_typed,
         # HAPLOGREP DISTANCE JC
         haplogrep_distance_jc_diff              = mean_haplogrep_distance_jc_imputed              - mean_haplogrep_distance_jc_typed,
         haplogrep_distance_jc_cutoff_diff       = mean_haplogrep_distance_jc_imputed_cutoff       - mean_haplogrep_distance_jc_typed)

mean_summary = combined_summary_maf %>%
  select(refpan_maf, contains("diff")) %>%
  group_by(refpan_maf) %>%
  summarise_all(funs(mean), na.rm = TRUE)

View(mean_summary)

median_summary = combined_summary_maf %>%
  select(refpan_maf, contains("diff")) %>%
  group_by(refpan_maf) %>%
  summarise_all(funs(median), na.rm = TRUE)

View(median_summary)

combined_summary %>%
  arrange(n_snps_array) %>%
  ggplot(aes(x = array, y = mean_haplogrep_concordance_imputed)) +
  geom_boxplot(na.rm = T) +
  coord_flip()

ggplot(combined_summary, aes(x = refpan_maf, y = mean_haplogrep_concordance_typed)) +
  geom_boxplot(na.rm = T)

ggplot(combined_summary, aes(x = refpan_maf, y = mean_haplogrep_concordance_imputed)) +
  geom_boxplot(na.rm = T)

ggplot(combined_summary, aes(x = refpan_maf, y = mean_haplogrep_concordance_imputed_cutoff)) +
  geom_boxplot(na.rm = T)

ggplot(combined_summary, aes(x = refpan_maf, y = mean_haplogrep_concordance_imputed_macro)) +
  geom_boxplot(na.rm = T)

ggplot(combined_summary, aes(x = refpan_maf, y = mean_haplogrep_concordance_imputed_macro_cutoff)) +
  geom_boxplot(na.rm = T)

ggplot(combined_summary, aes(x = refpan_maf, y = mean_haplogrep_distance_dl_typed)) +
  geom_boxplot(na.rm = T)

ggplot(combined_summary, aes(x = refpan_maf, y = mean_haplogrep_distance_dl_imputed)) +
  geom_boxplot(na.rm = T)

ggplot(combined_summary, aes(x = refpan_maf, y = mean_haplogrep_distance_dl_imputed_cutoff) +
  geom_boxplot(na.rm = T)
  
ggplot(combined_summary, aes(x = refpan_maf, y = combined_summary$mean_haplogrep_quality_truth)) +
  geom_boxplot(na.rm = T)

ggplot(combined_summary, aes(x = refpan_maf, y = combined_summary$mean_haplogrep_quality_typed)) +
  geom_boxplot(na.rm = T)

ggplot(combined_summary, aes(x = refpan_maf, y = combined_summary$mean_haplogrep_quality_imputed)) +
  geom_boxplot(na.rm = T)

ggplot(combined_summary, aes(x = refpan_maf, y = combined_summary$mean_haplogrep_quality_imputed_cutoff)) +
  geom_boxplot(na.rm = T)


# HiMC DIFFERENCES
ggplot(combined_summary_maf, aes(x = refpan_maf, y = himc_diff)) +
  geom_boxplot(na.rm = T) +
  scale_y_continuous(limits = c(-1, 1))

ggplot(combined_summary_maf, aes(x = refpan_maf, y = himc_cutoff_diff)) +
  geom_boxplot(na.rm = T) +
  scale_y_continuous(limits = c(-1, 1))

ggplot(combined_summary_maf, aes(x = refpan_maf, y = himc_macro_diff)) +
  geom_boxplot(na.rm = T) +
  scale_y_continuous(limits = c(-1, 1))

ggplot(combined_summary_maf, aes(x = refpan_maf, y = himc_macro_cutoff_diff)) +
  geom_boxplot(na.rm = T) +
  scale_y_continuous(limits = c(-1, 1))

# HAPLOGREP DIFFERENCES
ggplot(combined_summary_maf, aes(x = refpan_maf, y = haplogrep_diff)) +
  geom_boxplot(na.rm = T) +
  scale_y_continuous(limits = c(-1, 1))

ggplot(combined_summary_maf, aes(x = refpan_maf, y = haplogrep_cutoff_diff)) +
  geom_boxplot(na.rm = T) +
  scale_y_continuous(limits = c(-1, 1))

ggplot(combined_summary_maf, aes(x = refpan_maf, y = haplogrep_macro_diff)) +
  geom_boxplot(na.rm = T) +
  scale_y_continuous(limits = c(-1, 1))

ggplot(combined_summary_maf, aes(x = refpan_maf, y = haplogrep_macro_cutoff_diff)) +
  geom_boxplot(na.rm = T) +
  scale_y_continuous(limits = c(-1, 1))

# HAPLOGREP QUALITY DIFFERENCES
ggplot(combined_summary_maf, aes(x = refpan_maf, y = haplogrep_quality_diff)) +
  geom_boxplot(na.rm = T) +
  scale_y_continuous(limits = c(-1, 1))

ggplot(combined_summary_maf, aes(x = refpan_maf, y = haplogrep_quality_cutoff_diff)) +
  geom_boxplot(na.rm = T) +
  scale_y_continuous(limits = c(-1, 1))

ggplot(combined_summary_maf, aes(x = refpan_maf, y = haplogrep_quality_diff_truth_typed)) +
  geom_boxplot(na.rm = T) +
  scale_y_continuous(limits = c(-1, 1))

ggplot(combined_summary_maf, aes(x = refpan_maf, y = haplogrep_quality_diff_truth_imputed)) +
  geom_boxplot(na.rm = T) +
  scale_y_continuous(limits = c(-1, 1))

ggplot(combined_summary_maf, aes(x = refpan_maf, y = haplogrep_quality_diff_truth_imputed_cutoff)) +
  geom_boxplot(na.rm = T) +
  scale_y_continuous(limits = c(-1, 1))


# HAPLOGREP DISTANCE DIFFERENCES DL
ggplot(combined_summary_maf, aes(x = refpan_maf, y = haplogrep_distance_dl_diff)) +
  geom_boxplot(na.rm = T) +
  scale_y_continuous(limits = c(-1, 1))

ggplot(combined_summary_maf, aes(x = refpan_maf, y = haplogrep_distance_dl_cutoff_diff)) +
  geom_boxplot(na.rm = T) +
  scale_y_continuous(limits = c(-1, 1))


# HAPLOGREP DISTANCE DIFFERENCES LV
ggplot(combined_summary_maf, aes(x = refpan_maf, y = haplogrep_distance_lv_diff)) +
  geom_boxplot(na.rm = T) +
  scale_y_continuous(limits = c(-1, 1))

ggplot(combined_summary_maf, aes(x = refpan_maf, y = haplogrep_distance_lv_cutoff_diff)) +
  geom_boxplot(na.rm = T) +
  scale_y_continuous(limits = c(-1, 1))


# HAPLOGREP DISTANCE DIFFERENCES JC
ggplot(combined_summary_maf, aes(x = refpan_maf, y = haplogrep_distance_jc_diff)) +
  geom_boxplot(na.rm = T) +
  scale_y_continuous(limits = c(-1, 1))

ggplot(combined_summary_maf, aes(x = refpan_maf, y = haplogrep_distance_jc_cutoff_diff)) +
  geom_boxplot(na.rm = T) +
  scale_y_continuous(limits = c(-1, 1))





# END!