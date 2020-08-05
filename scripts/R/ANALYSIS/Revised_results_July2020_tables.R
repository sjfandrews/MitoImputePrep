library(tidyverse)
library(RColorBrewer)
library(emmeans)

add_signf_code = function(df) {
  df = df %>% 
    mutate(signf = if_else(p.value < 0.001, "***",
                           if_else(p.value >= 0.001 & p.value < 0.01, "**",
                                   if_else(p.value >= 0.01 & p.value < 0.05, "*",
                                           if_else(p.value >= 0.05 & p.value < 0.1, ".", "")))))
  return(df)
}

confidence.interval = function(data, interval = 0.95, na.rm = T) {
  m = mean(data, na.rm = na.rm)
  s = sd(data, na.rm = na.rm)
  if (na.rm == T) {
    n = length(!is.na(data)[!is.na(data) == "TRUE"])
  } else {
    n = length(data)
  }
  e = qnorm(interval + ((1 - interval) / 2)) * s / sqrt(n)
  l = m - e
  u = m + e
  return(c("mean" = m, "lower" = l, "upper" = u, "df" = round(n - 1, 0)))
  # suggestion: use t-value for alpha
}


pal_cols = brewer.pal(12, "Dark2")

mcmc.dir = "MCMC_Experiments"
mcmc.var = c("MCMC1", "MCMC5", "MCMC10", "MCMC20", "MCMC30")
khap.dir = "kHAP_Experiments"
khap.var = c("kHAP100", "kHAP250", "kHAP500", "kHAP1000", "kHAP2500", "kHAP5000", "kHAP10000", "kHAP20000", "kHAP30000")
khap.lab = c("khap_100", "khap_250", "khap_500", "khap_1000", "khap_2500", "khap_5000", "khap_10000", "khap_20000", "khap_30000")
maf.dir = "MAF_Experiments"
maf.panel = c("ReferencePanel_v2", "ReferencePanel_v4", "ReferencePanel_v3")
maf.var = c("MAF1%", "MAF0.5%", "MAF0.1%")
maf.lab = c("MAF>1%", "MAF>0.5%", "MAF>0.1%")

dir("/Volumes/TimMcInerney/MitoImpute/analyses/combined_summaries/")

maf_combined_summary_file  = "/Volumes/TimMcInerney/MitoImpute/analyses/combined_summaries/MAF_combined.csv"
khap_combined_summary_file = "/Volumes/TimMcInerney/MitoImpute/analyses/combined_summaries/kHAP_combined.csv"
maf_info_combined_file     = "/Volumes/TimMcInerney/MitoImpute/analyses/combined_summaries/MAF_MCC_concatenated.csv"
khap_info_combined_file    = "/Volumes/TimMcInerney/MitoImpute/analyses/combined_summaries/kHAP_MCC_concatenated.csv"
khap_haplogrep_file        = "/Volumes/TimMcInerney/MitoImpute/analyses/combined_summaries/kHAP_HaploGrep_haplogroups_concatenated.csv"
khap_himc_file             = "/Volumes/TimMcInerney/MitoImpute/analyses/combined_summaries/kHAP_HiMC_haplogroups_concatenated.csv"

out_prefix = "~/GitCode/MitoImputePrep/supplementary_information/R_tables/"

maf_combined_summary  = read_csv(maf_combined_summary_file)
khap_combined_summary = read_csv(khap_combined_summary_file)

maf_info_combined     = read_csv(maf_info_combined_file)
khap_info_combined    = read_csv(khap_info_combined_file)

khap_haplogrep         = read_csv(khap_haplogrep_file)
khap_himc              = read_csv(khap_himc_file)

maf_info_combined = maf_info_combined %>%
  mutate(array = as.factor(array),
         refpan_maf = factor(refpan_maf,
                             levels = maf.var,
                             labels = maf.lab,
                             ordered = T),
         k_hap      = factor(k_hap,
                             ordered = T))

khap_info_combined = khap_info_combined %>%
  mutate(array = as.factor(array),
         refpan_maf = factor(refpan_maf,
                             levels = maf.var,
                             labels = maf.lab,
                             ordered = T),
         k_hap = factor(k_hap,
                        levels = khap.var,
                        labels = khap.lab,
                        ordered = T))

khap_haplogrep = khap_haplogrep %>%
  mutate(array = as.factor(array),
         refpan_maf = factor(refpan_maf,
                             levels = maf.var,
                             labels = maf.lab,
                             ordered = T),
         k_hap = factor(k_hap,
                        levels = khap.var,
                        labels = khap.lab,
                        ordered = T)) %>%
  mutate_at(vars(contains("Haplogroup")), .funs = as.factor)

khap_himc = khap_himc %>%
  mutate(array = as.factor(array),
         refpan_maf = factor(refpan_maf,
                             levels = maf.var,
                             labels = maf.lab,
                             ordered = T),
         k_hap = factor(k_hap,
                        levels = khap.var,
                        labels = khap.lab,
                        ordered = T)) %>%
  mutate_at(vars(contains("Haplogroup")), .funs = as.factor)

maf_info_combined_summary = maf_info_combined %>%
  filter(type == 0) %>%
  group_by(array, refpan_maf, k_hap) %>%
  summarise(mean_info = mean(info, na.rm = T),
            mean_mcc  = mean(mcc, na.rm = T))

khap_info_combined_summary = khap_info_combined %>%
  filter(type == 0) %>%
  group_by(array, refpan_maf, k_hap) %>%
  summarise(mean_info = mean(info, na.rm = T),
            mean_mcc  = mean(mcc, na.rm = T))

khap_info_combined_recommendedOnly = khap_info_combined %>%
  filter(refpan_maf == "MAF>1%", k_hap == "khap_500")

khap_haplogrep_recommendedOnly = khap_haplogrep %>%
  filter(refpan_maf == "MAF>1%", k_hap == "khap_500")

khap_himc_recommendedOnly = khap_himc %>%
  filter(refpan_maf == "MAF>1%", k_hap == "khap_500")

maf_combined_summary = maf_combined_summary %>% 
  mutate(array = as.factor(array),
         refpan_maf = factor(refpan_maf,
                             levels = maf.var,
                             labels = maf.lab,
                             ordered = T),
         k_hap      = factor(k_hap,
                             ordered = T))

maf_info_combined_summary
khap_info_combined_summary
maf_combined_summary
khap_combined_summary


maf_combined_summary_means = maf_combined_summary %>%
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
glimpse(maf_combined_summary_means)

khap_combined_summary = khap_combined_summary %>% 
  mutate(array = as.factor(array),
         refpan_maf = factor(refpan_maf,
                             levels = maf.var,
                             labels = maf.lab,
                             ordered = T),
         k_hap = factor(k_hap,
                        levels = khap.var,
                        labels = khap.lab,
                        ordered = T))

khap_combined_summary_means = khap_combined_summary %>%
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
glimpse(khap_combined_summary_means)

khap_haplogrep_recommendedOnly_long = khap_haplogrep_recommendedOnly %>%
  gather(key = "key", value = "val", c(HaploGrep_imputed_cutoff_macro_match, HaploGrep_typed_macro_match)) %>%
  select(1:7, HaploGrep_Macrohaplogroup_truth, key, val) 

haplogrep_recommended_settings_out = paste0(out_prefix,     "haplogrep_recommended_settings.csv")
haplogrep_perhaplogroup = khap_haplogrep_recommendedOnly %>%
  group_by(HaploGrep_Macrohaplogroup_truth) %>%
  summarise(prop_match_imputed = mean(HaploGrep_imputed_cutoff_macro_match, na.rm = T),
            prop_match_typed = mean(HaploGrep_typed_macro_match, na.rm = T)) %>%
  mutate(diff = prop_match_imputed - prop_match_typed) %>%
  arrange(desc(diff)) %>%
  print(n = Inf) %>%
  write_csv(path = haplogrep_recommended_settings_out)

haplogrep_perhaplogroup %>%
  arrange(prop_match_typed) %>%
  print(n=Inf)

khap_himc_recommendedOnly_long = khap_himc_recommendedOnly %>%
  gather(key = "key", value = "val", c(HiMC_imputed_cutoff_macro_match, HiMC_typed_macro_match)) %>%
  select(1:7, HiMC_macrohaplogroup_truth, key, val) 

haplogrep_recommended_settings_out = paste0(out_prefix,     "haplogrep_recommended_settings.csv")

himc_recommended_settings_out = paste0(out_prefix,     "himc_recommended_settings.csv")
himc_perhaplogroup = khap_himc_recommendedOnly %>%
  group_by(HiMC_macrohaplogroup_truth) %>%
  summarise(prop_match_imputed = mean(HiMC_imputed_cutoff_macro_match, na.rm = T),
            prop_match_typed = mean(HiMC_typed_macro_match, na.rm = T)) %>%
  mutate(diff = prop_match_imputed - prop_match_typed) %>%
  arrange(desc(diff)) %>%
  print(n = Inf) %>%
  write_csv(path = himc_recommended_settings_out)

himc_perhaplogroup %>%
  arrange(prop_match_typed) %>%
  print(n=Inf)

maf_info_combined_summary
khap_info_combined_summary
maf_combined_summary_means
khap_combined_summary_means


# Matthew's Correlation Coefficient
maf_mcc_table_anova_out = paste0(out_prefix,     "maf_mcc_anova.csv")
maf_mcc_table_emmeans_out = paste0(out_prefix,   "maf_mcc_emmeans.csv")
maf_mcc_table_contrasts_out = paste0(out_prefix, "maf_mcc_contrasts.csv")

maf_mcc_cutoff_lm = lm(mean_mcc_cutoff ~ refpan_maf, data = maf_combined_summary_means)
maf_mcc_cutoff_lm_sum = summary(emmeans(maf_mcc_cutoff_lm, pairwise ~ refpan_maf))

write_csv(anova(maf_mcc_cutoff_lm),                        maf_mcc_table_anova_out)
write_csv(maf_mcc_cutoff_lm_sum$emmeans,                   maf_mcc_table_emmeans_out)
write_csv(add_signf_code(maf_mcc_cutoff_lm_sum$contrasts), maf_mcc_table_contrasts_out)

khap_mcc_table_anova_out = paste0(out_prefix,     "khap_mcc_anova.csv")
khap_mcc_table_emmeans_out = paste0(out_prefix,   "khap_mcc_emmeans.csv")
khap_mcc_table_contrasts_out = paste0(out_prefix, "khap_mcc_contrasts.csv")

khap_mcc_cutoff_lm = lm(mean_mcc_cutoff ~ k_hap, data = khap_combined_summary_means)
khap_mcc_cutoff_lm_sum = summary(emmeans(khap_mcc_cutoff_lm, pairwise ~ k_hap))

write_csv(anova(khap_mcc_cutoff_lm),                        khap_mcc_table_anova_out)
write_csv(khap_mcc_cutoff_lm_sum$emmeans,                   khap_mcc_table_emmeans_out)
write_csv(add_signf_code(khap_mcc_cutoff_lm_sum$contrasts), khap_mcc_table_contrasts_out)

# Matthew's Correlation Coefficient (TYPE 0 SNVs ONLY ie IMPUTED SNVs ONLY)
maf_mcc_t0_table_anova_out = paste0(out_prefix,     "maf_mcc_t0_anova.csv")
maf_mcc_t0_table_emmeans_out = paste0(out_prefix,   "maf_mcc_t0_emmeans.csv")
maf_mcc_t0_table_contrasts_out = paste0(out_prefix, "maf_mcc_t0_contrasts.csv")

maf_mcc_t0_cutoff_lm = lm(mean_mcc ~ refpan_maf, data = maf_info_combined_summary)
maf_mcc_t0_cutoff_lm_sum = summary(emmeans(maf_mcc_t0_cutoff_lm, pairwise ~ refpan_maf))

write_csv(anova(maf_mcc_t0_cutoff_lm),                        maf_mcc_t0_table_anova_out)
write_csv(maf_mcc_t0_cutoff_lm_sum$emmeans,                   maf_mcc_t0_table_emmeans_out)
write_csv(add_signf_code(maf_mcc_t0_cutoff_lm_sum$contrasts), maf_mcc_t0_table_contrasts_out)

khap_mcc_t0_table_anova_out = paste0(out_prefix,     "khap_info_t0_anova.csv")
khap_mcc_t0_table_emmeans_out = paste0(out_prefix,   "khap_info_t0_emmeans.csv")
khap_mcc_t0_table_contrasts_out = paste0(out_prefix, "khap_info_t0_contrasts.csv")

khap_mcc_t0_cutoff_lm = lm(mean_mcc ~ k_hap, data = khap_info_combined_summary)
khap_mcc_t0_cutoff_lm_sum = summary(emmeans(khap_mcc_t0_cutoff_lm, pairwise ~ k_hap))

write_csv(anova(khap_mcc_t0_cutoff_lm),                        khap_mcc_t0_table_anova_out)
write_csv(khap_mcc_t0_cutoff_lm_sum$emmeans,                   khap_mcc_t0_table_emmeans_out)
write_csv(add_signf_code(khap_mcc_t0_cutoff_lm_sum$contrasts), khap_mcc_t0_table_contrasts_out)


# IMPUTE2 INFO SCORE (TYPE 0 SNVs ONLY ie IMPUTED SNVs ONLY)
maf_info_t0_table_anova_out = paste0(out_prefix,     "maf_info_t0_anova.csv")
maf_info_t0_table_emmeans_out = paste0(out_prefix,   "maf_info_t0_emmeans.csv")
maf_info_t0_table_contrasts_out = paste0(out_prefix, "maf_info_t0_contrasts.csv")

maf_info_t0_cutoff_lm = lm(mean_info ~ refpan_maf, data = maf_info_combined_summary)
maf_info_t0_cutoff_lm_sum = summary(emmeans(maf_info_t0_cutoff_lm, pairwise ~ refpan_maf))

write_csv(anova(maf_info_t0_cutoff_lm),                        maf_info_t0_table_anova_out)
write_csv(maf_info_t0_cutoff_lm_sum$emmeans,                   maf_info_t0_table_emmeans_out)
write_csv(add_signf_code(maf_info_t0_cutoff_lm_sum$contrasts), maf_info_t0_table_contrasts_out)

khap_info_t0_table_anova_out = paste0(out_prefix,     "khap_info_t0_anova.csv")
khap_info_t0_table_emmeans_out = paste0(out_prefix,   "khap_info_t0_emmeans.csv")
khap_info_t0_table_contrasts_out = paste0(out_prefix, "khap_info_t0_contrasts.csv")

khap_info_t0_cutoff_lm = lm(mean_info ~ k_hap, data = khap_info_combined_summary)
khap_info_t0_cutoff_lm_sum = summary(emmeans(khap_info_t0_cutoff_lm, pairwise ~ k_hap))

write_csv(anova(khap_info_t0_cutoff_lm),                        khap_info_t0_table_anova_out)
write_csv(khap_info_t0_cutoff_lm_sum$emmeans,                   khap_info_t0_table_emmeans_out)
write_csv(add_signf_code(khap_info_t0_cutoff_lm_sum$contrasts), khap_info_t0_table_contrasts_out)


# HaploGrep2
maf_haplogrep_table_anova_out = paste0(out_prefix,     "maf_haplogrep_anova.csv")
maf_haplogrep_table_emmeans_out = paste0(out_prefix,   "maf_haplogrep_emmeans.csv")
maf_haplogrep_table_contrasts_out = paste0(out_prefix, "maf_haplogrep_contrasts.csv")

maf_haplogrep_cutoff_lm = lm(haplogrep_cutoff_diff ~ refpan_maf, data = maf_combined_summary_means)
maf_haplogrep_cutoff_lm_sum = summary(emmeans(maf_haplogrep_cutoff_lm, pairwise ~ refpan_maf))

write_csv(anova(maf_haplogrep_cutoff_lm),                        maf_haplogrep_table_anova_out)
write_csv(maf_haplogrep_cutoff_lm_sum$emmeans,                   maf_haplogrep_table_emmeans_out)
write_csv(add_signf_code(maf_haplogrep_cutoff_lm_sum$contrasts), maf_haplogrep_table_contrasts_out)

khap_haplogrep_table_anova_out = paste0(out_prefix,     "khap_haplogrep_anova.csv")
khap_haplogrep_table_emmeans_out = paste0(out_prefix,   "khap_haplogrep_emmeans.csv")
khap_haplogrep_table_contrasts_out = paste0(out_prefix, "khap_haplogrep_contrasts.csv")

khap_haplogrep_cutoff_lm = lm(haplogrep_cutoff_diff ~ k_hap, data = khap_combined_summary_means)
khap_haplogrep_cutoff_lm_sum = summary(emmeans(khap_haplogrep_cutoff_lm, pairwise ~ k_hap))

write_csv(anova(khap_haplogrep_cutoff_lm),                        khap_haplogrep_table_anova_out)
write_csv(khap_haplogrep_cutoff_lm_sum$emmeans,                   khap_haplogrep_table_emmeans_out)
write_csv(add_signf_code(khap_haplogrep_cutoff_lm_sum$contrasts), khap_haplogrep_table_contrasts_out)

# HaploGrep2 - MACRO
maf_macrohaplogrep_table_anova_out = paste0(out_prefix,     "maf_macrohaplogrep_anova.csv")
maf_macrohaplogrep_table_emmeans_out = paste0(out_prefix,   "maf_macrohaplogrep_emmeans.csv")
maf_macrohaplogrep_table_contrasts_out = paste0(out_prefix, "maf_macrohaplogrep_contrasts.csv")

maf_macrohaplogrep_cutoff_lm = lm(haplogrep_macro_cutoff_diff ~ refpan_maf, data = maf_combined_summary_means)
maf_macrohaplogrep_cutoff_lm_sum = summary(emmeans(maf_macrohaplogrep_cutoff_lm, pairwise ~ refpan_maf))

write_csv(anova(maf_macrohaplogrep_cutoff_lm),                        maf_macrohaplogrep_table_anova_out)
write_csv(maf_macrohaplogrep_cutoff_lm_sum$emmeans,                   maf_macrohaplogrep_table_emmeans_out)
write_csv(add_signf_code(maf_macrohaplogrep_cutoff_lm_sum$contrasts), maf_macrohaplogrep_table_contrasts_out)

khap_macrohaplogrep_table_anova_out = paste0(out_prefix,     "khap_macrohaplogrep_anova.csv")
khap_macrohaplogrep_table_emmeans_out = paste0(out_prefix,   "khap_macrohaplogrep_emmeans.csv")
khap_macrohaplogrep_table_contrasts_out = paste0(out_prefix, "khap_macrohaplogrep_contrasts.csv")

khap_macrohaplogrep_cutoff_lm = lm(haplogrep_macro_cutoff_diff ~ k_hap, data = khap_combined_summary_means)
khap_macrohaplogrep_cutoff_lm_sum = summary(emmeans(khap_macrohaplogrep_cutoff_lm, pairwise ~ k_hap))

write_csv(anova(khap_macrohaplogrep_cutoff_lm),                        khap_macrohaplogrep_table_anova_out)
write_csv(khap_macrohaplogrep_cutoff_lm_sum$emmeans,                   khap_macrohaplogrep_table_emmeans_out)
write_csv(add_signf_code(khap_macrohaplogrep_cutoff_lm_sum$contrasts), khap_macrohaplogrep_table_contrasts_out)

# HaploGrep2 - Quality Scores
maf_haplogrep_quality_table_anova_out = paste0(out_prefix,     "maf_haplogrep_quality_anova.csv")
maf_haplogrep_quality_table_emmeans_out = paste0(out_prefix,   "maf_haplogrep_quality_emmeans.csv")
maf_haplogrep_quality_table_contrasts_out = paste0(out_prefix, "maf_haplogrep_quality_contrasts.csv")

maf_haplogrep_quality_cutoff_lm = lm(haplogrep_quality_cutoff_diff ~ refpan_maf, data = maf_combined_summary_means)
maf_haplogrep_quality_cutoff_lm_sum = summary(emmeans(maf_haplogrep_quality_cutoff_lm, pairwise ~ refpan_maf))

write_csv(anova(maf_haplogrep_quality_cutoff_lm),                        maf_haplogrep_quality_table_anova_out)
write_csv(maf_haplogrep_quality_cutoff_lm_sum$emmeans,                   maf_haplogrep_quality_table_emmeans_out)
write_csv(add_signf_code(maf_haplogrep_quality_cutoff_lm_sum$contrasts), maf_haplogrep_quality_table_contrasts_out)

khap_haplogrep_quality_table_anova_out = paste0(out_prefix,     "khap_haplogrep_quality_anova.csv")
khap_haplogrep_quality_table_emmeans_out = paste0(out_prefix,   "khap_haplogrep_quality_emmeans.csv")
khap_haplogrep_quality_table_contrasts_out = paste0(out_prefix, "khap_haplogrep_quality_contrasts.csv")

khap_haplogrep_quality_cutoff_lm = lm(haplogrep_quality_cutoff_diff ~ k_hap, data = khap_combined_summary_means)
khap_haplogrep_quality_cutoff_lm_sum = summary(emmeans(khap_haplogrep_quality_cutoff_lm, pairwise ~ k_hap))

write_csv(anova(khap_haplogrep_quality_cutoff_lm),                        khap_haplogrep_quality_table_anova_out)
write_csv(khap_haplogrep_quality_cutoff_lm_sum$emmeans,                   khap_haplogrep_quality_table_emmeans_out)
write_csv(add_signf_code(khap_haplogrep_quality_cutoff_lm_sum$contrasts), khap_haplogrep_quality_table_contrasts_out)


# HiMC
maf_himc_table_anova_out = paste0(out_prefix,     "maf_himc_anova.csv")
maf_himc_table_emmeans_out = paste0(out_prefix,   "maf_himc_emmeans.csv")
maf_himc_table_contrasts_out = paste0(out_prefix, "maf_himc_contrasts.csv")

maf_himc_cutoff_lm = lm(himc_cutoff_diff ~ refpan_maf, data = maf_combined_summary_means)
maf_himc_cutoff_lm_sum = summary(emmeans(maf_himc_cutoff_lm, pairwise ~ refpan_maf))

write_csv(anova(maf_himc_cutoff_lm),                        maf_himc_table_anova_out)
write_csv(maf_himc_cutoff_lm_sum$emmeans,                   maf_himc_table_emmeans_out)
write_csv(add_signf_code(maf_himc_cutoff_lm_sum$contrasts), maf_himc_table_contrasts_out)

khap_himc_table_anova_out = paste0(out_prefix,     "khap_himc_anova.csv")
khap_himc_table_emmeans_out = paste0(out_prefix,   "khap_himc_emmeans.csv")
khap_himc_table_contrasts_out = paste0(out_prefix, "khap_himc_contrasts.csv")

khap_himc_cutoff_lm = lm(himc_cutoff_diff ~ k_hap, data = khap_combined_summary_means)
khap_himc_cutoff_lm_sum = summary(emmeans(khap_himc_cutoff_lm, pairwise ~ k_hap))

write_csv(anova(khap_himc_cutoff_lm),                        khap_himc_table_anova_out)
write_csv(khap_himc_cutoff_lm_sum$emmeans,                   khap_himc_table_emmeans_out)
write_csv(add_signf_code(khap_himc_cutoff_lm_sum$contrasts), khap_himc_table_contrasts_out)

# HiMC - MACRO
maf_macrohimc_table_anova_out = paste0(out_prefix,     "maf_macrohimc_anova.csv")
maf_macrohimc_table_emmeans_out = paste0(out_prefix,   "maf_macrohimc_emmeans.csv")
maf_macrohimc_table_contrasts_out = paste0(out_prefix, "maf_macrohimc_contrasts.csv")

maf_macrohimc_cutoff_lm = lm(himc_macro_cutoff_diff ~ refpan_maf, data = maf_combined_summary_means)
maf_macrohimc_cutoff_lm_sum = summary(emmeans(maf_macrohimc_cutoff_lm, pairwise ~ refpan_maf))

write_csv(anova(maf_macrohimc_cutoff_lm),                        maf_macrohimc_table_anova_out)
write_csv(maf_macrohimc_cutoff_lm_sum$emmeans,                   maf_macrohimc_table_emmeans_out)
write_csv(add_signf_code(maf_macrohimc_cutoff_lm_sum$contrasts), maf_macrohimc_table_contrasts_out)

khap_macrohimc_table_anova_out = paste0(out_prefix,     "khap_macrohimc_anova.csv")
khap_macrohimc_table_emmeans_out = paste0(out_prefix,   "khap_macrohimc_emmeans.csv")
khap_macrohimc_table_contrasts_out = paste0(out_prefix, "khap_macrohimc_contrasts.csv")

khap_macrohimc_cutoff_lm = lm(himc_macro_cutoff_diff ~ k_hap, data = khap_combined_summary_means)
khap_macrohimc_cutoff_lm_sum = summary(emmeans(khap_macrohimc_cutoff_lm, pairwise ~ k_hap))

write_csv(anova(khap_macrohimc_cutoff_lm),                        khap_macrohimc_table_anova_out)
write_csv(khap_macrohimc_cutoff_lm_sum$emmeans,                   khap_macrohimc_table_emmeans_out)
write_csv(add_signf_code(khap_macrohimc_cutoff_lm_sum$contrasts), khap_macrohimc_table_contrasts_out)

# HAPLOGREP HAPLOGROUPS DIFFERENCE
haplogrep_recommendedOnly_anova_out     = paste0(out_prefix,     "haplogrep_recommended_settings_anova.csv")
haplogrep_recommendedOnly_emmeans_out   = paste0(out_prefix,     "haplogrep_recommended_settings_emmeans.csv")
haplogrep_recommendedOnly_contrasts_out = paste0(out_prefix,     "haplogrep_recommended_settings_contrasts.csv")

haplogrep_recommendedOnly_lm = lm(val ~ key + HaploGrep_Macrohaplogroup_truth, data = khap_haplogrep_recommendedOnly_long)
haplogrep_recommendedOnly_lm_sum = summary(emmeans(haplogrep_recommendedOnly_lm, pairwise ~ key|HaploGrep_Macrohaplogroup_truth))

write_csv(anova(haplogrep_recommendedOnly_lm),                        haplogrep_recommendedOnly_anova_out)
write_csv(haplogrep_recommendedOnly_lm_sum$emmeans,                   haplogrep_recommendedOnly_emmeans_out)
write_csv(add_signf_code(haplogrep_recommendedOnly_lm_sum$contrasts), haplogrep_recommendedOnly_contrasts_out)

# HAPLOGREP HAPLOGROUPS DIFFERENCE
himc_recommendedOnly_anova_out     = paste0(out_prefix,     "himc_recommended_settings_anova.csv")
himc_recommendedOnly_emmeans_out   = paste0(out_prefix,     "himc_recommended_settings_emmeans.csv")
himc_recommendedOnly_contrasts_out = paste0(out_prefix,     "himc_recommended_settings_contrasts.csv")

himc_recommendedOnly_lm = lm(val ~ key + HiMC_macrohaplogroup_truth, data = khap_himc_recommendedOnly_long)
himc_recommendedOnly_lm_sum = summary(emmeans(himc_recommendedOnly_lm, pairwise ~ key|HiMC_macrohaplogroup_truth))

write_csv(anova(himc_recommendedOnly_lm),                        himc_recommendedOnly_anova_out)
write_csv(himc_recommendedOnly_lm_sum$emmeans,                   himc_recommendedOnly_emmeans_out)
write_csv(add_signf_code(himc_recommendedOnly_lm_sum$contrasts), himc_recommendedOnly_contrasts_out)

# BEST PERFORMING 95%CIs PER CHIP
haplogrep_recommendedOnly_mcc_out           = paste0(out_prefix,     "recommended_settings_mcc_rank.csv")
haplogrep_recommendedOnly_haplogrep_out     = paste0(out_prefix,     "recommended_settings_HaploGrep_rank.csv")
haplogrep_recommendedOnly_himc_out          = paste0(out_prefix,     "recommended_settings_HiMC_rank.csv")
haplogrep_recommendedOnly_summary_out       = paste0(out_prefix,     "recommended_settings_summary.csv")

n_snp_summary = khap_info_combined_recommendedOnly %>%
  group_by(array) %>%
  summarise(n_snps = n(),
            n_snps_t0 = length(type[type == 0])) %>%
  mutate(array = as.character(array))

khap_info_perchip_summaries = khap_info_combined_recommendedOnly %>%
  #filter(type == 0) %>%
  group_by(array) %>%
  summarise(mcc_mean  = confidence.interval(mcc)[1],
            mcc_lower = confidence.interval(mcc)[2],
            mcc_upper = confidence.interval(mcc)[3],
            mcc_df    = confidence.interval(mcc)[4]) %>%
  arrange(desc(mcc_mean)) %>%
  mutate(mcc_rank = 1:101,
         array = as.factor(array)) %>%
  print(n=Inf) %>%
  write_csv(path = haplogrep_recommendedOnly_mcc_out)

khap_haplogrep_perchip_summaries = khap_haplogrep_recommendedOnly %>%
  group_by(array) %>%
  summarise(haplogrep_imputed_match_mean  = confidence.interval(HaploGrep_imputed_cutoff_macro_match)[1],
            haplogrep_imputed_match_lower = confidence.interval(HaploGrep_imputed_cutoff_macro_match)[2],
            haplogrep_imputed_match_upper = confidence.interval(HaploGrep_imputed_cutoff_macro_match)[3],
            haplogrep_imputed_match_df    = confidence.interval(HaploGrep_imputed_cutoff_macro_match)[4],
            haplogrep_typed_match_mean    = confidence.interval(HaploGrep_typed_macro_match)[1],
            haplogrep_typed_match_lower   = confidence.interval(HaploGrep_typed_macro_match)[2],
            haplogrep_typed_match_upper   = confidence.interval(HaploGrep_typed_macro_match)[3],
            haplogrep_typed_match_df      = confidence.interval(HaploGrep_typed_macro_match)[4]) %>%
  arrange(desc(haplogrep_imputed_match_mean)) %>%
  mutate(haplogrep_diff               = haplogrep_imputed_match_mean - haplogrep_typed_match_mean,
         haplogrep_imputed_match_rank = 1:101) %>%
  arrange(desc(haplogrep_typed_match_mean)) %>%
  mutate(haplogrep_typed_match_rank = 1:101) %>%
  arrange(desc(haplogrep_diff)) %>%
  mutate(haplogrep_diff_rank = 1:101,
         array = as.factor(array)) %>%
  print(n=Inf) %>%
  write_csv(path = haplogrep_recommendedOnly_haplogrep_out)

khap_himc_perchip_summaries = khap_himc_recommendedOnly %>%
  group_by(array) %>%
  summarise(himc_imputed_match_mean  = confidence.interval(HiMC_imputed_cutoff_macro_match)[1],
            himc_imputed_match_lower = confidence.interval(HiMC_imputed_cutoff_macro_match)[2],
            himc_imputed_match_upper = confidence.interval(HiMC_imputed_cutoff_macro_match)[3],
            himc_imputed_match_df    = confidence.interval(HiMC_imputed_cutoff_macro_match)[4],
            himc_typed_match_mean    = confidence.interval(HiMC_typed_macro_match)[1],
            himc_typed_match_lower   = confidence.interval(HiMC_typed_macro_match)[2],
            himc_typed_match_upper   = confidence.interval(HiMC_typed_macro_match)[3],
            himc_typed_match_df      = confidence.interval(HiMC_typed_macro_match)[4]) %>%
  arrange(desc(himc_imputed_match_mean)) %>%
  mutate(himc_diff               = himc_imputed_match_mean - himc_typed_match_mean,
         himc_imputed_match_rank = 1:101) %>%
  arrange(desc(himc_typed_match_mean)) %>%
  mutate(himc_typed_match_rank = 1:101) %>%
  arrange(desc(himc_diff)) %>%
  mutate(himc_himc_diff_rank = 1:101,
         array = as.factor(array)) %>%
  print(n=Inf) %>%
  write_csv(path = haplogrep_recommendedOnly_himc_out)

combined_perchip_summaries = n_snp_summary %>%
  left_join(khap_info_perchip_summaries, by = "array") %>%
  left_join(khap_haplogrep_perchip_summaries, by = "array") %>%
  left_join(khap_himc_perchip_summaries, by = "array") %>%
  arrange(desc(mcc_mean, haplogrep_diff, himc_diff)) %>%
  write_csv(path = haplogrep_recommendedOnly_summary_out)

# BEST PERFORMING 95%CIs OVERALL 
khap_info_combined_recommendedOnly %>%
  #filter(type == 0) %>%
  summarise(mcc_mean  = confidence.interval(mcc)[1],
            mcc_lower = confidence.interval(mcc)[2],
            mcc_upper = confidence.interval(mcc)[3],
            mcc_df    = confidence.interval(mcc)[4]) %>%
  arrange(desc(mcc_mean)) %>%
  print(n=Inf)


khap_haplogrep_recommendedOnly %>%
  summarise(haplogrep_imputed_match_mean  = confidence.interval(HaploGrep_imputed_cutoff_macro_match)[1],
            haplogrep_imputed_match_lower = confidence.interval(HaploGrep_imputed_cutoff_macro_match)[2],
            haplogrep_imputed_match_upper = confidence.interval(HaploGrep_imputed_cutoff_macro_match)[3],
            haplogrep_imputed_match_df    = confidence.interval(HaploGrep_imputed_cutoff_macro_match)[4],
            haplogrep_typed_match_mean    = confidence.interval(HaploGrep_typed_macro_match)[1],
            haplogrep_typed_match_lower   = confidence.interval(HaploGrep_typed_macro_match)[2],
            haplogrep_typed_match_upper   = confidence.interval(HaploGrep_typed_macro_match)[3],
            haplogrep_typed_match_df      = confidence.interval(HaploGrep_typed_macro_match)[4]) %>%
  mutate(haplogrep_diff = haplogrep_imputed_match_mean - haplogrep_typed_match_mean) %>%
  print(n=Inf)

khap_himc_perchip_summaries = khap_himc_recommendedOnly %>%
  summarise(himc_imputed_match_mean  = confidence.interval(HiMC_imputed_cutoff_macro_match)[1],
            himc_imputed_match_lower = confidence.interval(HiMC_imputed_cutoff_macro_match)[2],
            himc_imputed_match_upper = confidence.interval(HiMC_imputed_cutoff_macro_match)[3],
            himc_imputed_match_df    = confidence.interval(HiMC_imputed_cutoff_macro_match)[4],
            himc_typed_match_mean    = confidence.interval(HiMC_typed_macro_match)[1],
            himc_typed_match_lower   = confidence.interval(HiMC_typed_macro_match)[2],
            himc_typed_match_upper   = confidence.interval(HiMC_typed_macro_match)[3],
            himc_typed_match_df      = confidence.interval(HiMC_typed_macro_match)[4]) %>%
  arrange(desc(himc_imputed_match_mean)) %>%
  mutate(himc_diff = himc_imputed_match_mean - himc_typed_match_mean) %>%
  print(n=Inf)



###############################################

# ADNI STUFF!

adni_haplogrep_file = "/Volumes/TimMcInerney/MitoImpute/data/ADNI_REDO/IMPUTED/ReferencePanel_v1_0.01/IMPUTE2/MitoImpute_ReferencePanel_v1_0.01_imputed_HaploGrep_haplogroups.csv"
adni_himc_file      = "/Volumes/TimMcInerney/MitoImpute/data/ADNI_REDO/IMPUTED/ReferencePanel_v1_0.01/IMPUTE2/MitoImpute_ReferencePanel_v1_0.01_imputed_HiMC_haplogroups.csv"
adni_info_file      = "/Volumes/TimMcInerney/MitoImpute/data/ADNI_REDO/IMPUTED/ReferencePanel_v1_0.01/IMPUTE2/MitoImpute_ReferencePanel_v1_0.01_imputed_info"
adni_mcc_file       = "/Volumes/TimMcInerney/MitoImpute/data/ADNI_REDO/IMPUTED/ReferencePanel_v1_0.01/IMPUTE2/MitoImpute_ReferencePanel_v1_0.01_imputed_imputed_MCC.csv"

adni_haplogrep = read_csv(adni_haplogrep_file)
adni_himc      = read_csv(adni_himc_file)
adni_info      = read_delim(adni_info_file, delim = " ")
adni_mcc       = read_csv(adni_mcc_file)



adni_mcc %>%
  #filter(type == 0) %>%
  summarise(mcc_mean  = confidence.interval(mcc)[1],
            mcc_lower = confidence.interval(mcc)[2],
            mcc_upper = confidence.interval(mcc)[3],
            mcc_df    = confidence.interval(mcc)[4]) %>%
  arrange(desc(mcc_mean)) %>%
  print(n=Inf)







# END!


