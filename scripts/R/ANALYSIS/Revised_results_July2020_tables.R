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

maf_combined_summary_file  = "/Volumes/TimMcInerney/MitoImpute/analyses/combined_summaries/MAF_combined.csv"
khap_combined_summary_file = "/Volumes/TimMcInerney/MitoImpute/analyses/combined_summaries/kHAP_combined.csv"

out_prefix = "~/GitCode/MitoImputePrep/supplementary_information/R_tables/"

maf_combined_summary  = read_csv(maf_combined_summary_file)
khap_combined_summary = read_csv(khap_combined_summary_file)

mcc_maf_combined_summary
mac_khap_combined_summary

maf_combined_summary = maf_combined_summary %>% 
  mutate(array = as.factor(array),
         refpan_maf = factor(refpan_maf,
                             levels = maf.var,
                             labels = maf.lab,
                             ordered = T),
         k_hap      = factor(k_hap,
                             ordered = T))

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
                             labels = maf.var,
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

# Matthew's Correlation Coefficient
maf_mcc_table_anova_out = paste0(out_prefix,     "maf_mcc_anova.csv")
maf_mcc_table_emmeans_out = paste0(out_prefix,   "maf_mcc_emmeans.csv")
maf_mcc_table_contrasts_out = paste0(out_prefix, "maf_mcc_contrasts.csv")

maf_mcc_cutoff_lm = lm(mean_mcc_cutoff ~ refpan_maf, data = maf_combined_summary_means)
maf_mcc_cutoff_lm_sum = summary(emmeans(maf_mcc_cutoff_lm, pairwise ~ refpan_maf))

write_csv(anova(maf_mcc_cutoff_lm),        maf_mcc_table_anova_out)
write_csv(maf_mcc_cutoff_lm_sum$emmeans,   maf_mcc_table_emmeans_out)
write_csv(add_signf_code(maf_mcc_cutoff_lm_sum$contrasts), maf_mcc_table_contrasts_out)

khap_mcc_table_anova_out = paste0(out_prefix,     "khap_mcc_anova.csv")
khap_mcc_table_emmeans_out = paste0(out_prefix,   "khap_mcc_emmeans.csv")
khap_mcc_table_contrasts_out = paste0(out_prefix, "khap_mcc_contrasts.csv")

khap_mcc_cutoff_lm = lm(mean_mcc_cutoff ~ k_hap, data = khap_combined_summary_means)
khap_mcc_cutoff_lm_sum = summary(emmeans(khap_mcc_cutoff_lm, pairwise ~ k_hap))

write_csv(anova(khap_mcc_cutoff_lm),        khap_mcc_table_anova_out)
write_csv(khap_mcc_cutoff_lm_sum$emmeans,   khap_mcc_table_emmeans_out)
write_csv(add_signf_code(khap_mcc_cutoff_lm_sum$contrasts), khap_mcc_table_contrasts_out)


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


# END!