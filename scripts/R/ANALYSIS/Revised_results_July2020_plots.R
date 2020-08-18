library(tidyverse)
require(scales)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)

pal_cols = brewer.pal(12, "Dark2")

'%!in%' = function(x,y)!('%in%'(x,y))

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

afr_mthg_1000g = paste0("L", 0:6)
amr_mthg_1000g = c("A", "B", "C")
eas_mthg_1000g = c("D", "F", "G", "N", "Y", "Z")
asn_mthg_1000g = c("M", "R")
eur_sas_mthg_1000g = c("H", "I", "J", "K", "T", "V", "W", "X", "U")
aus_oceania = c("O", "P", "Q", "S", "E")

cty = read_excel("~/GitCode/MitoImputePrep/supplementary_information/MitoImpute_SupplementaryTables_abridged2.xlsx", sheet = 1, col_names = T, skip = 1)

cty = cty %>%
  mutate(haplogroup_color = if_else(`HaploGrep Macrohaplogroup` %in% afr_mthg_1000g, "Africa", 
                                    if_else(`HaploGrep Macrohaplogroup` %in% amr_mthg_1000g, "Americas",
                                            if_else(`HaploGrep Macrohaplogroup` %in% eas_mthg_1000g, "East Asia", 
                                                    if_else(`HaploGrep Macrohaplogroup` %in% eur_sas_mthg_1000g, "Europe and India",
                                                            if_else(`HaploGrep Macrohaplogroup` %in% asn_mthg_1000g, "Broadly Asian",
                                                                    if_else(`HaploGrep Macrohaplogroup` %in% aus_oceania, "Australiasia and Oceania", "unassigned")))))))

cty %>%
  filter(haplogroup_color == "unassigned") %>%
  group_by(`HaploGrep Macrohaplogroup`) %>%
  summarise(n = n()) %>%
  print(n=Inf)

pal_cols = brewer.pal(12, "Dark2")

mcmc.dir = "MCMC_Experiments"
mcmc.var = c("MCMC1", "MCMC5", "MCMC10", "MCMC20", "MCMC30")
khap.dir = "khap=Experiments"
khap.var = c("kHAP100", "kHAP250", "kHAP500", "kHAP1000", "kHAP2500", "kHAP5000", "kHAP10000", "kHAP20000", "kHAP30000")
khap.lab = c("khap=100", "khap=250", "khap=500", "khap=1000", "khap=2500", "khap=5000", "khap=10000", "khap=20000", "khap=30000")
maf.dir = "MAF_Experiments"
maf.panel = c("ReferencePanel_v2", "ReferencePanel_v4", "ReferencePanel_v3")
maf.var = c("MAF1%", "MAF0.5%", "MAF0.1%")
maf.lab = c("MAF>1%", "MAF>0.5%", "MAF>0.1%")

dir("/Volumes/TimMcInerney/MitoImpute/analyses/combined_summaries/")

maf_combined_summary_file  = "/Volumes/TimMcInerney/MitoImpute/analyses/combined_summaries/MAF_combined.csv"
khap_combined_summary_file = "/Volumes/TimMcInerney/MitoImpute/analyses/combined_summaries/kHAP_combined.csv"
maf_info_combined_file     = "/Volumes/TimMcInerney/MitoImpute/analyses/combined_summaries/MAF_MCC_concatenated.csv"
khap_info_combined_file    = "/Volumes/TimMcInerney/MitoImpute/analyses/combined_summaries/kHAP_MCC_concatenated.csv"
maf_haplogrep_file         = "/Volumes/TimMcInerney/MitoImpute/analyses/combined_summaries/MAF_HaploGrep_haplogroups_concatenated.tsv"
maf_himc_file              = "/Volumes/TimMcInerney/MitoImpute/analyses/combined_summaries/MAF_HiMC_haplogroups_concatenated.csv"
khap_haplogrep_file        = "/Volumes/TimMcInerney/MitoImpute/analyses/combined_summaries/kHAP_HaploGrep_haplogroups_concatenated.csv"
khap_himc_file             = "/Volumes/TimMcInerney/MitoImpute/analyses/combined_summaries/kHAP_HiMC_haplogroups_concatenated.csv"

out_prefix = "~/GitCode/MitoImputePrep/supplementary_information/R_plots/"

maf_combined_summary  = read_csv(maf_combined_summary_file)
khap_combined_summary = read_csv(khap_combined_summary_file)

maf_info_combined     = read_csv(maf_info_combined_file)
khap_info_combined    = read_csv(khap_info_combined_file)

maf_haplogrep         = read_tsv(maf_haplogrep_file)
maf_himc              = read_csv(maf_himc_file)

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

maf_haplogrep = maf_haplogrep %>%
  mutate(array = as.factor(array),
         refpan_maf = factor(refpan_maf,
                             levels = maf.var,
                             labels = maf.lab,
                             ordered = T),
         k_hap      = factor(k_hap,
                             ordered = T)) %>%
  mutate_at(vars(contains("Haplogroup")), .funs = as.factor)

maf_himc = maf_himc %>%
  mutate(array = as.factor(array),
         refpan_maf = factor(refpan_maf,
                             levels = maf.var,
                             labels = maf.lab,
                             ordered = T),
         k_hap      = factor(k_hap,
                             ordered = T)) %>%
  mutate_at(vars(contains("Haplogroup")), .funs = as.factor)

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

maf_info_combined_recommendedOnly = maf_info_combined %>%
  filter(refpan_maf == "MAF>0.1%", k_hap == "kHAP500")

maf_haplogrep_recommendedOnly = maf_haplogrep %>%
  filter(refpan_maf == "MAF>0.1%", k_hap == "kHAP500")

maf_himc_recommendedOnly = maf_himc %>%
  filter(refpan_maf == "MAF>0.1%", k_hap == "kHAP500")

khap_info_combined_recommendedOnly = khap_info_combined %>%
  filter(refpan_maf == "MAF>1%", k_hap == "khap=500")

khap_haplogrep_recommendedOnly = khap_haplogrep %>%
  filter(refpan_maf == "MAF>1%", k_hap == "khap=500")

khap_himc_recommendedOnly = khap_himc %>%
  filter(refpan_maf == "MAF>1%", k_hap == "khap=500")

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

#haplogrep_recommendedOnly_long = khap_haplogrep_recommendedOnly %>%
haplogrep_recommendedOnly_long = maf_haplogrep_recommendedOnly %>%
  gather(key = "key", value = "val", c(HaploGrep_imputed_cutoff_macro_match, HaploGrep_typed_macro_match)) %>%
  select(1:7, HaploGrep_Macrohaplogroup_truth, key, val) 

haplogrep_recommended_settings_out = paste0(out_prefix,     "haplogrep_recommended_settings.csv")
#haplogrep_perhaplogroup = khap_haplogrep_recommendedOnly %>%
haplogrep_perhaplogroup = maf_haplogrep_recommendedOnly %>%
  group_by(HaploGrep_Macrohaplogroup_truth) %>%
  summarise(prop_match_imputed = mean(HaploGrep_imputed_cutoff_macro_match, na.rm = T),
            prop_match_typed = mean(HaploGrep_typed_macro_match, na.rm = T)) %>%
  mutate(diff = prop_match_imputed - prop_match_typed) %>%
  arrange(desc(diff)) %>%
  print(n = Inf) %>%
  #write_csv(path = haplogrep_recommended_settings_out)

haplogrep_perhaplogroup %>%
  arrange(prop_match_typed) %>%
  print(n=Inf)

haplogrep_perhaplogroup %>%
  arrange(diff) %>%
  print(n=Inf)

haplogrep_perhaplogroup %>%
  filter(diff < 0, HaploGrep_Macrohaplogroup_truth %in% eur_sas_mthg_1000g) %>%
  print(n=Inf)

#himc_recommendedOnly_long = khap_himc_recommendedOnly %>%
himc_recommendedOnly_long = maf_himc_recommendedOnly %>%
  gather(key = "key", value = "val", c(HiMC_imputed_cutoff_macro_match, HiMC_typed_macro_match)) %>%
  select(1:7, HiMC_macrohaplogroup_truth, key, val) 

himc_recommended_settings_out = paste0(out_prefix,     "himc_recommended_settings.csv")
#himc_perhaplogroup = khap_himc_recommendedOnly %>%
himc_perhaplogroup = maf_himc_recommendedOnly %>%
  group_by(HiMC_macrohaplogroup_truth) %>%
  filter(!is.na(HiMC_macrohaplogroup_truth)) %>%
  summarise(prop_match_imputed = mean(HiMC_imputed_cutoff_macro_match, na.rm = T),
            prop_match_typed = mean(HiMC_typed_macro_match, na.rm = T)) %>%
  mutate(diff = prop_match_imputed - prop_match_typed) %>%
  arrange(desc(diff)) %>%
  print(n = Inf) %>%
  #write_csv(path = himc_recommended_settings_out)

himc_perhaplogroup %>%
  arrange(prop_match_typed) %>%
  print(n=Inf)

himc_perhaplogroup %>%
  filter(startsWith(HiMC_macrohaplogroup_truth, "L")) %>%
  arrange(diff) %>%
  print(n=Inf)

himc_perhaplogroup %>%
  filter(diff < 0, HiMC_macrohaplogroup_truth %in% eur_sas_mthg_1000g) %>%
  print(n=Inf)

maf_info_combined_summary
khap_info_combined_summary
maf_combined_summary_means
khap_combined_summary_means


########## PLOTTING
#plot_dir = " ~/GitCode/MitoImputePrep/supplementary_information/R_plots/"
plot_dir = "~/GitCode/MitoImputePrep/Plots/August2020/"

plot_multiplier = 1
wd = 254 * plot_multiplier
ht = 117 * plot_multiplier
DPI = 600

# MCC
maf_mcc_plot = ggplot(maf_info_combined_summary, aes(x = refpan_maf, y = mean_mcc)) +
  geom_violin(na.rm = T, fill = pal_cols[2], alpha = (1/4), lwd = (1/4)) +
  geom_boxplot(na.rm = T, fill = pal_cols[2], alpha = (1/2), lwd = (1/4), width = (1/5), outlier.alpha = (1/1)) +
  labs(x = "Reference panel\nminor allele frequency",
       #y = "Matthew's correlation coefficient (MCC)",
       y = "") +
  scale_y_continuous(limits = c(0, 1)) +
  theme(axis.text.x = element_text(colour = "black", size = rel(1.125)),
        axis.text.y = element_text(colour = "black", size = rel(1.125)),
        axis.title.x = element_text(colour = "black", size = rel(1.5)),
        axis.title.y = element_text(colour = "black", size = rel(1.5)),
        legend.text = element_text(size = rel(1.0)),
        legend.title = element_text(size = rel(1.25)),
        legend.position= "top",
        plot.title = element_blank(),
        strip.text.x = element_text(size = rel(5/5)),
        strip.text.y = element_text(size = rel(5/5),),
        panel.grid.major.x = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.grid.minor.x = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = "black"),
        strip.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  coord_flip()

khap_mcc_plot = ggplot(khap_info_combined_summary, aes(x = reorder(k_hap, desc(k_hap)), y = mean_mcc)) +
  geom_violin(na.rm = T, fill = pal_cols[2], alpha = (1/4), lwd = (1/4)) +
  geom_boxplot(na.rm = T, fill = pal_cols[2], alpha = (1/2), lwd = (1/4), width = (1/5), outlier.alpha = (1/1)) +
  labs(x = "Number of included\nreference haplotypes",
       y = "Matthew's correlation coefficient (MCC)") +
  #scale_x_discrete(limits = rev(levels(k_hap))) +
  scale_y_continuous(limits = c(0, 1)) +
  theme(axis.text.x = element_text(colour = "black", size = rel(1.125)),
        axis.text.y = element_text(colour = "black", size = rel(1.125)),
        axis.title.x = element_text(colour = "black", size = rel(1.5)),
        axis.title.y = element_text(colour = "black", size = rel(1.5)),
        legend.text = element_text(size = rel(1.0)),
        legend.title = element_text(size = rel(1.25)),
        legend.position= "top",
        plot.title = element_blank(),
        strip.text.x = element_text(size = rel(5/5)),
        strip.text.y = element_text(size = rel(5/5),),
        panel.grid.major.x = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.grid.minor.x = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = "black"),
        strip.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  coord_flip()
ggarrange(maf_mcc_plot, khap_mcc_plot, align = "hv", nrow = 2)
ggsave(filename = paste0(plot_dir, "mean_mcc.png"), bg = "transparent", width = 250, height = 250, units = "mm", dpi = DPI)

# HAPLOGREP
maf_haplogrep_typed_plot = ggplot(maf_combined_summary_means, aes(x = refpan_maf, y = mean_haplogrep_concordance_typed_macro)) +
  geom_violin(na.rm = T, fill = pal_cols[2], alpha = (1/4), lwd = (1/4)) +
  geom_boxplot(na.rm = T, fill = pal_cols[2], alpha = (1/2), lwd = (1/4), width = (1/5), outlier.alpha = (1/1)) +
  labs(x = "Reference panel\nminor allele frequency",
       #y = "Mean macrohaplogroup concordance\n(HaploGrep2.0)",
       y = "") +
  scale_y_continuous(limits = c(0, 1)) +
  theme(axis.text.x = element_text(colour = "black", size = rel(1.125)),
        axis.text.y = element_text(colour = "black", size = rel(1.125)),
        axis.title.x = element_text(colour = "black", size = rel(1.5)),
        axis.title.y = element_text(colour = "black", size = rel(1.5)),
        legend.text = element_text(size = rel(1.0)),
        legend.title = element_text(size = rel(1.25)),
        legend.position= "top",
        plot.title = element_blank(),
        strip.text.x = element_text(size = rel(5/5)),
        strip.text.y = element_text(size = rel(5/5),),
        panel.grid.major.x = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.grid.minor.x = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = "black"),
        strip.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  coord_flip()

maf_haplogrep_imputed_plot = ggplot(maf_combined_summary_means, aes(x = refpan_maf, y = mean_haplogrep_concordance_imputed_macro)) +
  geom_violin(na.rm = T, fill = pal_cols[2], alpha = (1/4), lwd = (1/4)) +
  geom_boxplot(na.rm = T, fill = pal_cols[2], alpha = (1/2), lwd = (1/4), width = (1/5), outlier.alpha = (1/1)) +
  labs(x = "",
       #x = "Reference panel\nminor allele frequency",
       #y = "Mean macrohaplogroup concordance\n(HaploGrep2.0)",
       y = "") +
  scale_y_continuous(limits = c(0, 1)) +
  theme(axis.text.x = element_text(colour = "black", size = rel(1.125)),
        axis.text.y = element_text(colour = "black", size = rel(1.125)),
        axis.title.x = element_text(colour = "black", size = rel(1.5)),
        axis.title.y = element_text(colour = "black", size = rel(1.5)),
        legend.text = element_text(size = rel(1.0)),
        legend.title = element_text(size = rel(1.25)),
        legend.position= "top",
        plot.title = element_blank(),
        strip.text.x = element_text(size = rel(5/5)),
        strip.text.y = element_text(size = rel(5/5),),
        panel.grid.major.x = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.grid.minor.x = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = "black"),
        strip.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  coord_flip()

maf_haplogrep_diff_plot = ggplot(maf_combined_summary_means, aes(x = refpan_maf, y = haplogrep_macro_cutoff_diff)) +
  geom_violin(na.rm = T, fill = pal_cols[2], alpha = (1/4), lwd = (1/4)) +
  geom_boxplot(na.rm = T, fill = pal_cols[2], alpha = (1/2), lwd = (1/4), width = (1/5), outlier.alpha = (1/1)) +
  labs(x = "",
       #x = "Reference panel\nminor allele frequency",
       #y = "Mean difference in macrohaplogroup concordance\n(HaploGrep2.0)",
       y = "") +
  scale_y_continuous(limits = c(-0.25, 0.25)) +
  theme(axis.text.x = element_text(colour = "black", size = rel(1.125)),
        axis.text.y = element_text(colour = "black", size = rel(1.125)),
        axis.title.x = element_text(colour = "black", size = rel(1.5)),
        axis.title.y = element_text(colour = "black", size = rel(1.5)),
        legend.text = element_text(size = rel(1.0)),
        legend.title = element_text(size = rel(1.25)),
        legend.position= "top",
        plot.title = element_blank(),
        strip.text.x = element_text(size = rel(5/5)),
        strip.text.y = element_text(size = rel(5/5),),
        panel.grid.major.x = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.grid.minor.x = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = "black"),
        strip.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  coord_flip()

khap_haplogrep_typed_plot = ggplot(khap_combined_summary_means, aes(x = reorder(k_hap, desc(k_hap)), y = mean_haplogrep_concordance_typed_macro)) +
  geom_violin(na.rm = T, fill = pal_cols[2], alpha = (1/4), lwd = (1/4)) +
  geom_boxplot(na.rm = T, fill = pal_cols[2], alpha = (1/2), lwd = (1/4), width = (1/5), outlier.alpha = (1/1)) +
  labs(x = "Number of included\nreference haplotypes",
       y = "Mean macrohaplogroup concordance\n(HaploGrep2.0)") +
  scale_y_continuous(limits = c(0, 1)) +
  theme(axis.text.x = element_text(colour = "black", size = rel(1.125)),
        axis.text.y = element_text(colour = "black", size = rel(1.125)),
        axis.title.x = element_text(colour = "black", size = rel(1.5)),
        axis.title.y = element_text(colour = "black", size = rel(1.5)),
        legend.text = element_text(size = rel(1.0)),
        legend.title = element_text(size = rel(1.25)),
        legend.position= "top",
        plot.title = element_blank(),
        strip.text.x = element_text(size = rel(5/5)),
        strip.text.y = element_text(size = rel(5/5),),
        panel.grid.major.x = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.grid.minor.x = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = "black"),
        strip.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  coord_flip()

khap_haplogrep_imputed_plot = ggplot(khap_combined_summary_means, aes(x = reorder(k_hap, desc(k_hap)), y = mean_haplogrep_concordance_imputed_macro)) +
  geom_violin(na.rm = T, fill = pal_cols[2], alpha = (1/4), lwd = (1/4)) +
  geom_boxplot(na.rm = T, fill = pal_cols[2], alpha = (1/2), lwd = (1/4), width = (1/5), outlier.alpha = (1/1)) +
  labs(x = "",
       #x = "Number of included\nreference haplotypes",
       y = "Mean macrohaplogroup concordance\n(HaploGrep2.0)") +
  scale_y_continuous(limits = c(0, 1)) +
  theme(axis.text.x = element_text(colour = "black", size = rel(1.125)),
        axis.text.y = element_text(colour = "black", size = rel(1.125)),
        axis.title.x = element_text(colour = "black", size = rel(1.5)),
        axis.title.y = element_text(colour = "black", size = rel(1.5)),
        legend.text = element_text(size = rel(1.0)),
        legend.title = element_text(size = rel(1.25)),
        legend.position= "top",
        plot.title = element_blank(),
        strip.text.x = element_text(size = rel(5/5)),
        strip.text.y = element_text(size = rel(5/5),),
        panel.grid.major.x = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.grid.minor.x = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = "black"),
        strip.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  coord_flip()

khap_haplogrep_diff_plot = ggplot(khap_combined_summary_means, aes(x = reorder(k_hap, desc(k_hap)), y = haplogrep_macro_cutoff_diff)) +
  geom_violin(na.rm = T, fill = pal_cols[2], alpha = (1/4), lwd = (1/4)) +
  geom_boxplot(na.rm = T, fill = pal_cols[2], alpha = (1/2), lwd = (1/4), width = (1/5), outlier.alpha = (1/1)) +
  labs(x = "",
       #x = "Number of included\nreference haplotypes",
       y = "Mean difference in macrohaplogroup concordance\n(HaploGrep2.0)") +
  scale_y_continuous(limits = c(-0.25, 0.25)) +
  theme(axis.text.x = element_text(colour = "black", size = rel(1.125)),
        axis.text.y = element_text(colour = "black", size = rel(1.125)),
        axis.title.x = element_text(colour = "black", size = rel(1.5)),
        axis.title.y = element_text(colour = "black", size = rel(1.5)),
        legend.text = element_text(size = rel(1.0)),
        legend.title = element_text(size = rel(1.25)),
        legend.position= "top",
        plot.title = element_blank(),
        strip.text.x = element_text(size = rel(5/5)),
        strip.text.y = element_text(size = rel(5/5),),
        panel.grid.major.x = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.grid.minor.x = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = "black"),
        strip.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  coord_flip()

ggarrange(maf_haplogrep_typed_plot, maf_haplogrep_imputed_plot, maf_haplogrep_diff_plot,
          khap_haplogrep_typed_plot, khap_haplogrep_imputed_plot, khap_haplogrep_diff_plot,
          align = "hv", nrow = 2, ncol = 3)
ggsave(filename = paste0(plot_dir, "maf_haplgrep_comparison.png"), bg = "transparent", width = wd*2, height = ht*2, units = "mm", dpi = DPI)

# HIMC
maf_himc_typed_plot = ggplot(maf_combined_summary_means, aes(x = refpan_maf, y = mean_himc_concordance_typed_macro)) +
  geom_violin(na.rm = T, fill = pal_cols[2], alpha = (1/4), lwd = (1/4)) +
  geom_boxplot(na.rm = T, fill = pal_cols[2], alpha = (1/2), lwd = (1/4), width = (1/5), outlier.alpha = (1/1)) +
  labs(x = "Reference panel\nminor allele frequency",
       #y = "Mean macrohaplogroup concordance\n(Hi-MC)",
       y = "") +
  scale_y_continuous(limits = c(0, 1)) +
  theme(axis.text.x = element_text(colour = "black", size = rel(1.125)),
        axis.text.y = element_text(colour = "black", size = rel(1.125)),
        axis.title.x = element_text(colour = "black", size = rel(1.5)),
        axis.title.y = element_text(colour = "black", size = rel(1.5)),
        legend.text = element_text(size = rel(1.0)),
        legend.title = element_text(size = rel(1.25)),
        legend.position= "top",
        plot.title = element_blank(),
        strip.text.x = element_text(size = rel(5/5)),
        strip.text.y = element_text(size = rel(5/5),),
        panel.grid.major.x = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.grid.minor.x = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = "black"),
        strip.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  coord_flip()

maf_himc_imputed_plot = ggplot(maf_combined_summary_means, aes(x = refpan_maf, y = mean_himc_concordance_imputed_macro)) +
  geom_violin(na.rm = T, fill = pal_cols[2], alpha = (1/4), lwd = (1/4)) +
  geom_boxplot(na.rm = T, fill = pal_cols[2], alpha = (1/2), lwd = (1/4), width = (1/5), outlier.alpha = (1/1)) +
  labs(x = "",
       #x = "Reference panel\nminor allele frequency",
       #y = "Mean macrohaplogroup concordance\n(Hi-MC)",
       y = "") +
  scale_y_continuous(limits = c(0, 1)) +
  theme(axis.text.x = element_text(colour = "black", size = rel(1.125)),
        axis.text.y = element_text(colour = "black", size = rel(1.125)),
        axis.title.x = element_text(colour = "black", size = rel(1.5)),
        axis.title.y = element_text(colour = "black", size = rel(1.5)),
        legend.text = element_text(size = rel(1.0)),
        legend.title = element_text(size = rel(1.25)),
        legend.position= "top",
        plot.title = element_blank(),
        strip.text.x = element_text(size = rel(5/5)),
        strip.text.y = element_text(size = rel(5/5),),
        panel.grid.major.x = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.grid.minor.x = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = "black"),
        strip.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  coord_flip()

maf_himc_diff_plot = ggplot(maf_combined_summary_means, aes(x = refpan_maf, y = himc_macro_cutoff_diff)) +
  geom_violin(na.rm = T, fill = pal_cols[2], alpha = (1/4), lwd = (1/4)) +
  geom_boxplot(na.rm = T, fill = pal_cols[2], alpha = (1/2), lwd = (1/4), width = (1/5), outlier.alpha = (1/1)) +
  labs(x = "",
       #x = "Reference panel\nminor allele frequency",
       #y = "Mean macrohaplogroup concordance\n(Hi-MC)",
       y = "") +
  scale_y_continuous(limits = c(-0.65, 0.55)) +
  theme(axis.text.x = element_text(colour = "black", size = rel(1.125)),
        axis.text.y = element_text(colour = "black", size = rel(1.125)),
        axis.title.x = element_text(colour = "black", size = rel(1.5)),
        axis.title.y = element_text(colour = "black", size = rel(1.5)),
        legend.text = element_text(size = rel(1.0)),
        legend.title = element_text(size = rel(1.25)),
        legend.position= "top",
        plot.title = element_blank(),
        strip.text.x = element_text(size = rel(5/5)),
        strip.text.y = element_text(size = rel(5/5),),
        panel.grid.major.x = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.grid.minor.x = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = "black"),
        strip.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  coord_flip()

khap_himc_typed_plot = ggplot(khap_combined_summary_means, aes(x = reorder(k_hap, desc(k_hap)), y = mean_himc_concordance_typed_macro)) +
  geom_violin(na.rm = T, fill = pal_cols[2], alpha = (1/4), lwd = (1/4)) +
  geom_boxplot(na.rm = T, fill = pal_cols[2], alpha = (1/2), lwd = (1/4), width = (1/5), outlier.alpha = (1/1)) +
  labs(x = "Number of included\nreference haplotypes",
       y = "Mean macrohaplogroup concordance\n(Hi-MC)") +
  scale_y_continuous(limits = c(0, 1)) +
  theme(axis.text.x = element_text(colour = "black", size = rel(1.125)),
        axis.text.y = element_text(colour = "black", size = rel(1.125)),
        axis.title.x = element_text(colour = "black", size = rel(1.5)),
        axis.title.y = element_text(colour = "black", size = rel(1.5)),
        legend.text = element_text(size = rel(1.0)),
        legend.title = element_text(size = rel(1.25)),
        legend.position= "top",
        plot.title = element_blank(),
        strip.text.x = element_text(size = rel(5/5)),
        strip.text.y = element_text(size = rel(5/5),),
        panel.grid.major.x = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.grid.minor.x = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = "black"),
        strip.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  coord_flip()

khap_himc_imputed_plot = ggplot(khap_combined_summary_means, aes(x = reorder(k_hap, desc(k_hap)), y = mean_himc_concordance_imputed_macro)) +
  geom_violin(na.rm = T, fill = pal_cols[2], alpha = (1/4), lwd = (1/4)) +
  geom_boxplot(na.rm = T, fill = pal_cols[2], alpha = (1/2), lwd = (1/4), width = (1/5), outlier.alpha = (1/1)) +
  labs(x = "",
       #x = "Number of included\nreference haplotypes",
       y = "Mean macrohaplogroup concordance\n(Hi-MC)") +
  scale_y_continuous(limits = c(0, 1)) +
  theme(axis.text.x = element_text(colour = "black", size = rel(1.125)),
        axis.text.y = element_text(colour = "black", size = rel(1.125)),
        axis.title.x = element_text(colour = "black", size = rel(1.5)),
        axis.title.y = element_text(colour = "black", size = rel(1.5)),
        legend.text = element_text(size = rel(1.0)),
        legend.title = element_text(size = rel(1.25)),
        legend.position= "top",
        plot.title = element_blank(),
        strip.text.x = element_text(size = rel(5/5)),
        strip.text.y = element_text(size = rel(5/5),),
        panel.grid.major.x = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.grid.minor.x = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = "black"),
        strip.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  coord_flip()

khap_himc_diff_plot = ggplot(khap_combined_summary_means, aes(x = reorder(k_hap, desc(k_hap)), y = himc_macro_cutoff_diff)) +
  geom_violin(na.rm = T, fill = pal_cols[2], alpha = (1/4), lwd = (1/4)) +
  geom_boxplot(na.rm = T, fill = pal_cols[2], alpha = (1/2), lwd = (1/4), width = (1/5), outlier.alpha = (1/1)) +
  labs(x = "",
       #x = "Number of included\nreference haplotypes",x = "Number of included\nreference haplotypes",
       y = "Mean difference in macrohaplogroup concordance\n(Hi-MC)") +
  scale_y_continuous(limits = c(-0.65, 0.55)) +
  theme(axis.text.x = element_text(colour = "black", size = rel(1.125)),
        axis.text.y = element_text(colour = "black", size = rel(1.125)),
        axis.title.x = element_text(colour = "black", size = rel(1.5)),
        axis.title.y = element_text(colour = "black", size = rel(1.5)),
        legend.text = element_text(size = rel(1.0)),
        legend.title = element_text(size = rel(1.25)),
        legend.position= "top",
        plot.title = element_blank(),
        strip.text.x = element_text(size = rel(5/5)),
        strip.text.y = element_text(size = rel(5/5),),
        panel.grid.major.x = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.grid.minor.x = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = "black"),
        strip.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  coord_flip()

ggarrange(maf_himc_typed_plot, maf_himc_imputed_plot, maf_himc_diff_plot,
          khap_himc_typed_plot, khap_himc_imputed_plot, khap_himc_diff_plot,
          align = "hv", nrow = 2, ncol = 3)
ggsave(filename = paste0(plot_dir, "maf_himc_comparison.png"), bg = "transparent", width = wd*2, height = ht*2, units = "mm", dpi = DPI)

# HAPLOGROUPINGS PER 
ggplot(cty, aes(x = `HaploGrep Macrohaplogroup`, fill = haplogroup_color)) +
  geom_bar(na.rm = T) +
  labs(x = "Macrohaplogroup",
       y = "Frequency",
       fill = "Continental association") +
  guides(colour = F,
         fill = guide_legend(title.position="top", 
                             title.hjust = 0.5)) +
  scale_y_continuous(breaks = seq(0, 1e4, by = 2e3),
                     labels = comma) +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.text.x = element_text(colour = "black", size = rel(1.125)),
        axis.text.y = element_text(colour = "black", size = rel(1.125)),
        axis.title.x = element_text(colour = "black", size = rel(1.5)),
        axis.title.y = element_text(colour = "black", size = rel(1.5)),
        legend.text = element_text(size = rel(1.0)),
        legend.title = element_text(size = rel(1.25)),
        legend.position= "top",
        plot.title = element_blank(),
        strip.text.x = element_text(size = rel(5/5)),
        strip.text.y = element_text(size = rel(5/5),),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.grid.minor.y = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.background = element_rect(fill = "transparent", colour = "black"),
        strip.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))
ggsave(filename = paste0(plot_dir, "haplogroup_count.png"), bg = "transparent", width = wd, height = ht, units = "mm", dpi = DPI)




















# END!
