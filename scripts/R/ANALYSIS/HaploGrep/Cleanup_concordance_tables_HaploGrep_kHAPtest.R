require(ggplot2)
require(gridExtra)
require(tidyr)
require(emmeans)
require(dplyr)

#########################################################################
chip.table = read.table("~/GitCode/MitoImputePrep/scripts/INFORMATION_LISTS/b37_platforms.txt", header = F)
truth.table = read.table("~/GitCode/MitoImputePrep/metadata/HaploGrep_concordance/chrMT_1kg_diploid_haplogrep.txt", header = T)
container = "/Volumes/TimMcInerney/MitoImpute/data/HAPLOGROUPS/chips/" #BDCHP-1X10-HUMANHAP240S_11216501_A-b37/MCMC_Experiments/MCMC1"

truth.A = subset(truth.table, substr(truth.table$Haplogroup, 1, 1) == "L")
truth.A$Macrohaplogroup = substr(truth.A$Haplogroup, 1, 2)
truth.N = subset(truth.table, substr(truth.table$Haplogroup, 1, 1) != "L")
truth.N$Macrohaplogroup = substr(truth.N$Haplogroup, 1, 1)

truth.table = rbind(truth.A, truth.N)
truth.table = arrange(truth.table, truth.table$SampleID)

###############################################################################
## KHAP 

exp.dir = "kHAP_Experiments"
exp.var = c("kHAP100", "kHAP250", "kHAP500", "kHAP1000", "kHAP2500", "kHAP5000", "kHAP10000", "kHAP20000", "kHAP30000")
tmp_khap_df = data.frame(matrix(ncol = 1, nrow = nrow(chip.table)))
names(tmp_khap_df) = c("chip")
tmp_khap_df$chip = chip.table$V1
tmp_khap_df$experiment = exp.dir

main_khap_df = data.frame()

for (exp in 1:length(exp.var)) {
  out.file = paste0("/Volumes/TimMcInerney/MitoImpute/data/HAPLOGROUPS/combined/kHAP_MAF01/ConcordanceTables_", exp.var[exp], ".csv")
  tmp_khap_df$sub_experiment = exp.var[exp]
  for (chip in 1:nrow(chip.table)) {
    tmp.file = paste0(container, chip.table$V1[chip], "/", "ReferencePanel_v3/", exp.dir, "/", exp.var[exp], "/", "chrMT_1kg_", chip.table$V1[chip], "_imputed_", exp.var[exp], "_haplogrep.txt")
    if (file.exists(tmp.file) == T) {
      tmp_khap_df$imputed[chip] = T
      chip.table$imputed[chip] = T
    } else {
      tmp_khap_df$imputed[chip] = F
      chip.table$imputed[chip] = F
    }
  }
  
  tmp_khap_df$info_score[chip] = NA
  tmp_khap_df$typed_match[chip] = NA
  tmp_khap_df$typed_macro_match[chip] = NA
  tmp_khap_df$imputed_match[chip] = NA
  tmp_khap_df$imputed_macro_match[chip] = NA
  
  total = nrow(chip.table)
  total_imputed = nrow(subset(chip.table, chip.table$imputed == T))
  
  for (chip in 1:nrow(chip.table)) {
    message(paste0("WORKING ON CHIP FOR ", exp.var[exp], ":\t", chip, " / ", total))
    # TYPED FILE
    tmp1.file = paste0(container, chip.table$V1[chip], "/", exp.dir, "/", "chrMT_1kg_", chip.table$V1[chip], "_diploid.txt")
    if (file.exists(tmp1.file) == T) {
      tmp1.hg.table = read.table(tmp1.file, header = T)
      tmp1.hg.table$Range = NULL
      
      ## CREATE COLUMN FOR MACRO HAPLOGROUPS
      # USE FIRST LETTER AND NUMBER FOR AFRICAN CLADES
      tmp1.hg.table.A = subset(tmp1.hg.table, substr(tmp1.hg.table$Haplogroup, 1, 1) == "L")
      tmp1.hg.table.A$Macrohaplogroup = substr(tmp1.hg.table.A$Haplogroup, 1, 2)
      # USE FIRST LETTER FOR NON-AFRICAN CLADES
      tmp1.hg.table.N = subset(tmp1.hg.table, substr(tmp1.hg.table$Haplogroup, 1, 1) != "L")
      tmp1.hg.table.N$Macrohaplogroup = substr(tmp1.hg.table.N$Haplogroup, 1, 1)
      # COMBINED AND REARRANGE
      tmp1.hg.table = rbind(tmp1.hg.table.A, tmp1.hg.table.N)
      tmp1.hg.table = arrange(tmp1.hg.table, tmp1.hg.table$SampleID)
      
      tmp1.hg.table$typed_match = as.character(truth.table$Haplogroup) == as.character(tmp1.hg.table$Haplogroup)
      tmp1.hg.table$typed_macro_match = as.character(truth.table$Macrohaplogroup) == as.character(tmp1.hg.table$Macrohaplogroup)
      
      tmp_khap_df$typed_match[chip] = nrow(subset(tmp1.hg.table, tmp1.hg.table$typed_match == T)) /  nrow(tmp1.hg.table)
      tmp_khap_df$typed_macro_match[chip] = nrow(subset(tmp1.hg.table, tmp1.hg.table$typed_macro_match == T)) /  nrow(tmp1.hg.table)
    }
    
    # IMPUTED FILE
    tmp2.file = paste0(container, chip.table$V1[chip], "/", exp.dir, "/", exp.var[exp], "/", "chrMT_1kg_", chip.table$V1[chip], "_imputed_", exp.var[exp], "_haplogrep.txt")
    tmp3.file = paste0(container, chip.table$V1[chip], "/", exp.dir, "/", exp.var[exp], "/", "chrMT_1kg_", chip.table$V1[chip], "_imputed_", exp.var[exp], "_info")
    if (file.exists(tmp2.file) == T) {
      tmp2.hg.table = read.table(tmp2.file, header = T)
      tmp2.hg.table$Range = NULL
      
      ## GATHER INFO SCORE
      info.table = read.table(tmp3.file, header = T)
      tmp_khap_df$info_score[chip] = mean(info.table$info, na.rm = T)
      
      ## CREATE COLUMN FOR MACRO HAPLOGROUPS
      # USE FIRST LETTER AND NUMBER FOR AFRICAN CLADES
      tmp2.hg.table.A = subset(tmp2.hg.table, substr(tmp2.hg.table$Haplogroup, 1, 1) == "L")
      tmp2.hg.table.A$Macrohaplogroup = substr(tmp2.hg.table.A$Haplogroup, 1, 2)
      # USE FIRST LETTER FOR NON-AFRICAN CLADES
      tmp2.hg.table.N = subset(tmp2.hg.table, substr(tmp2.hg.table$Haplogroup, 1, 1) != "L")
      tmp2.hg.table.N$Macrohaplogroup = substr(tmp2.hg.table.N$Haplogroup, 1, 1)
      # COMBINED AND REARRANGE
      tmp2.hg.table = rbind(tmp2.hg.table.A, tmp2.hg.table.N)
      tmp2.hg.table = arrange(tmp2.hg.table, tmp2.hg.table$SampleID)
      
      tmp2.hg.table$imputed_match = as.character(truth.table$Haplogroup) == as.character(tmp2.hg.table$Haplogroup)
      tmp2.hg.table$imputed_macro_match = as.character(truth.table$Macrohaplogroup) == as.character(tmp2.hg.table$Macrohaplogroup)
      
      tmp_khap_df$imputed_match[chip] = nrow(subset(tmp2.hg.table, tmp2.hg.table$imputed_match == T)) / nrow(tmp2.hg.table)
      tmp_khap_df$imputed_macro_match[chip] = nrow(subset(tmp2.hg.table, tmp2.hg.table$imputed_macro_match == T)) / nrow(tmp2.hg.table)
    }
  }
  write.csv(tmp_khap_df, out.file, row.names = F, quote = F)
  message(paste0("WROTE ", out.file, " TO DISK"))
  main_khap_df = rbind(main_khap_df, tmp_khap_df)
}
main_khap_df$diff = main_khap_df$imputed_match - main_khap_df$typed_match
main_khap_df$diff_macro = main_khap_df$imputed_macro_match - main_khap_df$typed_macro_match
write.csv(main_khap_df, paste0("/Volumes/TimMcInerney/MitoImpute/data/HAPLOGROUPS/combined/kHAP_MAF01/ConcordanceTables_", exp.dir,"_COMBINED.csv"), row.names = F, quote = F)

exp.dir = "kHAP_Experiments"
exp.var = c("kHAP100", "kHAP250", "kHAP500", "kHAP1000", "kHAP2500", "kHAP5000", "kHAP10000", "kHAP20000", "kHAP30000")
main_khap_df$sub_experiment = factor(main_khap_df$sub_experiment, levels = exp.var)
k_hap_box = ggplot(main_khap_df, aes(x = sub_experiment, y = imputed_match)) +
  geom_violin(fill = "#feb600", na.rm = T) +
  geom_boxplot(width = 0.125, notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428") +
  theme_bw() +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                     limits = c(0, 1)) +
  labs(x = "Number of reference haplotypes used",
       y = "% concordance with resequenced dataset",
       title = expression(paste(bold("B."), " Number of included reference haplotypes (k_hap) variations")))
k_hap_box
ggsave(filename = paste0("/Volumes/TimMcInerney/MitoImpute/data/HAPLOGROUPS/plots/kHAP_MAF01/ConcordanceTables_", exp.dir,"_HaploGrep.png"), plot = k_hap_box, units = "mm", width = 297, height = 210, dpi = 300)

macro_k_hap_box = ggplot(main_khap_df, aes(x = sub_experiment, y = imputed_macro_match)) +
  geom_violin(fill = "#feb600", na.rm = T) +
  geom_boxplot(width = 0.125, notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428") +
  theme_bw() +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                     limits = c(0, 1)) +
  labs(x = "Number of reference haplotypes used",
       y = "% concordance with resequenced dataset",
       title = expression(paste(bold("B."), " Number of included reference haplotypes (k_hap) variations")))
macro_k_hap_box
ggsave(filename = paste0("/Volumes/TimMcInerney/MitoImpute/data/HAPLOGROUPS/plots/kHAP_MAF01/ConcordanceTables_macro_", exp.dir,"_HaploGrep.png"), plot = macro_k_hap_box, units = "mm", width = 297, height = 210, dpi = 300)

info_k_hap_box = ggplot(main_khap_df, aes(x = sub_experiment, y = info_score)) +
  geom_violin(fill = "#feb600", na.rm = T) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428") +
  theme_bw() +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                     limits = c(0, 1)) +
  labs(x = "Length of MCMC chain",
       y = expression(paste("IMPUTE2 info score (", italic("r")^"2", ")")),
       title = expression(paste(bold("B."), " Number of included reference haplotypes (k_hap) variations")))
info_k_hap_box
ggsave(filename = paste0("/Volumes/TimMcInerney/MitoImpute/data/HAPLOGROUPS/plots/kHAP_MAF01/ConcordanceTables_info_", exp.dir,"_HaploGrep.png"), plot = info_k_hap_box, units = "mm", width = 297, height = 210, dpi = 300)

## KHAP
# IMPUTED
stat.out.dir = "/Volumes/TimMcInerney/MitoImpute/data/HAPLOGROUPS/stat_tests/kHAP_MAF01/"
l_khap_imp = lm(imputed_match ~ sub_experiment, data = main_khap_df)
#anova(l_khap_imp)
#summary(l_khap_imp)
#emmeans(l_khap_imp, pairwise ~ sub_experiment)
l_khap_imp_s = summary(emmeans(l_khap_imp, pairwise ~ sub_experiment))
write.csv(data.frame(l_khap_imp_s$emmeans), paste0(stat.out.dir, "khap_imputed_emmeans.csv"), quote = F, row.names = F)
write.csv(data.frame(l_khap_imp_s$contrasts), paste0(stat.out.dir, "khap_imputed_contrasts.csv"), quote = F, row.names = F)

# IMPUTED MACRO
l_khap_imp_macro = lm(imputed_macro_match ~ sub_experiment, data = main_khap_df)
#anova(l_khap_imp_macro)
#summary(l_khap_imp_macro)
#emmeans(l_khap_imp_macro, pairwise ~ sub_experiment) # ADD type="response" IF IN LOG SCALE
l_khap_imp_macro_s = summary(emmeans(l_khap_imp_macro, pairwise ~ sub_experiment))
write.csv(data.frame(l_khap_imp_macro_s$emmeans), paste0(stat.out.dir, "khap_imputed_macro_emmeans.csv"), quote = F, row.names = F)
write.csv(data.frame(l_khap_imp_macro_s$contrasts), paste0(stat.out.dir, "khap_imputed_macro_contrasts.csv"), quote = F, row.names = F)

# DIFFERENCE
l_khap_imp_diff = lm(diff ~ sub_experiment, data = main_khap_df)
#anova(l_khap_imp_diff)
#summary(l_khap_imp_diff)
#emmeans(l_khap_imp_diff, pairwise ~ sub_experiment) # ADD type="response" IF IN LOG SCALE
l_khap_imp_diff_s = summary(emmeans(l_khap_imp_diff, pairwise ~ sub_experiment))
write.csv(data.frame(l_khap_imp_diff_s$emmeans), paste0(stat.out.dir, "khap_imputed_diff_emmeans.csv"), quote = F, row.names = F)
write.csv(data.frame(l_khap_imp_diff_s$contrasts), paste0(stat.out.dir, "khap_imputed_diff_contrasts.csv"), quote = F, row.names = F)

# DIFFERENCE MACRO
l_khap_imp_diff_macro = lm(diff_macro ~ sub_experiment, data = main_khap_df)
#anova(l_khap_imp_diff_macro)
#summary(l_khap_imp_diff_macro)
#emmeans(l_khap_imp_diff_macro, pairwise ~ sub_experiment) # ADD type="response" IF IN LOG SCALE
l_khap_imp_diff_macro_s = summary(emmeans(l_khap_imp_diff_macro, pairwise ~ sub_experiment))
write.csv(data.frame(l_khap_imp_diff_macro_s$emmeans), paste0(stat.out.dir, "khap_imputed_macro_diff_emmeans.csv"), quote = F, row.names = F)
write.csv(data.frame(l_khap_imp_diff_macro_s$contrasts), paste0(stat.out.dir, "khap_imputed_macro_diff_contrasts.csv"), quote = F, row.names = F)
