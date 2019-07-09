require(ggplot2)
require(gridExtra)
require(tidyr)
require(emmeans)
require(dplyr)
require(scales)
require(RColorBrewer)

mcmc_tables = read.csv("~/GitCode/MitoImputePrep/metadata/Concordance_tables/HiMC_HaploGrep/combined/ConcordanceTables_MCMC_Experiments_HiMC_COMBINED.csv", header = T)
maf_tables = read.csv("~/GitCode/MitoImputePrep/metadata/Concordance_tables/HiMC_HaploGrep/combined/ConcordanceTables_MAF_Experiments_HiMC_COMBINED.csv", header = T)
khap_tables = read.csv("~/GitCode/MitoImputePrep/metadata/Concordance_tables/HiMC_HaploGrep/combined/ConcordanceTables_kHAP_Experiments_HiMC_COMBINED.csv", header = T)

#plot_dir = "~/Desktop/MitoImpute_plots/HiMC/"
#plot_dir = "/Users/TimMcInerney/Dropbox/University/2019/ASMR_June2019/MCC/"
plot_dir = "/Volumes/TimMcInerney/MitoImpute/figures/HiMC/"

mcmc.var = c("MCMC1", "MCMC5", "MCMC10", "MCMC20", "MCMC30")
mcmc_tables$sub_experiment = factor(mcmc_tables$sub_experiment, levels = mcmc.var)

khap.var = c("kHAP100", "kHAP250", "kHAP500", "kHAP1000", "kHAP2500", "kHAP5000", "kHAP10000", "kHAP20000", "kHAP30000")
khap_tables$sub_experiment = factor(khap_tables$sub_experiment, levels = khap.var)

maf.var = c("MAF1%", "MAF0.5%", "MAF0.1%")
maf_tables$sub_experiment = factor(maf_tables$sub_experiment, levels = maf.var)

## MAF

ggplot(maf_tables, aes(x = sub_experiment, y = typed_match)) +
  geom_violin(fill = "#feb600", na.rm = T, lwd = rel(1/2)) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428", lwd = rel(1/2), fatten = rel(2.5)) +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                     limits = c(0, 1)) +
  labs(x = "Minor allele frequency (MAF)",
       y = "Proportion of haplogroups matching truth set",
       title = expression(paste(bold("Genotyped mtSNVs only"), ""))) +
  theme(axis.text.x = element_text(hjust = 1.0, vjust = 1.0, angle = 45, size = rel(1.0)),
        axis.text.y = element_text(size = rel(1.25)),
        axis.title.y = element_text(size = rel(1.0)),
        plot.title = element_text(size = rel(1.0)),
        #plot.title = element_blank(),
        panel.grid.major = element_line(colour = "black", size = rel(1/2)), 
        panel.grid.minor = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        plot.background = element_rect(fill = "transparent",colour = NA))
ggsave(paste0(plot_dir, "HiMC_MAF_typedMatch_boxplot_tp.png"), bg = "transparent", width = 24.5*2, height = 11*2, units = "cm", dpi = 300)
ggsave(paste0(plot_dir, "HiMC_MAF_typedMatch_boxplot.png"),    width = 24.5*2, height = 11*2, units = "cm", dpi = 300)

ggplot(maf_tables, aes(x = sub_experiment, y = imputed_match)) +
  geom_violin(fill = "#feb600", na.rm = T, lwd = rel(1/2)) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428", lwd = rel(1/2), fatten = rel(2.5)) +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                     limits = c(0, 1)) +
  labs(x = "Minor allele frequency (MAF)",
       y = "Proportion of haplogroups matching truth set",
       title = expression(paste(bold("Genotyped + imputed mtSNVs"), ""))) +
  theme(axis.text.x = element_text(hjust = 1.0, vjust = 1.0, angle = 45, size = rel(1.0)),
        axis.text.y = element_text(size = rel(1.25)),
        axis.title.y = element_text(size = rel(1.0)),
        plot.title = element_text(size = rel(1.0)),
        #plot.title = element_blank(),
        panel.grid.major = element_line(colour = "black", size = rel(1/2)), 
        panel.grid.minor = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        plot.background = element_rect(fill = "transparent",colour = NA))
ggsave(paste0(plot_dir, "HiMC_MAF_imputedMatch_boxplot_tp.png"), bg = "transparent", width = 297, height = 210, units = "mm", dpi = 300)
ggsave(paste0(plot_dir, "HiMC_MAF_imputedMatch_boxplot.png"), width = 297, height = 210, units = "mm", dpi = 300)

## KHAP

ggplot(khap_tables, aes(x = sub_experiment, y = typed_match)) +
  geom_violin(fill = "#feb600", na.rm = T, lwd = rel(1/2)) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428", lwd = rel(1/2), fatten = rel(2.5)) +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                     limits = c(0, 1)) +
  labs(x = bquote('Number of reference haplotypes (k'[HAP]~')'),
       y = "Proportion of haplogroups matching truth set",
       title = expression(paste(bold("Genotyped mtSNVs only"), ""))) +
  theme(axis.text.x = element_text(hjust = 1.0, vjust = 1.0, angle = 45, size = rel(1.0)),
        axis.text.y = element_text(size = rel(1.0)),
        axis.title.y = element_text(size = rel(1.0)),
        plot.title = element_text(size = rel(1.0)),
        #plot.title = element_blank(),
        panel.grid.major = element_line(colour = "black", size = rel(1/2)), 
        panel.grid.minor = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        plot.background = element_rect(fill = "transparent",colour = NA))
ggsave(paste0(plot_dir, "HiMC_kHAP_typedMatch_boxplot_tp.png"), bg = "transparent", width = 297, height = 210, units = "mm", dpi = 300)
ggsave(paste0(plot_dir, "HiMC_kHAP_typedMatch_boxplot.png"), width = 297, height = 210, units = "mm", dpi = 300)

ggplot(khap_tables, aes(x = sub_experiment, y = imputed_match)) +
  geom_violin(fill = "#feb600", na.rm = T, lwd = rel(1/2)) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428", lwd = rel(1/2), fatten = rel(2.5)) +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                     limits = c(0, 1)) +
  labs(x = bquote('Number of reference haplotypes (k'[HAP]~')'),
       y = "Proportion of haplogroups matching truth set",
       title = expression(paste(bold("Genotyped + imputed mtSNVs"), ""))) +
  theme(axis.text.x = element_text(hjust = 1.0, vjust = 1.0, angle = 45, size = rel(1.0)),
        axis.text.y = element_text(size = rel(1.0)),
        axis.title.y = element_text(size = rel(1.0)),
        plot.title = element_text(size = rel(1.0)),
        #plot.title = element_blank(),
        panel.grid.major = element_line(colour = "black", size = rel(1/2)), 
        panel.grid.minor = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        plot.background = element_rect(fill = "transparent",colour = NA))
ggsave(paste0(plot_dir, "HiMC_kHAP_imputedMatch_boxplot_tp.png"), bg = "transparent", width = 297, height = 210, units = "mm", dpi = 300)
ggsave(paste0(plot_dir, "HiMC_kHAP_imputedMatch_boxplot.png"), width = 297, height = 210, units = "mm", dpi = 300)

## MCMC 

ggplot(mcmc_tables, aes(x = sub_experiment, y = typed_match)) +
  geom_violin(fill = "#feb600", na.rm = T, lwd = rel(1/2)) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428", lwd = rel(1/2), fatten = rel(2.5)) +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                     limits = c(0, 1)) +
  labs(x = "Length of Markov chain Monte Carlo (MCMC)",
       y = "Proportion of haplogroups matching truth set",
       title = expression(paste(bold("Genotyped mtSNVs only"), ""))) +
  theme(axis.text.x = element_text(hjust = 1.0, vjust = 1.0, angle = 45, size = rel(1.0)),
        axis.text.y = element_text(size = rel(1.0)),
        axis.title.y = element_text(size = rel(1.0)),
        plot.title = element_text(size = rel(1.0)),
        #plot.title = element_blank(),
        panel.grid.major = element_line(colour = "black", size = rel(1/2)), 
        panel.grid.minor = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        plot.background = element_rect(fill = "transparent",colour = NA))
ggsave(paste0(plot_dir, "HiMC_MCMC_typedMatch_boxplot_tp.png"), bg = "transparent", width = 297, height = 210, units = "mm", dpi = 300)
ggsave(paste0(plot_dir, "HiMC_MCMC_typedMatch_boxplot.png"), width = 297, height = 210, units = "mm", dpi = 300)

ggplot(mcmc_tables, aes(x = sub_experiment, y = imputed_match)) +
  geom_violin(fill = "#feb600", na.rm = T, lwd = rel(1/2)) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428", lwd = rel(1/2), fatten = rel(2.5)) +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                     limits = c(0, 1)) +
  labs(x = "Length of Markov chain Monte Carlo (MCMC)",
       y = "Proportion of haplogroups matching truth set",
       title = expression(paste(bold("Genotyped + imputed mtSNVs"), ""))) +
  theme(axis.text.x = element_text(hjust = 1.0, vjust = 1.0, angle = 45, size = rel(1.0)),
        axis.text.y = element_text(size = rel(1.0)),
        axis.title.y = element_text(size = rel(1.0)),
        plot.title = element_text(size = rel(1.0)),
        #plot.title = element_blank(),
        panel.grid.major = element_line(colour = "black", size = rel(1/2)), 
        panel.grid.minor = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        plot.background = element_rect(fill = "transparent",colour = NA))
ggsave(paste0(plot_dir, "HiMC_MCMC_imputedMatch_boxplot_tp.png"), bg = "transparent", width = 297, height = 210, units = "mm", dpi = 300)
ggsave(paste0(plot_dir, "HiMC_MCMC_imputedMatch_boxplot.png"), width = 297, height = 210, units = "mm", dpi = 300)
















## FOR PUBLICATION:
# MAF
MAF_imp = read.csv("/Volumes/TimMcInerney/MitoImpute/analyses/concordance/MCC/MAF/ConcordanceTables_MAF_Experiments_MCC_imputed_genotype.csv", header = T) # LOAD INFO FILE

# Fix the data structure!
maf_tmp_typ = maf_tables[,1:4]
maf_tmp_typ$n_snps = MAF_imp$n.snps
maf_tmp_typ$hg.conc = maf_tables$typed_match
maf_tmp_typ$label = "Genotyped_only"

maf_tmp_imp = maf_tables[,1:4]
maf_tmp_imp$n_snps = MAF_imp$n.snps
maf_tmp_imp$hg.conc = maf_tables$imputed_match
maf_tmp_imp$label = "Genotyped_Imputed"

maf_tables_changed = rbind(maf_tmp_typ, maf_tmp_imp)

geno.var = c("Genotyped_only", "Genotyped_Imputed")
maf_tables_changed$label = factor(maf_tables_changed$label, levels = geno.var)
pd = position_dodge(width = 1.5)
maf_hg_changed = ggplot(maf_tables_changed, aes(y = hg.conc, fill = label)) +
  #geom_violin(fill = "#feb600", na.rm = T, lwd = rel(1/2)) +
  #geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428", lwd = rel(1/2), fatten = rel(2.5)) +
  stat_boxplot(geom = "errorbar", na.rm = T, lwd = (3/4), position=pd) +
  geom_boxplot(width = rel(1.5), notch = T, na.rm = T, lwd = rel(1/2), fatten = rel(2.5), outlier.colour = "#802428", position=pd, outlier.size = 2) +
  #geom_jitter(position=position_jitter(0.25), aes(size = n_snps, colour = label), na.rm = T, shape = 21, fill = NA, stroke = rel(0.5)) + #colour = "#f3e5b1"
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                     limits = c(0, 1)) +
  scale_x_discrete(breaks = NULL) +
  scale_fill_manual(values = c("#feb600", "#ea4e3c"),
                    name = "",
                    #name = expression(paste(italic("in silico"), " microarray")),
                    breaks = c("Genotyped_only", "Genotyped_Imputed"),
                    labels = c("Genotyped only", "Genotyped + Imputed")) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Minor allele frequency (MAF)",
       y = "Haplogroup concordance",
       title = expression(paste(bold("A."), " Genotyped mtSNVs only", "")),
       fill = "") +
  theme(axis.text.x = element_blank(), #axis.text.x = element_text(hjust = 1.0, vjust = 1.0, angle = 45, size = rel(1.0)),
        axis.text.y = element_text(size = rel(1.25)),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = rel(1.25), vjust = rel(2.5)),
        #plot.title = element_text(size = rel(1.0)),
        plot.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = rel(1.0)),
        legend.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.major = element_line(colour = alpha("black", 0.5), linetype = 2, size = rel(1/2)),
        panel.grid.minor = element_line(colour = alpha("black", 0.5), linetype = 2, size = rel(1/2)),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        plot.background = element_rect(fill = "transparent",colour = NA)) +
  facet_wrap(~sub_experiment, nrow = 1)
  #facet_grid(cols = vars(sub_experiment), space = "fixed")
maf_hg_changed
ggsave(paste0(plot_dir, "ForPublication/HiMC_MAF_HaplogroupConcordance_tp.png"), plot = maf_hg_changed, bg = "transparent", width = 24.5*2, height = 11*2, units = "cm", dpi = 300)
ggsave(paste0(plot_dir, "ForPublication/HiMC_MAF_HaplogroupConcordance.png"), plot = maf_hg_changed, width = 24.5*2, height = 11*2, units = "cm", dpi = 300)

# KHAP
KHAP_imp = read.csv("/Volumes/TimMcInerney/MitoImpute/analyses/concordance/MCC/KHAP/ConcordanceTables_KHAP_Experiments_MCC_imputed_genotype.csv", header = T) # LOAD INFO FILE

# Fix the data structure!
khap_tmp_typ = khap_tables[,1:4]
khap_tmp_typ$n_snps = KHAP_imp$n.snps
khap_tmp_typ$hg.conc = khap_tables$typed_match
khap_tmp_typ$label = "Genotyped_only"

khap_tmp_imp = khap_tables[,1:4]
khap_tmp_imp$n_snps = KHAP_imp$n.snps
khap_tmp_imp$hg.conc = khap_tables$imputed_match
khap_tmp_imp$label = "Genotyped_Imputed"

khap_tables_changed = rbind(khap_tmp_typ, khap_tmp_imp)

geno.var = c("Genotyped_only", "Genotyped_Imputed")
khap_tables_changed$label = factor(khap_tables_changed$label, levels = geno.var)

khap_hg_changed = ggplot(khap_tables_changed, aes(y = hg.conc, fill = label)) +
  #geom_violin(fill = "#feb600", na.rm = T, lwd = rel(1/2)) +
  #geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428", lwd = rel(1/2), fatten = rel(2.5)) +
  stat_boxplot(geom = "errorbar", na.rm = T, lwd = (3/4)) +
  geom_boxplot(notch = T, na.rm = T, lwd = rel(1/2), fatten = rel(2.5), outlier.colour = "#802428", outlier.size = 2) +
  #geom_jitter(position=position_jitter(0.25), aes(size = n_snps, colour = label), na.rm = T, shape = 21, fill = NA, stroke = rel(0.5)) + #colour = "#f3e5b1"
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                     limits = c(0, 1)) +
  scale_x_continuous(breaks = NULL) +
  scale_fill_manual(values = c("#feb600", "#ea4e3c"),
                    name = "",
                    #name = expression(paste(italic("in silico"), " microarray")),
                    breaks = c("Genotyped_only", "Genotyped_Imputed"),
                    labels = c("Genotyped only", "Genotyped + Imputed")) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Minor allele frequency (MAF)",
       y = "Haplogroup concordance",
       title = expression(paste(bold("A."), " Genotyped mtSNVs only", "")),
       fill = "") +
  theme(axis.text.x = element_blank(), #axis.text.x = element_text(hjust = 1.0, vjust = 1.0, angle = 45, size = rel(1.0)),
        axis.text.y = element_text(size = rel(1.25)),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = rel(1.25), vjust = rel(2.5)),
        #plot.title = element_text(size = rel(1.0)),
        plot.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = rel(1.0)),
        legend.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.major = element_line(colour = alpha("black", 0.5), linetype = 2, size = rel(1/2)),
        panel.grid.minor = element_line(colour = alpha("black", 0.5), linetype = 2, size = rel(1/2)),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        plot.background = element_rect(fill = "transparent",colour = NA)) +
  facet_grid(cols = vars(sub_experiment))
khap_hg_changed
ggsave(paste0(plot_dir, "ForPublication/HiMC_KHAP_HaplogroupConcordance_tp.png"), plot = khap_hg_changed, bg = "transparent", width = 24.5*2, height = 11*2, units = "cm", dpi = 300)
ggsave(paste0(plot_dir, "ForPublication/HiMC_KHAP_HaplogroupConcordance.png"), plot = khap_hg_changed, width = 24.5*2, height = 11*2, units = "cm", dpi = 300)

# MCMC
MCMC_imp = read.csv("/Volumes/TimMcInerney/MitoImpute/analyses/concordance/MCC/MCMC/ConcordanceTables_MCMC_Experiments_MCC_imputed_genotype.csv", header = T) # LOAD INFO FILE

# Fix the data structure!
mcmc_tmp_typ = mcmc_tables[,1:4]
mcmc_tmp_typ$n_snps = MCMC_imp$n.snps
mcmc_tmp_typ$hg.conc = mcmc_tables$typed_match
mcmc_tmp_typ$label = "Genotyped_only"

mcmc_tmp_imp = mcmc_tables[,1:4]
mcmc_tmp_imp$n_snps = MCMC_imp$n.snps
mcmc_tmp_imp$hg.conc = mcmc_tables$imputed_match
mcmc_tmp_imp$label = "Genotyped_Imputed"

mcmc_tables_changed = rbind(mcmc_tmp_typ, mcmc_tmp_imp)

geno.var = c("Genotyped_only", "Genotyped_Imputed")
mcmc_tables_changed$label = factor(mcmc_tables_changed$label, levels = geno.var)

mcmc_hg_changed = ggplot(mcmc_tables_changed, aes(y = hg.conc, fill = label)) +
  #geom_violin(fill = "#feb600", na.rm = T, lwd = rel(1/2)) +
  #geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428", lwd = rel(1/2), fatten = rel(2.5)) +
  stat_boxplot(geom = "errorbar", na.rm = T, lwd = (3/4)) +
  geom_boxplot(notch = T, na.rm = T, lwd = rel(1/2), fatten = rel(2.5), outlier.colour = "#802428", outlier.size = 2) +
  #geom_jitter(position=position_jitter(0.25), aes(size = n_snps, colour = label), na.rm = T, shape = 21, fill = NA, stroke = rel(0.5)) + #colour = "#f3e5b1"
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                     limits = c(0, 1)) +
  scale_x_continuous(breaks = NULL) +
  scale_fill_manual(values = c("#feb600", "#ea4e3c"),
                    name = "",
                    #name = expression(paste(italic("in silico"), " microarray")),
                    breaks = c("Genotyped_only", "Genotyped_Imputed"),
                    labels = c("Genotyped only", "Genotyped + Imputed")) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Minor allele frequency (MAF)",
       y = "Haplogroup concordance",
       title = expression(paste(bold("A."), " Genotyped mtSNVs only", "")),
       fill = "") +
  theme(axis.text.x = element_blank(), #axis.text.x = element_text(hjust = 1.0, vjust = 1.0, angle = 45, size = rel(1.0)),
        axis.text.y = element_text(size = rel(1.25)),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = rel(1.25), vjust = rel(2.5)),
        #plot.title = element_text(size = rel(1.0)),
        plot.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = rel(1.0)),
        legend.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.major = element_line(colour = alpha("black", 0.5), linetype = 2, size = rel(1/2)),
        panel.grid.minor = element_line(colour = alpha("black", 0.5), linetype = 2, size = rel(1/2)),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        plot.background = element_rect(fill = "transparent",colour = NA)) +
  facet_grid(cols = vars(sub_experiment))
mcmc_hg_changed
ggsave(paste0(plot_dir, "ForPublication/HiMC_MCMC_HaplogroupConcordance_tp.png"), plot = mcmc_hg_changed, bg = "transparent", width = 24.5*2, height = 11*2, units = "cm", dpi = 300)
ggsave(paste0(plot_dir, "ForPublication/HiMC_MCMC_HaplogroupConcordance.png"), plot = mcmc_hg_changed, width = 24.5*2, height = 11*2, units = "cm", dpi = 300)

