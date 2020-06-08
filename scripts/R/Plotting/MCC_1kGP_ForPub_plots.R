require(ggplot2)
require(gridExtra)
require(tidyr)
require(emmeans)
require(dplyr)
require(scales)
require(RColorBrewer)

container = "/Volumes/TimMcInerney/MitoImpute/analyses/concordance/MCC/"
plot_dir = "/Volumes/TimMcInerney/MitoImpute/figures/MCC/"

mcmc.var = c("MCMC1", "MCMC5", "MCMC10", "MCMC20", "MCMC30")
mcmc_tables$sub_experiment = factor(mcmc_tables$sub_experiment, levels = mcmc.var)

khap.var = c("kHAP100", "kHAP250", "kHAP500", "kHAP1000", "kHAP2500", "kHAP5000", "kHAP10000", "kHAP20000", "kHAP30000")
khap_tables$sub_experiment = factor(khap_tables$sub_experiment, levels = khap.var)

maf.var = c("MAF1%", "MAF0.5%", "MAF0.1%")
maf_tables$sub_experiment = factor(maf_tables$sub_experiment, levels = maf.var)

## FOR PUBLICATION:
# MAF
# MAF

MAF_imp = read.csv(paste0(container, "MAF/ConcordanceTables_MAF_Experiments_MCC_imputed_genotype.csv"), header = T)
MAF_typ = read.csv(paste0(container, "MAF/ConcordanceTables_MAF_Experiments_MCC_typed_genotype.csv"), header = T)
MAF_imp$version = "imputed"
MAF_typ$version = "genotyped"
#main_maf_df = rbind(MAF_imp, MAF_typ)
main_maf_df = MAF_imp

exp.dir = "MAF_Experiments"
exp.var = c("MAF1%", "MAF0.5%", "MAF0.1%")
main_maf_df$sub_experiment = factor(main_maf_df$sub_experiment, levels = exp.var)

maf_mcc_pub = ggplot(main_maf_df, aes(y = mcc)) +
  #geom_violin(fill = "#feb600", na.rm = T, lwd = rel(1/2)) +
  stat_boxplot(geom = "errorbar", na.rm = T, width = (1/2), lwd = (3/4)) +
  geom_boxplot(notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428", lwd = rel(1/2), fatten = rel(2.5), outlier.size = 2) +
  #geom_jitter(position=position_jitter(0.25), aes(size = n.snps), na.rm = T, shape = 21, fill = NA, colour = "#802428", stroke = rel(0.5)) + #colour = "#f3e5b1"
  #theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = rel(1.5)),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = rel(1.5), vjust = rel(2.5)),
        strip.text.x = element_text(size = rel(1.5)),
        #plot.title = element_text(size = rel(1.0)),
        plot.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = rel(1.0)),
        legend.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.major = element_line(colour = alpha("black", 0.5), linetype = 2, size = rel(1/2)),
        panel.grid.minor = element_line(colour = alpha("black", 0.5), linetype = 2, size = rel(1/2)),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        plot.background = element_rect(fill = "transparent",colour = NA)) +
  scale_y_continuous(breaks = seq(-1.0, 1.0, by = 0.2),
                     limits = c(-1, 1)) +
  scale_x_discrete(breaks = NULL) +
  labs(x = bquote('Number of reference haplotypes (k'[HAP]~')'),
       #x = "Length of KHAP chain",
       y = "Matthew's Correlation Coefficient",
       title = "") +
  facet_wrap(~sub_experiment, nrow = 1)
maf_mcc_pub
ggsave(paste0(plot_dir, "ForPublication/MCC_MAF_HaplogroupConcordance_tp.png"), plot = maf_mcc_pub, bg = "transparent", width = 24.5*2, height = 11*2, units = "cm", dpi = 300)
ggsave(paste0(plot_dir, "ForPublication/MCC_MAF_HaplogroupConcordance.png"), plot = maf_mcc_pub, width = 24.5*2, height = 11*2, units = "cm", dpi = 300)

# KHAP
KHAP_imp = read.csv(paste0(container, "KHAP/ConcordanceTables_KHAP_Experiments_MCC_imputed_genotype.csv"), header = T)
KHAP_typ = read.csv(paste0(container, "KHAP/ConcordanceTables_KHAP_Experiments_MCC_typed_genotype.csv"), header = T)
KHAP_imp$version = "imputed"
KHAP_typ$version = "genotyped"
#main_khap_df = rbind(KHAP_imp, KHAP_typ)
main_khap_df = KHAP_imp

exp.dir = "kHAP_Experiments"
exp.var = c("kHAP100", "kHAP250", "kHAP500", "kHAP1000", "kHAP2500", "kHAP5000", "kHAP10000", "kHAP20000", "kHAP30000")
main_khap_df$sub_experiment = factor(main_khap_df$sub_experiment, levels = exp.var)

khap_mcc_pub = ggplot(main_khap_df, aes(y = mcc)) +
  #geom_violin(fill = "#feb600", na.rm = T, lwd = rel(1/2)) +
  stat_boxplot(geom = "errorbar", na.rm = T, width = (1/2), lwd = (3/4)) +
  geom_boxplot(notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428", lwd = rel(1/2), fatten = rel(2.5), outlier.size = 2) +
  #geom_jitter(position=position_jitter(0.25), aes(size = n.snps), na.rm = T, shape = 21, fill = NA, colour = "#802428", stroke = rel(0.5)) + #colour = "#f3e5b1"
  #theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = rel(1.5)),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = rel(1.5), vjust = rel(2.5)),
        strip.text.x = element_text(size = rel(1.5)),
        #plot.title = element_text(size = rel(1.0)),
        plot.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = rel(1.0)),
        legend.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.major = element_line(colour = "black", size = rel(1/2)), 
        panel.grid.minor = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        plot.background = element_rect(fill = "transparent",colour = NA)) +
  scale_y_continuous(breaks = seq(-1.0, 1.0, by = 0.2),
                     limits = c(-1, 1)) +
  scale_x_discrete(breaks = NULL) +
  labs(x = bquote('Number of reference haplotypes (k'[HAP]~')'),
       #x = "Length of KHAP chain",
       y = "Matthew's Correlation Coefficient",
       title = "") +
  facet_wrap(~sub_experiment, nrow = 1)
khap_mcc_pub
ggsave(paste0(plot_dir, "ForPublication/MCC_KHAP_HaplogroupConcordance_tp.png"), plot = khap_mcc_pub, bg = "transparent", width = 24.5*2, height = 11*2, units = "cm", dpi = 300)
ggsave(paste0(plot_dir, "ForPublication/MCC_KHAP_HaplogroupConcordance.png"), plot = khap_mcc_pub, width = 24.5*2, height = 11*2, units = "cm", dpi = 300)

# MCMC
MCMC_imp = read.csv(paste0(container, "MCMC/ConcordanceTables_MCMC_Experiments_MCC_imputed_genotype.csv"), header = T)
MCMC_typ = read.csv(paste0(container, "MCMC/ConcordanceTables_MCMC_Experiments_MCC_typed_genotype.csv"), header = T)
MCMC_imp$version = "imputed"
MCMC_typ$version = "genotyped"
#main_mcmc_df = rbind(MCMC_imp, MCMC_typ)
main_mcmc_df = MCMC_imp

exp.dir = "MCMC"
exp.var = c("MCMC1", "MCMC5", "MCMC10", "MCMC20", "MCMC30")
main_mcmc_df$sub_experiment = factor(mcmc_tables$sub_experiment, levels = exp.var)

mcmc_mcc_pub = ggplot(main_mcmc_df, aes(y = mcc)) +
  #geom_violin(fill = "#feb600", na.rm = T, lwd = rel(1/2)) +
  stat_boxplot(geom = "errorbar", na.rm = T, width = (1/2), lwd = (3/4)) +
  geom_boxplot(notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428", lwd = rel(1/2), fatten = rel(2.5), outlier.size = 2) +
  #geom_jitter(position=position_jitter(0.25), aes(size = n.snps), na.rm = T, shape = 21, fill = NA, colour = "#802428", stroke = rel(0.5)) + #colour = "#f3e5b1"
  #theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = rel(1.5)),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = rel(1.5), vjust = rel(2.5)),
        strip.text.x = element_text(size = rel(1.5)),
        #plot.title = element_text(size = rel(1.0)),
        plot.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = rel(1.0)),
        legend.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.major = element_line(colour = "black", size = rel(1/2)), 
        panel.grid.minor = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        plot.background = element_rect(fill = "transparent",colour = NA)) +
  scale_y_continuous(breaks = seq(-1.0, 1.0, by = 0.2),
                     limits = c(-1, 1)) +
  scale_x_discrete(breaks = NULL) +
  labs(x = bquote('Number of reference haplotypes (k'[HAP]~')'),
       #x = "Length of MCMC chain",
       y = "Matthew's Correlation Coefficient",
       title = "") +
  facet_wrap(~sub_experiment, nrow = 1)
mcmc_mcc_pub
ggsave(paste0(plot_dir, "ForPublication/MCC_MCMC_HaplogroupConcordance_tp.png"), plot = mcmc_mcc_pub, bg = "transparent", width = 24.5*2, height = 11*2, units = "cm", dpi = 300)
ggsave(paste0(plot_dir, "ForPublication/MCC_MCMC_HaplogroupConcordance.png"), plot = mcmc_mcc_pub, width = 24.5*2, height = 11*2, units = "cm", dpi = 300)
