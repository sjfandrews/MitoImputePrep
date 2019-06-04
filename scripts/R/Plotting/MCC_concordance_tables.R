require(ggplot2)
require(gridExtra)
require(tidyr)
require(emmeans)
require(dplyr)
require(scales)

container = "/Volumes/TimMcInerney/MitoImpute/analyses/concordance/MCC/"

# MCMC

MCMC_imp = read.csv(paste0(container, "MCMC/ConcordanceTables_MCMC_Experiments_MCC_imputed_genotype.csv"), header = T)
MCMC_typ = read.csv(paste0(container, "MCMC/ConcordanceTables_MCMC_Experiments_MCC_typed_genotype.csv"), header = T)
MCMC_imp$version = "imputed"
MCMC_typ$version = "genotyped"
main_mcmc_df = rbind(MCMC_imp, MCMC_typ)
main_mcmc_df = MCMC_imp

exp.dir = "MCMC_Experiments"
exp.var = c("MCMC1", "MCMC5", "MCMC10", "MCMC20", "MCMC30")
main_mcmc_df$sub_experiment = factor(main_mcmc_df$sub_experiment, levels = exp.var)

mcmc_mcc_info_regression = ggplot(main_mcmc_df, aes(x = mcc, y = info, colour = sub_experiment)) +
  geom_point() +
  theme(axis.text.x = element_text(hjust = 1.0, angle = 45),
        axis.title.x = element_text(vjust = -1.0),
        legend.position = "bottom",
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.major = element_line(colour = "black", size = 0.25), 
        panel.grid.minor = element_line(colour = "black", linetype = 2, size = 0.25),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        plot.background = element_rect(fill = "transparent",colour = NA)) +
  labs(y = expression(paste("IMPUTE2 info score (", italic("r")^"2", ")")),
       x = "Matthew's Correlation Coefficient")
ggsave(filename = paste0(container, "Plots/MCMC_MCCvINFO_regression.png"), plot = mcmc_mcc_info_regression, width = 297, height = 210, units = "mm", dpi = 300)
ggsave(filename = paste0(container, "Plots/MCMC_MCCvINFO_regression_tp.png"), bg = "transparent", plot = mcmc_mcc_info_regression, width = 297, height = 210, units = "mm", dpi = 300)

mcmc_mcc_box = ggplot(main_mcmc_df, aes(x = sub_experiment, y = mcc)) +
  geom_violin(fill = "#feb600", na.rm = T, lwd = rel(1/2)) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428", lwd = rel(1/2), fatten = rel(2.5)) +
  ##geom_jitter(position=position_jitter(0.25), aes(size = n.snps), na.rm = T, shape = 21, fill = NA, colour = "#802428", stroke = rel(0.5)) + #colour = "#f3e5b1"
  #theme_bw() +
  theme(axis.text.x = element_text(hjust = 1.0, angle = 45),
        axis.title.x = element_text(vjust = -1.0),
        legend.position = "bottom",
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.major = element_line(colour = "black", size = 0.25), 
        panel.grid.minor = element_line(colour = "black", linetype = 2, size = 0.25),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        plot.background = element_rect(fill = "transparent",colour = NA)) +
  scale_y_continuous(breaks = seq(-1.0, 1.0, by = 0.2),
                     limits = c(-1, 1)) +
  labs(x = "",
       #x = "Length of MCMC chain",
       y = "Matthew's Correlation Coefficient",
       title = expression(paste(bold("A."), "")))
mcmc_mcc_box
ggsave(filename = paste0(container, "Plots/MCMC_MCC.png"), plot = mcmc_mcc_box, width = 297, height = 210, units = "mm", dpi = 300)
ggsave(filename = paste0(container, "Plots/MCMC_MCC_tp.png"), bg = "transparent", plot = mcmc_mcc_box, width = 297, height = 210, units = "mm", dpi = 300)

mcmc_info_box = ggplot(main_mcmc_df, aes(x = sub_experiment, y = info)) +
  geom_violin(fill = "#feb600", na.rm = T, lwd = rel(1/2)) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428", lwd = rel(1/2), fatten = rel(2.5)) +
  #theme_bw() +
  theme(axis.text.x = element_text(hjust = 1.0, angle = 45),
        axis.title.x = element_text(vjust = -1.0),
        legend.position = "bottom",
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.major = element_line(colour = "black", size = 0.25), 
        panel.grid.minor = element_line(colour = "black", linetype = 2, size = 0.25),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        plot.background = element_rect(fill = "transparent",colour = NA)) +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                     limits = c(0, 1)) +
  labs(x = "",
       #x = "Length of MCMC chain",
       y = expression(paste("IMPUTE2 info score (", italic("r")^"2", ")")),
       #title = expression(paste(bold("B."), " Markov chain Monte Carlo (MCMC) length variations")))
       title = expression(paste(bold("B."), "")))
mcmc_info_box
ggsave(filename = paste0(container, "Plots/MCMC_info.png"), plot = mcmc_info_box, width = 297, height = 210, units = "mm", dpi = 300)
ggsave(filename = paste0(container, "Plots/MCMC_info_tp.png"), bg = "transparent", plot = mcmc_info_box, width = 297, height = 210, units = "mm", dpi = 300)

mcmc_infoComb_box = ggplot(main_mcmc_df, aes(x = sub_experiment, y = info_comb)) +
  geom_violin(fill = "#feb600", na.rm = T, lwd = rel(1/2)) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428", lwd = rel(1/2), fatten = rel(2.5)) +
  #theme_bw() +
  theme(axis.text.x = element_text(hjust = 1.0, angle = 45),
        axis.title.x = element_text(vjust = -1.0),
        legend.position = "bottom",
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.major = element_line(colour = "black", size = 0.25), 
        panel.grid.minor = element_line(colour = "black", linetype = 2, size = 0.25),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        plot.background = element_rect(fill = "transparent",colour = NA)) +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                     limits = c(0, 1)) +
  labs(x = "",
       #x = "Length of MCMC chain",
       y = expression(paste("IMPUTE2 info score (", italic("r")^"2", ")")),
       title = expression(paste(bold("A."), " Markov chain Monte Carlo (MCMC) length variations")))
mcmc_infoComb_box
ggsave(filename = paste0(container, "Plots/MCMC_infoComb.png"), plot = mcmc_infoComb_box, width = 297, height = 210, units = "mm", dpi = 300)
ggsave(filename = paste0(container, "Plots/MCMC_infoComb_tp.png"), bg = "transparent", plot = mcmc_infoComb_box, width = 297, height = 210, units = "mm", dpi = 300)

mcmc_cert_box = ggplot(main_mcmc_df, aes(x = sub_experiment, y = certainty)) +
  geom_violin(fill = "#feb600", na.rm = T, lwd = rel(1/2)) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428", lwd = rel(1/2), fatten = rel(2.5)) +
  #theme_bw() +
  theme(axis.text.x = element_text(hjust = 1.0, angle = 45),
        axis.title.x = element_text(vjust = -1.0),
        legend.position = "bottom",
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.major = element_line(colour = "black", size = 0.25), 
        panel.grid.minor = element_line(colour = "black", linetype = 2, size = 0.25),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        plot.background = element_rect(fill = "transparent",colour = NA)) +
  labs(x = "",
       #x = "Length of MCMC chain",
       y = "Certainty",
       title = expression(paste(bold("C."), "")))
mcmc_cert_box
ggsave(filename = paste0(container, "Plots/MCMC_cert.png"), plot = mcmc_cert_box, width = 297, height = 210, units = "mm", dpi = 300)
ggsave(filename = paste0(container, "Plots/MCMC_cert_tp.png"), bg = "transparent", plot = mcmc_cert_box, width = 297, height = 210, units = "mm", dpi = 300)

mcmc_conc_box = ggplot(main_mcmc_df, aes(x = sub_experiment, y = concordance)) +
  geom_violin(fill = "#feb600", na.rm = T, lwd = rel(1/2)) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428", lwd = rel(1/2), fatten = rel(2.5)) +
  #geom_jitter(position=position_jitter(0.25), aes(size = n.snps), na.rm = T, shape = 21, fill = NA, colour = "#802428", stroke = rel(0.5)) + #colour = "#f3e5b1"
  #theme_bw() +
  theme(axis.text.x = element_text(hjust = 1.0, angle = 45),
        axis.title.x = element_text(vjust = -1.0),
        legend.position = "bottom",
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.major = element_line(colour = "black", size = 0.25), 
        panel.grid.minor = element_line(colour = "black", linetype = 2, size = 0.25),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        plot.background = element_rect(fill = "transparent",colour = NA)) +
  labs(x = "",
       #x = "Length of MCMC chain",
       y = "Concordance",
       title = expression(paste(bold("D."), "")))
mcmc_conc_box
ggsave(filename = paste0(container, "Plots/MCMC_conc.png"), plot = mcmc_conc_box, width = 297, height = 210, units = "mm", dpi = 300)
ggsave(filename = paste0(container, "Plots/MCMC_conc_tp.png"), bg = "transparent", plot = mcmc_conc_box, width = 297, height = 210, units = "mm", dpi = 300)

mcmc_plots = grid.arrange(top = "Markov chain Monte Carlo (MCMC) length variations",
                          arrangeGrob(mcmc_mcc_box, mcmc_info_box, mcmc_cert_box, mcmc_conc_box, ncol = 2, nrow = 2)
                          )
ggsave(filename = paste0(container, "Plots/MCMC.png"), plot = mcmc_plots, width = 297, height = 210, units = "mm", dpi = 300)
ggsave(filename = paste0(container, "Plots/MCMC_tp.png"), bg = "transparent", plot = mcmc_plots, width = 297, height = 210, units = "mm", dpi = 300)

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

khap_mcc_info_regression = ggplot(main_khap_df, aes(x = mcc, y = info, colour = sub_experiment)) +
  geom_point() +
  theme(axis.text.x = element_text(hjust = 1.0, angle = 45),
        axis.title.x = element_text(vjust = -1.0),
        legend.position = "bottom",
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.major = element_line(colour = "black", size = 0.25), 
        panel.grid.minor = element_line(colour = "black", linetype = 2, size = 0.25),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        plot.background = element_rect(fill = "transparent",colour = NA)) +
  labs(y = expression(paste("IMPUTE2 info score (", italic("r")^"2", ")")),
       x = "Matthew's Correlation Coefficient")
#ggsave(filename = paste0(container, "Plots/MCMC_MCCvINFO_regression.png"), plot = mcmc_mcc_info_regression, width = 297, height = 210, units = "mm", dpi = 300)

khap_mcc_box = ggplot(main_khap_df, aes(x = sub_experiment, y = mcc)) +
  geom_violin(fill = "#feb600", na.rm = T, lwd = rel(1/2)) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428", lwd = rel(1/2), fatten = rel(2.5)) +
  #geom_jitter(position=position_jitter(0.25), aes(size = n.snps), na.rm = T, shape = 21, fill = NA, colour = "#802428", stroke = rel(0.5)) + #colour = "#f3e5b1"
  #theme_bw() +
  theme(axis.text.x = element_text(hjust = 1.0, angle = 45),
        axis.title.x = element_text(vjust = -1.0),
        legend.position = "bottom",
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.major = element_line(colour = "black", size = 0.25), 
        panel.grid.minor = element_line(colour = "black", linetype = 2, size = 0.25),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        plot.background = element_rect(fill = "transparent",colour = NA)) +
  scale_y_continuous(breaks = seq(-1.0, 1.0, by = 0.2),
                     limits = c(-1, 1)) +
  labs(x = "",
       #x = "Length of KHAP chain",
       y = "Matthew's Correlation Coefficient",
       title = expression(paste(bold("A."), "")))
khap_mcc_box
ggsave(filename = paste0(container, "Plots/kHAP_MCC.png"), plot = khap_mcc_box, width = 297, height = 210, units = "mm", dpi = 300)
ggsave(filename = paste0(container, "Plots/kHAP_MCC_tp.png"), bg = "transparent", plot = khap_mcc_box, width = 297, height = 210, units = "mm", dpi = 300)

khap_info_box = ggplot(main_khap_df, aes(x = sub_experiment, y = info)) +
  geom_violin(fill = "#feb600", na.rm = T, lwd = rel(1/2)) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428", lwd = rel(1/2), fatten = rel(2.5)) +
  #theme_bw() +
  theme(axis.text.x = element_text(hjust = 1.0, angle = 45),
        axis.title.x = element_text(vjust = -1.0),
        legend.position = "bottom",
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.major = element_line(colour = "black", size = 0.25), 
        panel.grid.minor = element_line(colour = "black", linetype = 2, size = 0.25),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        plot.background = element_rect(fill = "transparent",colour = NA)) +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                     limits = c(0, 1)) +
  labs(x = "",
       #x = "Length of KHAP chain",
       y = expression(paste("IMPUTE2 info score (", italic("r")^"2", ")")),
       title = expression(paste(bold("B."), "")))
khap_info_box
ggsave(filename = paste0(container, "Plots/kHAP_info.png"), plot = khap_info_box, width = 297, height = 210, units = "mm", dpi = 300)
ggsave(filename = paste0(container, "Plots/kHAP_info_tp.png"), bg = "transparent", plot = khap_info_box, width = 297, height = 210, units = "mm", dpi = 300)

khap_infoComb_box = ggplot(main_khap_df, aes(x = sub_experiment, y = info_comb)) +
  geom_violin(fill = "#feb600", na.rm = T, lwd = rel(1/2)) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428", lwd = rel(1/2), fatten = rel(2.5)) +
  #theme_bw() +
  theme(axis.text.x = element_text(hjust = 1.0, angle = 45),
        axis.title.x = element_text(vjust = -1.0),
        legend.position = "bottom",
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.major = element_line(colour = "black", size = 0.25), 
        panel.grid.minor = element_line(colour = "black", linetype = 2, size = 0.25),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        plot.background = element_rect(fill = "transparent",colour = NA)) +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                     limits = c(0, 1)) +
  labs(x = "",
       #x = "Length of KHAP chain",
       y = expression(paste("IMPUTE2 info score (", italic("r")^"2", ")")),
       title = expression(paste(bold("A."), " Markov chain Monte Carlo (KHAP) length variations")))
khap_infoComb_box
ggsave(filename = paste0(container, "Plots/kHAP_infoComb.png"), plot = khap_infoComb_box, width = 297, height = 210, units = "mm", dpi = 300)
ggsave(filename = paste0(container, "Plots/kHAP_infoComb_tp.png"), bg = "transparent", plot = khap_infoComb_box, width = 297, height = 210, units = "mm", dpi = 300)

khap_cert_box = ggplot(main_khap_df, aes(x = sub_experiment, y = certainty)) +
  geom_violin(fill = "#feb600", na.rm = T, lwd = rel(1/2)) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428", lwd = rel(1/2), fatten = rel(2.5)) +
  #theme_bw() +
  theme(axis.text.x = element_text(hjust = 1.0, angle = 45),
        axis.title.x = element_text(vjust = -1.0),
        legend.position = "bottom",
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.major = element_line(colour = "black", size = 0.25), 
        panel.grid.minor = element_line(colour = "black", linetype = 2, size = 0.25),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        plot.background = element_rect(fill = "transparent",colour = NA)) +
  labs(x = "",
       #x = "Length of KHAP chain",
       y = "Certainty",
       title = expression(paste(bold("C."), "")))
khap_cert_box
ggsave(filename = paste0(container, "Plots/kHAP_cert.png"), plot = khap_cert_box, width = 297, height = 210, units = "mm", dpi = 300)
ggsave(filename = paste0(container, "Plots/kHAP_cert_tp.png"), bg = "transparent", plot = khap_cert_box, width = 297, height = 210, units = "mm", dpi = 300)

khap_conc_box = ggplot(main_khap_df, aes(x = sub_experiment, y = concordance)) +
  geom_violin(fill = "#feb600", na.rm = T, lwd = rel(1/2)) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428", lwd = rel(1/2), fatten = rel(2.5)) +
  #geom_jitter(position=position_jitter(0.25), aes(size = n.snps), na.rm = T, shape = 21, fill = NA, colour = "#802428", stroke = rel(0.5)) + #colour = "#f3e5b1"
  #theme_bw() +
  theme(axis.text.x = element_text(hjust = 1.0, angle = 45),
        axis.title.x = element_text(vjust = -1.0),
        legend.position = "bottom",
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.major = element_line(colour = "black", size = 0.25), 
        panel.grid.minor = element_line(colour = "black", linetype = 2, size = 0.25),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        plot.background = element_rect(fill = "transparent",colour = NA)) +
  labs(x = "",
       #x = "Length of KHAP chain",
       y = "Concordance",
       title = expression(paste(bold("D."), "")))
khap_conc_box
ggsave(filename = paste0(container, "Plots/kHAP_conc.png"), plot = khap_conc_box, width = 297, height = 210, units = "mm", dpi = 300)
ggsave(filename = paste0(container, "Plots/kHAP_conc_tp.png"), bg = "transparent", plot = khap_conc_box, width = 297, height = 210, units = "mm", dpi = 300)

khap_plots = grid.arrange(top = "Number of included reference haplotypes (k_hap) variations",
                          arrangeGrob(khap_mcc_box, khap_info_box, khap_cert_box, khap_conc_box, ncol = 2, nrow = 2)
                          )
ggsave(filename = paste0(container, "Plots/KHAP.png"), plot = khap_plots, width = 297, height = 210, units = "mm", dpi = 300)
ggsave(filename = paste0(container, "Plots/KHAP_tp.png"), bg = "transparent", plot = khap_plots, width = 297, height = 210, units = "mm", dpi = 300)

## MAF

MAF_imp = read.csv(paste0(container, "MAF/ConcordanceTables_MAF_Experiments_MCC_imputed_genotype.csv"), header = T)
MAF_typ = read.csv(paste0(container, "MAF/ConcordanceTables_MAF_Experiments_MCC_typed_genotype.csv"), header = T)
MAF_imp$version = "imputed"
MAF_typ$version = "genotyped"
#main_maf_df = rbind(MAF_imp, MAF_typ)
main_maf_df = MAF_imp

exp.dir = "MAF_Experiments"
exp.var = c("MAF1%", "MAF0.5%", "MAF0.1%")
main_maf_df$sub_experiment = factor(main_maf_df$sub_experiment, levels = exp.var)

maf_mcc_info_regression = ggplot(main_maf_df, aes(x = mcc, y = info, colour = sub_experiment)) +
  geom_point() +
  #stat_smooth(method = "lm") +
  theme(axis.text.x = element_text(hjust = 1.0, angle = 45),
        axis.title.x = element_text(vjust = -1.0),
        legend.position = "bottom",
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.major = element_line(colour = "black", size = 0.25), 
        panel.grid.minor = element_line(colour = "black", linetype = 2, size = 0.25),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        plot.background = element_rect(fill = "transparent",colour = NA)) +
  labs(y = expression(paste("IMPUTE2 info score (", italic("r")^"2", ")")),
       x = "Matthew's Correlation Coefficient")
#ggsave(filename = paste0(container, "Plots/MCMC_MCCvINFO_regression.png"), plot = mcmc_mcc_info_regression, width = 297, height = 210, units = "mm", dpi = 300)


maf_mcc_box = ggplot(main_maf_df, aes(x = sub_experiment, y = mcc)) +
  geom_violin(fill = "#feb600", na.rm = T, lwd = rel(1/2)) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428", lwd = rel(1/2), fatten = rel(2.5)) +
  #geom_jitter(position=position_jitter(0.25), aes(size = n.snps), na.rm = T, shape = 21, fill = NA, colour = "#802428", stroke = rel(0.5)) + #colour = "#f3e5b1"
  #theme_bw() +
  theme(axis.text.x = element_text(hjust = 1.0, angle = 45),
        axis.title.x = element_text(vjust = -1.0),
        legend.position = "bottom",
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.major = element_line(colour = "black", size = 0.25), 
        panel.grid.minor = element_line(colour = "black", linetype = 2, size = 0.25),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        plot.background = element_rect(fill = "transparent",colour = NA)) +
  scale_y_continuous(breaks = seq(-1.0, 1.0, by = 0.2),
                     limits = c(-1, 1)) +
  labs(x = "",
       #x = "Length of MAF chain",
       y = "Matthew's Correlation Coefficient",
       title = expression(paste(bold("A."), "")))
maf_mcc_box
ggsave(filename = paste0(container, "Plots/MAF_MCC.png"), plot = maf_mcc_box, width = 297, height = 210, units = "mm", dpi = 300)
ggsave(filename = paste0(container, "Plots/MAF_MCC_tp.png"), bg = "transparent", plot = maf_mcc_box, width = 297, height = 210, units = "mm", dpi = 300)

maf_info_box = ggplot(main_maf_df, aes(x = sub_experiment, y = info)) +
  geom_violin(fill = "#feb600", na.rm = T, lwd = rel(1/2)) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428", lwd = rel(1/2), fatten = rel(2.5)) +
  #theme_bw() +
  theme(axis.text.x = element_text(hjust = 1.0, angle = 45),
        axis.title.x = element_text(vjust = -1.0),
        legend.position = "bottom",
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.major = element_line(colour = "black", size = 0.25), 
        panel.grid.minor = element_line(colour = "black", linetype = 2, size = 0.25),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        plot.background = element_rect(fill = "transparent",colour = NA)) +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                     limits = c(0, 1)) +
  labs(x = "",
       #x = "Length of MAF chain",
       y = expression(paste("IMPUTE2 info score (", italic("r")^"2", ")")),
       title = expression(paste(bold("B."), "")))
maf_info_box
ggsave(filename = paste0(container, "Plots/MAF_info.png"), plot = maf_info_box, width = 297, height = 210, units = "mm", dpi = 300)
ggsave(filename = paste0(container, "Plots/MAF_info_tp.png"), bg = "transparent", plot = maf_info_box, width = 297, height = 210, units = "mm", dpi = 300)

maf_infoComb_box = ggplot(main_maf_df, aes(x = sub_experiment, y = info_comb)) +
  geom_violin(fill = "#feb600", na.rm = T, lwd = rel(1/2)) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428", lwd = rel(1/2), fatten = rel(2.5)) +
  #theme_bw() +
  theme(axis.text.x = element_text(hjust = 1.0, angle = 45),
        axis.title.x = element_text(vjust = -1.0),
        legend.position = "bottom",
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.major = element_line(colour = "black", size = 0.25), 
        panel.grid.minor = element_line(colour = "black", linetype = 2, size = 0.25),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        plot.background = element_rect(fill = "transparent",colour = NA)) +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                     limits = c(0, 1)) +
  labs(x = "",
       #x = "Length of MAF chain",
       y = expression(paste("IMPUTE2 info score (", italic("r")^"2", ")")),
       title = expression(paste(bold("A."), " Markov chain Monte Carlo (MAF) length variations")))
maf_infoComb_box
ggsave(filename = paste0(container, "Plots/MAF_infoComb.png"), plot = maf_infoComb_box, width = 297, height = 210, units = "mm", dpi = 300)
ggsave(filename = paste0(container, "Plots/MAF_infoComb_tp.png"), bg = "transparent", plot = maf_infoComb_box, width = 297, height = 210, units = "mm", dpi = 300)

maf_cert_box = ggplot(main_maf_df, aes(x = sub_experiment, y = certainty)) +
  geom_violin(fill = "#feb600", na.rm = T, lwd = rel(1/2)) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428", lwd = rel(1/2), fatten = rel(2.5)) +
  #theme_bw() +
  theme(axis.text.x = element_text(hjust = 1.0, angle = 45),
        axis.title.x = element_text(vjust = -1.0),
        legend.position = "bottom",
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.major = element_line(colour = "black", size = 0.25), 
        panel.grid.minor = element_line(colour = "black", linetype = 2, size = 0.25),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        plot.background = element_rect(fill = "transparent",colour = NA)) +
  labs(x = "",
       #x = "Length of MAF chain",
       y = "Certainty",
       title = expression(paste(bold("C."), "")))
maf_cert_box
ggsave(filename = paste0(container, "Plots/MAF_cert.png"), plot = maf_cert_box, width = 297, height = 210, units = "mm", dpi = 300)
ggsave(filename = paste0(container, "Plots/MAF_cert_tp.png"), bg = "transparent", plot = maf_cert_box, width = 297, height = 210, units = "mm", dpi = 300)

maf_conc_box = ggplot(main_maf_df, aes(x = sub_experiment, y = concordance)) +
  geom_violin(fill = "#feb600", na.rm = T, lwd = rel(1/2)) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428", lwd = rel(1/2), fatten = rel(2.5)) +
  #geom_jitter(position=position_jitter(0.25), aes(size = n.snps), na.rm = T, shape = 21, fill = NA, colour = "#802428", stroke = rel(0.5)) + #colour = "#f3e5b1"
  #theme_bw() +
  theme(axis.text.x = element_text(hjust = 1.0, angle = 45),
        axis.title.x = element_text(vjust = -1.0),
        legend.position = "bottom",
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.major = element_line(colour = "black", size = 0.25), 
        panel.grid.minor = element_line(colour = "black", linetype = 2, size = 0.25),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        plot.background = element_rect(fill = "transparent",colour = NA)) +
  labs(x = "",
       #x = "Length of MAF chain",
       y = "Concordance",
       title = expression(paste(bold("D."), "")))
maf_conc_box
ggsave(filename = paste0(container, "Plots/MAF_conc.png"), plot = maf_conc_box, width = 297, height = 210, units = "mm", dpi = 300)
ggsave(filename = paste0(container, "Plots/MAF_conc_tp.png"), bg = "transparent", plot = maf_conc_box, width = 297, height = 210, units = "mm", dpi = 300)

maf_plots = grid.arrange(top = "Minor allele frequencies (MAF) variations",
                         arrangeGrob(maf_mcc_box, maf_info_box, maf_cert_box, maf_conc_box, ncol = 2, nrow = 2)
                         )
ggsave(filename = paste0(container, "Plots/MAF.png"), plot = maf_plots, width = 297, height = 210, units = "mm", dpi = 300)
ggsave(filename = paste0(container, "Plots/MAF_tp.png"), bg = "transparent", plot = maf_plots, width = 297, height = 210, units = "mm", dpi = 300)

regressionPlots = grid.arrange(arrangeGrob(mcmc_mcc_info_regression, khap_mcc_info_regression, maf_mcc_info_regression,
                                           ncol = 1, nrow = 3))

ggsave(filename = paste0(container, "Plots/MCC_v_INFO_regression.png"), plot = regressionPlots, width = 297, height = 297, units = "mm", dpi = 300)
ggsave(filename = paste0(container, "Plots/MCC_v_INFO_regression_tp.png"), bg = "transparent", plot = regressionPlots, width = 297, height = 297, units = "mm", dpi = 300)

