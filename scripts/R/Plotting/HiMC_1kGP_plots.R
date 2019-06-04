mcmc_tables = read.csv("~/GitCode/MitoImputePrep/metadata/Concordance_tables/HiMC_HaploGrep/combined/ConcordanceTables_MCMC_Experiments_HiMC_COMBINED.csv", header = T)
maf_tables = read.csv("~/GitCode/MitoImputePrep/metadata/Concordance_tables/HiMC_HaploGrep/combined/ConcordanceTables_MAF_Experiments_HiMC_COMBINED.csv", header = T)
khap_tables = read.csv("~/GitCode/MitoImputePrep/metadata/Concordance_tables/HiMC_HaploGrep/combined/ConcordanceTables_kHAP_Experiments_HiMC_COMBINED.csv", header = T)

plot_dir = "~/Desktop/MitoImpute_plots/HiMC/"

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
        axis.text.y = element_text(size = rel(1.0)),
        axis.title.y = element_text(size = rel(1.0)),
        plot.title = element_text(size = rel(1.0)),
        #plot.title = element_blank(),
        panel.grid.major = element_line(colour = "black", size = rel(1/2)), 
        panel.grid.minor = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.background = element_rect(fill = "transparent",colour = "black"),
        plot.background = element_rect(fill = "transparent",colour = NA))
ggsave(paste0(plot_dir, "HiMC_MAF_typedMatch_boxplot_tp.png"), bg = "transparent", width = 297, height = 210, units = "mm", dpi = 300)
ggsave(paste0(plot_dir, "HiMC_MAF_typedMatch_boxplot.png"), width = 297, height = 210, units = "mm", dpi = 300)

ggplot(maf_tables, aes(x = sub_experiment, y = imputed_match)) +
  geom_violin(fill = "#feb600", na.rm = T, lwd = rel(1/2)) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428", lwd = rel(1/2), fatten = rel(2.5)) +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                     limits = c(0, 1)) +
  labs(x = "Minor allele frequency (MAF)",
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
