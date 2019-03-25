require(ggplot2)
require(scales)

# MCMC

MCMC_imp = read.csv("~/Desktop/ConcordanceTables_MCMC_Experiments_MCC_imputed_genotype.csv", header = T)
MCMC_typ = read.csv("~/Desktop/ConcordanceTables_MCMC_Experiments_MCC_imputed_genotype.csv", header = T)
MCMC_imp$version = "imputed"
MCMC_typ$version = "genotyped"
#main_mcmc_df = rbind(MCMC_imp, MCMC_typ)
main_mcmc_df = MCMC_imp

exp.dir = "MCMC_Experiments"
exp.var = c("MCMC1", "MCMC5", "MCMC10", "MCMC20", "MCMC30")
main_mcmc_df$sub_experiment = factor(main_mcmc_df$sub_experiment, levels = exp.var)

mcmc_mcc_box = ggplot(main_mcmc_df, aes(x = sub_experiment, y = abs(mcc))) +
  geom_violin(fill = "#feb600", na.rm = T) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428") +
  theme_bw() +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                     limits = c(0, 1)) +
  labs(x = "Length of MCMC chain",
       y = "Matthew's Correlation Coefficient",
       title = expression(paste(bold("A."), " Markov chain Monte Carlo (MCMC) length variations")))
mcmc_mcc_box

mcmc_info_box = ggplot(main_mcmc_df, aes(x = sub_experiment, y = info)) +
  geom_violin(fill = "#feb600", na.rm = T) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428") +
  theme_bw() +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                     limits = c(0, 1)) +
  labs(x = "Length of MCMC chain",
       y = expression(paste("IMPUTE2 info score (", italic("r")^"2", ")")),
       title = expression(paste(bold("A."), " Markov chain Monte Carlo (MCMC) length variations")))
mcmc_info_box

mcmc_infoComb_box = ggplot(main_mcmc_df, aes(x = sub_experiment, y = info_comb)) +
  geom_violin(fill = "#feb600", na.rm = T) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428") +
  theme_bw() +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                     limits = c(0, 1)) +
  labs(x = "Length of MCMC chain",
       y = expression(paste("IMPUTE2 info score (", italic("r")^"2", ")")),
       title = expression(paste(bold("A."), " Markov chain Monte Carlo (MCMC) length variations")))
mcmc_infoComb_box

mcmc_cert_box = ggplot(main_mcmc_df, aes(x = sub_experiment, y = certainty)) +
  geom_violin(fill = "#feb600", na.rm = T) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428") +
  theme_bw() +
  labs(x = "Length of MCMC chain",
       y = expression(paste("IMPUTE2 info score (", italic("r")^"2", ")")),
       title = expression(paste(bold("A."), " Markov chain Monte Carlo (MCMC) length variations")))
mcmc_cert_box

mcmc_conc_box = ggplot(main_mcmc_df, aes(x = sub_experiment, y = concordance)) +
  geom_violin(fill = "#feb600", na.rm = T) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428") +
  theme_bw() +
  labs(x = "Length of MCMC chain",
       y = expression(paste("IMPUTE2 info score (", italic("r")^"2", ")")),
       title = expression(paste(bold("A."), " Markov chain Monte Carlo (MCMC) length variations")))
mcmc_conc_box

# KHAP
KHAP_imp = read.csv("~/Desktop/ConcordanceTables_KHAP_Experiments_MCC_imputed_genotype.csv", header = T)
KHAP_typ = read.csv("~/Desktop/ConcordanceTables_KHAP_Experiments_MCC_imputed_genotype.csv", header = T)
KHAP_imp$version = "imputed"
KHAP_typ$version = "genotyped"
#main_khap_df = rbind(KHAP_imp, KHAP_typ)
main_khap_df = KHAP_imp

exp.dir = "kHAP_Experiments"
exp.var = c("kHAP100", "kHAP250", "kHAP500", "kHAP1000", "kHAP2500", "kHAP5000", "kHAP10000", "kHAP20000", "kHAP30000")
main_khap_df$sub_experiment = factor(main_khap_df$sub_experiment, levels = exp.var)

khap_mcc_box = ggplot(main_khap_df, aes(x = sub_experiment, y = abs(mcc))) +
  geom_violin(fill = "#feb600", na.rm = T) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428") +
  theme_bw() +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                     limits = c(0, 1)) +
  labs(x = "Length of KHAP chain",
       y = "Matthew's Correlation Coefficient",
       title = expression(paste(bold("A."), " Markov chain Monte Carlo (KHAP) length variations")))
khap_mcc_box

khap_info_box = ggplot(main_khap_df, aes(x = sub_experiment, y = info)) +
  geom_violin(fill = "#feb600", na.rm = T) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428") +
  theme_bw() +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                     limits = c(0, 1)) +
  labs(x = "Length of KHAP chain",
       y = expression(paste("IMPUTE2 info score (", italic("r")^"2", ")")),
       title = expression(paste(bold("A."), " Markov chain Monte Carlo (KHAP) length variations")))
khap_info_box

khap_infoComb_box = ggplot(main_khap_df, aes(x = sub_experiment, y = info_comb)) +
  geom_violin(fill = "#feb600", na.rm = T) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428") +
  theme_bw() +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                     limits = c(0, 1)) +
  labs(x = "Length of KHAP chain",
       y = expression(paste("IMPUTE2 info score (", italic("r")^"2", ")")),
       title = expression(paste(bold("A."), " Markov chain Monte Carlo (KHAP) length variations")))
khap_infoComb_box

khap_cert_box = ggplot(main_khap_df, aes(x = sub_experiment, y = certainty)) +
  geom_violin(fill = "#feb600", na.rm = T) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428") +
  theme_bw() +
  labs(x = "Length of KHAP chain",
       y = expression(paste("IMPUTE2 info score (", italic("r")^"2", ")")),
       title = expression(paste(bold("A."), " Markov chain Monte Carlo (KHAP) length variations")))
khap_cert_box

khap_conc_box = ggplot(main_khap_df, aes(x = sub_experiment, y = concordance)) +
  geom_violin(fill = "#feb600", na.rm = T) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428") +
  theme_bw() +
  labs(x = "Length of KHAP chain",
       y = expression(paste("IMPUTE2 info score (", italic("r")^"2", ")")),
       title = expression(paste(bold("A."), " Markov chain Monte Carlo (KHAP) length variations")))
khap_conc_box

## MAF

MAF_imp = read.csv("~/Desktop/ConcordanceTables_MAF_Experiments_MCC_imputed_genotype.csv", header = T)
MAF_typ = read.csv("~/Desktop/ConcordanceTables_MAF_Experiments_MCC_imputed_genotype.csv", header = T)
MAF_imp$version = "imputed"
MAF_typ$version = "genotyped"
main_maf_df = rbind(MAF_imp, MAF_typ)
main_maf_df = MAF_imp

exp.dir = "MAF_Experiments"
exp.var = c("MAF1%", "MAF0.5%", "MAF0.1%")
main_maf_df$sub_experiment = factor(main_maf_df$sub_experiment, levels = exp.var)

maf_mcc_box = ggplot(main_maf_df, aes(x = sub_experiment, y = abs(mcc))) +
  geom_violin(fill = "#feb600", na.rm = T) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428") +
  theme_bw() +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                     limits = c(0, 1)) +
  labs(x = "Length of MAF chain",
       y = "Matthew's Correlation Coefficient",
       title = expression(paste(bold("A."), " Markov chain Monte Carlo (MAF) length variations")))
maf_mcc_box

maf_info_box = ggplot(main_maf_df, aes(x = sub_experiment, y = info)) +
  geom_violin(fill = "#feb600", na.rm = T) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428") +
  theme_bw() +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                     limits = c(0, 1)) +
  labs(x = "Length of MAF chain",
       y = expression(paste("IMPUTE2 info score (", italic("r")^"2", ")")),
       title = expression(paste(bold("A."), " Markov chain Monte Carlo (MAF) length variations")))
maf_info_box

maf_infoComb_box = ggplot(main_maf_df, aes(x = sub_experiment, y = info_comb)) +
  geom_violin(fill = "#feb600", na.rm = T) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428") +
  theme_bw() +
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                     limits = c(0, 1)) +
  labs(x = "Length of MAF chain",
       y = expression(paste("IMPUTE2 info score (", italic("r")^"2", ")")),
       title = expression(paste(bold("A."), " Markov chain Monte Carlo (MAF) length variations")))
maf_infoComb_box

maf_cert_box = ggplot(main_maf_df, aes(x = sub_experiment, y = certainty)) +
  geom_violin(fill = "#feb600", na.rm = T) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428") +
  theme_bw() +
  labs(x = "Length of MAF chain",
       y = expression(paste("IMPUTE2 info score (", italic("r")^"2", ")")),
       title = expression(paste(bold("A."), " Markov chain Monte Carlo (MAF) length variations")))
maf_cert_box

maf_conc_box = ggplot(main_maf_df, aes(x = sub_experiment, y = concordance)) +
  geom_violin(fill = "#feb600", na.rm = T) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428") +
  theme_bw() +
  labs(x = "Length of MAF chain",
       y = expression(paste("IMPUTE2 info score (", italic("r")^"2", ")")),
       title = expression(paste(bold("A."), " Markov chain Monte Carlo (MAF) length variations")))
maf_conc_box
