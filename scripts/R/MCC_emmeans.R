require(emmeans)
require(dplyr)
require(scales)

container = "/Volumes/TimMcInerney/MitoImpute/analyses/concordance/MCC/"
stat.out.dir = "/Volumes/TimMcInerney/MitoImpute/analyses/concordance/MCC/statistical/"

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

lm_MCMC = lm(mcc ~ sub_experiment, data = main_mcmc_df)
#anova(lm_MCMC)
#summary(lm_MCMC)
#emmeans(lm_MCMC, pairwise ~ sub_experiment)
lm_MCMC.df = summary(emmeans(lm_MCMC, pairwise ~ sub_experiment))
lm_MCMC.df$emmeans
lm_MCMC.df$contrasts

write.csv(lm_MCMC.df$contrasts, paste0(stat.out.dir, "MCMC_MCC_contrasts.csv"), row.names = F, quote = F)
write.csv(lm_MCMC.df$emmeans, paste0(stat.out.dir, "MCMC_MCC_emmeans.csv"), row.names = F, quote = F)

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

lm_KHAP = lm(mcc ~ sub_experiment, data = main_khap_df)
#anova(lm_MCMC)
#summary(lm_MCMC)
#emmeans(lm_MCMC, pairwise ~ sub_experiment)
lm_KHAP.df = summary(emmeans(lm_MCMC, pairwise ~ sub_experiment))
lm_KHAP.df$emmeans
lm_KHAP.df$contrasts

write.csv(lm_KHAP.df$contrasts, paste0(stat.out.dir, "MCMC_MCC_contrasts.csv"), row.names = F, quote = F)
write.csv(lm_KHAP.df$emmeans, paste0(stat.out.dir, "MCMC_MCC_emmeans.csv"), row.names = F, quote = F)