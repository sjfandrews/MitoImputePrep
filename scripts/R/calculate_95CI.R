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
}

out.dir = "~/GitCode/MitoImputePrep/metadata/Concordance_tables/MCC/95CI_tables/"

kHAP = read.csv("~/GitCode/MitoImputePrep/metadata/Concordance_tables/ConcordanceTables_kHAP_Combined.csv", header = T)
MCMC = read.csv("~/GitCode/MitoImputePrep/metadata/Concordance_tables/ConcordanceTables_MCMC_Combined.csv", header = T)
MAF = read.csv("~/GitCode/MitoImputePrep/metadata/Concordance_tables/ConcordanceTables_MAF_Combined.csv", header = T)

kHAP$diff = kHAP$Imputed.hg.Conc - kHAP$Typed.hg.Conc
MCMC$diff = MCMC$Imputed.hg.Conc - MCMC$Typed.hg.Conc
MAF$diff = MAF$Imputed.hg.Conc - MAF$Typed.hg.Conc

write.csv(aggregate(kHAP$Imputed.hg.Conc~kHAP$kHAP, FUN=confidence.interval), paste0(out.dir, "kHAP_imputed_HiMC_95.csv"), row.names = F, quote = F)
write.csv(aggregate(MCMC$Imputed.hg.Conc~MCMC$MCMC, FUN=confidence.interval), paste0(out.dir, "MCMC_imputed_HiMC_95.csv"), row.names = F, quote = F)
write.csv(aggregate(MAF$Imputed.hg.Conc~MAF$MAF, FUN=confidence.interval), paste0(out.dir, "MAF_imputed_HiMC_95.csv"), row.names = F, quote = F)

write.csv(aggregate(kHAP$Typed.hg.Conc~kHAP$kHAP, FUN=confidence.interval), paste0(out.dir, "kHAP_typed_HiMC_95.csv"), row.names = F, quote = F)
write.csv(aggregate(MCMC$Typed.hg.Conc~MCMC$MCMC, FUN=confidence.interval), paste0(out.dir, "MCMC_typed_HiMC_95.csv"), row.names = F, quote = F)
write.csv(aggregate(MAF$Typed.hg.Conc~MAF$MAF, FUN=confidence.interval), paste0(out.dir, "MAF_typed_HiMC_95.csv"), row.names = F, quote = F)

write.csv(aggregate(kHAP$diff~kHAP$kHAP, FUN=confidence.interval), paste0(out.dir, "kHAP_diff_HiMC_95.csv"), row.names = F, quote = F)
write.csv(aggregate(MCMC$diff~MCMC$MCMC, FUN=confidence.interval), paste0(out.dir, "MCMC_diff_HiMC_95.csv"), row.names = F, quote = F)
write.csv(aggregate(MAF$diff~MAF$MAF, FUN=confidence.interval), paste0(out.dir, "MAF_diff_HiMC_95.csv"), row.names = F, quote = F)

kHAP = read.csv("/Volumes/TimMcInerney/MitoImpute/analyses/concordance/MCC/KHAP/ConcordanceTables_kHAP_Experiments_MCC_imputed_genotype.csv", header = T)
MCMC = read.csv("/Volumes/TimMcInerney/MitoImpute/analyses/concordance/MCC/MCMC/ConcordanceTables_MCMC_Experiments_MCC_imputed_genotype.csv", header = T)
MAF = read.csv("/Volumes/TimMcInerney/MitoImpute/analyses/concordance/MCC/MAF/ConcordanceTables_MAF_Experiments_MCC_imputed_genotype.csv", header = T)

write.csv(aggregate(kHAP$mcc~kHAP$sub_experiment, FUN=confidence.interval), paste0(out.dir, "kHAP_imputed_mcc_95.csv"), row.names = F, quote = F)
write.csv(aggregate(MCMC$mcc~MCMC$sub_experiment, FUN=confidence.interval), paste0(out.dir, "MCMC_imputed_mcc_95.csv"), row.names = F, quote = F)
write.csv(aggregate(MAF$mcc~MAF$sub_experiment, FUN=confidence.interval), paste0(out.dir, "MAF_imputed_mcc_95.csv"), row.names = F, quote = F)

write.csv(aggregate(kHAP$info~kHAP$sub_experiment, FUN=confidence.interval), paste0(out.dir, "kHAP_imputed_info_95.csv"), row.names = F, quote = F)
write.csv(aggregate(MCMC$info~MCMC$sub_experiment, FUN=confidence.interval), paste0(out.dir, "MCMC_imputed_info_95.csv"), row.names = F, quote = F)
write.csv(aggregate(MAF$info~MAF$sub_experiment, FUN=confidence.interval), paste0(out.dir, "MCMC_imputed_info_95.csv"), row.names = F, quote = F)

write.csv(aggregate(kHAP$concordance~kHAP$sub_experiment, FUN=confidence.interval), paste0(out.dir, "kHAP_imputed_concordance_95.csv"), row.names = F, quote = F)
write.csv(aggregate(MCMC$concordance~MCMC$sub_experiment, FUN=confidence.interval), paste0(out.dir, "MCMC_imputed_concordance_95.csv"), row.names = F, quote = F)
write.csv(aggregate(MAF$concordance~MAF$sub_experiment, FUN=confidence.interval), paste0(out.dir, "MCMC_imputed_concordance_95.csv"), row.names = F, quote = F)

ADNI = read.csv("/Volumes/TimMcInerney/MitoImpute/data/ADNI_REDO/IMPUTED/ReferencePanel_v2/IMPUTE2/MitoImpute_ReferencePanel_v2_imputed_typed_MCC.csv", header = T)
confidence.interval(ADNI$mcc)
confidence.interval(ADNI$info)
confidence.interval(ADNI$concodance)
