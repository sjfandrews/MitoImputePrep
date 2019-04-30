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
  return(c("mean" = m, "lower" = l, "upper" = u, "n" = as.integer(n)))
}

kHAP = read.csv("~/GitCode/MitoImputePrep/metadata/Concordance_tables/ConcordanceTables_kHAP_Combined.csv", header = T)
MCMC = read.csv("~/GitCode/MitoImputePrep/metadata/Concordance_tables/ConcordanceTables_MCMC_Combined.csv", header = T)
MAF = read.csv("~/GitCode/MitoImputePrep/metadata/Concordance_tables/ConcordanceTables_MAF_Combined.csv", header = T)

kHAP$diff = kHAP$Imputed.hg.Conc - kHAP$Typed.hg.Conc
MCMC$diff = MCMC$Imputed.hg.Conc - MCMC$Typed.hg.Conc
MAF$diff = MAF$Imputed.hg.Conc - MAF$Typed.hg.Conc

aggregate(kHAP$Imputed.hg.Conc~kHAP$kHAP, FUN=confidence.interval)
aggregate(MCMC$Imputed.hg.Conc~MCMC$MCMC, FUN=confidence.interval)
aggregate(MAF$Imputed.hg.Conc~MAF$MAF, FUN=confidence.interval)

aggregate(kHAP$diff~kHAP$kHAP, FUN=confidence.interval)
aggregate(MCMC$diff~MCMC$MCMC, FUN=confidence.interval)
aggregate(MAF$diff~MAF$MAF, FUN=confidence.interval)

kHAP = read.csv("/Volumes/TimMcInerney/MitoImpute/analyses/concordance/MCC/KHAP/ConcordanceTables_kHAP_Experiments_MCC_imputed_genotype.csv", header = T)
MCMC = read.csv("/Volumes/TimMcInerney/MitoImpute/analyses/concordance/MCC/MCMC/ConcordanceTables_MCMC_Experiments_MCC_imputed_genotype.csv", header = T)
MAF = read.csv("/Volumes/TimMcInerney/MitoImpute/analyses/concordance/MCC/MAF/ConcordanceTables_MAF_Experiments_MCC_imputed_genotype.csv", header = T)

aggregate(kHAP$mcc~kHAP$sub_experiment, FUN=confidence.interval)
aggregate(MCMC$mcc~MCMC$sub_experiment, FUN=confidence.interval)
aggregate(MAF$mcc~MAF$sub_experiment, FUN=confidence.interval)

aggregate(kHAP$info~kHAP$sub_experiment, FUN=confidence.interval)
aggregate(MCMC$info~MCMC$sub_experiment, FUN=confidence.interval)
aggregate(MAF$info~MAF$sub_experiment, FUN=confidence.interval)

aggregate(kHAP$concordance~kHAP$sub_experiment, FUN=confidence.interval)
aggregate(MCMC$concordance~MCMC$sub_experiment, FUN=confidence.interval)
aggregate(MAF$concordance~MAF$sub_experiment, FUN=confidence.interval)

ADNI = read.csv("/Volumes/TimMcInerney/MitoImpute/data/ADNI_REDO/IMPUTED/ReferencePanel_v2/IMPUTE2/MitoImpute_ReferencePanel_v2_imputed_typed_MCC.csv", header = T)
confidence.interval(ADNI$mcc)
confidence.interval(ADNI$info)
confidence.interval(ADNI$concodance)
