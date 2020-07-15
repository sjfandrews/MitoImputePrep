library(tidyverse)

haplogrep_file = paste0("/Volumes/TimMcInerney/MitoImpute/data/HAPLOGROUPS/combined/ConcordanceTables_MAF_Experiments_COMBINED.csv")
mcc_file       = paste0("/Volumes/TimMcInerney/MitoImpute/analyses/concordance/MCC/MAF/ConcordanceTables_MAF_Experiments_MCC_imputed_genotype.csv")

haplogrep_df = read_csv(haplogrep_file)
mcc_df       = read_csv(mcc_file)

haplogrep_df = haplogrep_df %>% filter(sub_experiment == "MAF1%")
mcc_df       = mcc_df %>% filter(sub_experiment == "MAF1%")

shared_names = names(haplogrep_df)[names(haplogrep_df) %in% names(mcc_df)]

joined_df = haplogrep_df %>%
  inner_join(mcc_df, by = shared_names) %>%
  write_csv(path = "~/GitCode/MitoImputePrep/metadata/Concordance_tables/MCC/HaploGrep_MCC_joined.csv")


joined_df

needed_names = c("array", "MCMC", "MAF", "k_hap", "imputed", "mean.info", "SNP.Ref.Only", "SNP.Ref.Samp", "SNP.Samp.Only", "TOTAL", "Retained.After.Filt", "Typed.hg.Conc", "Imputed.hg.Conc", "n.snps", "allele_freq", "mcc", "concordance", "certainty", "info_type0", "concord_type0", "r2_type0")

shared_named_excel = names(joined_df) %in% needed_names

names(joined_df)[shared_named_excel]
