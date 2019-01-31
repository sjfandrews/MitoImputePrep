library(tidyverse)
library(readxl)
library(HiMC); data(nodes)
library(dplyr)

'%!in%' = function(x,y)!('%in%'(x,y))

## READ IN FILE TO RECONCILE SAMPLE IDs
rect_csv = read.csv("/Users/u5015730/GitCode/MitoImputePrep/metadata/ADNI_samples_BOTH.csv", header = T) 

## READ IN THE LIST OF CHIPS
chips = read.table("~/GitCode/MitoImputePrep/scripts/INFORMATION_LISTS/b37_platforms.txt")
names(chips) = c("array")

info.cutoff = 0.3

## READ IN .ped AND .map FILES FOR TRUTH SET 
full_1kGP = generate_snp_data("/Volumes/TimMcInerney/MitoImpute/data/ADNI/Timpute/ADNI_reseq/DECOMPOSED/adni_mito_genomes_180214_norm_decomposed_firstAlt.map",
                              "/Volumes/TimMcInerney/MitoImpute/data/ADNI/Timpute/ADNI_reseq/DECOMPOSED/adni_mito_genomes_180214_norm_decomposed_firstAlt.ped")
full_1kGP_hg_FULL = HiMC::getClassifications(full_1kGP)


### WORK ON THE THE ith IMPUTED SET
## MAKE DATA FRAMES FOR EACH MINOR ALLELE FREQ

COLS = c("array", "MCMC", "imputed", "mean.info", "SNP.Ref.Only", "SNP.Ref.Samp", "SNP.Samp.Only", "TOTAL", "Retained.After.Filt", "Typed.hg.Conc", "Imputed.hg.Conc", "Typed.Mhg.Conc", "Imputed.Mhg.Conc")

#CONC_TABLE = data.frame(matrix(ncol = 6, nrow = ))
MCMC1 = data.frame(matrix(ncol = length(COLS), nrow = 1))

names(MCMC1) = COLS

# SET A NEW COLUMN DESIGNATING THE MAF AND ONE WHETHER THEY HAVE BEEN IMPUTED
MCMC1$MCMC = "1"
MCMC1$imputed = FALSE
mcmc = "1"

# LOOK AT ADNI RESULTS
i = 1
DIR = "/Volumes/TimMcInerney/MitoImpute/data/ADNI/Timpute/IMPUTE2/MCMC1/"
EXT= "mitoimpute_adni_imputed_MCMC1"
MCMC1$array[i] = as.character("ADNI")
MCMC1$MCMC[i] = mcmc
MCMC1$imputed[i] = T

## READ IN THE INFO SCORE
info.score <- read_delim(paste0(DIR, EXT, "_info"), delim = " ")
mean.info = mean(info.score$info)
exclude.info.score = subset(info.score, !(info.score$info >= info.cutoff))
exclude.pos = as.character(exclude.info.score$position)

MCMC1$mean.info[i] = mean.info
MCMC1$SNP.Ref.Only[i] = nrow(subset(info.score, info.score$type == 0))
MCMC1$SNP.Ref.Samp[i] = nrow(subset(info.score, info.score$type == 2))
MCMC1$SNP.Samp.Only[i] = nrow(subset(info.score, info.score$type == 3))
MCMC1$TOTAL[i] = sum(MCMC1$SNP.Ref.Only[i], MCMC1$SNP.Ref.Samp[i], MCMC1$SNP.Samp.Only[i])
MCMC1$Retained.After.Filt[i] = nrow(subset(info.score, info.score$info >= info.cutoff))

## FILTER VARIANTS IN TRUTH SET
full_1kGP_cut = full_1kGP[, -(which(names(full_1kGP) %in% exclude.pos))]
full_1kGP_hg = HiMC::getClassifications(full_1kGP_cut)
full_index_exclude = which(names(full_1kGP) %in% exclude.pos)
if (length(full_index_exclude) > 0) {
  full_1kGP_cut = full_1kGP[, -full_index_exclude]
  full_1kGP_hg = HiMC::getClassifications(full_1kGP_cut)
  message("TRUTH HAPLOTYPES CLASSIFIED")
} else {
  message("CAUTION! NOT FILTERING")
  message("PULLING FROM FULL LIST")
  full_1kGP_hg = full_1kGP_hg_FULL
}


## READ IN .ped AND .map FILES FOR IMPUTED SET 
message("READING IN .ped AND .map FILES FOR IMPUTED ADNI")
imputed_1kGP <- generate_snp_data(paste0(DIR, EXT, ".map"),
                                  paste0(DIR, EXT, ".ped"))
imputed_index_exclude = which(names(imputed_1kGP) %in% exclude.pos)
if (length(imputed_index_exclude) > 0) {
  imputed_1kGP_cut = imputed_1kGP[, -imputed_index_exclude]
  imputed_1kGP_hg = HiMC::getClassifications(imputed_1kGP_cut)
} else {
  imputed_1kGP_hg = HiMC::getClassifications(imputed_1kGP)
}
message("IMPUTED HAPLOTYPES CLASSIFIED")

DIR = "/Volumes/TimMcInerney/MitoImpute/data/ADNI/Timpute/PLINK/"
EXT = "mitoimpute_adni"

typed_1kGP <- generate_snp_data(paste0(DIR, EXT, ".map"),
                                paste0(DIR, EXT, ".ped"))
typed_index_exclude = which(names(typed_1kGP) %in% exclude.pos)
if (length(typed_index_exclude) > 0) {
  typed_1kGP_cut = typed_1kGP[, -typed_index_exclude]
  typed_1kGP_hg = HiMC::getClassifications(typed_1kGP_cut)
} else {
  typed_1kGP_hg = HiMC::getClassifications(typed_1kGP)
}
message("TYPED HAPLOTYPES CLASSIFIED")

## RECONCILE SAMPLE IDs
full_1kGP_hg$Individual = as.character(full_1kGP_hg$Individual)

for (i in 1:nrow(full_1kGP_hg)) {
  full_1kGP_hg$Individual[i] = as.character(rect_csv$TYPED_LABEL[match(full_1kGP_hg$Individual[i], rect_csv$SAMPLE_NUMBER)])
}

full_1kGP_hg = arrange(full_1kGP_hg, full_1kGP_hg$Individual)

## DEFINE META HAPLOGROUPS
difficult_mhgs = c("HV", "JT", "unclassified")

full_1kGP_hg_D = subset(full_1kGP_hg, full_1kGP_hg$haplogroup %in% difficult_mhgs)
full_1kGP_hg_D$metahaplogroup = full_1kGP_hg_D$haplogroup
full_1kGP_hg = subset(full_1kGP_hg, full_1kGP_hg$haplogroup %!in% difficult_mhgs)
full_1kGP_hg_A = subset(full_1kGP_hg, substr(full_1kGP_hg$haplogroup, 1, 1) == "L")
full_1kGP_hg_A$metahaplogroup = substr(full_1kGP_hg_A$haplogroup, 1, 2)
full_1kGP_hg_N = subset(full_1kGP_hg, substr(full_1kGP_hg$haplogroup, 1, 1) != "L")
full_1kGP_hg_N$metahaplogroup = substr(full_1kGP_hg_N$haplogroup, 1, 1)
full_1kGP_hg = rbind(full_1kGP_hg_A, full_1kGP_hg_N, full_1kGP_hg_D)
full_1kGP_hg = arrange(full_1kGP_hg, full_1kGP_hg$Individual)

typed_1kGP_hg_D = subset(typed_1kGP_hg, typed_1kGP_hg$haplogroup %in% difficult_mhgs)
typed_1kGP_hg_D$metahaplogroup = typed_1kGP_hg_D$haplogroup
typed_1kGP_hg = subset(typed_1kGP_hg, typed_1kGP_hg$haplogroup %!in% difficult_mhgs)
typed_1kGP_hg_A = subset(typed_1kGP_hg, substr(typed_1kGP_hg$haplogroup, 1, 1) == "L")
typed_1kGP_hg_A$metahaplogroup = substr(typed_1kGP_hg_A$haplogroup, 1, 2)
typed_1kGP_hg_N = subset(typed_1kGP_hg, substr(typed_1kGP_hg$haplogroup, 1, 1) != "L")
typed_1kGP_hg_N$metahaplogroup = substr(typed_1kGP_hg_N$haplogroup, 1, 1)
typed_1kGP_hg = rbind(typed_1kGP_hg_A, typed_1kGP_hg_N, typed_1kGP_hg_D)
typed_1kGP_hg = arrange(typed_1kGP_hg, typed_1kGP_hg$Individual)

imputed_1kGP_hg_D = subset(imputed_1kGP_hg, imputed_1kGP_hg$haplogroup %in% difficult_mhgs)
imputed_1kGP_hg_D$metahaplogroup = imputed_1kGP_hg_D$haplogroup
imputed_1kGP_hg = subset(imputed_1kGP_hg, imputed_1kGP_hg$haplogroup %!in% difficult_mhgs)
imputed_1kGP_hg_A = subset(imputed_1kGP_hg, substr(imputed_1kGP_hg$haplogroup, 1, 1) == "L")
imputed_1kGP_hg_A$metahaplogroup = substr(imputed_1kGP_hg_A$haplogroup, 1, 2)
imputed_1kGP_hg_N = subset(imputed_1kGP_hg, substr(imputed_1kGP_hg$haplogroup, 1, 1) != "L")
imputed_1kGP_hg_N$metahaplogroup = substr(imputed_1kGP_hg_N$haplogroup, 1, 1)
imputed_1kGP_hg = rbind(imputed_1kGP_hg_A, imputed_1kGP_hg_N, imputed_1kGP_hg_D)
imputed_1kGP_hg = arrange(imputed_1kGP_hg, imputed_1kGP_hg$Individual)

## CALCULATE HAPLOGROUP CONCORDANCE
conc_typed = full_1kGP_hg$haplogroup == typed_1kGP_hg$haplogroup
conc_imputed = full_1kGP_hg$haplogroup == imputed_1kGP_hg$haplogroup
conc_typed_pc = length(conc_typed[conc_typed == T]) / length(conc_typed)
conc_imputed_pc = length(conc_imputed[conc_imputed == T]) / length(conc_imputed)
MCMC1$Typed.hg.Conc[1] = conc_typed_pc
MCMC1$Imputed.hg.Conc[1] = conc_imputed_pc

M_conc_typed = full_1kGP_hg$metahaplogroup == typed_1kGP_hg$metahaplogroup
M_conc_imputed = full_1kGP_hg$metahaplogroup == imputed_1kGP_hg$metahaplogroup
M_conc_typed_pc = length(M_conc_typed[M_conc_typed == T]) / length(M_conc_typed)
M_conc_imputed_pc = length(M_conc_imputed[M_conc_imputed == T]) / length(M_conc_imputed)
MCMC1$Typed.Mhg.Conc[1] = M_conc_typed_pc
MCMC1$Imputed.Mhg.Conc[1] = M_conc_imputed_pc

write.csv(MCMC1, "~/GitCode/MitoImputePrep/metadata/ADNI_concordance_HiMC.csv", quote = F, row.names = F)

## UNIFY TABLE
unity.table = data.frame(matrix(ncol = 7, nrow = nrow(full_1kGP_hg)))
names(unity.table) = c("sample", "wgs", "wgs_meta", "typed", "imputed", "typed_meta", "imputed_meta")
unity.table$sample = full_1kGP_hg$Individual
unity.table$wgs = full_1kGP_hg$haplogroup
unity.table$wgs_meta = full_1kGP_hg$metahaplogroup
unity.table$typed = typed_1kGP_hg$haplogroup
unity.table$typed_meta = typed_1kGP_hg$metahaplogroup
unity.table$imputed = imputed_1kGP_hg$haplogroup
unity.table$imputed_meta = imputed_1kGP_hg$metahaplogroup

unity.table$typed_match = unity.table$wgs == unity.table$typed
unity.table$typed_meta_match = unity.table$wgs_meta == unity.table$typed_meta
unity.table$imputed_match = unity.table$wgs == unity.table$imputed
unity.table$imputed_meta_match = unity.table$wgs_meta == unity.table$imputed_meta

write.csv(unity.table, "~/GitCode/MitoImputePrep/metadata/ADNI_concordance_HiMC_FULL.csv", quote = F, row.names = F)

## PER HAPLOGROUP TABLE
hgs = names(table(unity.table$wgs))

hgs_df = data.frame(matrix(nrow = 2, ncol = length(hgs) + 1))
names(hgs_df) = c("chip", hgs)

hgs_df$chip[1] = "typed"
hgs_df$chip[2] = "imputed"

for (hg in 1:length(hgs)) { # FOR EACH HAPLOGROUP IN THE WGS (TRUTH) COLUMN ...
  t_t = subset(unity.table, unity.table$typed_match) # SUBSET TO ONLY TYPED MATCHES
  t_i = subset(unity.table, unity.table$imputed_match) # SUBSET TO ONLY IMPUTED MATCHES
  #v = nrow(t) / nrow(s)
  hgs_df[1, hgs[hg]] = nrow(t_t) / nrow(unity.table) # WRITE PROPORTION MATCHING TO CELL FOR TYPED
  hgs_df[2, hgs[hg]] = nrow(t_t) / nrow(unity.table) # WRITE PROPORTION MATCHING TO CELL FOR IMPUTED
}

###################################### HAPLOGREP ###################################### 

HaploGrep_WGS = read.table("/Volumes/TimMcInerney/MitoImpute/data/ADNI/Timpute/ADNI_reseq/adni_mito_genomes_180214_haplogrep2.txt", header = T)
HaploGrep_TYP = read.table("/Volumes/TimMcInerney/MitoImpute/data/ADNI/Timpute/VCF/mitoimpute_adni/mitoimpute_adni_diploid.txt", header = T)
HaploGrep_IMP = read.table("/Volumes/TimMcInerney/MitoImpute/data/ADNI/Timpute/IMPUTE2/MCMC1/mitoimpute_adni_imputed_MCMC1_haplogrep.txt", header = T)

HaploGrep_WGS$SampleID = as.character(HaploGrep_WGS$SampleID)
HaploGrep_WGS$Haplogroup = as.character(HaploGrep_WGS$Haplogroup)
HaploGrep_TYP$Haplogroup = as.character(HaploGrep_TYP$Haplogroup)
HaploGrep_IMP$Haplogroup = as.character(HaploGrep_IMP$Haplogroup)

HaploGrep_WGS$ORIGINAL = HaploGrep_WGS$SampleID
HaploGrep_WGS$BOTH = NA
for (i in 1:nrow(HaploGrep_WGS)) {
  if (is.na(match(HaploGrep_WGS$ORIGINAL[i], rect_csv$SAMPLE_NUMBER))) {
    HaploGrep_WGS$BOTH[i] = F
  } else {
    HaploGrep_WGS$SampleID[i] = as.character(rect_csv$TYPED_LABEL[match(HaploGrep_WGS$ORIGINAL[i], rect_csv$SAMPLE_NUMBER)])
    HaploGrep_WGS$BOTH[i] = T
  }
}

HaploGrep_WGS = subset(HaploGrep_WGS, HaploGrep_WGS$BOTH == T)
HaploGrep_WGS = arrange(HaploGrep_WGS, HaploGrep_WGS$SampleID)

df = data.frame(cbind(HaploGrep_WGS$SampleID, HaploGrep_WGS$Haplogroup))
names(df) = c("SampleID", "WGS")
df$TYP_HG = NA
df$IMP_HG = NA

for (i in 1:nrow(df)) {
  df$TYP_HG[i] = HaploGrep_TYP$Haplogroup[match(df$SampleID[i], HaploGrep_TYP$SampleID)]
  df$IMP_HG[i] = HaploGrep_IMP$Haplogroup[match(paste0(df$SampleID[i], "_", df$SampleID[i]), HaploGrep_IMP$SampleID)]
}
df$TYP_MATCH = df$WGS == df$TYP_HG
df$IMP_MATCH = df$WGS == df$IMP_HG

A_df = subset(df, substr(df$WGS, 1, 1) == "L")
A_df$META_WGS = substr(A_df$WGS, 1, 2)
A_df$META_TYP_HG = substr(A_df$TYP_HG, 1, 2)
A_df$META_IMP_HG = substr(A_df$IMP_HG, 1, 2)

N_df = subset(df, substr(df$WGS, 1, 1) != "L")
N_df$META_WGS = substr(N_df$WGS, 1, 1)
N_df$META_TYP_HG = substr(N_df$TYP_HG, 1, 1)
N_df$META_IMP_HG = substr(N_df$IMP_HG, 1, 1)

df = rbind(N_df, A_df)
df = arrange(df, df$SampleID)

df$META_TYP_MATCH = df$META_WGS == df$META_TYP_HG
df$META_IMP_MATCH = df$META_WGS == df$META_IMP_HG

out_table = data.frame(matrix(ncol = 3, nrow = 2))
names(out_table) = c("HG", "TYPED", "IMPUTED")

out_table$HG[1] = "Haplogroup"
out_table$TYPED[1] = nrow(subset(df, df$TYP_MATCH == T)) / nrow(df)
out_table$IMPUTED[1] = nrow(subset(df, df$IMP_MATCH == T)) / nrow(df)

out_table$HG[2] = "Meta-Haplogroup"
out_table$TYPED[2] = nrow(subset(df, df$META_TYP_MATCH == T)) / nrow(df)
out_table$IMPUTED[2] = nrow(subset(df, df$META_IMP_MATCH == T)) / nrow(df)

write.csv(out_table, "~/GitCode/MitoImputePrep/metadata/ADNI_concordance_HaploGrep.csv", quote = F, row.names = F)
