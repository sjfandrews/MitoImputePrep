library(tidyverse)
library(readxl)
library(HiMC); data(nodes)
library(dplyr)

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

COLS = c("array", "MCMC", "imputed", "mean.info", "SNP.Ref.Only", "SNP.Ref.Samp", "SNP.Samp.Only", "TOTAL", "Retained.After.Filt", "Typed.hg.Conc", "Imputed.hg.Conc")

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

## CALCULATE HAPLOGROUP CONCORDANCE
conc_typed = full_1kGP_hg$haplogroup == typed_1kGP_hg$haplogroup
conc_imputed = full_1kGP_hg$haplogroup == imputed_1kGP_hg$haplogroup
conc_typed_pc = length(conc_typed[conc_typed == T]) / length(conc_typed)
conc_imputed_pc = length(conc_imputed[conc_imputed == T]) / length(conc_imputed)
MCMC1$Typed.hg.Conc[i] = conc_typed_pc
MCMC1$Imputed.hg.Conc[i] = conc_imputed_pc

write.csv(MCMC1, "~/GitCode/MitoImputePrep/metadata/ADNI_concordance_HiMC.csv", quote = F, row.names = F)

rect_csv = read.csv("/Users/u5015730/GitCode/MitoImputePrep/metadata/ADNI_samples_BOTH.csv", header = T)
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

out_table = data.frame(matrix(ncol = 2, nrow = 1))
names(out_table) = c("TYPED", "IMPUTED")

out_table$TYPED = nrow(subset(df, df$TYP_MATCH == T)) / nrow(df)
out_table$IMPUTED = nrow(subset(df, df$IMP_MATCH == T)) / nrow(df)

