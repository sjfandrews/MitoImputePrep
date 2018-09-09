library(tidyverse)
library(readxl)
library(HiMC); data(nodes)

## READ IN THE LIST OF CHIPS
chips = read.table("~/GitCode/MitoImputePrep/scripts/INFORMATION_LISTS/b37_platforms.txt")
names(chips) = c("array")

info.cutoff = 0.3

## MAKE DATA FRAMES FOR EACH MINOR ALLELE FREQ
MAF1pc = chips
MAF0.5pc = chips
MAF0.1pc = chips

# SET A NEW COLUMN DESIGNATING THE MAF AND ONE WHETHER THEY HAVE BEEN IMPUTED
MAF1pc$MAF = "0.1"
MAF1pc$imputed = FALSE
MAF0.5pc$MAF = "0.05"
MAF0.5pc$imputed = FALSE
MAF0.1pc$MAF = "0.01"
MAF0.1pc$imputed = FALSE

## IF AN ARRAY HAD IMPUTATION PERFORMED ON IT, ASSIGN VALUE TRUE TO imputed COLUMN
for (i in 1:length(chips$array)) {
  DIR = "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/" #BDCHP-1X10-HUMANHAP240S_11216501_A-b37/ReferencePanel_v2/"
  if (file.exists(paste0(DIR, chips$array[i], "/ReferencePanel_v2/chrMT_1kg_", chips$array[i], "_imputed_info"))) {
    MAF1pc$imputed[i] = T
  }
}

i = 1

## READ IN THE INFO SCORE
info.score <- read_delim(paste0(DIR, chips$array[i], "/ReferencePanel_v2/chrMT_1kg_", chips$array[i], "_imputed_info"), delim = " ")
exclude.info.score  = subset(info.score, !(info.score$info >= info.cutoff))
exclude.pos = as.character(exclude.info.score$position)

##  Readin .ped files 
full_1kGP = generate_snp_data("/Volumes/TimMcInerney/MitoImpute/data/PLINK/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.map",
                              "/Volumes/TimMcInerney/MitoImpute/data/PLINK/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.ped")
full_1kGP_cut = full_1kGP[, -which(names(full_1kGP) %in% exclude.pos)]
full_1kGP_hg = HiMC::getClassifications(full_1kGP_cut)

imputed_1kGP <- generate_snp_data(paste0(DIR, chips$array[i], "/ReferencePanel_v2/chrMT_1kg_", chips$array[i], "_imputed.map"),
                                  paste0(DIR, chips$array[i], "/ReferencePanel_v2/chrMT_1kg_", chips$array[i], "_imputed.ped"))
imputed_1kGP_cut = imputed_1kGP[, -which(names(imputed_1kGP) %in% exclude.pos)]
imputed_1kGP_hg = HiMC::getClassifications(imputed_1kGP_cut)

typed_1kGP <- generate_snp_data(paste0(DIR, chips$array[i], "/ReferencePanel_v2/chrMT_1kg_", chips$array[i], ".map"),
                                paste0(DIR, chips$array[i], "/ReferencePanel_v2/chrMT_1kg_", chips$array[i], ".ped"))
typed_1kGP_cut = typed_1kGP[, -which(names(typed_1kGP) %in% exclude.pos)]
typed_1kGP_hg = HiMC::getClassifications(typed_1kGP)

conc_typed = full_1kGP_hg$haplogroup == typed_1kGP_hg$haplogroup
conc_imputed = full_1kGP_hg$haplogroup == imputed_1kGP_hg$haplogroup

