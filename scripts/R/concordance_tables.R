library(tidyverse)
library(readxl)
library(HiMC); data(nodes)

## READ IN THE LIST OF CHIPS
chips = read.table("~/GitCode/MitoImputePrep/scripts/INFORMATION_LISTS/b37_platforms.txt")
names(chips) = c("array")

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

for (i in 1:length(chips$array)) {
  DIR = "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/" #BDCHP-1X10-HUMANHAP240S_11216501_A-b37/ReferencePanel_v2/"
  if (file.exists(paste0(DIR, chips$array[i], "/ReferencePanel_v2/chrMT_1kg_", chips$array[i], "_imputed_info"))) {
    MAF1pc$imputed[i] = T
  }
}

info.score <- read_delim(paste0(DIR, chips$array[1], "/ReferencePanel_v2/chrMT_1kg_", chips$array[1], "_imputed_info"), delim = " ")

##  Readin .ped files 
full_1kGP = generate_snp_data("/Volumes/TimMcInerney/MitoImpute/data/PLINK/1kGP_chrMT_SNPonly.map",
                              "/Volumes/TimMcInerney/MitoImpute/data/PLINK/1kGP_chrMT_SNPonly.ped")
full_1kGP_hg = HiMC::getClassifications(full_1kGP)

imputed_1kGP <- generate_snp_data(paste0(DIR, chips$array[1], "/ReferencePanel_v2/chrMT_1kg_", chips$array[1], "_imputed.map"),
                                  paste0(DIR, chips$array[1], "/ReferencePanel_v2/chrMT_1kg_", chips$array[1], "_imputed.ped"))
imputed_1kGP_hg = HiMC::getClassifications(imputed_1kGP)

typed_1kGP <- generate_snp_data(paste0(DIR, chips$array[1], "/ReferencePanel_v2/chrMT_1kg_", chips$array[1], ".map"),
                                paste0(DIR, chips$array[1], "/ReferencePanel_v2/chrMT_1kg_", chips$array[1], ".ped"))
typed_1kGP_hg = HiMC::getClassifications(typed_1kGP)

conc_typed = full_1kGP_hg$haplogroup == typed_1kGP_hg$haplogroup
conc_imputed = full_1kGP_hg$haplogroup == imputed_1kGP_hg$haplogroup

