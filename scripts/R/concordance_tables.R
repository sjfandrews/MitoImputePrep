library(tidyverse)
library(readxl)
library(HiMC); data(nodes)

## READ IN THE LIST OF CHIPS
chips = read.table("~/GitCode/MitoImputePrep/scripts/INFORMATION_LISTS/b37_platforms.txt")
names(chips) = c("array")

info.cutoff = 0.3

## READ IN .ped AND .map FILES FOR TRUTH SET 
full_1kGP = generate_snp_data("/Volumes/TimMcInerney/MitoImpute/data/PLINK/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.map",
                              "/Volumes/TimMcInerney/MitoImpute/data/PLINK/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.ped")
full_1kGP_cut = full_1kGP[, -(which(names(full_1kGP) %in% exclude.pos))]
full_1kGP_hg = HiMC::getClassifications(full_1kGP_cut)
full_index_exclude = which(names(full_1kGP) %in% exclude.pos)
if (length(full_index_exclude) > 0) {
  full_1kGP_cut = full_1kGP[, -full_index_exclude]
  full_1kGP_hg = HiMC::getClassifications(full_1kGP_cut)
} else {
  full_1kGP_hg = HiMC::getClassifications(full_1kGP)
}


### WORK ON THE THE ith IMPUTED SET
## MAKE DATA FRAMES FOR EACH MINOR ALLELE FREQ
MAF1pc = data.frame(matrix(ncol = 6, nrow = 0))
MAF0.5pc = data.frame(matrix(ncol = 6, nrow = 0))
MAF0.1pc = data.frame(matrix(ncol = 6, nrow = 0))

#MAF1pc = chips
#MAF0.5pc = chips
#MAF0.1pc = chips

# SET A NEW COLUMN DESIGNATING THE MAF AND ONE WHETHER THEY HAVE BEEN IMPUTED
MAF1pc$MAF = "0.1"
MAF1pc$imputed = FALSE
MAF0.5pc$MAF = "0.05"
MAF0.5pc$imputed = FALSE
MAF0.1pc$MAF = "0.01"
MAF0.1pc$imputed = FALSE

## IF AN ARRAY HAD IMPUTATION PERFORMED ON IT, ASSIGN VALUE TRUE TO imputed COLUMN
for (i in 1:length(chips$array)) {
  print(paste0(i, " / ", length(chips$array)))
  DIR = "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/"
  
  # REF PANEL v2 (MAF >= 1%)
  MAF1pc$array[i] = chips$array[i]
  MAF1pc$MAF[i] = "0.1"
  if (file.exists(paste0(DIR, chips$array[i], "/ReferencePanel_v2/chrMT_1kg_", chips$array[i], "_imputed_info"))) {
    
    MAF1pc$imputed[i] = T
    
    ## READ IN THE INFO SCORE
    info.score <- read_delim(paste0(DIR, chips$array[i], "/ReferencePanel_v2/chrMT_1kg_", chips$array[i], "_imputed_info"), delim = " ")
    exclude.info.score  = subset(info.score, !(info.score$info >= info.cutoff))
    exclude.pos = as.character(exclude.info.score$position)
    
    MAF1pc$SNP.Ref.Only[i] = nrow(subset(info.score, info.score$type == 0))
    MAF1pc$SNP.Ref.Samp[i] = nrow(subset(info.score, info.score$type == 2))
    MAF1pc$SNP.Samp.Only[i] = nrow(subset(info.score, info.score$type == 3))
    MAF1pc$Retained.After.Filt[i] = nrow(subset(info.score, info.score$info >= info.cutoff))
    
    ## READ IN .ped AND .map FILES FOR IMPUTED SET 
    imputed_1kGP <- generate_snp_data(paste0(DIR, chips$array[i], "/ReferencePanel_v2/chrMT_1kg_", chips$array[i], "_imputed.map"),
                                      paste0(DIR, chips$array[i], "/ReferencePanel_v2/chrMT_1kg_", chips$array[i], "_imputed.ped"))
    imputed_1kGP_cut = imputed_1kGP[, -which(names(imputed_1kGP) %in% exclude.pos)]
    if (length(imputed_index_exclude) > 0) {
      imputed_1kGP_cut = imputed_1kGP[, -imputed_index_exclude]
      imputed_1kGP_hg = HiMC::getClassifications(imputed_1kGP_cut)
    } else {
      imputed_1kGP_hg = HiMC::getClassifications(imputed_1kGP)
    }
    
    typed_1kGP <- generate_snp_data(paste0(DIR, chips$array[i], "/ReferencePanel_v2/chrMT_1kg_", chips$array[i], ".map"),
                                    paste0(DIR, chips$array[i], "/ReferencePanel_v2/chrMT_1kg_", chips$array[i], ".ped"))
    typed_index_exclude = which(names(typed_1kGP) %in% exclude.pos)
    if (length(typed_index_exclude) > 0) {
      typed_1kGP_cut = typed_1kGP[, -typed_index_exclude]
      typed_1kGP_hg = HiMC::getClassifications(typed_1kGP_cut)
    } else {
      typed_1kGP_hg = HiMC::getClassifications(typed_1kGP)
    }
    
    ## CALCULATE HAPLOGROUP CONCORDANCE
    conc_typed = full_1kGP_hg$haplogroup == typed_1kGP_hg$haplogroup
    conc_imputed = full_1kGP_hg$haplogroup == imputed_1kGP_hg$haplogroup
    conc_typed_pc = length(conc_typed[conc_typed == T]) / length(conc_typed)
    conc_imputed_pc = length(conc_imputed[conc_imputed == T]) / length(conc_imputed)
    MAF1pc$Typed.hg.Conc[i] = conc_typed_pc
    MAF1pc$Imputed.hg.Conc[i] = conc_imputed_pc
  } else {
    MAF1pc$imputed[i] = F
    MAF1pc$SNP.Ref.Only[i] = NA
    MAF1pc$SNP.Ref.Samp[i] = NA
    MAF1pc$SNP.Samp.Only[i] = NA
    MAF1pc$Retained.After.Filt[i] = NA
    MAF1pc$Typed.hg.Conc[i] = NA
    MAF1pc$Imputed.hg.Conc[i] = NA
  }
}

head(MAF1pc)

