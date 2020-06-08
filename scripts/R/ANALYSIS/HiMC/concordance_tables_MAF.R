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
full_1kGP_hg_FULL = HiMC::getClassifications(full_1kGP)


### WORK ON THE THE ith IMPUTED SET
## MAKE DATA FRAMES FOR EACH MINOR ALLELE FREQ

COLS = c("array", "MAF", "imputed", "mean.info", "SNP.Ref.Only", "SNP.Ref.Samp", "SNP.Samp.Only", "TOTAL", "Retained.After.Filt", "Typed.hg.Conc", "Imputed.hg.Conc")

#CONC_TABLE = data.frame(matrix(ncol = 6, nrow = ))
MAF1pc = data.frame(matrix(ncol = length(COLS), nrow = nrow(chips)))
MAF0.5pc = data.frame(matrix(ncol = length(COLS), nrow = nrow(chips)))
MAF0.1pc = data.frame(matrix(ncol = length(COLS), nrow = nrow(chips)))

names(MAF1pc) = COLS
names(MAF0.5pc) = COLS
names(MAF0.1pc) = COLS


#MAF1pc = chips
#MAF0.5pc = chips
#MAF0.1pc = chips

# SET A NEW COLUMN DESIGNATING THE MAF AND ONE WHETHER THEY HAVE BEEN IMPUTED
MAF1pc$MAF = "0.01"
MAF1pc$imputed = FALSE
MAF0.5pc$MAF = "0.005"
MAF0.5pc$imputed = FALSE
MAF0.1pc$MAF = "0.001"
MAF0.1pc$imputed = FALSE

# REF PANEL v2 (MAF >= 1%)
for (i in 1:length(chips$array)) {
  print(paste0(i, " / ", length(chips$array)))
  DIR = "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/"
  maf = "0.01"
  
  # REF PANEL v2 (MAF >= 1%)
  MAF1pc$array[i] = as.character(chips$array[i])
  MAF1pc$MAF[i] = maf
  if (file.exists(paste0(DIR, chips$array[i], "/ReferencePanel_v2/chrMT_1kg_", chips$array[i], "_imputed_info"))) {
    ## IF AN ARRAY HAD IMPUTATION PERFORMED ON IT, ASSIGN VALUE TRUE TO imputed COLUMN
    MAF1pc$imputed[i] = T
    
    ## READ IN THE INFO SCORE
    info.score <- read_delim(paste0(DIR, chips$array[i], "/ReferencePanel_v2/chrMT_1kg_", chips$array[i], "_imputed_info"), delim = " ")
    mean.info = mean(info.score$info)
    exclude.info.score = subset(info.score, !(info.score$info >= info.cutoff))
    exclude.pos = as.character(exclude.info.score$position)
  
  	MAF1pc$mean.info[i] = mean.info
    MAF1pc$SNP.Ref.Only[i] = nrow(subset(info.score, info.score$type == 0))
    MAF1pc$SNP.Ref.Samp[i] = nrow(subset(info.score, info.score$type == 2))
    MAF1pc$SNP.Samp.Only[i] = nrow(subset(info.score, info.score$type == 3))
    MAF1pc$TOTAL[i] = sum(MAF1pc$SNP.Ref.Only[i], MAF1pc$SNP.Ref.Samp[i], MAF1pc$SNP.Samp.Only[i])
    MAF1pc$Retained.After.Filt[i] = nrow(subset(info.score, info.score$info >= info.cutoff))
    
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
    imputed_1kGP <- generate_snp_data(paste0(DIR, chips$array[i], "/ReferencePanel_v2/chrMT_1kg_", chips$array[i], "_imputed.map"),
                                      paste0(DIR, chips$array[i], "/ReferencePanel_v2/chrMT_1kg_", chips$array[i], "_imputed.ped"))
    imputed_index_exclude = which(names(imputed_1kGP) %in% exclude.pos)
    if (length(imputed_index_exclude) > 0) {
      imputed_1kGP_cut = imputed_1kGP[, -imputed_index_exclude]
      imputed_1kGP_hg = HiMC::getClassifications(imputed_1kGP_cut)
    } else {
      imputed_1kGP_hg = HiMC::getClassifications(imputed_1kGP)
    }
    message("IMPUTED HAPLOTYPES CLASSIFIED")
    
    typed_1kGP <- generate_snp_data(paste0(DIR, chips$array[i], "/ReferencePanel_v2/chrMT_1kg_", chips$array[i], ".map"),
                                    paste0(DIR, chips$array[i], "/ReferencePanel_v2/chrMT_1kg_", chips$array[i], ".ped"))
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
write.csv(MAF1pc, "/Users/u5015730/GitCode/MitoImputePrep/metadata/ConcordanceTables_MAF1pc.csv", quote = F, row.names = F)

# REF PANEL v3 (MAF >= 0.1%)
refPanel = "ReferencePanel_v3"
for (i in 1:length(chips$array)) {
  print(paste0(i, " / ", length(chips$array)))
  DIR = "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/"
  maf = "0.001"
  
  # REF PANEL v3 (MAF >= 0.1%)
  MAF0.1pc$array[i] = as.character(chips$array[i])
  MAF0.1pc$MAF[i] = maf
  if (file.exists(paste0(DIR, chips$array[i], "/ReferencePanel_v3/chrMT_1kg_", chips$array[i], "_imputed_info"))) {
    ## IF AN ARRAY HAD IMPUTATION PERFORMED ON IT, ASSIGN VALUE TRUE TO imputed COLUMN
    MAF0.1pc$imputed[i] = T
    
    ## READ IN THE INFO SCORE
    info.score <- read_delim(paste0(DIR, chips$array[i], "/", refPanel, "/chrMT_1kg_", chips$array[i], "_imputed_info"), delim = " ")
    mean.info = mean(info.score$info)
    exclude.info.score = subset(info.score, !(info.score$info >= info.cutoff))
    exclude.pos = as.character(exclude.info.score$position)
    
    MAF0.1pc$mean.info[i] = mean.info
    MAF0.1pc$SNP.Ref.Only[i] = nrow(subset(info.score, info.score$type == 0))
    MAF0.1pc$SNP.Ref.Samp[i] = nrow(subset(info.score, info.score$type == 2))
    MAF0.1pc$SNP.Samp.Only[i] = nrow(subset(info.score, info.score$type == 3))
    MAF0.1pc$TOTAL[i] = sum(MAF0.1pc$SNP.Ref.Only[i], MAF0.1pc$SNP.Ref.Samp[i], MAF0.1pc$SNP.Samp.Only[i])
    MAF0.1pc$Retained.After.Filt[i] = nrow(subset(info.score, info.score$info >= info.cutoff))
    
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
    imputed_1kGP <- generate_snp_data(paste0(DIR, chips$array[i], "/", refPanel, "/chrMT_1kg_", chips$array[i], "_imputed.map"),
                                      paste0(DIR, chips$array[i], "/", refPanel, "/chrMT_1kg_", chips$array[i], "_imputed.ped"))
    imputed_index_exclude = which(names(imputed_1kGP) %in% exclude.pos)
    if (length(imputed_index_exclude) > 0) {
      imputed_1kGP_cut = imputed_1kGP[, -imputed_index_exclude]
      imputed_1kGP_hg = HiMC::getClassifications(imputed_1kGP_cut)
    } else {
      imputed_1kGP_hg = HiMC::getClassifications(imputed_1kGP)
    }
    message("IMPUTED HAPLOTYPES CLASSIFIED")
    
    typed_1kGP <- generate_snp_data(paste0(DIR, chips$array[i], "/", refPanel, "/chrMT_1kg_", chips$array[i], ".map"),
                                    paste0(DIR, chips$array[i], "/", refPanel, "/chrMT_1kg_", chips$array[i], ".ped"))
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
    MAF0.1pc$Typed.hg.Conc[i] = conc_typed_pc
    MAF0.1pc$Imputed.hg.Conc[i] = conc_imputed_pc
  } else {
    MAF0.1pc$imputed[i] = F
    MAF0.1pc$SNP.Ref.Only[i] = NA
    MAF0.1pc$SNP.Ref.Samp[i] = NA
    MAF0.1pc$SNP.Samp.Only[i] = NA
    MAF0.1pc$Retained.After.Filt[i] = NA
    MAF0.1pc$Typed.hg.Conc[i] = NA
    MAF0.1pc$Imputed.hg.Conc[i] = NA
  }
}
write.csv(MAF0.1pc, "/Users/u5015730/GitCode/MitoImputePrep/metadata/ConcordanceTables_MAF0-1pc.csv", quote = F, row.names = F)

# REF PANEL v4 (MAF >= 0.5%)
refPanel = "ReferencePanel_v4"
for (i in 1:length(chips$array)) {
  print(paste0(i, " / ", length(chips$array)))
  DIR = "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/"
  maf = "0.005"
  
  # REF PANEL v4 (MAF >= 0.5%)
  MAF0.5pc$array[i] = as.character(chips$array[i])
  MAF0.5pc$MAF[i] = maf
  if (file.exists(paste0(DIR, chips$array[i], "/", refPanel, "/chrMT_1kg_", chips$array[i], "_imputed_info"))) {
    ## IF AN ARRAY HAD IMPUTATION PERFORMED ON IT, ASSIGN VALUE TRUE TO imputed COLUMN
    MAF0.5pc$imputed[i] = T
    
    ## READ IN THE INFO SCORE
    info.score <- read_delim(paste0(DIR, chips$array[i], "/", refPanel, "/chrMT_1kg_", chips$array[i], "_imputed_info"), delim = " ")
    mean.info = mean(info.score$info)
    exclude.info.score = subset(info.score, !(info.score$info >= info.cutoff))
    exclude.pos = as.character(exclude.info.score$position)
    
    MAF0.5pc$mean.info[i] = mean.info
    MAF0.5pc$SNP.Ref.Only[i] = nrow(subset(info.score, info.score$type == 0))
    MAF0.5pc$SNP.Ref.Samp[i] = nrow(subset(info.score, info.score$type == 2))
    MAF0.5pc$SNP.Samp.Only[i] = nrow(subset(info.score, info.score$type == 3))
    MAF0.5pc$TOTAL[i] = sum(MAF0.5pc$SNP.Ref.Only[i], MAF0.5pc$SNP.Ref.Samp[i], MAF0.5pc$SNP.Samp.Only[i])
    MAF0.5pc$Retained.After.Filt[i] = nrow(subset(info.score, info.score$info >= info.cutoff))
    
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
    imputed_1kGP <- generate_snp_data(paste0(DIR, chips$array[i], "/", refPanel, "/chrMT_1kg_", chips$array[i], "_imputed.map"),
                                      paste0(DIR, chips$array[i], "/", refPanel, "/chrMT_1kg_", chips$array[i], "_imputed.ped"))
    imputed_index_exclude = which(names(imputed_1kGP) %in% exclude.pos)
    if (length(imputed_index_exclude) > 0) {
      imputed_1kGP_cut = imputed_1kGP[, -imputed_index_exclude]
      imputed_1kGP_hg = HiMC::getClassifications(imputed_1kGP_cut)
    } else {
      imputed_1kGP_hg = HiMC::getClassifications(imputed_1kGP)
    }
    message("IMPUTED HAPLOTYPES CLASSIFIED")
    
    typed_1kGP <- generate_snp_data(paste0(DIR, chips$array[i], "/", refPanel, "/chrMT_1kg_", chips$array[i], ".map"),
                                    paste0(DIR, chips$array[i], "/", refPanel, "/chrMT_1kg_", chips$array[i], ".ped"))
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
    MAF0.5pc$Typed.hg.Conc[i] = conc_typed_pc
    MAF0.5pc$Imputed.hg.Conc[i] = conc_imputed_pc
  } else {
    MAF0.5pc$imputed[i] = F
    MAF0.5pc$SNP.Ref.Only[i] = NA
    MAF0.5pc$SNP.Ref.Samp[i] = NA
    MAF0.5pc$SNP.Samp.Only[i] = NA
    MAF0.5pc$Retained.After.Filt[i] = NA
    MAF0.5pc$Typed.hg.Conc[i] = NA
    MAF0.5pc$Imputed.hg.Conc[i] = NA
  }
}
write.csv(MAF0.5pc, "/Users/u5015730/GitCode/MitoImputePrep/metadata/ConcordanceTables_MAF0-5pc.csv", quote = F, row.names = F)

COMB = rbind(MAF1pc, MAF0.5pc, MAF0.1pc)
write.csv(COMB, "/Users/u5015730/GitCode/MitoImputePrep/metadata/ConcordanceTables_MAF_Combined.csv", quote = F, row.names = F)
