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

COLS = c("array", "MCMC", "imputed", "mean.info", "SNP.Ref.Only", "SNP.Ref.Samp", "SNP.Samp.Only", "TOTAL", "Retained.After.Filt", "Typed.hg.Conc", "Imputed.hg.Conc")

#CONC_TABLE = data.frame(matrix(ncol = 6, nrow = ))
MCMC1 = data.frame(matrix(ncol = length(COLS), nrow = nrow(chips)))
MCMC5 = data.frame(matrix(ncol = length(COLS), nrow = nrow(chips)))
MCMC10 = data.frame(matrix(ncol = length(COLS), nrow = nrow(chips)))
MCMC20 = data.frame(matrix(ncol = length(COLS), nrow = nrow(chips)))

names(MCMC1) = COLS
names(MCMC5) = COLS
names(MCMC10) = COLS
names(MCMC20) = COLS


#MAF1pc = chips
#MAF0.5pc = chips
#MAF0.1pc = chips

# SET A NEW COLUMN DESIGNATING THE MAF AND ONE WHETHER THEY HAVE BEEN IMPUTED
MCMC1$MCMC = "1"
MCMC1$imputed = FALSE
MCMC5$MCMC = "5"
MCMC5$imputed = FALSE
MCMC10$MCMC = "10"
MCMC10$imputed = FALSE
MCMC20$MCMC = "20"
MCMC20$imputed = FALSE

# MCMC = 1
for (i in 1:length(chips$array)) { 
  print(paste0(i, " / ", length(chips$array)))
  DIR = "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/"
  mcmc = "1"
  
  # REF PANEL v2 (MAF >= 1%)
  MCMC1$array[i] = as.character(chips$array[i])
  MCMC1$MCMC[i] = mcmc
  if (file.exists(paste0(DIR, chips$array[i], "/MCMC_Experiments/MCMC1/chrMT_1kg_", chips$array[i], "_imputed_MCMC", mcmc,"_info"))) {
    ## IF AN ARRAY HAD IMPUTATION PERFORMED ON IT, ASSIGN VALUE TRUE TO imputed COLUMN
    MCMC1$imputed[i] = T
    
    ## READ IN THE INFO SCORE
    info.score <- read_delim(paste0(DIR, chips$array[i], "/MCMC_Experiments/MCMC1/chrMT_1kg_", chips$array[i], "_imputed_MCMC", mcmc,"_info"), delim = " ")
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
    message(paste("READING IN .ped AND .map FILES FOR ", chips$array[i], " MCMC = ", mcmc))
    imputed_1kGP <- generate_snp_data(paste0(DIR, chips$array[i], "/MCMC_Experiments/MCMC1/chrMT_1kg_", chips$array[i], "_imputed_MCMC", mcmc, ".map"),
                                      paste0(DIR, chips$array[i], "/MCMC_Experiments/MCMC1/chrMT_1kg_", chips$array[i], "_imputed_MCMC", mcmc, ".ped"))
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
    MCMC1$Typed.hg.Conc[i] = conc_typed_pc
    MCMC1$Imputed.hg.Conc[i] = conc_imputed_pc
  } else {
    MCMC1$imputed[i] = F
    MCMC1$SNP.Ref.Only[i] = NA
    MCMC1$SNP.Ref.Samp[i] = NA
    MCMC1$SNP.Samp.Only[i] = NA
    MCMC1$Retained.After.Filt[i] = NA
    MCMC1$Typed.hg.Conc[i] = NA
    MCMC1$Imputed.hg.Conc[i] = NA
  }
}
write.csv(MCMC1, "~/GitCode/MitoImputePrep/metadata/ConcordanceTables_MCMC1.csv", quote = F, row.names = F)

# MCMC = 5
for (i in 1:length(chips$array)) { 
  print(paste0(i, " / ", length(chips$array)))
  DIR = "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/"
  mcmc = "5"
  
  # REF PANEL v2 (MAF >= 1%)
  MCMC5$array[i] = as.character(chips$array[i])
  MCMC5$MCMC[i] = mcmc
  if (file.exists(paste0(DIR, chips$array[i], "/MCMC_Experiments/MCMC5/chrMT_1kg_", chips$array[i], "_imputed_MCMC", mcmc,"_info"))) {
    ## IF AN ARRAY HAD IMPUTATION PERFORMED ON IT, ASSIGN VALUE TRUE TO imputed COLUMN
    MCMC5$imputed[i] = T
    
    ## READ IN THE INFO SCORE
    info.score <- read_delim(paste0(DIR, chips$array[i], "/MCMC_Experiments/MCMC5/chrMT_1kg_", chips$array[i], "_imputed_MCMC", mcmc,"_info"), delim = " ")
    mean.info = mean(info.score$info)
    exclude.info.score = subset(info.score, !(info.score$info >= info.cutoff))
    exclude.pos = as.character(exclude.info.score$position)
    
    MCMC5$mean.info[i] = mean.info
    MCMC5$SNP.Ref.Only[i] = nrow(subset(info.score, info.score$type == 0))
    MCMC5$SNP.Ref.Samp[i] = nrow(subset(info.score, info.score$type == 2))
    MCMC5$SNP.Samp.Only[i] = nrow(subset(info.score, info.score$type == 3))
    MCMC5$TOTAL[i] = sum(MCMC5$SNP.Ref.Only[i], MCMC5$SNP.Ref.Samp[i], MCMC5$SNP.Samp.Only[i])
    MCMC5$Retained.After.Filt[i] = nrow(subset(info.score, info.score$info >= info.cutoff))
    
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
    message(paste("READING IN .ped AND .map FILES FOR ", chips$array[i], " MCMC = ", mcmc))
    imputed_1kGP <- generate_snp_data(paste0(DIR, chips$array[i], "/MCMC_Experiments/MCMC5/chrMT_1kg_", chips$array[i], "_imputed_MCMC", mcmc, ".map"),
                                      paste0(DIR, chips$array[i], "/MCMC_Experiments/MCMC5/chrMT_1kg_", chips$array[i], "_imputed_MCMC", mcmc, ".ped"))
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
    MCMC5$Typed.hg.Conc[i] = conc_typed_pc
    MCMC5$Imputed.hg.Conc[i] = conc_imputed_pc
  } else {
    MCMC5$imputed[i] = F
    MCMC5$SNP.Ref.Only[i] = NA
    MCMC5$SNP.Ref.Samp[i] = NA
    MCMC5$SNP.Samp.Only[i] = NA
    MCMC5$Retained.After.Filt[i] = NA
    MCMC5$Typed.hg.Conc[i] = NA
    MCMC5$Imputed.hg.Conc[i] = NA
  }
}
write.csv(MCMC5, "~/GitCode/MitoImputePrep/metadata/ConcordanceTables_MCMC5.csv", quote = F, row.names = F)

# MCMC = 10
for (i in 1:length(chips$array)) { 
  print(paste0(i, " / ", length(chips$array)))
  DIR = "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/"
  mcmc = mcmc
  
  # REF PANEL v2 (MAF >= 1%)
  MCMC10$array[i] = as.character(chips$array[i])
  MCMC10$MCMC[i] = "1"
  if (file.exists(paste0(DIR, chips$array[i], "/MCMC_Experiments/MCMC10/chrMT_1kg_", chips$array[i], "_imputed_MCMC", mcmc,"_info"))) {
    ## IF AN ARRAY HAD IMPUTATION PERFORMED ON IT, ASSIGN VALUE TRUE TO imputed COLUMN
    MCMC10$imputed[i] = T
    
    ## READ IN THE INFO SCORE
    info.score <- read_delim(paste0(DIR, chips$array[i], "/MCMC_Experiments/MCMC10/chrMT_1kg_", chips$array[i], "_imputed_MCMC", mcmc,"_info"), delim = " ")
    mean.info = mean(info.score$info)
    exclude.info.score = subset(info.score, !(info.score$info >= info.cutoff))
    exclude.pos = as.character(exclude.info.score$position)
    
    MCMC10$mean.info[i] = mean.info
    MCMC10$SNP.Ref.Only[i] = nrow(subset(info.score, info.score$type == 0))
    MCMC10$SNP.Ref.Samp[i] = nrow(subset(info.score, info.score$type == 2))
    MCMC10$SNP.Samp.Only[i] = nrow(subset(info.score, info.score$type == 3))
    MCMC10$TOTAL[i] = sum(MCMC10$SNP.Ref.Only[i], MCMC10$SNP.Ref.Samp[i], MCMC10$SNP.Samp.Only[i])
    MCMC10$Retained.After.Filt[i] = nrow(subset(info.score, info.score$info >= info.cutoff))
    
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
    message(paste("READING IN .ped AND .map FILES FOR ", chips$array[i], " MCMC = ", mcmc))
    imputed_1kGP <- generate_snp_data(paste0(DIR, chips$array[i], "/MCMC_Experiments/MCMC10/chrMT_1kg_", chips$array[i], "_imputed_MCMC", mcmc, ".map"),
                                      paste0(DIR, chips$array[i], "/MCMC_Experiments/MCMC10/chrMT_1kg_", chips$array[i], "_imputed_MCMC", mcmc, ".ped"))
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
    MCMC10$Typed.hg.Conc[i] = conc_typed_pc
    MCMC10$Imputed.hg.Conc[i] = conc_imputed_pc
  } else {
    MCMC10$imputed[i] = F
    MCMC10$SNP.Ref.Only[i] = NA
    MCMC10$SNP.Ref.Samp[i] = NA
    MCMC10$SNP.Samp.Only[i] = NA
    MCMC10$Retained.After.Filt[i] = NA
    MCMC10$Typed.hg.Conc[i] = NA
    MCMC10$Imputed.hg.Conc[i] = NA
  }
}
write.csv(MCMC10, "~/GitCode/MitoImputePrep/metadata/ConcordanceTables_MCMC10.csv", quote = F, row.names = F)

# MCMC = 20
for (i in 1:length(chips$array)) { 
  print(paste0(i, " / ", length(chips$array)))
  DIR = "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/"
  mcmc = mcmc
  
  # REF PANEL v2 (MAF >= 1%)
  MCMC20$array[i] = as.character(chips$array[i])
  MCMC20$MCMC[i] = "1"
  if (file.exists(paste0(DIR, chips$array[i], "/MCMC_Experiments/MCMC20/chrMT_1kg_", chips$array[i], "_imputed_MCMC", mcmc,"_info"))) {
    ## IF AN ARRAY HAD IMPUTATION PERFORMED ON IT, ASSIGN VALUE TRUE TO imputed COLUMN
    MCMC20$imputed[i] = T
    
    ## READ IN THE INFO SCORE
    info.score <- read_delim(paste0(DIR, chips$array[i], "/MCMC_Experiments/MCMC20/chrMT_1kg_", chips$array[i], "_imputed_MCMC", mcmc,"_info"), delim = " ")
    mean.info = mean(info.score$info)
    exclude.info.score = subset(info.score, !(info.score$info >= info.cutoff))
    exclude.pos = as.character(exclude.info.score$position)
    
    MCMC20$mean.info[i] = mean.info
    MCMC20$SNP.Ref.Only[i] = nrow(subset(info.score, info.score$type == 0))
    MCMC20$SNP.Ref.Samp[i] = nrow(subset(info.score, info.score$type == 2))
    MCMC20$SNP.Samp.Only[i] = nrow(subset(info.score, info.score$type == 3))
    MCMC20$TOTAL[i] = sum(MCMC20$SNP.Ref.Only[i], MCMC20$SNP.Ref.Samp[i], MCMC20$SNP.Samp.Only[i])
    MCMC20$Retained.After.Filt[i] = nrow(subset(info.score, info.score$info >= info.cutoff))
    
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
    message(paste("READING IN .ped AND .map FILES FOR ", chips$array[i], " MCMC = ", mcmc))
    imputed_1kGP <- generate_snp_data(paste0(DIR, chips$array[i], "/MCMC_Experiments/MCMC20/chrMT_1kg_", chips$array[i], "_imputed_MCMC", mcmc, ".map"),
                                      paste0(DIR, chips$array[i], "/MCMC_Experiments/MCMC20/chrMT_1kg_", chips$array[i], "_imputed_MCMC", mcmc, ".ped"))
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
    MCMC20$Typed.hg.Conc[i] = conc_typed_pc
    MCMC20$Imputed.hg.Conc[i] = conc_imputed_pc
  } else {
    MCMC20$imputed[i] = F
    MCMC20$SNP.Ref.Only[i] = NA
    MCMC20$SNP.Ref.Samp[i] = NA
    MCMC20$SNP.Samp.Only[i] = NA
    MCMC20$Retained.After.Filt[i] = NA
    MCMC20$Typed.hg.Conc[i] = NA
    MCMC20$Imputed.hg.Conc[i] = NA
  }
}
write.csv(MCMC20, "~/GitCode/MitoImputePrep/metadata/ConcordanceTables_MCMC20.csv", quote = F, row.names = F)

COMB = rbind(MCMC1, MCMC5, MCMC10, MCMC20)
write.csv(COMB, "/Users/u5015730/GitCode/MitoImputePrep/metadata/ConcordanceTables_MCMC_Combined.csv", quote = F, row.names = F)
