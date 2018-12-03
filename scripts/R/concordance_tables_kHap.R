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

COLS = c("array", "kHAP", "imputed", "SNP.Ref.Only", "SNP.Ref.Samp", "SNP.Samp.Only", "TOTAL", "Retained.After.Filt", "Typed.hg.Conc", "Imputed.hg.Conc")

#CONC_TABLE = data.frame(matrix(ncol = 6, nrow = ))
kHAP100 = data.frame(matrix(ncol = length(COLS), nrow = nrow(chips)))
kHAP250 = data.frame(matrix(ncol = length(COLS), nrow = nrow(chips)))
kHAP500 = data.frame(matrix(ncol = length(COLS), nrow = nrow(chips)))
kHAP1000 = data.frame(matrix(ncol = length(COLS), nrow = nrow(chips)))
kHAP2500 = data.frame(matrix(ncol = length(COLS), nrow = nrow(chips)))
kHAP5000 = data.frame(matrix(ncol = length(COLS), nrow = nrow(chips)))
kHAP10000 = data.frame(matrix(ncol = length(COLS), nrow = nrow(chips)))
kHAP20000 = data.frame(matrix(ncol = length(COLS), nrow = nrow(chips)))
kHAP30000 = data.frame(matrix(ncol = length(COLS), nrow = nrow(chips)))

names(kHAP100) = COLS
names(kHAP250) = COLS
names(kHAP500) = COLS
names(kHAP1000) = COLS
names(kHAP2500) = COLS
names(kHAP5000) = COLS
names(kHAP10000) = COLS
names(kHAP20000) = COLS
names(kHAP30000) = COLS


#MAF1pc = chips
#MAF0.5pc = chips
#MAF0.1pc = chips

# SET A NEW COLUMN DESIGNATING THE MAF AND ONE WHETHER THEY HAVE BEEN IMPUTED
kHAP100$kHAP = "100"
kHAP100$imputed = FALSE
kHAP250$kHAP = "250"
kHAP250$imputed = FALSE
kHAP500$kHAP = "500"
kHAP500$imputed = FALSE
kHAP1000$kHAP = "1000"
kHAP1000$imputed = FALSE
kHAP2500$kHAP = "2500"
kHAP2500$imputed = FALSE
kHAP5000$kHAP = "5000"
kHAP5000$imputed = FALSE
kHAP10000$kHAP = "10000"
kHAP10000$imputed = FALSE
kHAP20000$kHAP = "20000"
kHAP20000$imputed = FALSE
kHAP30000$kHAP = "30000"
kHAP30000$imputed = FALSE


# kHAP100 = 1
for (i in 1:length(chips$array)) { 
  print(paste0(i, " / ", length(chips$array)))
  DIR = "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/"
  k.hap = "100"
  
  # REF PANEL v2 (MAF >= 1%)
  kHAP100$array[i] = as.character(chips$array[i])
  kHAP100$kHAP[i] = k.hap
  if (file.exists(paste0(DIR, chips$array[i], "/kHAP_Experiments/kHAP100/chrMT_1kg_", chips$array[i], "_imputed_kHAP", k.hap,"_info"))) {
    ## IF AN ARRAY HAD IMPUTATION PERFORMED ON IT, ASSIGN VALUE TRUE TO imputed COLUMN
    kHAP100$imputed[i] = T
    
    ## READ IN THE INFO SCORE
    info.score <- read_delim(paste0(DIR, chips$array[i], "/kHAP_Experiments/kHAP100/chrMT_1kg_", chips$array[i], "_imputed_kHAP", k.hap,"_info"), delim = " ")
    exclude.info.score  = subset(info.score, !(info.score$info >= info.cutoff))
    exclude.pos = as.character(exclude.info.score$position)
    
    kHAP100$SNP.Ref.Only[i] = nrow(subset(info.score, info.score$type == 0))
    kHAP100$SNP.Ref.Samp[i] = nrow(subset(info.score, info.score$type == 2))
    kHAP100$SNP.Samp.Only[i] = nrow(subset(info.score, info.score$type == 3))
    kHAP100$TOTAL[i] = sum(kHAP100$SNP.Ref.Only[i], kHAP100$SNP.Ref.Samp[i], kHAP100$SNP.Samp.Only[i])
    kHAP100$Retained.After.Filt[i] = nrow(subset(info.score, info.score$info >= info.cutoff))
    
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
    message(paste("READING IN .ped AND .map FILES FOR ", chips$array[i], " k_hap = ", k.hap))
    imputed_1kGP <- generate_snp_data(paste0(DIR, chips$array[i], "/kHAP_Experiments/kHAP100/chrMT_1kg_", chips$array[i], "_imputed_kHAP", k.hap, ".map"),
                                      paste0(DIR, chips$array[i], "/kHAP_Experiments/kHAP100/chrMT_1kg_", chips$array[i], "_imputed_kHAP", k.hap, ".ped"))
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
    kHAP100$Typed.hg.Conc[i] = conc_typed_pc
    kHAP100$Imputed.hg.Conc[i] = conc_imputed_pc
  } else {
    kHAP100$imputed[i] = F
    kHAP100$SNP.Ref.Only[i] = NA
    kHAP100$SNP.Ref.Samp[i] = NA
    kHAP100$SNP.Samp.Only[i] = NA
    kHAP100$Retained.After.Filt[i] = NA
    kHAP100$Typed.hg.Conc[i] = NA
    kHAP100$Imputed.hg.Conc[i] = NA
  }
}
write.csv(kHAP100, "~/GitCode/MitoImputePrep/metadata/ConcordanceTables_kHAP100.csv", quote = F, row.names = F)



COMB = rbind(kHAP100,
             kHAP250,
             kHAP500,
             kHAP1000,
             kHAP2500,
             kHAP5000,
             kHAP10000,
             kHAP20000,
             kHAP30000)

write.csv(COMB, "/Users/u5015730/GitCode/MitoImputePrep/metadata/ConcordanceTables_kHAP_Combined.csv", quote = F, row.names = F)
