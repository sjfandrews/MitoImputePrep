library(tidyverse)
library(readxl)
library(HiMC); data(nodes)
library(dplyr)

'%!in%' = function(x,y)!('%in%'(x,y))
difficult_mhgs = c("HV", "JT", "unclassified")

## READ IN THE LIST OF CHIPS
chips = read.table("~/GitCode/MitoImputePrep/scripts/INFORMATION_LISTS/b37_platforms.txt")
names(chips) = c("array")

info.cutoff = 0.3

## READ IN .ped AND .map FILES FOR TRUTH SET 
full_1kGP = generate_snp_data("/Volumes/TimMcInerney/MitoImpute/data/PLINK/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.map",
                              "/Volumes/TimMcInerney/MitoImpute/data/PLINK/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.ped")
full_1kGP_hg_FULL = HiMC::getClassifications(full_1kGP)

full_1kGP_meta = full_1kGP_hg_FULL
full_1kGP_meta_D = subset(full_1kGP_meta, full_1kGP_meta$haplogroup %in% difficult_mhgs)
full_1kGP_meta_D$metahaplogroup = full_1kGP_meta_D$haplogroup
full_1kGP_meta = subset(full_1kGP_meta, full_1kGP_meta$haplogroup %!in% difficult_mhgs)
full_1kGP_meta_A = subset(full_1kGP_meta, substr(full_1kGP_meta$haplogroup, 1 ,1) == "L")
full_1kGP_meta_A$metahaplogroup = substr(full_1kGP_meta_A$haplogroup, 1, 2)
full_1kGP_meta_N = subset(full_1kGP_meta, substr(full_1kGP_meta$haplogroup, 1 ,1) != "L")
full_1kGP_meta_N$metahaplogroup = substr(full_1kGP_meta_N$haplogroup, 1, 1)

full_1kGP_meta = rbind(full_1kGP_meta_D, full_1kGP_meta_N, full_1kGP_meta_A)

hg.df = data.frame(cbind(full_1kGP_hg_FULL$Individual, full_1kGP_hg_FULL$haplogroup))
names(hg.df) = c("Individual", "WGS")
hg.df$chip = NA
hg.df$typed = NA
hg.df$imputed = NA
hg.df$match.t = NA
hg.df$match.i = NA

for (i in 1:length(chips$array)) { # length(chips$array)
  print(paste0(i, " / ", length(chips$array)))
  DIR = "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/"
  mcmc = "1"
  hg.df$chip = chips$array[i]
  
  # REF PANEL v2 (MAF >= 1%)
  if (file.exists(paste0(DIR, chips$array[i], "/MCMC_Experiments/MCMC1/chrMT_1kg_", chips$array[i], "_imputed_MCMC", mcmc,"_info"))) {
    ## READ IN THE INFO SCORE
    info.score <- read_delim(paste0(DIR, chips$array[i], "/MCMC_Experiments/MCMC1/chrMT_1kg_", chips$array[i], "_imputed_MCMC", mcmc,"_info"), delim = " ")
    mean.info = mean(info.score$info)
    exclude.info.score = subset(info.score, !(info.score$info >= info.cutoff))
    exclude.pos = as.character(exclude.info.score$position)
    
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
    
    ## MERGE TABLES
    for (j in 1:nrow(hg.df)) {
      hg.df$typed[j] = typed_1kGP_hg$haplogroup[match(hg.df$Individual[j], typed_1kGP_hg$Individual)]
      hg.df$imputed[j] = imputed_1kGP_hg$haplogroup[match(hg.df$Individual[j], imputed_1kGP_hg$Individual)]
      hg.df$match.t[j] = hg.df$WGS[j] == hg.df$typed[j]
      hg.df$match.i[j] = hg.df$WGS[j] == hg.df$imputed[j]
    }
    write.csv(hg.df, paste0("~/GitCode/MitoImputePrep/metadata/haplogroup_concordance/", chips$array[i], "_haplogroupings.csv"))
  }
}
