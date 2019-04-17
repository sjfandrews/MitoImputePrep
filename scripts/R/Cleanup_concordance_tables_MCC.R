require(ggplot2)
require(gridExtra)
require(tidyr)
require(emmeans)
require(dplyr)

#########################################################################
chip.table = read.table("~/GitCode/MitoImputePrep/scripts/INFORMATION_LISTS/b37_platforms.txt", header = F)
#truth.table = read.table("~/GitCode/MitoImputePrep/metadata/HaploGrep_concordance/chrMT_1kg_diploid_haplogrep.txt", header = T)
container = "/g/data1a/te53/MitoImpute/data/STRANDS/" #BDCHP-1X10-HUMANHAP240S_11216501_A-b37/MCMC_Experiments/MCMC1"
outward_dir = "/g/data1a/te53/MitoImpute/analyses/concordance/MCC/"
#container = "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/"
#outward_dir = "/Volumes/TimMcInerney/MitoImpute/analyses/concordance/MCC/"

###############################################################################
## MCMC 

exp.dir = "MCMC_Experiments"
exp.var = c("MCMC1", "MCMC5", "MCMC10", "MCMC20", "MCMC30")
tmp_mcmc_imputed_df = data.frame(matrix(ncol = 1, nrow = nrow(chip.table)))
names(tmp_mcmc_imputed_df) = c("chip")
tmp_mcmc_imputed_df$chip = chip.table$V1
tmp_mcmc_imputed_df$experiment = exp.dir

main_mcmc_imputed_df = data.frame()
main_mcmc_imputed_df_t0 = data.frame()
main_mcmc_imputed_df_t1 = data.frame()
main_mcmc_imputed_df_t2 = data.frame()
main_mcmc_imputed_df_t3 = data.frame()

main_mcmc_typed_df = data.frame()
main_mcmc_typed_df_t0 = data.frame()
main_mcmc_typed_df_t1 = data.frame()
main_mcmc_typed_df_t2 = data.frame()
main_mcmc_typed_df_t3 = data.frame()

for (exp in 1:length(exp.var)) {
  message(paste0(exp.var[exp]))
  imp.out.file = paste0(outward_dir, "MCMC/ConcordanceTables_", exp.var[exp], "_MCC_imputed_genotype.csv")
  imp.out.file.t0 = paste0(outward_dir, "MCMC/ConcordanceTables_", exp.var[exp], "_MCC_imputed_genotype_t0.csv")
  imp.out.file.t1 = paste0(outward_dir, "MCMC/ConcordanceTables_", exp.var[exp], "_MCC_imputed_genotype_t1.csv")
  imp.out.file.t2 = paste0(outward_dir, "MCMC/ConcordanceTables_", exp.var[exp], "_MCC_imputed_genotype_t2.csv")
  imp.out.file.t3 = paste0(outward_dir, "MCMC/ConcordanceTables_", exp.var[exp], "_MCC_imputed_genotype_t3.csv")
  
  typ.out.file = paste0(outward_dir, "MCMC/ConcordanceTables_", exp.var[exp], "_MCC_typed_genotype.csv")
  typ.out.file.t0 = paste0(outward_dir, "MCMC/ConcordanceTables_", exp.var[exp], "_MCC_typed_genotype_t0.csv")
  typ.out.file.t1 = paste0(outward_dir, "MCMC/ConcordanceTables_", exp.var[exp], "_MCC_typed_genotype_t1.csv")
  typ.out.file.t2 = paste0(outward_dir, "MCMC/ConcordanceTables_", exp.var[exp], "_MCC_typed_genotype_t2.csv")
  typ.out.file.t3 = paste0(outward_dir, "MCMC/ConcordanceTables_", exp.var[exp], "_MCC_typed_genotype_t3.csv")
  
  tmp_mcmc_imputed_df$sub_experiment = exp.var[exp]
  #tmp_mcmc_typed_df$sub_experiment = exp.var[exp]
  for (chip in 1:nrow(chip.table)) {
    tmp.imp.file = paste0(container, chip.table$V1[chip], "/", exp.dir, "/", exp.var[exp], "/", "chrMT_1kg_", chip.table$V1[chip], "_imputed_", exp.var[exp], "_imputed_MCC.csv")
    tmp.typ.file = paste0(container, chip.table$V1[chip], "/", exp.dir, "/", exp.var[exp], "/", "chrMT_1kg_", chip.table$V1[chip], "_imputed_", exp.var[exp], "_typed_MCC.csv")
    
    platform.snps = paste0(container, chip.table$V1[chip], "/", chip.table$V1[chip], "_MT_snps.txt")
    print(platform.snps)
    
    if (file.exists(tmp.imp.file) == T) {
      tmp_mcmc_imputed_df$imputed[chip] = T
      chip.table$imputed[chip] = T
      snps = read.table(platform.snps, header = F, sep = "\t")
      tmp_mcmc_imputed_df$n.snps[chip] = nrow(snps)
    } else {
      tmp_mcmc_imputed_df$imputed[chip] = F
      chip.table$imputed[chip] = F
      tmp_mcmc_imputed_df$n.snps[chip] = 0
    }
    
    #if (file.exists(platform.snps) == T) {
    #  snps = read.table(platform.snps, header = F, sep = "\t")
    #  tmp_mcmc_imputed_df$n.snps[chip] = nrow(snps)
    #} else {
    #  tmp_mcmc_imputed_df$n.snps[chip] = 0
    #}
    
    tmp_mcmc_imputed_df$allele_freq[chip] = NA
    tmp_mcmc_imputed_df$mcc[chip] = NA
    tmp_mcmc_imputed_df$concordance[chip] = NA
    tmp_mcmc_imputed_df$info[chip] = NA
    tmp_mcmc_imputed_df$certainty[chip] = NA
    tmp_mcmc_imputed_df$info_type0[chip] = NA
    tmp_mcmc_imputed_df$concord_type0[chip] = NA
    tmp_mcmc_imputed_df$r2_type0[chip] = NA
    tmp_mcmc_imputed_df$info_comb[chip] = NA
  }
  
  tmp_mcmc_typed_df = tmp_mcmc_imputed_df
  
  tmp_mcmc_imputed_df_t0 = tmp_mcmc_imputed_df
  tmp_mcmc_imputed_df_t1 = tmp_mcmc_imputed_df
  tmp_mcmc_imputed_df_t2 = tmp_mcmc_imputed_df
  tmp_mcmc_imputed_df_t3 = tmp_mcmc_imputed_df
  
  tmp_mcmc_typed_df_t0 = tmp_mcmc_imputed_df_t0
  tmp_mcmc_typed_df_t1 = tmp_mcmc_imputed_df_t1
  tmp_mcmc_typed_df_t2 = tmp_mcmc_imputed_df_t2
  tmp_mcmc_typed_df_t3 = tmp_mcmc_imputed_df_t3
  
  total = nrow(chip.table)
  total_imputed = nrow(subset(chip.table, chip.table$imputed == T))
  
  for (chip in 1:nrow(chip.table)) {
    message(paste0("WORKING ON CHIP FOR ", exp.var[exp], ":\t", chip, " / ", total))
    
    tmp.imp.file = paste0(container, chip.table$V1[chip], "/", exp.dir, "/", exp.var[exp], "/", "chrMT_1kg_", chip.table$V1[chip], "_imputed_", exp.var[exp], "_imputed_MCC.csv")
    tmp.typ.file = paste0(container, chip.table$V1[chip], "/", exp.dir, "/", exp.var[exp], "/", "chrMT_1kg_", chip.table$V1[chip], "_imputed_", exp.var[exp], "_typed_MCC.csv")
    # IMPUTED FILE
    if (file.exists(tmp.imp.file) == T) {
      tmp.imp.table = read.csv(tmp.imp.file, header = T)
      
      # THE FULL THING
      tmp_mcmc_imputed_df$allele_freq[chip]       = mean(tmp.imp.table$af, na.rm = T)
      tmp_mcmc_imputed_df$mcc[chip]               = mean(tmp.imp.table$mcc, na.rm = T)
      tmp_mcmc_imputed_df$concordance[chip]       = mean(tmp.imp.table$concodance, na.rm = T)
      tmp_mcmc_imputed_df$info[chip]              = mean(tmp.imp.table$info, na.rm = T)
      tmp_mcmc_imputed_df$certainty[chip]         = mean(tmp.imp.table$certainty, na.rm = T)
      tmp_mcmc_imputed_df$info_type0[chip]        = mean(tmp.imp.table$info_type0, na.rm = T)
      tmp_mcmc_imputed_df$concord_type0[chip]     = mean(tmp.imp.table$concord_type0, na.rm = T)
      tmp_mcmc_imputed_df$r2_type0[chip]          = mean(tmp.imp.table$r2_type0, na.rm = T)
      tmp_mcmc_imputed_df$info_comb[chip]         = mean(tmp.imp.table$info_comb, na.rm = T)
      
      # TYPE 0 ONLY!
      tmp.imp.table_t0 = subset(tmp.imp.table, tmp.imp.table$type == 0)
      
      tmp_mcmc_imputed_df_t0$allele_freq[chip]       = mean(tmp.imp.table_t0$af, na.rm = T)
      tmp_mcmc_imputed_df_t0$mcc[chip]               = mean(tmp.imp.table_t0$mcc, na.rm = T)
      tmp_mcmc_imputed_df_t0$concordance[chip]       = mean(tmp.imp.table_t0$concodance, na.rm = T)
      tmp_mcmc_imputed_df_t0$info[chip]              = mean(tmp.imp.table_t0$info, na.rm = T)
      tmp_mcmc_imputed_df_t0$certainty[chip]         = mean(tmp.imp.table_t0$certainty, na.rm = T)
      tmp_mcmc_imputed_df_t0$info_type0[chip]        = mean(tmp.imp.table_t0$info_type0, na.rm = T)
      tmp_mcmc_imputed_df_t0$concord_type0[chip]     = mean(tmp.imp.table_t0$concord_type0, na.rm = T)
      tmp_mcmc_imputed_df_t0$r2_type0[chip]          = mean(tmp.imp.table_t0$r2_type0, na.rm = T)
      tmp_mcmc_imputed_df_t0$info_comb[chip]         = mean(tmp.imp.table_t0$info_comb, na.rm = T)
      # TYPE 1 ONLY!
      tmp.imp.table_t1 = subset(tmp.imp.table, tmp.imp.table$type == 1)
      
      tmp_mcmc_imputed_df_t1$allele_freq[chip]       = mean(tmp.imp.table_t1$af, na.rm = T)
      tmp_mcmc_imputed_df_t1$mcc[chip]               = mean(tmp.imp.table_t1$mcc, na.rm = T)
      tmp_mcmc_imputed_df_t1$concordance[chip]       = mean(tmp.imp.table_t1$concodance, na.rm = T)
      tmp_mcmc_imputed_df_t1$info[chip]              = mean(tmp.imp.table_t1$info, na.rm = T)
      tmp_mcmc_imputed_df_t1$certainty[chip]         = mean(tmp.imp.table_t1$certainty, na.rm = T)
      tmp_mcmc_imputed_df_t1$info_type0[chip]        = mean(tmp.imp.table_t1$info_type0, na.rm = T)
      tmp_mcmc_imputed_df_t1$concord_type0[chip]     = mean(tmp.imp.table_t1$concord_type0, na.rm = T)
      tmp_mcmc_imputed_df_t1$r2_type0[chip]          = mean(tmp.imp.table_t1$r2_type0, na.rm = T)
      tmp_mcmc_imputed_df_t1$info_comb[chip]         = mean(tmp.imp.table_t1$info_comb, na.rm = T)
      # TYPE 2 ONLY!
      tmp.imp.table_t2 = subset(tmp.imp.table, tmp.imp.table$type == 2)
      
      tmp_mcmc_imputed_df_t2$allele_freq[chip]       = mean(tmp.imp.table_t2$af, na.rm = T)
      tmp_mcmc_imputed_df_t2$mcc[chip]               = mean(tmp.imp.table_t2$mcc, na.rm = T)
      tmp_mcmc_imputed_df_t2$concordance[chip]       = mean(tmp.imp.table_t2$concodance, na.rm = T)
      tmp_mcmc_imputed_df_t2$info[chip]              = mean(tmp.imp.table_t2$info, na.rm = T)
      tmp_mcmc_imputed_df_t2$certainty[chip]         = mean(tmp.imp.table_t2$certainty, na.rm = T)
      tmp_mcmc_imputed_df_t2$info_type0[chip]        = mean(tmp.imp.table_t2$info_type0, na.rm = T)
      tmp_mcmc_imputed_df_t2$concord_type0[chip]     = mean(tmp.imp.table_t2$concord_type0, na.rm = T)
      tmp_mcmc_imputed_df_t2$r2_type0[chip]          = mean(tmp.imp.table_t2$r2_type0, na.rm = T)
      tmp_mcmc_imputed_df_t2$info_comb[chip]         = mean(tmp.imp.table_t2$info_comb, na.rm = T)
      # TYPE 3 ONLY!
      tmp.imp.table_t3 = subset(tmp.imp.table, tmp.imp.table$type == 3)
      
      tmp_mcmc_imputed_df_t3$allele_freq[chip]       = mean(tmp.imp.table_t3$af, na.rm = T)
      tmp_mcmc_imputed_df_t3$mcc[chip]               = mean(tmp.imp.table_t3$mcc, na.rm = T)
      tmp_mcmc_imputed_df_t3$concordance[chip]       = mean(tmp.imp.table_t3$concodance, na.rm = T)
      tmp_mcmc_imputed_df_t3$info[chip]              = mean(tmp.imp.table_t3$info, na.rm = T)
      tmp_mcmc_imputed_df_t3$certainty[chip]         = mean(tmp.imp.table_t3$certainty, na.rm = T)
      tmp_mcmc_imputed_df_t3$info_type0[chip]        = mean(tmp.imp.table_t3$info_type0, na.rm = T)
      tmp_mcmc_imputed_df_t3$concord_type0[chip]     = mean(tmp.imp.table_t3$concord_type0, na.rm = T)
      tmp_mcmc_imputed_df_t3$r2_type0[chip]          = mean(tmp.imp.table_t3$r2_type0, na.rm = T)
      tmp_mcmc_imputed_df_t3$info_comb[chip]         = mean(tmp.imp.table_t3$info_comb, na.rm = T)
    }
    
    # TYPED FILE
    if (file.exists(tmp.typ.file) == T) {
      tmp.typ.table = read.csv(tmp.typ.file, header = T)
      
      # THE FULL THING
      tmp_mcmc_typed_df$allele_freq[chip]       = mean(tmp.typ.table$af, na.rm = T)
      tmp_mcmc_typed_df$mcc[chip]               = mean(tmp.typ.table$mcc, na.rm = T)
      tmp_mcmc_typed_df$concordance[chip]       = mean(tmp.typ.table$concodance, na.rm = T)
      tmp_mcmc_typed_df$info[chip]              = mean(tmp.typ.table$info, na.rm = T)
      tmp_mcmc_typed_df$certainty[chip]         = mean(tmp.typ.table$certainty, na.rm = T)
      tmp_mcmc_typed_df$info_type0[chip]        = mean(tmp.typ.table$info_type0, na.rm = T)
      tmp_mcmc_typed_df$concord_type0[chip]     = mean(tmp.typ.table$concord_type0, na.rm = T)
      tmp_mcmc_typed_df$r2_type0[chip]          = mean(tmp.typ.table$r2_type0, na.rm = T)
      tmp_mcmc_typed_df$info_comb[chip]         = mean(tmp.typ.table$info_comb, na.rm = T)
      
      # TYPE 0 ONLY!
      tmp.typ.table_t0 = subset(tmp.typ.table, tmp.typ.table$type == 0)
      
      tmp_mcmc_typed_df_t0$allele_freq[chip]       = mean(tmp.typ.table_t0$af, na.rm = T)
      tmp_mcmc_typed_df_t0$mcc[chip]               = mean(tmp.typ.table_t0$mcc, na.rm = T)
      tmp_mcmc_typed_df_t0$concordance[chip]       = mean(tmp.typ.table_t0$concodance, na.rm = T)
      tmp_mcmc_typed_df_t0$info[chip]              = mean(tmp.typ.table_t0$info, na.rm = T)
      tmp_mcmc_typed_df_t0$certainty[chip]         = mean(tmp.typ.table_t0$certainty, na.rm = T)
      tmp_mcmc_typed_df_t0$info_type0[chip]        = mean(tmp.typ.table_t0$info_type0, na.rm = T)
      tmp_mcmc_typed_df_t0$concord_type0[chip]     = mean(tmp.typ.table_t0$concord_type0, na.rm = T)
      tmp_mcmc_typed_df_t0$r2_type0[chip]          = mean(tmp.typ.table_t0$r2_type0, na.rm = T)
      tmp_mcmc_typed_df_t0$info_comb[chip]         = mean(tmp.typ.table_t0$info_comb, na.rm = T)
      # TYPE 1 ONLY!
      tmp.typ.table_t1 = subset(tmp.typ.table, tmp.typ.table$type == 1)
      
      tmp_mcmc_typed_df_t1$allele_freq[chip]       = mean(tmp.typ.table_t1$af, na.rm = T)
      tmp_mcmc_typed_df_t1$mcc[chip]               = mean(tmp.typ.table_t1$mcc, na.rm = T)
      tmp_mcmc_typed_df_t1$concordance[chip]       = mean(tmp.typ.table_t1$concodance, na.rm = T)
      tmp_mcmc_typed_df_t1$info[chip]              = mean(tmp.typ.table_t1$info, na.rm = T)
      tmp_mcmc_typed_df_t1$certainty[chip]         = mean(tmp.typ.table_t1$certainty, na.rm = T)
      tmp_mcmc_typed_df_t1$info_type0[chip]        = mean(tmp.typ.table_t1$info_type0, na.rm = T)
      tmp_mcmc_typed_df_t1$concord_type0[chip]     = mean(tmp.typ.table_t1$concord_type0, na.rm = T)
      tmp_mcmc_typed_df_t1$r2_type0[chip]          = mean(tmp.typ.table_t1$r2_type0, na.rm = T)
      tmp_mcmc_typed_df_t1$info_comb[chip]         = mean(tmp.typ.table_t1$info_comb, na.rm = T)
      # TYPE 2 ONLY!
      tmp.typ.table_t2 = subset(tmp.typ.table, tmp.typ.table$type == 2)
      
      tmp_mcmc_typed_df_t2$allele_freq[chip]       = mean(tmp.typ.table_t2$af, na.rm = T)
      tmp_mcmc_typed_df_t2$mcc[chip]               = mean(tmp.typ.table_t2$mcc, na.rm = T)
      tmp_mcmc_typed_df_t2$concordance[chip]       = mean(tmp.typ.table_t2$concodance, na.rm = T)
      tmp_mcmc_typed_df_t2$info[chip]              = mean(tmp.typ.table_t2$info, na.rm = T)
      tmp_mcmc_typed_df_t2$certainty[chip]         = mean(tmp.typ.table_t2$certainty, na.rm = T)
      tmp_mcmc_typed_df_t2$info_type0[chip]        = mean(tmp.typ.table_t2$info_type0, na.rm = T)
      tmp_mcmc_typed_df_t2$concord_type0[chip]     = mean(tmp.typ.table_t2$concord_type0, na.rm = T)
      tmp_mcmc_typed_df_t2$r2_type0[chip]          = mean(tmp.typ.table_t2$r2_type0, na.rm = T)
      tmp_mcmc_typed_df_t2$info_comb[chip]         = mean(tmp.typ.table_t2$info_comb, na.rm = T)
      # TYPE 3 ONLY!
      tmp.typ.table_t3 = subset(tmp.typ.table, tmp.typ.table$type == 3)
      
      tmp_mcmc_typed_df_t3$allele_freq[chip]       = mean(tmp.typ.table_t3$af, na.rm = T)
      tmp_mcmc_typed_df_t3$mcc[chip]               = mean(tmp.typ.table_t3$mcc, na.rm = T)
      tmp_mcmc_typed_df_t3$concordance[chip]       = mean(tmp.typ.table_t3$concodance, na.rm = T)
      tmp_mcmc_typed_df_t3$info[chip]              = mean(tmp.typ.table_t3$info, na.rm = T)
      tmp_mcmc_typed_df_t3$certainty[chip]         = mean(tmp.typ.table_t3$certainty, na.rm = T)
      tmp_mcmc_typed_df_t3$info_type0[chip]        = mean(tmp.typ.table_t3$info_type0, na.rm = T)
      tmp_mcmc_typed_df_t3$concord_type0[chip]     = mean(tmp.typ.table_t3$concord_type0, na.rm = T)
      tmp_mcmc_typed_df_t3$r2_type0[chip]          = mean(tmp.typ.table_t3$r2_type0, na.rm = T)
      tmp_mcmc_typed_df_t3$info_comb[chip]         = mean(tmp.typ.table_t3$info_comb, na.rm = T)
    }
  }
  write.csv(tmp_mcmc_imputed_df, imp.out.file, row.names = F, quote = F)
  #message(paste0("WROTE ", imp.out.file, " TO DISK"))
  write.csv(tmp_mcmc_imputed_df_t0, imp.out.file.t0, row.names = F, quote = F)
  #message(paste0("WROTE ", imp.out.file.t0, " TO DISK"))
  write.csv(tmp_mcmc_imputed_df_t1, imp.out.file.t1, row.names = F, quote = F)
  #message(paste0("WROTE ", imp.out.file.t1, " TO DISK"))
  write.csv(tmp_mcmc_imputed_df_t2, imp.out.file.t2, row.names = F, quote = F)
  #message(paste0("WROTE ", imp.out.file.t2, " TO DISK"))
  write.csv(tmp_mcmc_imputed_df_t3, imp.out.file.t3, row.names = F, quote = F)
  #message(paste0("WROTE ", imp.out.file.t3, " TO DISK"))
  #
  write.csv(tmp_mcmc_typed_df, typ.out.file, row.names = F, quote = F)
  #message(paste0("WROTE ", typ.out.file, " TO DISK"))
  write.csv(tmp_mcmc_typed_df_t0, typ.out.file.t0, row.names = F, quote = F)
  #message(paste0("WROTE ", typ.out.file.t0, " TO DISK"))
  write.csv(tmp_mcmc_typed_df_t1, typ.out.file.t1, row.names = F, quote = F)
  #message(paste0("WROTE ", typ.out.file.t1, " TO DISK"))
  write.csv(tmp_mcmc_typed_df_t2, typ.out.file.t2, row.names = F, quote = F)
  #message(paste0("WROTE ", typ.out.file.t2, " TO DISK"))
  write.csv(tmp_mcmc_typed_df_t3, typ.out.file.t3, row.names = F, quote = F)
  #message(paste0("WROTE ", typ.out.file.t3, " TO DISK"))
  
  main_mcmc_imputed_df = rbind(main_mcmc_imputed_df, tmp_mcmc_imputed_df)
  main_mcmc_imputed_df_t0 = rbind(main_mcmc_imputed_df_t0, tmp_mcmc_imputed_df_t0)
  main_mcmc_imputed_df_t1 = rbind(main_mcmc_imputed_df_t1, tmp_mcmc_imputed_df_t1)
  main_mcmc_imputed_df_t2 = rbind(main_mcmc_imputed_df_t2, tmp_mcmc_imputed_df_t2)
  main_mcmc_imputed_df_t3 = rbind(main_mcmc_imputed_df_t3, tmp_mcmc_imputed_df_t3)
  
  main_mcmc_typed_df = rbind(main_mcmc_typed_df, tmp_mcmc_typed_df)
  main_mcmc_typed_df_t0 = rbind(main_mcmc_typed_df_t0, tmp_mcmc_typed_df_t0)
  main_mcmc_typed_df_t1 = rbind(main_mcmc_typed_df_t1, tmp_mcmc_typed_df_t1)
  main_mcmc_typed_df_t2 = rbind(main_mcmc_typed_df_t2, tmp_mcmc_typed_df_t2)
  main_mcmc_typed_df_t3 = rbind(main_mcmc_typed_df_t3, tmp_mcmc_typed_df_t3)
}
write.csv(main_mcmc_imputed_df, paste0(outward_dir, "MCMC/ConcordanceTables_", exp.dir,"_MCC_imputed_genotype.csv"), row.names = F, quote = F)
write.csv(main_mcmc_imputed_df_t0, paste0(outward_dir, "MCMC/ConcordanceTables_", exp.dir,"_MCC_imputed_genotype_t0.csv"), row.names = F, quote = F)
write.csv(main_mcmc_imputed_df_t1, paste0(outward_dir, "MCMC/ConcordanceTables_", exp.dir,"_MCC_imputed_genotype_t1.csv"), row.names = F, quote = F)
write.csv(main_mcmc_imputed_df_t2, paste0(outward_dir, "MCMC/ConcordanceTables_", exp.dir,"_MCC_imputed_genotype_t2.csv"), row.names = F, quote = F)
write.csv(main_mcmc_imputed_df_t3, paste0(outward_dir, "MCMC/ConcordanceTables_", exp.dir,"_MCC_imputed_genotype_t3.csv"), row.names = F, quote = F)
write.csv(main_mcmc_typed_df, paste0(outward_dir, "MCMC/ConcordanceTables_", exp.dir,"_MCC_typed_genotype.csv"), row.names = F, quote = F)
write.csv(main_mcmc_typed_df_t0, paste0(outward_dir, "MCMC/ConcordanceTables_", exp.dir,"_MCC_typed_genotype_t0.csv"), row.names = F, quote = F)
write.csv(main_mcmc_typed_df_t1, paste0(outward_dir, "MCMC/ConcordanceTables_", exp.dir,"_MCC_typed_genotype_t1.csv"), row.names = F, quote = F)
write.csv(main_mcmc_typed_df_t2, paste0(outward_dir, "MCMC/ConcordanceTables_", exp.dir,"_MCC_typed_genotype_t2.csv"), row.names = F, quote = F)
write.csv(main_mcmc_typed_df_t3, paste0(outward_dir, "MCMC/ConcordanceTables_", exp.dir,"_MCC_typed_genotype_t3.csv"), row.names = F, quote = F)
message(paste0("WROTE ", paste0(outward_dir, "MCMC/ConcordanceTables_", exp.dir,"_MCC_imputed_genotype.csv"), " TO DISK"))
message(paste0("WROTE ", paste0(outward_dir, "MCMC/ConcordanceTables_", exp.dir,"_MCC_imputed_genotype_t0.csv"), " TO DISK"))
message(paste0("WROTE ", paste0(outward_dir, "MCMC/ConcordanceTables_", exp.dir,"_MCC_imputed_genotype_t1.csv"), " TO DISK"))
message(paste0("WROTE ", paste0(outward_dir, "MCMC/ConcordanceTables_", exp.dir,"_MCC_imputed_genotype_t2.csv"), " TO DISK"))
message(paste0("WROTE ", paste0(outward_dir, "MCMC/ConcordanceTables_", exp.dir,"_MCC_imputed_genotype_t3.csv"), " TO DISK"))
message(paste0("WROTE ", paste0(outward_dir, "MCMC/ConcordanceTables_", exp.dir,"_MCC_typed_genotype.csv"), " TO DISK"))
message(paste0("WROTE ", paste0(outward_dir, "MCMC/ConcordanceTables_", exp.dir,"_MCC_typed_genotype_t0.csv"), " TO DISK"))
message(paste0("WROTE ", paste0(outward_dir, "MCMC/ConcordanceTables_", exp.dir,"_MCC_typed_genotype_t1.csv"), " TO DISK"))
message(paste0("WROTE ", paste0(outward_dir, "MCMC/ConcordanceTables_", exp.dir,"_MCC_typed_genotype_t2.csv"), " TO DISK"))
message(paste0("WROTE ", paste0(outward_dir, "MCMC/ConcordanceTables_", exp.dir,"_MCC_typed_genotype_t3.csv"), " TO DISK"))
## KHAP 

exp.dir = "kHAP_Experiments"
exp.var = c("kHAP100", "kHAP250", "kHAP500", "kHAP1000", "kHAP2500", "kHAP5000", "kHAP10000", "kHAP20000", "kHAP30000")
tmp_khap_imputed_df = data.frame(matrix(ncol = 1, nrow = nrow(chip.table)))
names(tmp_khap_imputed_df) = c("chip")
tmp_khap_imputed_df$chip = chip.table$V1
tmp_khap_imputed_df$experiment = exp.dir

main_khap_imputed_df = data.frame()
main_khap_imputed_df_t0 = data.frame()
main_khap_imputed_df_t1 = data.frame()
main_khap_imputed_df_t2 = data.frame()
main_khap_imputed_df_t3 = data.frame()

main_khap_typed_df = data.frame()
main_khap_typed_df_t0 = data.frame()
main_khap_typed_df_t1 = data.frame()
main_khap_typed_df_t2 = data.frame()
main_khap_typed_df_t3 = data.frame()

for (exp in 1:length(exp.var)) {
  message(paste0(exp.var[exp]))
  imp.out.file = paste0(outward_dir, "KHAP/ConcordanceTables_", exp.var[exp], "_MCC_imputed_genotype.csv")
  imp.out.file.t0 = paste0(outward_dir, "KHAP/ConcordanceTables_", exp.var[exp], "_MCC_imputed_genotype_t0.csv")
  imp.out.file.t1 = paste0(outward_dir, "KHAP/ConcordanceTables_", exp.var[exp], "_MCC_imputed_genotype_t1.csv")
  imp.out.file.t2 = paste0(outward_dir, "KHAP/ConcordanceTables_", exp.var[exp], "_MCC_imputed_genotype_t2.csv")
  imp.out.file.t3 = paste0(outward_dir, "KHAP/ConcordanceTables_", exp.var[exp], "_MCC_imputed_genotype_t3.csv")
  
  typ.out.file = paste0(outward_dir, "KHAP/ConcordanceTables_", exp.var[exp], "_MCC_typed_genotype.csv")
  typ.out.file.t0 = paste0(outward_dir, "KHAP/ConcordanceTables_", exp.var[exp], "_MCC_typed_genotype_t0.csv")
  typ.out.file.t1 = paste0(outward_dir, "KHAP/ConcordanceTables_", exp.var[exp], "_MCC_typed_genotype_t1.csv")
  typ.out.file.t2 = paste0(outward_dir, "KHAP/ConcordanceTables_", exp.var[exp], "_MCC_typed_genotype_t2.csv")
  typ.out.file.t3 = paste0(outward_dir, "KHAP/ConcordanceTables_", exp.var[exp], "_MCC_typed_genotype_t3.csv")
  
  tmp_khap_imputed_df$sub_experiment = exp.var[exp]
  #tmp_khap_typed_df$sub_experiment = exp.var[exp]
  for (chip in 1:nrow(chip.table)) {
    tmp.imp.file = paste0(container, chip.table$V1[chip], "/", exp.dir, "/", exp.var[exp], "/", "chrMT_1kg_", chip.table$V1[chip], "_imputed_", exp.var[exp], "_imputed_MCC.csv")
    tmp.typ.file = paste0(container, chip.table$V1[chip], "/", exp.dir, "/", exp.var[exp], "/", "chrMT_1kg_", chip.table$V1[chip], "_imputed_", exp.var[exp], "_typed_MCC.csv")
    
    platform.snps = paste0(container, chip.table$V1[chip], "/", chip.table$V1[chip], "_MT_snps.txt")
    
    if (file.exists(tmp.imp.file) == T) {
      tmp_khap_imputed_df$imputed[chip] = T
      chip.table$imputed[chip] = T
    } else {
      tmp_khap_imputed_df$imputed[chip] = F
      chip.table$imputed[chip] = F
    }
    
    if (file.exists(platform.snps) == T) {
      snps = read.table(platform.snps, header = F, sep = "\t")
      tmp_khap_imputed_df$n.snps[chip] = nrow(snps)
    } else {
      tmp_khap_imputed_df$n.snps[chip] = 0
    }
    
    tmp_khap_imputed_df$n.snps[chip] = nrow(snps)
    tmp_khap_imputed_df$allele_freq[chip] = NA
    tmp_khap_imputed_df$mcc[chip] = NA
    tmp_khap_imputed_df$concordance[chip] = NA
    tmp_khap_imputed_df$info[chip] = NA
    tmp_khap_imputed_df$certainty[chip] = NA
    tmp_khap_imputed_df$info_type0[chip] = NA
    tmp_khap_imputed_df$concord_type0[chip] = NA
    tmp_khap_imputed_df$r2_type0[chip] = NA
    tmp_khap_imputed_df$info_comb[chip] = NA
  }
  
  tmp_khap_typed_df = tmp_khap_imputed_df
  
  tmp_khap_imputed_df_t0 = tmp_khap_imputed_df
  tmp_khap_imputed_df_t1 = tmp_khap_imputed_df
  tmp_khap_imputed_df_t2 = tmp_khap_imputed_df
  tmp_khap_imputed_df_t3 = tmp_khap_imputed_df
  
  tmp_khap_typed_df_t0 = tmp_khap_imputed_df_t0
  tmp_khap_typed_df_t1 = tmp_khap_imputed_df_t1
  tmp_khap_typed_df_t2 = tmp_khap_imputed_df_t2
  tmp_khap_typed_df_t3 = tmp_khap_imputed_df_t3
  
  total = nrow(chip.table)
  total_imputed = nrow(subset(chip.table, chip.table$imputed == T))
  
  for (chip in 1:nrow(chip.table)) {
    message(paste0("WORKING ON CHIP FOR ", exp.var[exp], ":\t", chip, " / ", total))
    
    tmp.imp.file = paste0(container, chip.table$V1[chip], "/", exp.dir, "/", exp.var[exp], "/", "chrMT_1kg_", chip.table$V1[chip], "_imputed_", exp.var[exp], "_imputed_MCC.csv")
    tmp.typ.file = paste0(container, chip.table$V1[chip], "/", exp.dir, "/", exp.var[exp], "/", "chrMT_1kg_", chip.table$V1[chip], "_imputed_", exp.var[exp], "_typed_MCC.csv")
    # IMPUTED FILE
    if (file.exists(tmp.imp.file) == T) {
      tmp.imp.table = read.csv(tmp.imp.file, header = T)
      
      # THE FULL THING
      tmp_khap_imputed_df$allele_freq[chip]       = mean(tmp.imp.table$af, na.rm = T)
      tmp_khap_imputed_df$mcc[chip]               = mean(tmp.imp.table$mcc, na.rm = T)
      tmp_khap_imputed_df$concordance[chip]       = mean(tmp.imp.table$concodance, na.rm = T)
      tmp_khap_imputed_df$info[chip]              = mean(tmp.imp.table$info, na.rm = T)
      tmp_khap_imputed_df$certainty[chip]         = mean(tmp.imp.table$certainty, na.rm = T)
      tmp_khap_imputed_df$info_type0[chip]        = mean(tmp.imp.table$info_type0, na.rm = T)
      tmp_khap_imputed_df$concord_type0[chip]     = mean(tmp.imp.table$concord_type0, na.rm = T)
      tmp_khap_imputed_df$r2_type0[chip]          = mean(tmp.imp.table$r2_type0, na.rm = T)
      tmp_khap_imputed_df$info_comb[chip]         = mean(tmp.imp.table$info_comb, na.rm = T)
      
      # TYPE 0 ONLY!
      tmp.imp.table_t0 = subset(tmp.imp.table, tmp.imp.table$type == 0)
      
      tmp_khap_imputed_df_t0$allele_freq[chip]       = mean(tmp.imp.table_t0$af, na.rm = T)
      tmp_khap_imputed_df_t0$mcc[chip]               = mean(tmp.imp.table_t0$mcc, na.rm = T)
      tmp_khap_imputed_df_t0$concordance[chip]       = mean(tmp.imp.table_t0$concodance, na.rm = T)
      tmp_khap_imputed_df_t0$info[chip]              = mean(tmp.imp.table_t0$info, na.rm = T)
      tmp_khap_imputed_df_t0$certainty[chip]         = mean(tmp.imp.table_t0$certainty, na.rm = T)
      tmp_khap_imputed_df_t0$info_type0[chip]        = mean(tmp.imp.table_t0$info_type0, na.rm = T)
      tmp_khap_imputed_df_t0$concord_type0[chip]     = mean(tmp.imp.table_t0$concord_type0, na.rm = T)
      tmp_khap_imputed_df_t0$r2_type0[chip]          = mean(tmp.imp.table_t0$r2_type0, na.rm = T)
      tmp_khap_imputed_df_t0$info_comb[chip]         = mean(tmp.imp.table_t0$info_comb, na.rm = T)
      # TYPE 1 ONLY!
      tmp.imp.table_t1 = subset(tmp.imp.table, tmp.imp.table$type == 1)
      
      tmp_khap_imputed_df_t1$allele_freq[chip]       = mean(tmp.imp.table_t1$af, na.rm = T)
      tmp_khap_imputed_df_t1$mcc[chip]               = mean(tmp.imp.table_t1$mcc, na.rm = T)
      tmp_khap_imputed_df_t1$concordance[chip]       = mean(tmp.imp.table_t1$concodance, na.rm = T)
      tmp_khap_imputed_df_t1$info[chip]              = mean(tmp.imp.table_t1$info, na.rm = T)
      tmp_khap_imputed_df_t1$certainty[chip]         = mean(tmp.imp.table_t1$certainty, na.rm = T)
      tmp_khap_imputed_df_t1$info_type0[chip]        = mean(tmp.imp.table_t1$info_type0, na.rm = T)
      tmp_khap_imputed_df_t1$concord_type0[chip]     = mean(tmp.imp.table_t1$concord_type0, na.rm = T)
      tmp_khap_imputed_df_t1$r2_type0[chip]          = mean(tmp.imp.table_t1$r2_type0, na.rm = T)
      tmp_khap_imputed_df_t1$info_comb[chip]         = mean(tmp.imp.table_t1$info_comb, na.rm = T)
      # TYPE 2 ONLY!
      tmp.imp.table_t2 = subset(tmp.imp.table, tmp.imp.table$type == 2)
      
      tmp_khap_imputed_df_t2$allele_freq[chip]       = mean(tmp.imp.table_t2$af, na.rm = T)
      tmp_khap_imputed_df_t2$mcc[chip]               = mean(tmp.imp.table_t2$mcc, na.rm = T)
      tmp_khap_imputed_df_t2$concordance[chip]       = mean(tmp.imp.table_t2$concodance, na.rm = T)
      tmp_khap_imputed_df_t2$info[chip]              = mean(tmp.imp.table_t2$info, na.rm = T)
      tmp_khap_imputed_df_t2$certainty[chip]         = mean(tmp.imp.table_t2$certainty, na.rm = T)
      tmp_khap_imputed_df_t2$info_type0[chip]        = mean(tmp.imp.table_t2$info_type0, na.rm = T)
      tmp_khap_imputed_df_t2$concord_type0[chip]     = mean(tmp.imp.table_t2$concord_type0, na.rm = T)
      tmp_khap_imputed_df_t2$r2_type0[chip]          = mean(tmp.imp.table_t2$r2_type0, na.rm = T)
      tmp_khap_imputed_df_t2$info_comb[chip]         = mean(tmp.imp.table_t2$info_comb, na.rm = T)
      # TYPE 3 ONLY!
      tmp.imp.table_t3 = subset(tmp.imp.table, tmp.imp.table$type == 3)
      
      tmp_khap_imputed_df_t3$allele_freq[chip]       = mean(tmp.imp.table_t3$af, na.rm = T)
      tmp_khap_imputed_df_t3$mcc[chip]               = mean(tmp.imp.table_t3$mcc, na.rm = T)
      tmp_khap_imputed_df_t3$concordance[chip]       = mean(tmp.imp.table_t3$concodance, na.rm = T)
      tmp_khap_imputed_df_t3$info[chip]              = mean(tmp.imp.table_t3$info, na.rm = T)
      tmp_khap_imputed_df_t3$certainty[chip]         = mean(tmp.imp.table_t3$certainty, na.rm = T)
      tmp_khap_imputed_df_t3$info_type0[chip]        = mean(tmp.imp.table_t3$info_type0, na.rm = T)
      tmp_khap_imputed_df_t3$concord_type0[chip]     = mean(tmp.imp.table_t3$concord_type0, na.rm = T)
      tmp_khap_imputed_df_t3$r2_type0[chip]          = mean(tmp.imp.table_t3$r2_type0, na.rm = T)
      tmp_khap_imputed_df_t3$info_comb[chip]         = mean(tmp.imp.table_t3$info_comb, na.rm = T)
    }
    
    # GENOTYPED FILE
    if (file.exists(tmp.typ.file) == T) {
      tmp.typ.table = read.csv(tmp.typ.file, header = T)
      
      # THE FULL THING
      tmp_khap_typed_df$allele_freq[chip]       = mean(tmp.typ.table$af, na.rm = T)
      tmp_khap_typed_df$mcc[chip]               = mean(tmp.typ.table$mcc, na.rm = T)
      tmp_khap_typed_df$concordance[chip]       = mean(tmp.typ.table$concodance, na.rm = T)
      tmp_khap_typed_df$info[chip]              = mean(tmp.typ.table$info, na.rm = T)
      tmp_khap_typed_df$certainty[chip]         = mean(tmp.typ.table$certainty, na.rm = T)
      tmp_khap_typed_df$info_type0[chip]        = mean(tmp.typ.table$info_type0, na.rm = T)
      tmp_khap_typed_df$concord_type0[chip]     = mean(tmp.typ.table$concord_type0, na.rm = T)
      tmp_khap_typed_df$r2_type0[chip]          = mean(tmp.typ.table$r2_type0, na.rm = T)
      tmp_khap_typed_df$info_comb[chip]         = mean(tmp.typ.table$info_comb, na.rm = T)
      
      # TYPE 0 ONLY!
      tmp.typ.table_t0 = subset(tmp.typ.table, tmp.typ.table$type == 0)
      
      tmp_khap_typed_df_t0$allele_freq[chip]       = mean(tmp.typ.table_t0$af, na.rm = T)
      tmp_khap_typed_df_t0$mcc[chip]               = mean(tmp.typ.table_t0$mcc, na.rm = T)
      tmp_khap_typed_df_t0$concordance[chip]       = mean(tmp.typ.table_t0$concodance, na.rm = T)
      tmp_khap_typed_df_t0$info[chip]              = mean(tmp.typ.table_t0$info, na.rm = T)
      tmp_khap_typed_df_t0$certainty[chip]         = mean(tmp.typ.table_t0$certainty, na.rm = T)
      tmp_khap_typed_df_t0$info_type0[chip]        = mean(tmp.typ.table_t0$info_type0, na.rm = T)
      tmp_khap_typed_df_t0$concord_type0[chip]     = mean(tmp.typ.table_t0$concord_type0, na.rm = T)
      tmp_khap_typed_df_t0$r2_type0[chip]          = mean(tmp.typ.table_t0$r2_type0, na.rm = T)
      tmp_khap_typed_df_t0$info_comb[chip]         = mean(tmp.typ.table_t0$info_comb, na.rm = T)
      # TYPE 1 ONLY!
      tmp.typ.table_t1 = subset(tmp.typ.table, tmp.typ.table$type == 1)
      
      tmp_khap_typed_df_t1$allele_freq[chip]       = mean(tmp.typ.table_t1$af, na.rm = T)
      tmp_khap_typed_df_t1$mcc[chip]               = mean(tmp.typ.table_t1$mcc, na.rm = T)
      tmp_khap_typed_df_t1$concordance[chip]       = mean(tmp.typ.table_t1$concodance, na.rm = T)
      tmp_khap_typed_df_t1$info[chip]              = mean(tmp.typ.table_t1$info, na.rm = T)
      tmp_khap_typed_df_t1$certainty[chip]         = mean(tmp.typ.table_t1$certainty, na.rm = T)
      tmp_khap_typed_df_t1$info_type0[chip]        = mean(tmp.typ.table_t1$info_type0, na.rm = T)
      tmp_khap_typed_df_t1$concord_type0[chip]     = mean(tmp.typ.table_t1$concord_type0, na.rm = T)
      tmp_khap_typed_df_t1$r2_type0[chip]          = mean(tmp.typ.table_t1$r2_type0, na.rm = T)
      tmp_khap_typed_df_t1$info_comb[chip]         = mean(tmp.typ.table_t1$info_comb, na.rm = T)
      # TYPE 2 ONLY!
      tmp.typ.table_t2 = subset(tmp.typ.table, tmp.typ.table$type == 2)
      
      tmp_khap_typed_df_t2$allele_freq[chip]       = mean(tmp.typ.table_t2$af, na.rm = T)
      tmp_khap_typed_df_t2$mcc[chip]               = mean(tmp.typ.table_t2$mcc, na.rm = T)
      tmp_khap_typed_df_t2$concordance[chip]       = mean(tmp.typ.table_t2$concodance, na.rm = T)
      tmp_khap_typed_df_t2$info[chip]              = mean(tmp.typ.table_t2$info, na.rm = T)
      tmp_khap_typed_df_t2$certainty[chip]         = mean(tmp.typ.table_t2$certainty, na.rm = T)
      tmp_khap_typed_df_t2$info_type0[chip]        = mean(tmp.typ.table_t2$info_type0, na.rm = T)
      tmp_khap_typed_df_t2$concord_type0[chip]     = mean(tmp.typ.table_t2$concord_type0, na.rm = T)
      tmp_khap_typed_df_t2$r2_type0[chip]          = mean(tmp.typ.table_t2$r2_type0, na.rm = T)
      tmp_khap_typed_df_t2$info_comb[chip]         = mean(tmp.typ.table_t2$info_comb, na.rm = T)
      # TYPE 3 ONLY!
      tmp.typ.table_t3 = subset(tmp.typ.table, tmp.typ.table$type == 3)
      
      tmp_khap_typed_df_t3$allele_freq[chip]       = mean(tmp.typ.table_t3$af, na.rm = T)
      tmp_khap_typed_df_t3$mcc[chip]               = mean(tmp.typ.table_t3$mcc, na.rm = T)
      tmp_khap_typed_df_t3$concordance[chip]       = mean(tmp.typ.table_t3$concodance, na.rm = T)
      tmp_khap_typed_df_t3$info[chip]              = mean(tmp.typ.table_t3$info, na.rm = T)
      tmp_khap_typed_df_t3$certainty[chip]         = mean(tmp.typ.table_t3$certainty, na.rm = T)
      tmp_khap_typed_df_t3$info_type0[chip]        = mean(tmp.typ.table_t3$info_type0, na.rm = T)
      tmp_khap_typed_df_t3$concord_type0[chip]     = mean(tmp.typ.table_t3$concord_type0, na.rm = T)
      tmp_khap_typed_df_t3$r2_type0[chip]          = mean(tmp.typ.table_t3$r2_type0, na.rm = T)
      tmp_khap_typed_df_t3$info_comb[chip]         = mean(tmp.typ.table_t3$info_comb, na.rm = T)
    }
  }
  write.csv(tmp_khap_imputed_df, imp.out.file, row.names = F, quote = F)
  #message(paste0("WROTE ", imp.out.file, " TO DISK"))
  write.csv(tmp_khap_imputed_df_t0, imp.out.file.t0, row.names = F, quote = F)
  #message(paste0("WROTE ", imp.out.file.t0, " TO DISK"))
  write.csv(tmp_khap_imputed_df_t1, imp.out.file.t1, row.names = F, quote = F)
  #message(paste0("WROTE ", imp.out.file.t1, " TO DISK"))
  write.csv(tmp_khap_imputed_df_t2, imp.out.file.t2, row.names = F, quote = F)
  #message(paste0("WROTE ", imp.out.file.t2, " TO DISK"))
  write.csv(tmp_khap_imputed_df_t3, imp.out.file.t3, row.names = F, quote = F)
  #message(paste0("WROTE ", imp.out.file.t3, " TO DISK"))
  #
  write.csv(tmp_khap_typed_df, typ.out.file, row.names = F, quote = F)
  #message(paste0("WROTE ", typ.out.file, " TO DISK"))
  write.csv(tmp_khap_typed_df_t0, typ.out.file.t0, row.names = F, quote = F)
  #message(paste0("WROTE ", typ.out.file.t0, " TO DISK"))
  write.csv(tmp_khap_typed_df_t1, typ.out.file.t1, row.names = F, quote = F)
  #message(paste0("WROTE ", typ.out.file.t1, " TO DISK"))
  write.csv(tmp_khap_typed_df_t2, typ.out.file.t2, row.names = F, quote = F)
  #message(paste0("WROTE ", typ.out.file.t2, " TO DISK"))
  write.csv(tmp_khap_typed_df_t3, typ.out.file.t3, row.names = F, quote = F)
  #message(paste0("WROTE ", typ.out.file.t3, " TO DISK"))
  
  main_khap_imputed_df = rbind(main_khap_imputed_df, tmp_khap_imputed_df)
  main_khap_imputed_df_t0 = rbind(main_khap_imputed_df_t0, tmp_khap_imputed_df_t0)
  main_khap_imputed_df_t1 = rbind(main_khap_imputed_df_t1, tmp_khap_imputed_df_t1)
  main_khap_imputed_df_t2 = rbind(main_khap_imputed_df_t2, tmp_khap_imputed_df_t2)
  main_khap_imputed_df_t3 = rbind(main_khap_imputed_df_t3, tmp_khap_imputed_df_t3)
  
  main_khap_typed_df = rbind(main_khap_typed_df, tmp_khap_typed_df)
  main_khap_typed_df_t0 = rbind(main_khap_typed_df_t0, tmp_khap_typed_df_t0)
  main_khap_typed_df_t1 = rbind(main_khap_typed_df_t1, tmp_khap_typed_df_t1)
  main_khap_typed_df_t2 = rbind(main_khap_typed_df_t2, tmp_khap_typed_df_t2)
  main_khap_typed_df_t3 = rbind(main_khap_typed_df_t3, tmp_khap_typed_df_t3)
}
write.csv(main_khap_imputed_df, paste0(outward_dir, "KHAP/ConcordanceTables_", exp.dir,"_MCC_imputed_genotype.csv"), row.names = F, quote = F)
write.csv(main_khap_imputed_df_t0, paste0(outward_dir, "KHAP/ConcordanceTables_", exp.dir,"_MCC_imputed_genotype_t0.csv"), row.names = F, quote = F)
write.csv(main_khap_imputed_df_t1, paste0(outward_dir, "KHAP/ConcordanceTables_", exp.dir,"_MCC_imputed_genotype_t1.csv"), row.names = F, quote = F)
write.csv(main_khap_imputed_df_t2, paste0(outward_dir, "KHAP/ConcordanceTables_", exp.dir,"_MCC_imputed_genotype_t2.csv"), row.names = F, quote = F)
write.csv(main_khap_imputed_df_t3, paste0(outward_dir, "KHAP/ConcordanceTables_", exp.dir,"_MCC_imputed_genotype_t3.csv"), row.names = F, quote = F)
write.csv(main_khap_typed_df, paste0(outward_dir, "KHAP/ConcordanceTables_", exp.dir,"_MCC_typed_genotype.csv"), row.names = F, quote = F)
write.csv(main_khap_typed_df_t0, paste0(outward_dir, "KHAP/ConcordanceTables_", exp.dir,"_MCC_typed_genotype_t0.csv"), row.names = F, quote = F)
write.csv(main_khap_typed_df_t1, paste0(outward_dir, "KHAP/ConcordanceTables_", exp.dir,"_MCC_typed_genotype_t1.csv"), row.names = F, quote = F)
write.csv(main_khap_typed_df_t2, paste0(outward_dir, "KHAP/ConcordanceTables_", exp.dir,"_MCC_typed_genotype_t2.csv"), row.names = F, quote = F)
write.csv(main_khap_typed_df_t3, paste0(outward_dir, "KHAP/ConcordanceTables_", exp.dir,"_MCC_typed_genotype_t3.csv"), row.names = F, quote = F)
message(paste0("WROTE ", paste0(outward_dir, "KHAP/ConcordanceTables_", exp.dir,"_MCC_imputed_genotype.csv"), " TO DISK"))
message(paste0("WROTE ", paste0(outward_dir, "KHAP/ConcordanceTables_", exp.dir,"_MCC_imputed_genotype_t0.csv"), " TO DISK"))
message(paste0("WROTE ", paste0(outward_dir, "KHAP/ConcordanceTables_", exp.dir,"_MCC_imputed_genotype_t1.csv"), " TO DISK"))
message(paste0("WROTE ", paste0(outward_dir, "KHAP/ConcordanceTables_", exp.dir,"_MCC_imputed_genotype_t2.csv"), " TO DISK"))
message(paste0("WROTE ", paste0(outward_dir, "KHAP/ConcordanceTables_", exp.dir,"_MCC_imputed_genotype_t3.csv"), " TO DISK"))
message(paste0("WROTE ", paste0(outward_dir, "KHAP/ConcordanceTables_", exp.dir,"_MCC_typed_genotype.csv"), " TO DISK"))
message(paste0("WROTE ", paste0(outward_dir, "KHAP/ConcordanceTables_", exp.dir,"_MCC_typed_genotype_t0.csv"), " TO DISK"))
message(paste0("WROTE ", paste0(outward_dir, "KHAP/ConcordanceTables_", exp.dir,"_MCC_typed_genotype_t1.csv"), " TO DISK"))
message(paste0("WROTE ", paste0(outward_dir, "KHAP/ConcordanceTables_", exp.dir,"_MCC_typed_genotype_t2.csv"), " TO DISK"))
message(paste0("WROTE ", paste0(outward_dir, "KHAP/ConcordanceTables_", exp.dir,"_MCC_typed_genotype_t3.csv"), " TO DISK"))

## MAF 

exp.dir = "MAF_Experiments"
ref.panel = c("ReferencePanel_v2", "ReferencePanel_v4", "ReferencePanel_v3")
exp.var = c("MAF1%", "MAF0.5%", "MAF0.1%")

tmp_maf_imputed_df = data.frame(matrix(ncol = 1, nrow = nrow(chip.table)))
names(tmp_maf_imputed_df) = c("chip")
tmp_maf_imputed_df$chip = chip.table$V1
tmp_maf_imputed_df$experiment = exp.dir

main_maf_imputed_df = data.frame()
main_maf_imputed_df_t0 = data.frame()
main_maf_imputed_df_t1 = data.frame()
main_maf_imputed_df_t2 = data.frame()
main_maf_imputed_df_t3 = data.frame()

main_maf_typed_df = data.frame()
main_maf_typed_df_t0 = data.frame()
main_maf_typed_df_t1 = data.frame()
main_maf_typed_df_t2 = data.frame()
main_maf_typed_df_t3 = data.frame()

for (exp in 1:length(exp.var)) {
  message(paste0(exp.var[exp]))
  # out.file = paste0("/Volumes/TimMcInerney/MitoImpute/data/HAPLOGROUPS/combined/ConcordanceTables_", ref.panel[exp], ".csv")
  imp.out.file = paste0(outward_dir, "MAF/ConcordanceTables_", ref.panel[exp], "_MCC_imputed_genotype.csv")
  imp.out.file.t0 = paste0(outward_dir, "MAF/ConcordanceTables_", ref.panel[exp], "_MCC_imputed_genotype_t0.csv")
  imp.out.file.t1 = paste0(outward_dir, "MAF/ConcordanceTables_", ref.panel[exp], "_MCC_imputed_genotype_t1.csv")
  imp.out.file.t2 = paste0(outward_dir, "MAF/ConcordanceTables_", ref.panel[exp], "_MCC_imputed_genotype_t2.csv")
  imp.out.file.t3 = paste0(outward_dir, "MAF/ConcordanceTables_", ref.panel[exp], "_MCC_imputed_genotype_t3.csv")
  
  typ.out.file = paste0(outward_dir, "MAF/ConcordanceTables_", ref.panel[exp], "_MCC_typed_genotype.csv")
  typ.out.file.t0 = paste0(outward_dir, "MAF/ConcordanceTables_", ref.panel[exp], "_MCC_typed_genotype_t0.csv")
  typ.out.file.t1 = paste0(outward_dir, "MAF/ConcordanceTables_", ref.panel[exp], "_MCC_typed_genotype_t1.csv")
  typ.out.file.t2 = paste0(outward_dir, "MAF/ConcordanceTables_", ref.panel[exp], "_MCC_typed_genotype_t2.csv")
  typ.out.file.t3 = paste0(outward_dir, "MAF/ConcordanceTables_", ref.panel[exp], "_MCC_typed_genotype_t3.csv")
  
  tmp_maf_imputed_df$sub_experiment = exp.var[exp]
  #tmp_maf_typed_df$sub_experiment = exp.var[exp]
  for (chip in 1:nrow(chip.table)) {
    # tmp1.file = paste0(container, chip.table$V1[chip], "/", ref.panel[exp], "/", "chrMT_1kg_", chip.table$V1[chip], "_diploid.txt")
    tmp.imp.file = paste0(container, chip.table$V1[chip], "/", ref.panel[exp], "/MCMC1/", "chrMT_1kg_", chip.table$V1[chip], "_imputed_MCMC1", "_imputed_MCC.csv")
    tmp.typ.file = paste0(container, chip.table$V1[chip], "/", ref.panel[exp], "/MCMC1/", "chrMT_1kg_", chip.table$V1[chip], "_imputed_MCMC1", "_typed_MCC.csv")
    
    platform.snps = paste0(container, chip.table$V1[chip], "/", chip.table$V1[chip], "_MT_snps.txt")
    
    if (file.exists(tmp.imp.file) == T) {
      tmp_maf_imputed_df$imputed[chip] = T
      chip.table$imputed[chip] = T
      #message(T)
    } else {
      tmp_maf_imputed_df$imputed[chip] = F
      chip.table$imputed[chip] = F
    }
    
    if (file.exists(platform.snps) == T) {
      snps = read.table(platform.snps, header = F, sep = "\t")
      tmp_maf_imputed_df$n.snps[chip] = nrow(snps)
    } else {
      tmp_maf_imputed_df$n.snps[chip] = 0
    }
    
    tmp_maf_imputed_df$n.snps[chip] = nrow(snps)
    tmp_maf_imputed_df$allele_freq[chip] = NA
    tmp_maf_imputed_df$mcc[chip] = NA
    tmp_maf_imputed_df$concordance[chip] = NA
    tmp_maf_imputed_df$info[chip] = NA
    tmp_maf_imputed_df$certainty[chip] = NA
    tmp_maf_imputed_df$info_type0[chip] = NA
    tmp_maf_imputed_df$concord_type0[chip] = NA
    tmp_maf_imputed_df$r2_type0[chip] = NA
    tmp_maf_imputed_df$info_comb[chip] = NA
  }
  
  tmp_maf_typed_df = tmp_maf_imputed_df
  
  tmp_maf_imputed_df_t0 = tmp_maf_imputed_df
  tmp_maf_imputed_df_t1 = tmp_maf_imputed_df
  tmp_maf_imputed_df_t2 = tmp_maf_imputed_df
  tmp_maf_imputed_df_t3 = tmp_maf_imputed_df
  
  tmp_maf_typed_df_t0 = tmp_maf_imputed_df_t0
  tmp_maf_typed_df_t1 = tmp_maf_imputed_df_t1
  tmp_maf_typed_df_t2 = tmp_maf_imputed_df_t2
  tmp_maf_typed_df_t3 = tmp_maf_imputed_df_t3
  
  total = nrow(chip.table)
  total_imputed = nrow(subset(chip.table, chip.table$imputed == T))
  
  for (chip in 1:nrow(chip.table)) {
    message(paste0("WORKING ON CHIP FOR ", exp.var[exp], ":\t", chip, " / ", total))
    
    # tmp1.file = paste0(container, chip.table$V1[chip], "/", ref.panel[exp], "/", "chrMT_1kg_", chip.table$V1[chip], "_diploid.txt")
    tmp.imp.file = paste0(container, chip.table$V1[chip], "/", ref.panel[exp], "/MCMC1/", "chrMT_1kg_", chip.table$V1[chip], "_imputed_MCMC1", "_imputed_MCC.csv")
    tmp.typ.file = paste0(container, chip.table$V1[chip], "/", ref.panel[exp], "/MCMC1/", "chrMT_1kg_", chip.table$V1[chip], "_imputed_MCMC1", "_typed_MCC.csv")
    
    # IMPUTED FILE
    if (file.exists(tmp.imp.file) == T) {
      tmp.imp.table = read.csv(tmp.imp.file, header = T)
      
      # THE FULL THING
      tmp_maf_imputed_df$allele_freq[chip]       = mean(tmp.imp.table$af, na.rm = T)
      tmp_maf_imputed_df$mcc[chip]               = mean(tmp.imp.table$mcc, na.rm = T)
      tmp_maf_imputed_df$concordance[chip]       = mean(tmp.imp.table$concodance, na.rm = T)
      tmp_maf_imputed_df$info[chip]              = mean(tmp.imp.table$info, na.rm = T)
      tmp_maf_imputed_df$certainty[chip]         = mean(tmp.imp.table$certainty, na.rm = T)
      tmp_maf_imputed_df$info_type0[chip]        = mean(tmp.imp.table$info_type0, na.rm = T)
      tmp_maf_imputed_df$concord_type0[chip]     = mean(tmp.imp.table$concord_type0, na.rm = T)
      tmp_maf_imputed_df$r2_type0[chip]          = mean(tmp.imp.table$r2_type0, na.rm = T)
      tmp_maf_imputed_df$info_comb[chip]         = mean(tmp.imp.table$info_comb, na.rm = T)
      
      # TYPE 0 ONLY!
      tmp.imp.table_t0 = subset(tmp.imp.table, tmp.imp.table$type == 0)
      
      tmp_maf_imputed_df_t0$allele_freq[chip]       = mean(tmp.imp.table_t0$af, na.rm = T)
      tmp_maf_imputed_df_t0$mcc[chip]               = mean(tmp.imp.table_t0$mcc, na.rm = T)
      tmp_maf_imputed_df_t0$concordance[chip]       = mean(tmp.imp.table_t0$concodance, na.rm = T)
      tmp_maf_imputed_df_t0$info[chip]              = mean(tmp.imp.table_t0$info, na.rm = T)
      tmp_maf_imputed_df_t0$certainty[chip]         = mean(tmp.imp.table_t0$certainty, na.rm = T)
      tmp_maf_imputed_df_t0$info_type0[chip]        = mean(tmp.imp.table_t0$info_type0, na.rm = T)
      tmp_maf_imputed_df_t0$concord_type0[chip]     = mean(tmp.imp.table_t0$concord_type0, na.rm = T)
      tmp_maf_imputed_df_t0$r2_type0[chip]          = mean(tmp.imp.table_t0$r2_type0, na.rm = T)
      tmp_maf_imputed_df_t0$info_comb[chip]         = mean(tmp.imp.table_t0$info_comb, na.rm = T)
      # TYPE 1 ONLY!
      tmp.imp.table_t1 = subset(tmp.imp.table, tmp.imp.table$type == 1)
      
      tmp_maf_imputed_df_t1$allele_freq[chip]       = mean(tmp.imp.table_t1$af, na.rm = T)
      tmp_maf_imputed_df_t1$mcc[chip]               = mean(tmp.imp.table_t1$mcc, na.rm = T)
      tmp_maf_imputed_df_t1$concordance[chip]       = mean(tmp.imp.table_t1$concodance, na.rm = T)
      tmp_maf_imputed_df_t1$info[chip]              = mean(tmp.imp.table_t1$info, na.rm = T)
      tmp_maf_imputed_df_t1$certainty[chip]         = mean(tmp.imp.table_t1$certainty, na.rm = T)
      tmp_maf_imputed_df_t1$info_type0[chip]        = mean(tmp.imp.table_t1$info_type0, na.rm = T)
      tmp_maf_imputed_df_t1$concord_type0[chip]     = mean(tmp.imp.table_t1$concord_type0, na.rm = T)
      tmp_maf_imputed_df_t1$r2_type0[chip]          = mean(tmp.imp.table_t1$r2_type0, na.rm = T)
      tmp_maf_imputed_df_t1$info_comb[chip]         = mean(tmp.imp.table_t1$info_comb, na.rm = T)
      # TYPE 2 ONLY!
      tmp.imp.table_t2 = subset(tmp.imp.table, tmp.imp.table$type == 2)
      
      tmp_maf_imputed_df_t2$allele_freq[chip]       = mean(tmp.imp.table_t2$af, na.rm = T)
      tmp_maf_imputed_df_t2$mcc[chip]               = mean(tmp.imp.table_t2$mcc, na.rm = T)
      tmp_maf_imputed_df_t2$concordance[chip]       = mean(tmp.imp.table_t2$concodance, na.rm = T)
      tmp_maf_imputed_df_t2$info[chip]              = mean(tmp.imp.table_t2$info, na.rm = T)
      tmp_maf_imputed_df_t2$certainty[chip]         = mean(tmp.imp.table_t2$certainty, na.rm = T)
      tmp_maf_imputed_df_t2$info_type0[chip]        = mean(tmp.imp.table_t2$info_type0, na.rm = T)
      tmp_maf_imputed_df_t2$concord_type0[chip]     = mean(tmp.imp.table_t2$concord_type0, na.rm = T)
      tmp_maf_imputed_df_t2$r2_type0[chip]          = mean(tmp.imp.table_t2$r2_type0, na.rm = T)
      tmp_maf_imputed_df_t2$info_comb[chip]         = mean(tmp.imp.table_t2$info_comb, na.rm = T)
      # TYPE 3 ONLY!
      tmp.imp.table_t3 = subset(tmp.imp.table, tmp.imp.table$type == 3)
      
      tmp_maf_imputed_df_t3$allele_freq[chip]       = mean(tmp.imp.table_t3$af, na.rm = T)
      tmp_maf_imputed_df_t3$mcc[chip]               = mean(tmp.imp.table_t3$mcc, na.rm = T)
      tmp_maf_imputed_df_t3$concordance[chip]       = mean(tmp.imp.table_t3$concodance, na.rm = T)
      tmp_maf_imputed_df_t3$info[chip]              = mean(tmp.imp.table_t3$info, na.rm = T)
      tmp_maf_imputed_df_t3$certainty[chip]         = mean(tmp.imp.table_t3$certainty, na.rm = T)
      tmp_maf_imputed_df_t3$info_type0[chip]        = mean(tmp.imp.table_t3$info_type0, na.rm = T)
      tmp_maf_imputed_df_t3$concord_type0[chip]     = mean(tmp.imp.table_t3$concord_type0, na.rm = T)
      tmp_maf_imputed_df_t3$r2_type0[chip]          = mean(tmp.imp.table_t3$r2_type0, na.rm = T)
      tmp_maf_imputed_df_t3$info_comb[chip]         = mean(tmp.imp.table_t3$info_comb, na.rm = T)
    }
    
    # GENOTYPED FILE
    if (file.exists(tmp.typ.file) == T) {
      tmp.typ.table = read.csv(tmp.typ.file, header = T)
      
      # THE FULL THING
      tmp_maf_typed_df$allele_freq[chip]       = mean(tmp.typ.table$af, na.rm = T)
      tmp_maf_typed_df$mcc[chip]               = mean(tmp.typ.table$mcc, na.rm = T)
      tmp_maf_typed_df$concordance[chip]       = mean(tmp.typ.table$concodance, na.rm = T)
      tmp_maf_typed_df$info[chip]              = mean(tmp.typ.table$info, na.rm = T)
      tmp_maf_typed_df$certainty[chip]         = mean(tmp.typ.table$certainty, na.rm = T)
      tmp_maf_typed_df$info_type0[chip]        = mean(tmp.typ.table$info_type0, na.rm = T)
      tmp_maf_typed_df$concord_type0[chip]     = mean(tmp.typ.table$concord_type0, na.rm = T)
      tmp_maf_typed_df$r2_type0[chip]          = mean(tmp.typ.table$r2_type0, na.rm = T)
      tmp_maf_typed_df$info_comb[chip]         = mean(tmp.typ.table$info_comb, na.rm = T)
      
      # TYPE 0 ONLY!
      tmp.typ.table_t0 = subset(tmp.typ.table, tmp.typ.table$type == 0)
      
      tmp_maf_typed_df_t0$allele_freq[chip]       = mean(tmp.typ.table_t0$af, na.rm = T)
      tmp_maf_typed_df_t0$mcc[chip]               = mean(tmp.typ.table_t0$mcc, na.rm = T)
      tmp_maf_typed_df_t0$concordance[chip]       = mean(tmp.typ.table_t0$concodance, na.rm = T)
      tmp_maf_typed_df_t0$info[chip]              = mean(tmp.typ.table_t0$info, na.rm = T)
      tmp_maf_typed_df_t0$certainty[chip]         = mean(tmp.typ.table_t0$certainty, na.rm = T)
      tmp_maf_typed_df_t0$info_type0[chip]        = mean(tmp.typ.table_t0$info_type0, na.rm = T)
      tmp_maf_typed_df_t0$concord_type0[chip]     = mean(tmp.typ.table_t0$concord_type0, na.rm = T)
      tmp_maf_typed_df_t0$r2_type0[chip]          = mean(tmp.typ.table_t0$r2_type0, na.rm = T)
      tmp_maf_typed_df_t0$info_comb[chip]         = mean(tmp.typ.table_t0$info_comb, na.rm = T)
      # TYPE 1 ONLY!
      tmp.typ.table_t1 = subset(tmp.typ.table, tmp.typ.table$type == 1)
      
      tmp_maf_typed_df_t1$allele_freq[chip]       = mean(tmp.typ.table_t1$af, na.rm = T)
      tmp_maf_typed_df_t1$mcc[chip]               = mean(tmp.typ.table_t1$mcc, na.rm = T)
      tmp_maf_typed_df_t1$concordance[chip]       = mean(tmp.typ.table_t1$concodance, na.rm = T)
      tmp_maf_typed_df_t1$info[chip]              = mean(tmp.typ.table_t1$info, na.rm = T)
      tmp_maf_typed_df_t1$certainty[chip]         = mean(tmp.typ.table_t1$certainty, na.rm = T)
      tmp_maf_typed_df_t1$info_type0[chip]        = mean(tmp.typ.table_t1$info_type0, na.rm = T)
      tmp_maf_typed_df_t1$concord_type0[chip]     = mean(tmp.typ.table_t1$concord_type0, na.rm = T)
      tmp_maf_typed_df_t1$r2_type0[chip]          = mean(tmp.typ.table_t1$r2_type0, na.rm = T)
      tmp_maf_typed_df_t1$info_comb[chip]         = mean(tmp.typ.table_t1$info_comb, na.rm = T)
      # TYPE 2 ONLY!
      tmp.typ.table_t2 = subset(tmp.typ.table, tmp.typ.table$type == 2)
      
      tmp_maf_typed_df_t2$allele_freq[chip]       = mean(tmp.typ.table_t2$af, na.rm = T)
      tmp_maf_typed_df_t2$mcc[chip]               = mean(tmp.typ.table_t2$mcc, na.rm = T)
      tmp_maf_typed_df_t2$concordance[chip]       = mean(tmp.typ.table_t2$concodance, na.rm = T)
      tmp_maf_typed_df_t2$info[chip]              = mean(tmp.typ.table_t2$info, na.rm = T)
      tmp_maf_typed_df_t2$certainty[chip]         = mean(tmp.typ.table_t2$certainty, na.rm = T)
      tmp_maf_typed_df_t2$info_type0[chip]        = mean(tmp.typ.table_t2$info_type0, na.rm = T)
      tmp_maf_typed_df_t2$concord_type0[chip]     = mean(tmp.typ.table_t2$concord_type0, na.rm = T)
      tmp_maf_typed_df_t2$r2_type0[chip]          = mean(tmp.typ.table_t2$r2_type0, na.rm = T)
      tmp_maf_typed_df_t2$info_comb[chip]         = mean(tmp.typ.table_t2$info_comb, na.rm = T)
      # TYPE 3 ONLY!
      tmp.typ.table_t3 = subset(tmp.typ.table, tmp.typ.table$type == 3)
      
      tmp_maf_typed_df_t3$allele_freq[chip]       = mean(tmp.typ.table_t3$af, na.rm = T)
      tmp_maf_typed_df_t3$mcc[chip]               = mean(tmp.typ.table_t3$mcc, na.rm = T)
      tmp_maf_typed_df_t3$concordance[chip]       = mean(tmp.typ.table_t3$concodance, na.rm = T)
      tmp_maf_typed_df_t3$info[chip]              = mean(tmp.typ.table_t3$info, na.rm = T)
      tmp_maf_typed_df_t3$certainty[chip]         = mean(tmp.typ.table_t3$certainty, na.rm = T)
      tmp_maf_typed_df_t3$info_type0[chip]        = mean(tmp.typ.table_t3$info_type0, na.rm = T)
      tmp_maf_typed_df_t3$concord_type0[chip]     = mean(tmp.typ.table_t3$concord_type0, na.rm = T)
      tmp_maf_typed_df_t3$r2_type0[chip]          = mean(tmp.typ.table_t3$r2_type0, na.rm = T)
      tmp_maf_typed_df_t3$info_comb[chip]         = mean(tmp.typ.table_t3$info_comb, na.rm = T)
    }
  }
  write.csv(tmp_maf_imputed_df, imp.out.file, row.names = F, quote = F)
  #message(paste0("WROTE ", imp.out.file, " TO DISK"))
  write.csv(tmp_maf_imputed_df_t0, imp.out.file.t0, row.names = F, quote = F)
  #message(paste0("WROTE ", imp.out.file.t0, " TO DISK"))
  write.csv(tmp_maf_imputed_df_t1, imp.out.file.t1, row.names = F, quote = F)
  #message(paste0("WROTE ", imp.out.file.t1, " TO DISK"))
  write.csv(tmp_maf_imputed_df_t2, imp.out.file.t2, row.names = F, quote = F)
  #message(paste0("WROTE ", imp.out.file.t2, " TO DISK"))
  write.csv(tmp_maf_imputed_df_t3, imp.out.file.t3, row.names = F, quote = F)
  #message(paste0("WROTE ", imp.out.file.t3, " TO DISK"))
  #
  write.csv(tmp_maf_typed_df, typ.out.file, row.names = F, quote = F)
  #message(paste0("WROTE ", typ.out.file, " TO DISK"))
  write.csv(tmp_maf_typed_df_t0, typ.out.file.t0, row.names = F, quote = F)
  #message(paste0("WROTE ", typ.out.file.t0, " TO DISK"))
  write.csv(tmp_maf_typed_df_t1, typ.out.file.t1, row.names = F, quote = F)
  #message(paste0("WROTE ", typ.out.file.t1, " TO DISK"))
  write.csv(tmp_maf_typed_df_t2, typ.out.file.t2, row.names = F, quote = F)
  #message(paste0("WROTE ", typ.out.file.t2, " TO DISK"))
  write.csv(tmp_maf_typed_df_t3, typ.out.file.t3, row.names = F, quote = F)
  #message(paste0("WROTE ", typ.out.file.t3, " TO DISK"))
  
  main_maf_imputed_df = rbind(main_maf_imputed_df, tmp_maf_imputed_df)
  main_maf_imputed_df_t0 = rbind(main_maf_imputed_df_t0, tmp_maf_imputed_df_t0)
  main_maf_imputed_df_t1 = rbind(main_maf_imputed_df_t1, tmp_maf_imputed_df_t1)
  main_maf_imputed_df_t2 = rbind(main_maf_imputed_df_t2, tmp_maf_imputed_df_t2)
  main_maf_imputed_df_t3 = rbind(main_maf_imputed_df_t3, tmp_maf_imputed_df_t3)
  
  main_maf_typed_df = rbind(main_maf_typed_df, tmp_maf_typed_df)
  main_maf_typed_df_t0 = rbind(main_maf_typed_df_t0, tmp_maf_typed_df_t0)
  main_maf_typed_df_t1 = rbind(main_maf_typed_df_t1, tmp_maf_typed_df_t1)
  main_maf_typed_df_t2 = rbind(main_maf_typed_df_t2, tmp_maf_typed_df_t2)
  main_maf_typed_df_t3 = rbind(main_maf_typed_df_t3, tmp_maf_typed_df_t3)
}
write.csv(main_maf_imputed_df, paste0(outward_dir, "MAF/ConcordanceTables_", exp.dir,"_MCC_imputed_genotype.csv"), row.names = F, quote = F)
write.csv(main_maf_imputed_df_t0, paste0(outward_dir, "MAF/ConcordanceTables_", exp.dir,"_MCC_imputed_genotype_t0.csv"), row.names = F, quote = F)
write.csv(main_maf_imputed_df_t1, paste0(outward_dir, "MAF/ConcordanceTables_", exp.dir,"_MCC_imputed_genotype_t1.csv"), row.names = F, quote = F)
write.csv(main_maf_imputed_df_t2, paste0(outward_dir, "MAF/ConcordanceTables_", exp.dir,"_MCC_imputed_genotype_t2.csv"), row.names = F, quote = F)
write.csv(main_maf_imputed_df_t3, paste0(outward_dir, "MAF/ConcordanceTables_", exp.dir,"_MCC_imputed_genotype_t3.csv"), row.names = F, quote = F)
write.csv(main_maf_typed_df, paste0(outward_dir, "MAF/ConcordanceTables_", exp.dir,"_MCC_typed_genotype.csv"), row.names = F, quote = F)
write.csv(main_maf_typed_df_t0, paste0(outward_dir, "MAF/ConcordanceTables_", exp.dir,"_MCC_typed_genotype_t0.csv"), row.names = F, quote = F)
write.csv(main_maf_typed_df_t1, paste0(outward_dir, "MAF/ConcordanceTables_", exp.dir,"_MCC_typed_genotype_t1.csv"), row.names = F, quote = F)
write.csv(main_maf_typed_df_t2, paste0(outward_dir, "MAF/ConcordanceTables_", exp.dir,"_MCC_typed_genotype_t2.csv"), row.names = F, quote = F)
write.csv(main_maf_typed_df_t3, paste0(outward_dir, "MAF/ConcordanceTables_", exp.dir,"_MCC_typed_genotype_t3.csv"), row.names = F, quote = F)
message(paste0("WROTE ", paste0(outward_dir, "MAF/ConcordanceTables_", exp.dir,"_MCC_imputed_genotype.csv"), " TO DISK"))
message(paste0("WROTE ", paste0(outward_dir, "MAF/ConcordanceTables_", exp.dir,"_MCC_imputed_genotype_t0.csv"), " TO DISK"))
message(paste0("WROTE ", paste0(outward_dir, "MAF/ConcordanceTables_", exp.dir,"_MCC_imputed_genotype_t1.csv"), " TO DISK"))
message(paste0("WROTE ", paste0(outward_dir, "MAF/ConcordanceTables_", exp.dir,"_MCC_imputed_genotype_t2.csv"), " TO DISK"))
message(paste0("WROTE ", paste0(outward_dir, "MAF/ConcordanceTables_", exp.dir,"_MCC_imputed_genotype_t3.csv"), " TO DISK"))
message(paste0("WROTE ", paste0(outward_dir, "MAF/ConcordanceTables_", exp.dir,"_MCC_typed_genotype.csv"), " TO DISK"))
message(paste0("WROTE ", paste0(outward_dir, "MAF/ConcordanceTables_", exp.dir,"_MCC_typed_genotype_t0.csv"), " TO DISK"))
message(paste0("WROTE ", paste0(outward_dir, "MAF/ConcordanceTables_", exp.dir,"_MCC_typed_genotype_t1.csv"), " TO DISK"))
message(paste0("WROTE ", paste0(outward_dir, "MAF/ConcordanceTables_", exp.dir,"_MCC_typed_genotype_t2.csv"), " TO DISK"))
message(paste0("WROTE ", paste0(outward_dir, "MAF/ConcordanceTables_", exp.dir,"_MCC_typed_genotype_t3.csv"), " TO DISK"))

message("")
message("END!")
message("")
## END