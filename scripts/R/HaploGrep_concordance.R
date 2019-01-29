require(dplyr)
'%!in%' = function(x,y)!('%in%'(x,y))

strands_file = "~/GitCode/MitoImputePrep/scripts/INFORMATION_LISTS/b37_platforms.txt"
WGS_file = "~/GitCode/MitoImputePrep/metadata/HaploGrep_concordance/chrMT_1kg_diploid.haplogrep.txt"
geno_dir = "/Volumes/TimMcInerney/MitoImpute/analyses/Haplogroups/HaploGrep/ReferencePanel_v5/genotyped/"
imp_dir = "/Volumes/TimMcInerney/MitoImpute/analyses/Haplogroups/HaploGrep/ReferencePanel_v5/imputed/"
comp_dir = "/Volumes/TimMcInerney/MitoImpute/analyses/Haplogroups/HaploGrep/ReferencePanel_v5/comparison/"

if (!endsWith(geno_dir, "/")) {
  geno_dir = paste0(geno_dir, "/")
}
if (!endsWith(imp_dir, "/")) {
  imp_dir = paste0(imp_dir, "/")
}
if (!endsWith(comp_dir, "/")) {
  comp_dir = paste0(comp_dir, "/")
}

analysed_chips = sub(x = dir(imp_dir), pattern = "_imputed_MCMC1_haplogrep.txt", replacement = "")
analysed_chips = sub(x = analysed_chips, pattern = "chrMT_1kg_", replacement = "")

strands = read.table(strands_file, header = F)
names(strands) = c("chip")
strands = subset(strands, strands$chip %in% analysed_chips)

# GATHER WGS HAPLOGROUPING
wgs_table = read.table(WGS_file, header = T)
afr_hgs = subset(wgs_table, substr(wgs_table$Haplogroup, 1, 1) == "L")
afr_hgs$MetaHaplogroup = substr(afr_hgs$Haplogroup, 1, 2)
eua_hgs = subset(wgs_table, substr(wgs_table$Haplogroup, 1, 1) != "L")
eua_hgs$MetaHaplogroup = substr(eua_hgs$Haplogroup, 1, 1)
wgs_table = rbind(afr_hgs, eua_hgs)

wgs_table = arrange(wgs_table, wgs_table$SampleID)

hgs = names(table(wgs_table$Haplogroup))
mhgs = names(table(wgs_table$MetaHaplogroup))


### META HAPLOGROUPS
# CREATE A DATA FRAME WITH THE NUMBER OF ROWS = #CHIPS and NUMBER OF COLUMNS #HAPLOGROUPS+1
imputed_df = data.frame(matrix(nrow = length(dir(analysed_chips)), ncol = length(mhgs) + 1))
names(imputed_df) = c("chip", mhgs)

typed_df = data.frame(matrix(nrow = length(dir(analysed_chips)), ncol = length(mhgs) + 1))
names(typed_df) = c("chip", mhgs)

# CALCULATE CONCORDANCE FOR EACH CHIP
for (chip in 1:length(analysed_chips)) {
  message(paste0(chip, " / ", length(analysed_chips)))
  if (file.exists(paste0(geno_dir, "chrMT_1kg_", analysed_chips[chip], "chrMT_1kg_", analysed_chips[chip], "_diploid.txt")) && file.exists(paste0(imp_dir, "chrMT_1kg_", analysed_chips[chip], "_imputed_MCMC1_haplogrep.txt"))) {
    tmp_table = data.frame(cbind(as.character(wgs_table$SampleID), wgs_table$MetaHaplogroup))
    names(tmp_table) = c("sample", "wgs")
    temp_gen = read.table(paste0(geno_dir, "chrMT_1kg_", analysed_chips[chip], "chrMT_1kg_", analysed_chips[chip], "_diploid.txt"), header = T)
    temp_imp = read.table(paste0(imp_dir, "chrMT_1kg_", analysed_chips[chip], "_imputed_MCMC1_haplogrep.txt"), header = T)
    
    # GENOTYPED
    tmp_gen_afr = subset(temp_gen, substr(temp_gen$Haplogroup, 1, 1) == "L")
    tmp_gen_afr$MetaHaplogroup = substr(tmp_gen_afr$Haplogroup, 1, 2)
    tmp_gen_eua = subset(temp_gen, substr(temp_gen$Haplogroup, 1, 1) != "L")
    tmp_gen_eua$MetaHaplogroup = substr(tmp_gen_eua$Haplogroup, 1, 1)
    temp_gen = rbind(tmp_gen_afr, tmp_gen_eua)
    temp_gen = arrange(temp_gen, temp_gen$SampleID)
    
    tmp_table$genotyped = temp_gen$MetaHaplogroup
    tmp_table$geno_match = tmp_table$wgs == tmp_table$genotyped
    
    # IMPUTED
    tmp_imp_afr = subset(temp_imp, substr(temp_imp$Haplogroup, 1, 1) == "L")
    tmp_imp_afr$MetaHaplogroup = substr(tmp_imp_afr$Haplogroup, 1, 2)
    tmp_imp_eua = subset(temp_imp, substr(temp_imp$Haplogroup, 1, 1) != "L")
    tmp_imp_eua$MetaHaplogroup = substr(tmp_imp_eua$Haplogroup, 1, 1)
    temp_imp = rbind(tmp_imp_afr, tmp_imp_eua)
    temp_imp = arrange(temp_imp, temp_imp$SampleID)
    
    tmp_table$imputed = temp_imp$MetaHaplogroup
    tmp_table$imp_match = tmp_table$wgs == tmp_table$imputed
    
    # SAVE THE FILE
    write.csv(tmp_table, paste0(comp_dir, "meta_haplogroup/", analysed_chips[chip], "_HaploGrep_concordance.csv"), row.names = F, quote = F)
  }
}




# CREATE A DATA FRAME WITH THE NUMBER OF ROWS = #CHIPS and NUMBER OF COLUMNS #HAPLOGROUPS+1
imputed_df = data.frame(matrix(nrow = length(dir(analysed_chips)), ncol = length(hgs) + 1))
names(imputed_df) = c("chip", hgs)

typed_df = data.frame(matrix(nrow = length(dir(analysed_chips)), ncol = length(hgs) + 1))
names(typed_df) = c("chip", hgs)

# RUN THE CALCULATIONS
for (chip in 1:length(dir(conc_dir))) { # FOR EACH CHIP
  print(dir(conc_dir)[chip])
  imputed_df$chip[chip] = sub(x = dir(conc_dir)[chip], pattern = "_haplogroupings.csv", replacement = "") # STRIP EXTENSION FROM FILE NAME AND PLACE CHIP NAME IN COLUMN
  typed_df$chip[chip] = sub(x = dir(conc_dir)[chip], pattern = "_haplogroupings.csv", replacement = "") # STRIP EXTENSION FROM FILE NAME AND PLACE CHIP NAME IN COLUMN
  
  tmp_chip_hgs = read.csv(paste0(conc_dir, dir(conc_dir)[chip]), header = T) # READ THE HAPLOGROUP CONCORDANCE CSV IN
  
  for (hg in 1:length(hgs)) { # FOR EACH HAPLOGROUP IN THE WGS (TRUTH) COLUMN ...
    s = subset(tmp_chip_hgs, tmp_chip_hgs$WGS == hgs[hg]) # SUBSET THE HAPLOGROUP CONCORDANCE CSV FOR ONLY THAT HAPLOGROUP
    t_i = subset(s, s$match.i) # SUBSET TO ONLY IMPUTED MATCHES
    t_t = subset(s, s$match.t) # SUBSET TO ONLY TYPED MATCHES
    #v = nrow(t) / nrow(s)
    imputed_df[chip, hgs[hg]] = nrow(t_i) / nrow(s) # WRITE PROPORTION MATCHING TO CELL FOR IMPUTED
    typed_df[chip, hgs[hg]] = nrow(t_t) / nrow(s) # WRITE PROPORTION MATCHING TO CELL FOR TYPED
  }
}

hgs_means = data.frame(hgs ,colMeans(typed_df[,2:ncol(typed_df)], na.rm = T), colMeans(imputed_df[,2:ncol(imputed_df)], na.rm = T)) # GET THE MEAN VALUE FOR EACH HAPLOGROUP
names(hgs_means) = c("haplogroup", "typed_true", "imputed_true")
hgs_means$diff = hgs_means$imputed_true - hgs_means$typed_true # CALCULATE THE DIFFERENCE

write.csv(imputed_df, paste0(out_dir, "imputed_haplogroup_concordance.csv"), row.names = F, quote = F)
write.csv(typed_df, paste0(out_dir, "typed_haplogroup_concordance.csv"), row.names = F, quote = F)
write.csv(hgs_means, paste0(out_dir, "haplogroup_means.csv"), row.names = F, quote = F)



### META HAPLOGROUP ANALYSIS
# SPECIFY LISTS OF HAPLOGROUPS
mhgs = c(LETTERS, "HV", "JT", "unclassified")
easy_mhgs = LETTERS
difficult_mhgs = c("HV", "JT", "unclassified")

# GATHER WGS HAPLOGROUPING AND TRUNCATE HAPLOGROUPING TO META-HAPLOGROUPS
Primary_table = read.csv(paste0(conc_dir, as.character(strands$chip[1]), "_haplogroupings.csv"), header = T)
difficult_df = subset(Primary_table, Primary_table$WGS %in% difficult_mhgs) # SUBSET FOR HGs JT, HV, & unclassified
easy_df = subset(Primary_table, Primary_table$WGS %!in% difficult_mhgs) # SUBSET FOR THE REST
easy_df$WGS = substr(easy_df$WGS, 1, 1) # TRUNCATE HAPLOGROUPS IN WGS TO META-HAPLOGROUP ONLY
Primary_table = rbind(easy_df, difficult_df) # STITCH THEM BACK TOGETHER
hgs = names(table(Primary_table$WGS))

# CREATE A DATA FRAME WITH THE NUMBER OF ROWS = #CHIPS and NUMBER OF COLUMNS #HAPLOGROUPS+1
imputed_df = data.frame(matrix(nrow = length(dir(conc_dir)), ncol = length(hgs) + 1))
names(imputed_df) = c("chip", hgs)

typed_df = data.frame(matrix(nrow = length(dir(conc_dir)), ncol = length(hgs) + 1))
names(typed_df) = c("chip", hgs)

# RUN THE CALCULATIONS
for (chip in 1:length(dir(conc_dir))) { # FOR EACH CHIP
  print(dir(conc_dir)[chip])
  imputed_df$chip[chip] = sub(x = dir(conc_dir)[chip], pattern = "_haplogroupings.csv", replacement = "") # STRIP EXTENSION FROM FILE NAME AND PLACE CHIP NAME IN COLUMN
  typed_df$chip[chip] = sub(x = dir(conc_dir)[chip], pattern = "_haplogroupings.csv", replacement = "") # STRIP EXTENSION FROM FILE NAME AND PLACE CHIP NAME IN COLUMN
  
  tmp_chip_hgs = read.csv(paste0(conc_dir, dir(conc_dir)[chip]), header = T) # READ THE HAPLOGROUP CONCORDANCE CSV IN
  
  # TRUNCATE HAPLOGROUPING TO META-HAPLOGROUPS
  difficult_df = subset(tmp_chip_hgs, tmp_chip_hgs$WGS %in% difficult_mhgs) # SUBSET FOR HGs JT, HV, & unclassified
  easy_df = subset(tmp_chip_hgs, tmp_chip_hgs$WGS %!in% difficult_mhgs) # SUBSET FOR THE REST
  easy_df$WGS = substr(easy_df$WGS, 1, 1) # TRUNCATE HAPLOGROUPS IN WGS TO META-HAPLOGROUP ONLY
  easy_df$typed = as.character(easy_df$typed) # CONVERT FROM FACTOR TO CHARACTER
  easy_df$imputed = as.character(easy_df$imputed) # CONVERT FROM FACTOR TO CHARACTER
  
  # TRUNCATE THE TYPED AND IMPUTED HAPLOGROUPS TO META-HAPLOGROUPS IF THEY AREN'T JT OR HV OR unclassified
  for (i in 1:nrow(easy_df)) {
    #print(i)
    if (!any(easy_df$typed[i] == difficult_mhgs)) {
      easy_df$typed[i] = substr(easy_df$typed[i], 1, 1)
    }
    if (!any(easy_df$imputed[i] == difficult_mhgs)) {
      easy_df$imputed[i] = substr(easy_df$imputed[i], 1, 1)
    }
  }
  tmp_chip_hgs = rbind(easy_df, difficult_df)
  # 
  
  for (hg in 1:length(hgs)) { # FOR EACH HAPLOGROUP IN THE WGS (TRUTH) COLUMN ...
    s = subset(tmp_chip_hgs, tmp_chip_hgs$WGS == hgs[hg]) # SUBSET THE HAPLOGROUP CONCORDANCE CSV FOR ONLY THAT HAPLOGROUP
    t_i = subset(s, s$match.i) # SUBSET TO ONLY IMPUTED MATCHES
    t_t = subset(s, s$match.t) # SUBSET TO ONLY TYPED MATCHES
    #v = nrow(t) / nrow(s)
    imputed_df[chip, hgs[hg]] = nrow(t_i) / nrow(s) # WRITE PROPORTION MATCHING TO CELL FOR IMPUTED
    typed_df[chip, hgs[hg]] = nrow(t_t) / nrow(s) # WRITE PROPORTION MATCHING TO CELL FOR TYPED
  }
}

hgs_means = data.frame(hgs ,colMeans(typed_df[,2:ncol(typed_df)], na.rm = T), colMeans(imputed_df[,2:ncol(imputed_df)], na.rm = T)) # GET THE MEAN VALUE FOR EACH HAPLOGROUP
names(hgs_means) = c("haplogroup", "typed_true", "imputed_true")
hgs_means$diff = hgs_means$imputed_true - hgs_means$typed_true # CALCULATE THE DIFFERENCE

write.csv(imputed_df, paste0(out_dir2, "imputed_haplogroup_concordance.csv"), row.names = F, quote = F)
write.csv(typed_df, paste0(out_dir2, "typed_haplogroup_concordance.csv"), row.names = F, quote = F)
write.csv(hgs_means, paste0(out_dir2, "haplogroup_means.csv"), row.names = F, quote = F)
