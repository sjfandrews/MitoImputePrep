'%!in%' = function(x,y)!('%in%'(x,y))

strands_file = "~/GitCode/MitoImputePrep/scripts/INFORMATION_LISTS/b37_platforms.txt"
conc_dir = "~/GitCode/MitoImputePrep/metadata/haplogroup_concordance"
out_dir = "~/GitCode/MitoImputePrep/metadata/haplogroup_concordance_allChips/"
out_dir2 = "~/GitCode/MitoImputePrep/metadata/haplogroup_concordance_allChips/meta_haplogroups/"

if (!endsWith(conc_dir, "/")) {
  conc_dir = paste0(conc_dir, "/")
}
if (!endsWith(out_dir, "/")) {
  out_dir = paste0(out_dir, "/")
}
if (!endsWith(out_dir2, "/")) {
  out_dir2 = paste0(out_dir2, "/")
}

analysed_chips = sub(x = dir(conc_dir), pattern = "_haplogroupings.csv", replacement = "")

strands = read.table(strands_file, header = F)
names(strands) = c("chip")

# GATHER WGS HAPLOGROUPING
Primary_table = read.csv(paste0(conc_dir, as.character(strands$chip[1]), "_haplogroupings.csv"), header = T)
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
mhgs = c(LETTERS, "HV", "JT", "unclassified")
easy_mhgs = LETTERS
difficult_mhgs = c("HV", "JT", "unclassified")
easy_df = subset(tmp_chip_hgs, tmp_chip_hgs$WGS %!in% difficult_mhgs)
easy_df$WGS = substr(easy_df$WGS, 1, 1)
easy_df$typed = as.character(easy_df$typed)
easy_df$imputed = as.character(easy_df$imputed)

for (i in 1:nrow(easy_df)) {
  print(i)
  if (!any(easy_df$typed[i] == difficult_mhgs)) {
    easy_df$typed[i] = substr(easy_df$typed[i], 1, 1)
  }
  if (!any(easy_df$imputed[i] == difficult_mhgs)) {
    easy_df$imputed[i] = substr(easy_df$imputed[i], 1, 1)
  }
}


easy_df$match.t = easy_df$typed == easy_df$WGS
easy_df$match.i = easy_df$imputed == easy_df$WGS
difficult_df = subset(tmp_chip_hgs, tmp_chip_hgs$WGS %in% difficult_mhgs)

fixed_df = rbind(easy_df, difficult_df)
