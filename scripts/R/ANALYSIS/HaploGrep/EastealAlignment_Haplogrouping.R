library(dplyr)

haplogrep_file = "/Volumes/TimMcInerney/Master_Alignment/VCF/hsapiensCRS7k_ambigANDgap2missing_HaploGrep.txt"
haplogrep_csv = "/Volumes/TimMcInerney/Master_Alignment/VCF/hsapiensCRS7k_ambigANDgap2missing_HaploGrep.csv"
haplogroup_table_csv = "/Volumes/TimMcInerney/Master_Alignment/VCF/hsapiensCRS7k_ambigANDgap2missing_HaploGrep_freTable.csv"

haplogroups = read.table(haplogrep_file, header = T, sep = "\t")

nonAFR = subset(haplogroups, substr(haplogroups$Haplogroup, 1, 1) != "L")
nonAFR$Macrohaplogroup = substr(nonAFR$Haplogroup, 1, 1)

AFR = subset(haplogroups, substr(haplogroups$Haplogroup, 1, 1) == "L")
AFR$Macrohaplogroup = substr(AFR$Haplogroup, 1, 2)

fixed_haplogroups = rbind(nonAFR, AFR)
fixed_haplogroups = arrange(fixed_haplogroups, fixed_haplogroups$SampleID)

write.csv(fixed_haplogroups, haplogrep_csv, quote = F, row.names = F)

freqTable = as.data.frame(table(fixed_haplogroups$Macrohaplogroup))
names(freqTable) = c("Haplogroup", "Freq")
freqTable

write.csv(freqTable, haplogroup_table_csv, quote = F, row.names = F)
