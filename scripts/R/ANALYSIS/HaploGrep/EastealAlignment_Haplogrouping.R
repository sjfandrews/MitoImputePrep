library(dplyr)

Easteal_haplogrep_file = "/Volumes/TimMcInerney/Master_Alignment/VCF/hsapiensCRS7k_ambigANDgap2missing_HaploGrep.txt"
Easteal_haplogrep_csv = "/Volumes/TimMcInerney/Master_Alignment/VCF/hsapiensCRS7k_ambigANDgap2missing_HaploGrep.csv"
Easteal_haplogroup_table_csv = "/Volumes/TimMcInerney/Master_Alignment/VCF/hsapiensCRS7k_ambigANDgap2missing_HaploGrep_freTable.csv"

Easteal_haplogroups = read.table(Easteal_haplogrep_file, header = T, sep = "\t")

Easteal_nonAFR = subset(Easteal_haplogroups, substr(Easteal_haplogroups$Haplogroup, 1, 1) != "L")
Easteal_nonAFR$Macrohaplogroup = substr(Easteal_nonAFR$Haplogroup, 1, 1)

Easteal_AFR = subset(haplogroups, substr(Easteal_haplogroups$Haplogroup, 1, 1) == "L")
Easteal_AFR$Macrohaplogroup = substr(Easteal_AFR$Haplogroup, 1, 2)

Easteal_fixed_haplogroups = rbind(Easteal_nonAFR, Easteal_AFR)
Easteal_fixed_haplogroups = arrange(Easteal_fixed_haplogroups, Easteal_fixed_haplogroups$SampleID)

write.csv(Easteal_fixed_haplogroups, Easteal_haplogrep_csv, quote = F, row.names = F)

freqTable = as.data.frame(table(Easteal_fixed_haplogroups$Macrohaplogroup))
names(freqTable) = c("Haplogroup", "Freq")
freqTable

write.csv(freqTable, Easteal_haplogroup_table_csv, quote = F, row.names = F)
