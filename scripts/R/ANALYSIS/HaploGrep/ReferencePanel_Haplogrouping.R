library(dplyr)

# FILTERED REFERENCE PANEL

haplogrep_file = "/Volumes/TimMcInerney/MitoImpute/data/VCF/McInerney_Master_Alignment_July18_2018_ambigANDgap2missing_tagFilled_filt_noInvar.txt"
haplogrep_csv = "/Volumes/TimMcInerney/MitoImpute/data/VCF/McInerney_Master_Alignment_July18_2018_ambigANDgap2missing_tagFilled_filt_noInvar.csv"
haplogroup_table_csv = "/Volumes/TimMcInerney/MitoImpute/data/VCF/McInerney_Master_Alignment_July18_2018_ambigANDgap2missing_tagFilled_filt_noInvar_freqTable.csv"

haplogroups = read.table(haplogrep_file, header = T, sep = "\t")

nonAFR = subset(haplogroups, substr(haplogroups$Haplogroup, 1, 1) != "L")
nonAFR$Macrohaplogroup = substr(nonAFR$Haplogroup, 1, 1)

AFR = subset(haplogroups, substr(haplogroups$Haplogroup, 1, 1) == "L")
AFR$Macrohaplogroup = substr(AFR$Haplogroup, 1, 2)

fixed_haplogroups = rbind(nonAFR, AFR)
fixed_haplogroups = arrange(fixed_haplogroups, fixed_haplogroups$SampleID)

write.csv(fixed_haplogroups, haplogrep_csv, quote = F, row.names = F)

## UNFILTERED REFERENCE PANEL

full_haplogrep_file = "/Volumes/TimMcInerney/MitoImpute/data/VCF/McInerney_Master_Alignment_July18_2018_ambigANDgap2missing_tagFilled_noInvar_HaploGrep.txt"
full_haplogrep_csv = "/Volumes/TimMcInerney/MitoImpute/data/VCF/McInerney_Master_Alignment_July18_2018_ambigANDgap2missing_tagFilled_noInvar_HaploGrep.csv"
haplogroup_table_csv = "/Volumes/TimMcInerney/MitoImpute/data/VCF/McInerney_Master_Alignment_July18_2018_ambigANDgap2missing_tagFilled_noInvar_HaploGrep_freqTable.csv"

full_haplogroups = read.table(full_haplogrep_file, header = T, sep = "\t")
full_nonAFR = subset(full_haplogroups, substr(full_haplogroups$Haplogroup, 1, 1) != "L")
full_nonAFR$Macrohaplogroup = substr(full_nonAFR$Haplogroup, 1, 1)

full_AFR = subset(full_haplogroups, substr(full_haplogroups$Haplogroup, 1, 1) == "L")
full_AFR$Macrohaplogroup = substr(full_AFR$Haplogroup, 1, 2)

full_fixed_haplogroups = rbind(full_nonAFR, full_AFR)
full_fixed_haplogroups = arrange(full_fixed_haplogroups, full_fixed_haplogroups$SampleID)

write.csv(full_fixed_haplogroups, full_haplogrep_csv, quote = F, row.names = F)

## EASTEAL MASTER REFERENCE ALIGNMENT
Easteal_haplogrep_file = "/Volumes/TimMcInerney/Master_Alignment/VCF/hsapiensCRS7k_ambigANDgap2missing_HaploGrep.txt"
Easteal_haplogroups = read.table(Easteal_haplogrep_file, header = T, sep = "\t")

Easteal_nonAFR = subset(Easteal_haplogroups, substr(Easteal_haplogroups$Haplogroup, 1, 1) != "L")
Easteal_nonAFR$Macrohaplogroup = substr(Easteal_nonAFR$Haplogroup, 1, 1)

Easteal_AFR = subset(haplogroups, substr(Easteal_haplogroups$Haplogroup, 1, 1) == "L")
Easteal_AFR$Macrohaplogroup = substr(Easteal_AFR$Haplogroup, 1, 2)

Easteal_fixed_haplogroups = rbind(Easteal_nonAFR, Easteal_AFR)
Easteal_fixed_haplogroups = arrange(Easteal_fixed_haplogroups, Easteal_fixed_haplogroups$SampleID)
######################

# FILTERED HAPLOGROUPS TABLE
table(fixed_haplogroups$Macrohaplogroup)

# FULL HAPLOGROUPS TABLE
table(full_fixed_haplogroups$Macrohaplogroup)

# DIFFERENCE HAPLOGROUP TABLE
(table(full_fixed_haplogroups$Macrohaplogroup) - table(fixed_haplogroups$Macrohaplogroup)) / sum(table(full_fixed_haplogroups$Macrohaplogroup) - table(fixed_haplogroups$Macrohaplogroup)) * 100
nrow(full_fixed_haplogroups) - nrow(fixed_haplogroups)

freqTable = as.data.frame(table(Easteal_fixed_haplogroups$Macrohaplogroup))
freqTable$Freq1 = as.data.frame(table(fixed_haplogroups$Macrohaplogroup))[,2]
freqTable$Freq2 = as.data.frame(table(full_fixed_haplogroups$Macrohaplogroup))[,2]
freqTable$Diff = as.data.frame((table(full_fixed_haplogroups$Macrohaplogroup) - table(fixed_haplogroups$Macrohaplogroup)))[,2]
freqTable$DiffPC = as.data.frame((table(full_fixed_haplogroups$Macrohaplogroup) - table(fixed_haplogroups$Macrohaplogroup)) / sum(table(full_fixed_haplogroups$Macrohaplogroup) - table(fixed_haplogroups$Macrohaplogroup)))[,2]
names(freqTable) = c("Haplogroup", "RefAln", "FullRefPanel", "FiltRefPanel", "Diff", "DiffPC")
freqTable

write.csv(freqTable, haplogroup_table_csv, quote = F, row.names = F)
