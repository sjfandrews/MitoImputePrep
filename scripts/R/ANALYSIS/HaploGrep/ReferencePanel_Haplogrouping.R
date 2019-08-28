library(dplyr)

haplogrep_file = "/Volumes/TimMcInerney/MitoImpute/data/VCF/McInerney_Master_Alignment_July18_2018_ambigANDgap2missing_tagFilled_filt_noInvar.txt"
haplogrep_csv = "/Volumes/TimMcInerney/MitoImpute/data/VCF/McInerney_Master_Alignment_July18_2018_ambigANDgap2missing_tagFilled_filt_noInvar.csv"

haplogroups = read.table(haplogrep_file, header = T, sep = "\t")

nonAFR = subset(haplogroups, substr(haplogroups$Haplogroup, 1, 1) != "L")
nonAFR$Macrohaplogroup = substr(nonAFR$Haplogroup, 1, 1)

AFR = subset(haplogroups, substr(haplogroups$Haplogroup, 1, 1) == "L")
AFR$Macrohaplogroup = substr(AFR$Haplogroup, 1, 2)

fixed_haplogroups = rbind(nonAFR, AFR)
fixed_haplogroups = arrange(fixed_haplogroups, fixed_haplogroups$SampleID)

write.csv(fixed_haplogroups, haplogrep_csv, quote = F, row.names = F)

full_haplogrep_file = "/Volumes/TimMcInerney/MitoImpute/data/VCF/McInerney_Master_Alignment_July18_2018_ambigANDgap2missing_tagFilled_noInvar_HaploGrep.txt"
full_haplogrep_csv = "/Volumes/TimMcInerney/MitoImpute/data/VCF/McInerney_Master_Alignment_July18_2018_ambigANDgap2missing_tagFilled_noInvar_HaploGrep.csv"

full_haplogroups = read.table(full_haplogrep_file, header = T, sep = "\t")
full_nonAFR = subset(full_haplogroups, substr(full_haplogroups$Haplogroup, 1, 1) != "L")
full_nonAFR$Macrohaplogroup = substr(full_nonAFR$Haplogroup, 1, 1)

full_AFR = subset(full_haplogroups, substr(full_haplogroups$Haplogroup, 1, 1) == "L")
full_AFR$Macrohaplogroup = substr(full_AFR$Haplogroup, 1, 2)

full_fixed_haplogroups = rbind(full_nonAFR, full_AFR)
full_fixed_haplogroups = arrange(full_fixed_haplogroups, full_fixed_haplogroups$SampleID)

write.csv(full_fixed_haplogroups, full_haplogrep_csv, quote = F, row.names = F)

# FILTERED HAPLOGROUPS TABLE
table(fixed_haplogroups$Macrohaplogroup)

# FULL HAPLOGROUPS TABLE
table(full_fixed_haplogroups$Macrohaplogroup)

# DIFFERENCE HAPLOGROUP TABLE
(table(full_fixed_haplogroups$Macrohaplogroup) - table(fixed_haplogroups$Macrohaplogroup)) / sum(table(full_fixed_haplogroups$Macrohaplogroup) - table(fixed_haplogroups$Macrohaplogroup)) * 100
nrow(full_fixed_haplogroups) - nrow(fixed_haplogroups)
