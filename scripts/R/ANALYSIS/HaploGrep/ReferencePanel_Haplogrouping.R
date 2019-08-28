library(dplyr)

haplogrep_file = "/Volumes/TimMcInerney/MitoImpute/data/VCF/McInerney_Master_Alignment_July18_2018_ambigANDgap2missing_tagFilled_filt_noInvar.txt"

haplogroups = read.table(haplogrep_file, header = T)

nonAFR = subset(haplogroups, substr(haplogroups$Haplogroup, 1, 1) != "L")
nonAFR$Macrohaplogroup = substr(nonAFR$Haplogroup, 1, 1)

AFR = subset(haplogroups, substr(haplogroups$Haplogroup, 1, 1) == "L")
AFR$Macrohaplogroup = substr(AFR$Haplogroup, 1, 2)

fixed_haplogroups = rbind(nonAFR, AFR)
fixed_haplogroups = arrange(fixed_haplogroups, fixed_haplogroups$SampleID)

write.csv(fixed_haplogroups, haplogrep_file, quote = F, row.names = F)

full_haplogrep_file = "/Volumes/TimMcInerney/MitoImpute/data/VCF/McInerney_Master_Alignment_July18_2018_ambigANDgap2missing_tagFilled_filt_noInvar.txt"
