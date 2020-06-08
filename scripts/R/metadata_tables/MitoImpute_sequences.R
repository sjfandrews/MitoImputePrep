full = read.csv("~/GitCode/MitoImputePrep/metadata/MitoImpute_sequences/MitoImpute_sequences_full_2018-07-18.csv", header = F)
QC = read.csv("~/GitCode/MitoImputePrep/metadata/MitoImpute_sequences/MitoImpute_sequences_QC_2018-07-18.csv", header = F)

table(full$V1 %in% QC$V1)
full$V2 = full$V1 %in% QC$V1
names(full) = c("Sequence", "QC_Pass")

write.csv(full, "~/GitCode/MitoImputePrep/metadata/MitoImpute_sequences_full-with-pass_2018-07-18.csv", row.names = F, quote = F)
