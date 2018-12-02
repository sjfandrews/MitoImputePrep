wgs = read.csv("~/GitCode/MitoImputePrep/metadata/ADNI_wgs_Mitochondrial_Haplotypes.csv", header = T)
typ = read.table("~/GitCode/MitoImputePrep/metadata/mito_snps_rcrs_Samples.txt", header = F)
names(typ) = "SAMPLE"

wgs$ID = NA
wgs$in.ADNI = NA

for (i in 1:nrow(wgs)) {
  x = unlist(strsplit(as.character(wgs$PATNO[i]), "_"))[3]
  wgs$ID[i] = x
  if (paste0("ADNI_", x) %in% typ$SAMPLE) {
    wgs$in.ADNI[i] = T
  } else {
    wgs$in.ADNI[i] = F
  }
}

inBoth = subset(wgs, wgs$in.ADNI == T)
write.csv(inBoth, "~/GitCode/MitoImputePrep/metadata/ADNI_samples_BOTH.csv", row.names = F)
write.csv(wgs, "~/GitCode/MitoImputePrep/metadata/ADNI_samples.csv", row.names = F)
