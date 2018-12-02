wgs = read.csv("/Users/u5015730/GitCode/MitoImputePrep/metadata/ADNI_wgs_Mitochondrial_Haplotypes.csv", header = T)
typ = read.table("/Users/u5015730/GitCode/MitoImputePrep/metadata/mito_snps_rcrs_Samples.txt", header = F)
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
