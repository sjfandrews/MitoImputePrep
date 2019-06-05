message("")
message("ARGUMENT HINTS!")
message("ARGUMENT 1:    VCF SNP POSITIONS FILE")
message("")

args = commandArgs(trailingOnly = TRUE) # Set arguments from the command line

snp.file = args[1] # Input EBSP .log.csv file

snps = read.table(snp.file)

df = data.frame(matrix(nrow = nrow(snps), ncol = 2))
df$X1 = "MT"
df$X2 = snps$V1
write.table(df, snp.file, col.names = F, row.names = F, quote = F, sep = "\t")