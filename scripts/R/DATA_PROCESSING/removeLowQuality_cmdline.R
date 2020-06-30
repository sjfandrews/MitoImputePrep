require(phangorn)

message("")
message("ARGUMENT HINTS!")
message("ARGUMENT 1:    .fasta FILE")
message("ARGUMENT 2:    OUTPUT HIGH QUALITY SEQUENCES FILE")
message("ARGUMENT 3:    MAXIMUM NUMBER OF MISSING ALLOWED")
message("ARGUMENT 4:    MAXIMUM NUMBER OF GAPS ALLOWED")
message("")

args = commandArgs(trailingOnly = TRUE) # Set arguments from the command line

fasta.file = args[1] # SPECIFY IN THE FASTA FILE
out.file = args[2] # SPECIFY THE OUTPUT FILE
missing.allowed = args[3]
gaps.allowed = args[4]

message("")
if (is.na(out.file) == TRUE | is.null(out.file) == TRUE) {
  message('OUTPUT FILE NOT SPECIFIED')
  out.file = paste0(sub("^([^.]*).*", "\\1", fasta.file), "_highQuality.txt")
  message(paste('OUTPUT FILE ASSIGNED TO ', out.file))
} else {
  message(paste('OUTPUT FILE ASSIGNED TO ', out.file))
}

if (is.na(missing.allowed) == TRUE | is.null(missing.allowed) == TRUE) {
  message("NUMBER OF MISSING SITES ALLOWED NOT SPECIFIED")
  message("DEFAULTING TO 5")
  missing.allowed = 5
} else {
  missing.allowed = as.numeric(missing.allowed)
}

if (is.na(gaps.allowed) == TRUE | is.null(gaps.allowed) == TRUE) {
  message("NUMBER OF GAPS ALLOWED NOT SPECIFIED")
  message("DEFAULTING TO 7")
  gaps.allowed = 7
} else {
  gaps.allowed = as.numeric(gaps.allowed)
}

aln = read.dna(fasta.file, format = 'fasta') # READ THE FASTA FILE
message(fasta.file, " SUCCESSFULLY PARSED")
message("")

Index = 1:nrow(aln) # BUILD AN INDEX FOR THE DATA FRAME
df = data.frame(Index) # BUILD DATA FRAME
df$seq = row.names(aln) # ASSIGN THE SEQUENCE LIST TO DATA FRAME

for (i in 1:nrow(aln)) {
  if ((i %% 1000) == 0) {
    cat(paste("CALCULATED BASE FREQUENCIES FOR SAMPLE ", i, "/", nrow(aln), "\r"))
  }
  tb = base.freq(aln[i,], freq = T, all = T)
  df$a[i] = as.numeric(tb[1])
  df$c[i] = as.numeric(tb[2])
  df$g[i] = as.numeric(tb[3])
  df$t[i] = as.numeric(tb[4])
  df$r[i] = as.numeric(tb[5])
  df$m[i] = as.numeric(tb[6])
  df$w[i] = as.numeric(tb[7])
  df$s[i] = as.numeric(tb[8])
  df$k[i] = as.numeric(tb[9])
  df$y[i] = as.numeric(tb[10])
  df$v[i] = as.numeric(tb[11])
  df$h[i] = as.numeric(tb[12])
  df$d[i] = as.numeric(tb[13])
  df$b[i] = as.numeric(tb[14])
  df$n[i] = as.numeric(tb[15])
  df$gap[i] = as.numeric(tb[16])
  df$missing[i] = as.numeric(tb[17])
}
#head(df)
#df.missing = subset(df, df$n > 5)
#df.gap = subset(df, df$gap > 7)

df.highQual = subset(df, df$n <= missing.allowed)
df.highQual = subset(df.highQual, df.highQual$gap <= gaps.allowed)
#df.lowQual = subset(df, df$n > 5)
#df.lowQual = subset(df.lowQual, df$gap > 7)

write.table(df.highQual$seq, out.file, sep = '\t', quote = F, col.names = F, row.names = F) 
message(out.file, " SUCCESSFULLY WRITTEN")
