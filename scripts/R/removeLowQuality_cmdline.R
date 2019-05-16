suppressPackageStartupMessages(require(phangorn))
suppressPackageStartupMessages(require(dplyr))

message("")
message("ARGUMENT HINTS!")
message("ARGUMENT 1:    .fasta FILE")
message("ARGUMENT 2:    OUTPUT HIGH QUALITY SEQUENCES FILE")
message("ARGUMENT 3:    MAXIMUM NUMBER OF MISSING ALLOWED")
message("ARGUMENT 4:    MAXIMUM NUMBER OF GAPS ALLOWED")
message("")

args <- commandArgs(trailingOnly = T) # Set arguments from the command line

fasta.file <- args[1] # SPECIFY IN THE FASTA FILE
out.file <- args[2] # SPECIFY THE OUTPUT FILE
missing.allowed <- args[3]
gaps.allowed <- args[4]

message("")
if (is.na(out.file) == TRUE | is.null(out.file) == TRUE) {
  message("OUTPUT FILE NOT SPECIFIED")
  out.file <- paste0(sub("^([^.]*).*", "\\1", fasta.file), "_highQuality.txt")
  message(paste("OUTPUT FILE ASSIGNED TO ", out.file))
} else {
  message(paste("OUTPUT FILE ASSIGNED TO ", out.file))
}

if (is.na(missing.allowed) == TRUE | is.null(missing.allowed) == TRUE) {
  message("NUMBER OF MISSING SITES ALLOWED NOT SPECIFIED")
  message("DEFAULTING TO 5")
  missing.allowed <- 5
} else {
  missing.allowed <- as.numeric(missing.allowed)
}

if (is.na(gaps.allowed) == TRUE | is.null(gaps.allowed) == TRUE) {
  message("NUMBER OF GAPS ALLOWED NOT SPECIFIED")
  message("DEFAULTING TO 7")
  gaps.allowed <- 7
} else {
  gaps.allowed <- as.numeric(gaps.allowed)
}

message("READING ", fasta.file)
aln <- read.dna(fasta.file, format = "fasta") # READ THE FASTA FILE
message(fasta.file, " SUCCESSFULLY PARSED")
message("")

dfr <- tibble(Index = 1:nrow(aln), # BUILD AN INDEX FOR THE DATA FRAME
             seq = row.names(aln)) # ASSIGN SEQUENCE LIST TO DATA FRAME

base_freq <- Vectorize(function (i) {
  base.freq(aln[i, ], freq = T, all = T)
})

message("COUNTING ALLELES")
freqs <- base_freq(dfr$Index) %>%
  t %>% # transpose so that columns are bases
  as_tibble %>% # turn into tibble
  rename(gap = "-", missing = "?") %>% # rename gap and missing
  bind_cols(dfr, .) # prepend sequence names

message("FILTERING AND WRITING OUTPUT")
freqs %>%
  filter(n <= missing.allowed & gap <= gaps.allowed) %>%
  pull(seq) %>% # get sequence names
  writeLines(con = out.file, sep = "\n") # output list of seq names

message(out.file, " SUCCESSFULLY WRITTEN")
