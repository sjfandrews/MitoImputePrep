generate_snp_data <- function (map_file, ped_file) {
  map <- read.csv(map_file, sep = "\t", header = F, stringsAsFactors = F)
  ped <- read.csv(ped_file, sep = " ", header = F,
                  stringsAsFactors = F, colClasses = "character")
  infocols <- ped[, 1:6]

  #start from the first genotype and iterate by two
  A_dex <- seq(7, ncol(ped), by = 2)
  A1 <- as.matrix(ped[, A_dex])
  A2 <- as.matrix(ped[, A_dex + 1]); rm(ped)
  #Join A1 and A2 into a single column
  A1_A2 <- matrix(paste(A1, A2), nrow = nrow(A1)); rm(A1, A2)
  noprint <- gc()
  #Concatenate the info and genotype columns
  out <- cbind(infocols, as.data.frame(A1_A2, stringsAsFactors = F)); rm(A1_A2)
  #make header with ped columns and snp locations from map
  names(out) <- c("Family", "Individual", "Father", "Mother",
                  "Sex", "Phenotype", map[, 4])
  return(out)
}

generate_snp_data_lowmem <- function (map_file, ped_file) {
  map <- read.csv(map_file, sep = "\t", header = F, stringsAsFactors = F)
  ped <- read.csv(ped_file, sep = " ", header = F,
                  stringsAsFactors = F, colClasses = "character")
  snp_data <- data.frame(seq(1, nrow(ped), by = 1))
  for (i in seq(1, 6, by = 1)) {
    snp_data[, i] <- ped[, i]
  }
  for (i in seq(7, ncol(ped), by = 2)) {
    index <- ( (i - 7) / 2) + 1
    snp_data[, index + 6] <- paste(ped[, i], ped[, i + 1])
  }
  #make header with ped columns and snp locations from map
  new_header <- c("Family", "Individual", "Father", "Mother",
                  "Sex", "Phenotype", map[, 4])
  names(snp_data) <- new_header
  return(snp_data)
}
