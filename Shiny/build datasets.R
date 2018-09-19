library(tidyverse)
library(pbapply)
library(HiMC); data(nodes)
library(taRifx)

##  Function
generate_snp_data_fixed <- function (map_file, ped_file) 
{
  map <- read.csv(map_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  header_row <- c("Family", "Individual", "Father", "Mother", 
                  "Sex", "Phenotype")
  snps = map[, 4]
  new_header = c(header_row, snps)
  ped <- read.csv(ped_file, sep = " ", header = FALSE, stringsAsFactors = FALSE, colClasses = 'character')
  range1 = seq(1, 6, by = 1)
  snp_data = data.frame(seq(1, nrow(ped), by = 1))
  for (i in range1) {
    snp_data[, i] = ped[, i]
  }
  range2 = seq(7, ncol(ped), by = 2)
  for (i in range2) {
    index = ((i - 7)/2) + 1
    snp_data[, index + 6] = paste(ped[, i], ped[, i + 1])
  }
  names(snp_data) <- new_header
  return(snp_data)
}


setwd("~/Dropbox/STRANDS")
setwd("~/Desktop/STRANDS")
setwd("~/Desktop/STRANDS_ref3")
##===============================##
##  WGS plink Files 
##===============================##
wgs.map <- '~/Dropbox/src/MitoImputePrep/DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.map'
wgs.ped <- '~/Dropbox/src/MitoImputePrep/DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.ped'
wgs.dat <- generate_snp_data_fixed(wgs.map, wgs.ped)  

# Assign haplogorups
MTwgs.classifications <- HiMC::getClassifications(wgs.dat) 
MTwgs.classifications <- as.tibble(MTwgs.classifications)
  
##===============================##
#  Typed only
##===============================##

# file names
typ.map <- list.files(path = "~/Dropbox/STRANDS", recursive = TRUE, pattern = "*b37.map")
typ.ped <- list.files(path = "~/Dropbox/STRANDS", recursive = TRUE, pattern = "*b37.ped")
typ.names <- typ.map %>% as.tibble() %>% separate(value, c('platform', 'file'), sep = '/')

typ.map <- list.files(path = "~/Desktop/STRANDS", recursive = TRUE, pattern = "*b37.map")
typ.ped <- list.files(path = "~/Desktop/STRANDS", recursive = TRUE, pattern = "*b37.ped")
typ.names <- typ.map %>% as.tibble() %>% separate(value, c('platform', 'reference', 'file'), sep = '/')

# read in files
typ.dat <- mapply(generate_snp_data_fixed, typ.map, typ.ped, SIMPLIFY = F)

# assign haplogroups
MTtyp.classifications <- pblapply(typ.dat, HiMC::getClassifications)
MTtyp.classifications <- lapply(MTtyp.classifications, as.tibble)

# Join Typed and WGS classifications 
MT_haps.out <- lapply(MTtyp.classifications, function(x){
  out <- x %>% 
    left_join(MTwgs.classifications, by = 'Individual', suffix = c("_typ", "_wgs")) %>% 
    as.tibble()
  out
})
names(MT_haps.out) <- typ.names$platform

MT_haps <- MT_haps.out[imp.names$platform] 

saveRDS(MT_haps, "~/Dropbox/Research/PostDoc-MSSM/3_mitoWAX/3_Scripts/ShinnyApp/MT_haps.rds")

##===============================##
##  Imputed plink files 
##===============================##

## file names
# 0.01 MAF Reference 
imp.map <- list.files(path = "~/Dropbox/STRANDS", recursive = TRUE, pattern = "*imputed.map")
imp.ped <- list.files(path = "~/Dropbox/STRANDS", recursive = TRUE, pattern = "*imputed.ped")
imp.names <- imp.map %>% as.tibble() %>% separate(value, c('platform', 'file'), sep = '/')

# 0.005 MAF Reference 
imp.map <- list.files(path = "~/Desktop/STRANDS", recursive = TRUE, pattern = "*imputed.map")
imp.ped <- list.files(path = "~/Desktop/STRANDS", recursive = TRUE, pattern = "*imputed.ped")
imp.names <- imp.map %>% as.tibble() %>% separate(value, c('platform', 'reference', 'file'), sep = '/')

# 0.001 MAF Reference 
imp.map <- list.files(path = "~/Desktop/STRANDS_ref3", recursive = TRUE, pattern = "*.map")
imp.map <- grep('Imputed_', imp.map, value = T)
imp.ped <- list.files(path = "~/Desktop/STRANDS_ref3", recursive = TRUE, pattern = "*.ped")
imp.ped <- grep('Imputed_', imp.ped, value = T)
imp.names <- imp.map %>% as.tibble() %>% separate(value, c('platform', 'reference', 'file'), sep = '/')

# read in files
imp.dat <- mapply(generate_snp_data_fixed, imp.map, imp.ped, SIMPLIFY = F)
imp.dat <- lapply(imp.dat, function(x){
  out <- x[,-c(grep("\\<57\\>", colnames(x)),
               grep("\\<62\\>", colnames(x)),
               grep("\\<63\\>", colnames(x)),
               grep("\\<72\\>", colnames(x)),
               grep("\\<185\\>", colnames(x)),
               grep("\\<186\\>", colnames(x)),
               grep("\\<189\\>", colnames(x)),
               grep("\\<195\\>", colnames(x)),
               grep("\\<228\\>", colnames(x)),
               grep("\\<295\\>", colnames(x)),
               grep("\\<723\\>", colnames(x)),
               grep("\\<750\\>", colnames(x)),
               grep("\\<930\\>", colnames(x)),
               grep("\\<961\\>", colnames(x)),
               grep("\\<1692\\>", colnames(x)),
               grep("\\<2831\\>", colnames(x)),
               grep("\\<3200\\>", colnames(x)),
               grep("\\<3552\\>", colnames(x)),
               grep("\\<3796\\>", colnames(x)),
               grep("\\<3921\\>", colnames(x)),
               grep("\\<4454\\>", colnames(x)),
               grep("\\<4562\\>", colnames(x)),
               grep("\\<4769\\>", colnames(x)),
               grep("\\<5774\\>", colnames(x)),
               grep("\\<5894\\>", colnames(x)),
               grep("\\<6221\\>", colnames(x)),
               grep("\\<7196\\>", colnames(x)),
               grep("\\<7624\\>", colnames(x)),
               grep("\\<8014\\>", colnames(x)),
               grep("\\<8080\\>", colnames(x)),
               grep("\\<8860\\>", colnames(x)),
               grep("\\<9824\\>", colnames(x)),
               grep("\\<10097\\>", colnames(x)),
               grep("\\<10410\\>", colnames(x)),
               grep("\\<10754\\>", colnames(x)),
               grep("\\<12633\\>", colnames(x)),
               grep("\\<12738\\>", colnames(x)),
               grep("\\<12930\\>", colnames(x)),
               grep("\\<12950\\>", colnames(x)),
               grep("\\<14470\\>", colnames(x)),
               grep("\\<15884\\>", colnames(x)),
               grep("\\<15954\\>", colnames(x)),
               grep("\\<16111\\>", colnames(x)),
               grep("\\<16114\\>", colnames(x)),
               grep("\\<16129\\>", colnames(x)),
               grep("\\<16147\\>", colnames(x)),
               grep("\\<16166\\>", colnames(x)),
               grep("\\<16176\\>", colnames(x)),
               grep("\\<16182\\>", colnames(x)),
               grep("\\<16183\\>", colnames(x)),
               grep("\\<16188\\>", colnames(x)),
               grep("\\<16232\\>", colnames(x)),
               grep("\\<16240\\>", colnames(x)),
               grep("\\<16257\\>", colnames(x)),
               grep("\\<16258\\>", colnames(x)),
               grep("\\<16265\\>", colnames(x)),
               grep("\\<16266\\>", colnames(x)),
               grep("\\<16286\\>", colnames(x)),
               grep("\\<16291\\>", colnames(x)),
               grep("\\<16293\\>", colnames(x)),
               grep("\\<16318\\>", colnames(x)),
               grep("\\<16327\\>", colnames(x)))]
  out
})
imp.dat <- lapply(test, as.tibble)
names(imp.dat) <- imp.names$platform

imp.dat01 <- imp.dat
saveRDS(imp.dat, "~/Dropbox/Research/PostDoc-MSSM/3_mitoWAX/3_Scripts/ShinnyApp/imp.dat.rds")

imp.dat005 <- imp.dat
saveRDS(imp.dat005, "~/Dropbox/Research/PostDoc-MSSM/3_mitoWAX/3_Scripts/ShinnyApp/imp.dat005.rds")

imp.dat001 <- imp.dat
saveRDS(imp.dat001, "~/Dropbox/Research/PostDoc-MSSM/3_mitoWAX/3_Scripts/ShinnyApp/imp.dat001.rds")

##===============================##
##  Info Score Files 
##===============================##
HiMC <- tibble(
  himc = 'yes',
  pos = as.numeric(c('10115', '1018', '10398', '10400', '10550', '11177', '11251', '11719', '11947', '12007', '12308', '12414', '12705', '13263', '13368', '13506', '13708', '13789', '14178', '14318', '1438', '14470', '14560', '14668', '14766', '14905', '15043', '15326', '15452', '15535', '16111', '16189', '16271', '16362', '16390', '16391', '16391', '1719', '1736', '2092', '3505', '3552', '3594', '4580', '4769', '4883', '4917', '4977', '5178', '5442', '6371', '7028', '825', '8251', '8414', '8468', '8703', '9042', '9055', '9347', '9950')),
  Haplogroup = c("L2", "L3", "K1", "M", "K", "B2", "JT", "R0", "W", "A2", "U", "N2", "R", "C", "T", "L2'3'4'6", "J", "L1", "L1", "C", "H2", "K", ".", "D4", "HV", "T", "N1a1b", "H2a2a", "JT", "B4b'd'e", "A2", "T1", "JT", "L4", "L2", "I", "I", "X2", "A", "D1", "W", "C", "L3'4", "V", "H2a", "M80'D", "T", "B2", "D ", "L0", "X ", "H", "L2'3'4'6", "N1a1b", "D4", "L2'3'4'6", "D2", "L0", "U8b", "L0", "B2"))

# 0.01 MAF Reference 
info <- list.files(path = "~/Dropbox/STRANDS", recursive = TRUE, pattern = "*_info")
info.names <- info %>% as.tibble() %>% separate(value, c('platform', 'file'), sep = '/')

# 0.005 MAF Reference 
info <- list.files(path = "~/Desktop/STRANDS", recursive = TRUE, pattern = "*_imputed_info")
info <- info[!grepl('by_sample', info)]

# 0.001 MAF Reference 
info <- list.files(path = "~/Desktop/STRANDS_ref3", recursive = TRUE, pattern = "*_imputed_info")
info <- info[!grepl('by_sample', info)]

info.names <- info %>% as.tibble() %>% 
  filter(!grepl('by_sample', value)) %>% 
  separate(value, c('platform', 'reference', 'file'), sep = '/')

info.dat <- lapply(info, read_delim, delim = " ")
names(info.dat) <- info.names$platform

imp.info <- lapply(info.dat, function(x){
  out <- x %>% mutate(info_comb = ifelse(info_type0 == -1, info,info_type0 )) %>% 
    left_join(HiMC, by = c('position' = 'pos')) %>%
    mutate(himc = ifelse(is.na(himc), 'no', himc))
  out
})

imp.info01 <- imp.info
saveRDS(imp.info01, "~/Dropbox/Research/PostDoc-MSSM/3_mitoWAX/3_Scripts/ShinnyApp/imp.info01.rds")

imp.info005 <- imp.info 
saveRDS(imp.info005, "~/Dropbox/Research/PostDoc-MSSM/3_mitoWAX/3_Scripts/ShinnyApp/imp.info005.rds")

imp.info001 <- imp.info 
saveRDS(imp.info001, "~/Dropbox/Research/PostDoc-MSSM/3_mitoWAX/3_Scripts/ShinnyApp/imp.info001.rds")


##===============================##
##   
##===============================##















