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

##===============================##
##  Read in Typed and WGS haplogroups 
##===============================##

mt.haps <- readRDS('~/Dropbox/src/MitoImputePrep/Shiny/MT_haps.rds')

wd.path <- '~/Dropbox/src/MitoImputePrep/data/STRANDS/'

##===============================##
##  Read in reference 
##===============================##

dat <- read_csv('~/Dropbox/src/MitoImputePrep/metadata/ConcordanceTables_MAF1pc.csv') 
dat2 <- read_csv('~/Dropbox/src/MitoImputePrep/metadata/ConcordanceTables_MAF0-1pc.csv')

out <- left_join(dat, dat2, by = 'array', suffix = c('.1pc', '.01pc')) %>% 
  mutate(improvment = Imputed.hg.Conc.1pc - Typed.hg.Conc.1pc)

select(out, array, TOTAL.1pc, TOTAL.01pc, Imputed.hg.Conc.1pc, Imputed.hg.Conc.01pc, Typed.hg.Conc.1pc, improvment) %>% 
  arrange(improvment) %>% 
  summarise(mean.improvment = mean(improvment, na.rm = T), 
            mean.typed = mean(Typed.hg.Conc.1pc, na.rm = T), 
            mean.imputed = mean(Imputed.hg.Conc.1pc, na.rm = T))

ggplot(out, aes(y = Imputed.hg.Conc.1pc, x = Typed.hg.Conc.1pc)) + 
  geom_point() + geom_abline(intercept = 0, slope = 1, colour = 'red', linetype = 2) +
  xlim(0,1) + ylim(0,1) + theme_bw() + theme(aspect.ratio=1) + labs(x = "Haplogroup Concordance \n Typed vs WGS", y = "Haplogroup Concordance \n Typed + Imputed vs WGS")

ggplot(out, aes(x = Imputed.hg.Conc.1pc, y = Imputed.hg.Conc.01pc)) + 
  geom_point() + geom_abline(intercept = 0, slope = 1, colour = 'red', linetype = 2) +
  xlim(0,1) + ylim(0,1) + theme_bw() + theme(aspect.ratio=1) + 
  labs(x = "Haplogroup Concordance \n MAF > 1%", y = "Haplogroup Concordance \n MAF > 0.1%")


##===============================##
##  Imputed plink files 
##===============================##

## file names
# MCMC 1 Imputation 
imp_MCMC1.map <- list.files(path = wd.path, recursive = TRUE, pattern = "*imputed_MCMC1.map")
imp_MCMC1.ped <- list.files(path = wd.path, recursive = TRUE, pattern = "*imputed_MCMC1.ped")
imp_MCMC1.names <- imp_MCMC1.map %>% as.tibble() %>% separate(value, c('platform', 'experiment', 'MCMC', 'file'), sep = '/')
imp.dat <- mapply(generate_snp_data_fixed, paste0(wd.path, imp_MCMC1.map), paste0(wd.path, imp_MCMC1.ped), SIMPLIFY = F)

# MCMC 5 Imputation 
imp_MCMC5.map <- list.files(path = wd.path, recursive = TRUE, pattern = "*imputed_MCMC5.map")
imp_MCMC5.ped <- list.files(path = wd.path, recursive = TRUE, pattern = "*imputed_MCMC5.ped")
imp_MCMC5.names <- imp_MCMC5.map %>% as.tibble() %>% separate(value, c('platform', 'experiment', 'MCMC', 'file'), sep = '/')
imp.dat <- mapply(generate_snp_data_fixed, paste0(wd.path, imp_MCMC5.map), paste0(wd.path, imp_MCMC5.ped), SIMPLIFY = F)
names(imp.dat) <- imp_MCMC5.names$platform

# MCMC 10 Imputation 
imp_MCMC10.map <- list.files(path = wd.path, recursive = TRUE, pattern = "*imputed_MCMC10.map")
imp_MCMC10.ped <- list.files(path = wd.path, recursive = TRUE, pattern = "*imputed_MCMC10.ped")
imp_MCMC10.names <- imp_MCMC10.map %>% as.tibble() %>% separate(value, c('platform', 'experiment', 'MCMC', 'file'), sep = '/')
imp.dat <- mapply(generate_snp_data_fixed, paste0(wd.path, imp_MCMC10.map), paste0(wd.path, imp_MCMC10.ped), SIMPLIFY = F)
names(imp.dat) <- imp_MCMC10.names$platform

# MCMC 20 Imputation 
imp_MCMC20.map <- list.files(path = wd.path, recursive = TRUE, pattern = "*imputed_MCMC20.map")
imp_MCMC20.ped <- list.files(path = wd.path, recursive = TRUE, pattern = "*imputed_MCMC20.ped")
imp_MCMC20.names <- imp_MCMC20.map %>% as.tibble() %>% separate(value, c('platform', 'experiment', 'MCMC', 'file'), sep = '/')
imp.dat <- mapply(generate_snp_data_fixed, paste0(wd.path, imp_MCMC20.map), paste0(wd.path, imp_MCMC20.ped), SIMPLIFY = F)
names(imp.dat) <- imp_MCMC20.names$platform

# read in files

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
imp.dat <- lapply(imp.dat, as.tibble)


##===============================##
##  Info Score Files 
##===============================##
HiMC <- tibble(
  himc = 'yes',
  pos = as.numeric(c('10115', '1018', '10398', '10400', '10550', '11177', '11251', '11719', '11947', '12007', '12308', '12414', '12705', '13263', '13368', '13506', '13708', '13789', '14178', '14318', '1438', '14470', '14560', '14668', '14766', '14905', '15043', '15326', '15452', '15535', '16111', '16189', '16271', '16362', '16390', '16391', '16391', '1719', '1736', '2092', '3505', '3552', '3594', '4580', '4769', '4883', '4917', '4977', '5178', '5442', '6371', '7028', '825', '8251', '8414', '8468', '8703', '9042', '9055', '9347', '9950')),
  Haplogroup = c("L2", "L3", "K1", "M", "K", "B2", "JT", "R0", "W", "A2", "U", "N2", "R", "C", "T", "L2'3'4'6", "J", "L1", "L1", "C", "H2", "K", ".", "D4", "HV", "T", "N1a1b", "H2a2a", "JT", "B4b'd'e", "A2", "T1", "JT", "L4", "L2", "I", "I", "X2", "A", "D1", "W", "C", "L3'4", "V", "H2a", "M80'D", "T", "B2", "D ", "L0", "X ", "H", "L2'3'4'6", "N1a1b", "D4", "L2'3'4'6", "D2", "L0", "U8b", "L0", "B2"))

# MCMC 5 Imputation 
MCMC5_info <- list.files(path = wd.path, recursive = TRUE, pattern = "*MCMC5_info") 
MCMC5_info <- MCMC5_info[grepl("sample", MCMC5_info) == F]
MCMC5_info.names <- MCMC5_info %>% as.tibble() %>% separate(value, c('platform', 'experiment', 'MCMC', 'file'), sep = '/')
info.dat <- lapply(paste0(wd.path, MCMC5_info), read_delim, delim = " ")
names(info.dat) <- MCMC5_info.names$platform

# MCMC 10 Imputation 
MCMC10_info <- list.files(path = wd.path, recursive = TRUE, pattern = "*MCMC10_info") 
MCMC10_info <- MCMC10_info[grepl("sample", MCMC10_info) == F]
MCMC10_info.names <- MCMC10_info %>% as.tibble() %>% separate(value, c('platform', 'experiment', 'MCMC', 'file'), sep = '/')
info.dat <- lapply(paste0(wd.path, MCMC10_info), read_delim, delim = " ")
names(info.dat) <- MCMC10_info.names$platform

# MCMC 20 Imputation 
MCMC20_info <- list.files(path = wd.path, recursive = TRUE, pattern = "*MCMC20_info") 
MCMC20_info <- MCMC20_info[grepl("sample", MCMC20_info) == F]
MCMC20_info.names <- MCMC20_info %>% as.tibble() %>% separate(value, c('platform', 'experiment', 'MCMC', 'file'), sep = '/')
info.dat <- lapply(paste0(wd.path, MCMC20_info), read_delim, delim = " ")
names(info.dat) <- MCMC20_info.names$platform


info.dat <- lapply(info.dat, function(x){
  out <- x %>% mutate(info_comb = ifelse(info_type0 == -1, info,info_type0 )) %>% 
    left_join(HiMC, by = c('position' = 'pos')) %>%
    mutate(himc = ifelse(is.na(himc), 'no', himc))
  out
})



##===============================##
##  Haplogroup assocations 
##===============================##
hapclassify <- function(x, y){
  ## Filter SNPs
  rm.info <- filter(y, info >= 0.3)
  imp.dat_filt <- x[ ,colnames(x) %in% c('Individual', rm.info$position)]
  
  ## Assign haplogroups
  MTimp.classifications <- HiMC::getClassifications(as.data.frame(imp.dat_filt))
  MTimp.classifications
}

imp.concordance <- function(x){
  out <- x %>%
    count(haplogroup, haplogroup_wgs) %>% 
    mutate(perc = (n/sum(n))) %>% 
    mutate(match = haplogroup == haplogroup_wgs) %>% 
    group_by(match) %>%
    summarise(sum = round(sum(perc), 2)) 
}

imp.haps <- mapply(hapclassify, imp.dat, imp.info, SIMPLIFY = FALSE)
imp.haps <- lapply(imp.haps, as.tibble)
hap.dat <- mapply(left_join,  mt.haps, imp.haps, SIMPLIFY = FALSE)

mcmc_conc20 <- map(hap.dat, imp.concordance) %>% 
  bind_rows() %>% 
  filter(match == TRUE) %>% 
  mutate(array = names(hap.dat)) %>%
  rename(mcmc20.hg.conc = sum) %>% 
  select(array, mcmc20.hg.conc) %>% 
  right_join(dat, by = 'array') %>%
  mutate(improvement = mcmc20.hg.conc - Imputed.hg.Conc) %>% 
  arrange(improvement) %>% 
  print(n = Inf)
  
ggplot(mcmc_conc20, aes(y =  mcmc20.hg.conc, x = Imputed.hg.Conc)) + 
  geom_point() + geom_abline(intercept = 0, slope = 1, colour = 'red', linetype = 2) +
  xlim(0,1) + ylim(0,1) + theme_bw() + theme(aspect.ratio=1) + 
  labs(x = "Haplogroup Concordance \n MCMC 30", y = "Haplogroup Concordance \n MCMC 20")
  
  
ggplot(mcmc_conc, aes(y = mcmc1.hg.conc, x = Typed.hg.Conc)) + 
  geom_point() + geom_abline(intercept = 0, slope = 1, colour = 'red', linetype = 2) +
  xlim(0,1) + ylim(0,1) + theme_bw() + theme(aspect.ratio=1) + 
  labs(x = "Haplogroup Concordance \n Typed vs WGS", y = "Haplogroup Concordance \n Typed + Imputed vs WGS")




mcmc5_conc
mcmc_conc10





























