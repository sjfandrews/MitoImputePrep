suppressPackageStartupMessages(library(tidyverse))
library(ggforce)
library(readxl)
library(HiMC); data(nodes)

# PACKAGES THAT WILL ACTUALLY BE NEEDED
# dplyr
# readr
# purrr

mccr <- function (act, pred) 
{
  TP <- sum(act %in% 1 & pred %in% 1)
  TN <- sum(act %in% 0 & pred %in% 0)
  FP <- sum(act %in% 0 & pred %in% 1)
  FN <- sum(act %in% 1 & pred %in% 0)
  denom <- as.double(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)
  if (any((TP + FP) == 0, (TP + FN) == 0, (TN + FP) == 0, (TN + FN) == 0)) 
    denom <- 1
  mcc <- ((TP * TN) - (FP * FN))/sqrt(denom)
  return(mcc)
}

##  Read in Thousand Genomes WGS 
#     Plink .ped files 
#     VCF Files
wgs.ped = "/Volumes/TimMcInerney/MitoImpute/data/PLINK/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.ped"
wgs.map = "/Volumes/TimMcInerney/MitoImpute/data/PLINK/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.map"
wgs.vcf = "/Volumes/TimMcInerney/MitoImpute/data/VCF/1kGP_chrMT_SNPonly.vcf.gz"

wgs_1kg.ped <- generate_snp_data(wgs.map, wgs.ped)
wgs_1kg.vcf <- read_tsv(wgs.vcf, comment = '##', na = c(".", "", "NA"))

##  Read in Thousand Genomes WGS - Typed only
#     Plink .ped files 
typ.map = "/Volumes/TimMcInerney/MitoImpute/data/HAPLOGROUPS/chips/BDCHP-1X10-HUMANHAP240S_11216501_A-b37/MCMC_Experiments/MCMC1/chrMT_1kg_BDCHP-1X10-HUMANHAP240S_11216501_A-b37_imputed_MCMC1.map"
typ.ped = "/Volumes/TimMcInerney/MitoImpute/data/HAPLOGROUPS/chips/BDCHP-1X10-HUMANHAP240S_11216501_A-b37/MCMC_Experiments/MCMC1/chrMT_1kg_BDCHP-1X10-HUMANHAP240S_11216501_A-b37_imputed_MCMC1.ped"
typ.vcf = "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/BDCHP-1X10-HUMANHAP240S_11216501_A-b37/ReferencePanel_v2/chrMT_1kg_BDCHP-1X10-HUMANHAP240S_11216501_A-b37.vcf.gz"

typ_1kg.ped <- generate_snp_data(typ.map, typ.ped)
typ_1kg.vcf <- read_tsv(typ.vcf, comment = '##', na = c(".", "", "NA"))

##  Read in Thousand Genomes WGS - Imputed only
#     Plink .ped files 
imp.map = "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/BDCHP-1X10-HUMANHAP240S_11216501_A-b37/ReferencePanel_v2/MCMC1/chrMT_1kg_BDCHP-1X10-HUMANHAP240S_11216501_A-b37_imputed_MCMC1.map"
imp.ped = "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/BDCHP-1X10-HUMANHAP240S_11216501_A-b37/ReferencePanel_v2/MCMC1/chrMT_1kg_BDCHP-1X10-HUMANHAP240S_11216501_A-b37_imputed_MCMC1.ped"
imp.vcf = "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/BDCHP-1X10-HUMANHAP240S_11216501_A-b37/ReferencePanel_v2/MCMC1/chrMT_1kg_BDCHP-1X10-HUMANHAP240S_11216501_A-b37_imputed_MCMC1_haplogrep.vcf.gz"
imp.info = "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/BDCHP-1X10-HUMANHAP240S_11216501_A-b37/ReferencePanel_v2/MCMC1/chrMT_1kg_BDCHP-1X10-HUMANHAP240S_11216501_A-b37_imputed_MCMC1_info"

imp_1kg.ped <- generate_snp_data(imp.map, imp.ped) 
imp_1kg.ped <- imp_1kg.ped[,-c(grep("\\<189\\>", colnames(imp_1kg.ped)), grep("\\<16183\\>", colnames(imp_1kg.ped)))]
imp_1kg.vcf <- read_tsv(imp.vcf, comment = '##', na = c(".", "", "NA"))
imp_1kg.info <- read_delim(imp.info, delim = " ")

## Data Munging 
##  obtain intersect of SNPs from imputed and wgs 
snp.intersect <- wgs_1kg.vcf %>% 
  select(POS) %>% 
  semi_join(imp_1kg.vcf, by = 'POS') %>%
  mutate(POS = paste0('mt', POS))

##  obtain intersect of SNPs from imputed and wgs 
# THIS IS FOR TYPED, BUT WE DONT NEED IT
snp.intersect.typ <- wgs_1kg.vcf %>% 
  select(POS) %>% 
  semi_join(typ_1kg.vcf, by = 'POS') %>%
  mutate(POS = paste0('mt', POS))

##  Munge imputed dataframe
# Drop all columns but position and genotypes
imp_1kg <- imp_1kg.vcf %>% 
  select(-`#CHROM`, -ID, -REF, -ALT, -QUAL, -FILTER, -INFO, -FORMAT) %>% 
  distinct(POS, .keep_all = T) %>% 
  mutate(POS = paste0('mt', POS))

# Transpose dataframe
imp_1kg <- imp_1kg %>%
  gather(key = var_name, value = value, 2:ncol(imp_1kg)) %>% 
  spread_(key = names(imp_1kg)[1],value = 'value')

# substitute allels for NA, 0, 1 calls, change formate to interger
imp_1kg <- as.tibble(sapply(imp_1kg, function(x) gsub('\\:.*', "", x)))
imp_1kg[imp_1kg == './.'] <- NA
imp_1kg[imp_1kg == '.'] <- NA
imp_1kg[imp_1kg == '0/0'] <- 0
imp_1kg[imp_1kg == '0'] <- 0
imp_1kg[imp_1kg == '1/1'] <- 1
imp_1kg[imp_1kg == '1'] <- 1
imp_1kg <- imp_1kg %>% 
  select(var_name, snp.intersect$POS) %>% 
  mutate_at(vars(contains('mt')), as.integer)

##  Munge ADNI wgs dataframe
# Drop all columns but position and genotypes
wgs_1kg <- wgs_1kg.vcf %>% 
  select(-`#CHROM`, -ID, -REF, -ALT, -QUAL, -FILTER, -INFO, -FORMAT) %>% 
  distinct(POS, .keep_all = T) %>% 
  mutate(POS = paste0('mt', POS))

# Transpose dataframe
wgs_1kg <- wgs_1kg %>%
  gather(key = var_name, value = value, 2:ncol(wgs_1kg)) %>% 
  spread_(key = names(wgs_1kg)[1],value = 'value')

# substitute allels for NA, 0, 1 calls, change formate to interger
wgs_1kg[wgs_1kg == './.'] <- NA
wgs_1kg[wgs_1kg == '.'] <- NA
wgs_1kg[wgs_1kg == '0/0'] <- 0
wgs_1kg[wgs_1kg == '0'] <- 0
wgs_1kg[wgs_1kg == '1/1'] <- 1
wgs_1kg[wgs_1kg == '1'] <- 1
wgs_1kg <- wgs_1kg %>% 
  select(var_name, snp.intersect$POS) %>% 
  mutate_at(vars(contains('mt')), as.integer) 

##  Munge ADNI typ dataframe
# Drop all columns but position and genotypes
typ_1kg <- typ_1kg.vcf %>% 
  select(-`#CHROM`, -ID, -REF, -ALT, -QUAL, -FILTER, -INFO, -FORMAT) %>% 
  distinct(POS, .keep_all = T) %>% 
  mutate(POS = paste0('mt', POS))

# Transpose dataframe
typ_1kg <- typ_1kg %>%
  gather(key = var_name, value = value, 2:ncol(typ_1kg)) %>% 
  spread_(key = names(typ_1kg)[1],value = 'value')

# substitute allels for NA, 0, 1 calls, change formate to interger
typ_1kg[typ_1kg == './.'] <- NA
typ_1kg[typ_1kg == '.'] <- NA
typ_1kg[typ_1kg == '0/0'] <- 0
typ_1kg[typ_1kg == '0'] <- 0
typ_1kg[typ_1kg == '1/1'] <- 1
typ_1kg[typ_1kg == '1'] <- 1
typ_1kg <- typ_1kg %>% 
  select(var_name, snp.intersect.typ$POS) %>% 
  mutate_at(vars(contains('mt')), as.integer) 

##
## JUST USED FOR HiMC, DONT NEED FOR MCC
imp_1kg.info <- mutate(imp_1kg.info, info_comb = ifelse(info_type0 == -1, info,info_type0 ))
imp_1kg.info <- mutate(imp_1kg.info, himc = ifelse(position %in% c(825, 1018, 1438, 1719, 1736, 2092, 3505, 3552, 3594, 4580, 4769, 4917, 4977, 5178, 5442, 6371, 7028, 8251, 8414, 8468, 8703, 9042, 9055, 9347, 9950, 10115, 10398, 10398, 10400, 10550, 11177, 11251, 11947, 12007, 12308, 12705, 13263, 13368, 13506, 13708, 13789, 14178, 14318, 14470, 14560, 14668, 14766, 15043, 15326, 15452, 15535, 16111, 16189, 16391), 'yes', 'no'))

## SNP Concordance Statisitics
# NOTE: SHOULD ONLY REALLY BE LOOKING AT TYPE 1 and 2 SNPS (in REF ONLY and SAMPLE ONLY)
# TYPE 0 SHOULD HAVE A MCC = 1
# TYPE 2 MAY BE PERFECT, SO LOOK AT THIS
# ONLY WANT TO LOOK AT MCC FOR IMPUTED SNPs

## calculate MCC between imputed and typed SNPs
mccr.geno <- unlist(map2(wgs_1kg[,2:ncol(wgs_1kg)], imp_1kg[,2:ncol(imp_1kg)], function(a,b) mccr(a,b)))
mccr.geno.typ <- unlist(map2(wgs_1kg[,2:ncol(wgs_1kg)], typ_1kg[,2:ncol(typ_1kg)], function(a,b) mccr(a,b)))

##  Calculate concordance between imputed and typed SNPs
concodance.geno <- unlist(map2(imp_1kg[,2:ncol(imp_1kg)], wgs_1kg[,2:ncol(wgs_1kg)], function(a, b) sum(a == b, na.rm = T)/length(a)))

##  Calculate allele frequency of typed SNPs from WGS data
af <- unlist(map(wgs_1kg[,2:ncol(wgs_1kg)], function(a) sum(a, na.rm = T)/length(a)))

##  Dataframe for summary stats
summary.stats <- tibble(mtSNP = colnames(wgs_1kg[,2:ncol(wgs_1kg)]), 
                        pos = as.integer(gsub('mt', '', colnames(wgs_1kg[,2:ncol(wgs_1kg)]))),
                        af = af,
                        mcc = mccr.geno,
                        concodance = concodance.geno)

##  merge on info.score file
summary.stats <- left_join(summary.stats, select(imp_1kg.info, c(-snp_id, -rs_id)), by = c('pos' = 'position')) 
summary.stats <- mutate(summary.stats, info.cat = cut_width(summary.stats$info, 0.25, boundary = 0))
##  basic summary stats
summary.stats %>% count(af > 0.01)
summary.stats %>% count(mcc > 0.4)
summary.stats %>% count(concodance < 0.9)
summary.stats %>% count(info > 0.3); summary.stats %>% count(info > 0.5)

##  Plots by bp position - all SNPs
# MCC
ggplot(summary.stats, aes(x = pos, y = mcc, colour = as.factor(type), size = af)) + geom_point() + theme_bw() + 
  labs(x = 'mtDNA position', y = 'MCC', title = 'Imputed mtSNP allele vs sequenced mtSNP allele') + 
  guides(colour=guide_legend(title="Impute2 Info")) + 
  guides(size=guide_legend(title="Allele Frequency"))

# CONCORDANCE
ggplot(summary.stats, aes(x = pos, y = concodance, colour = info.cat, size = af)) + geom_point() + theme_bw() + 
  labs(x = 'mtDNA position', y = 'Concordance', title = 'Imputed mtSNP allele vs sequenced mtSNP allele') + 
  guides(colour=guide_legend(title="Impute2 Info")) + 
  guides(size=guide_legend(title="Allele Frequency"))

# INFO SCORE
ggplot(summary.stats, aes(x = pos, y = info_comb, colour = as.factor(type), size = af)) + geom_point() + theme_bw() + 
  labs(x = 'mtDNA position', y = 'Info Score', title = 'Info Score of imputed mtSNPs') + 
  guides(colour=guide_legend(title="SNP Type")) 
