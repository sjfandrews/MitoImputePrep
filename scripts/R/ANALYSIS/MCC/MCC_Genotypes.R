library(tidyr)
library(dplyr)
library(readr)
library(purrr)
library(tibble)
library(ggplot2)

printTime = function(end.time) {
  if (end.time[3] < 60) {
    out.time = end.time[3]
    return(paste0('Process completed in ',round(out.time,2),' seconds'))
  }
  if (end.time[3] >= 60 & end.time[3] < 3600) {
    out.time = end.time[3] / 60
    return(paste0('Process completed in ',round(out.time,2),' minutes'))
  }
  if (end.time[3] >= 3600 & end.time[3] < 86400) {
    out.time = (end.time[3] / 60) / 60
    return(paste0('Process completed in ',round(out.time,2),' hours'))
  }
  if (end.time[3] >= 86400) {
    out.time = ((end.time[3] / 60) / 60) / 24
    return(paste0('Process completed in ',round(out.time,2),' days'))
  }
}

mccr <- function (act, pred) {
  TP <- sum(act %in% 1 & pred %in% 1)
  TN <- sum(act %in% 0 & pred %in% 0)
  FP <- sum(act %in% 0 & pred %in% 1)
  FN <- sum(act %in% 1 & pred %in% 0)
  denom <- as.double((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  if (any((TP + FP) == 0, (TP + FN) == 0, (TN + FP) == 0, (TN + FN) == 0)) 
    denom <- 1
  mcc <- ((TP * TN) - (FP * FN))/sqrt(denom)
  return(mcc)
}

start.time = proc.time() # Start the timer!

### SPECIFY ARGUMENTS!
print("USAGE HINTS")
print("ARGUMENT 1:  WGS VCF FILE")
print("ARGUMENT 2:  Genotyped VCF FILE")
print("ARGUMENT 3:  Imputed VCF FILE")
print("ARGUMENT 4:  IMPUTE2 INFO FILE")
print("ARGUMENT 5:  OUTPUT FILE PREFIX")

args = commandArgs(trailingOnly = TRUE) # Set arguments from the command line

wgs.vcf = args[1] # WGS VCF file
typ.vcf = args[2] # Genotyped VCF file
imp.vcf = args[3] # Imputed VCF file
imp.info = args[4] # IMPUTE2 INFO file
out.file = args[5] # Outfile prefix

# FIX OUTFILE
out_file_imp = paste0(out.file, "_imputed_MCC.csv")
out_file_typ = paste0(out.file, "_typed_MCC.csv")
#out_file_imp = "~/GitCode/MitoImputePrep/metadata/MCC_files/BDCHP-1X10-HUMANHAP240S_11216501_A-b37/BDCHP-1X10-HUMANHAP240S_11216501_A-b37_imputed_MCC.txt"
#out_file_typ = "~/GitCode/MitoImputePrep/metadata/MCC_files/BDCHP-1X10-HUMANHAP240S_11216501_A-b37/BDCHP-1X10-HUMANHAP240S_11216501_A-b37_typed_MCC.txt"

## DATA READING
# READ IN WGS VCF FILE
#wgs.vcf = "/Volumes/TimMcInerney/MitoImpute/data/VCF/1kGP_chrMT_SNPonly.vcf.gz"
wgs_1kg.vcf = read_tsv(wgs.vcf, comment = '##', na = c(".", "", "NA"))

# READ IN THE TYPED VCF FILE
#typ.vcf = "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/BDCHP-1X10-HUMANHAP240S_11216501_A-b37/ReferencePanel_v2/chrMT_1kg_BDCHP-1X10-HUMANHAP240S_11216501_A-b37.vcf.gz"
typ_1kg.vcf = read_tsv(typ.vcf, comment = '##', na = c(".", "", "NA"))

# READ IN THE IMPUTED VCF FILE AND INFO FILE
#imp.vcf = "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/BDCHP-1X10-HUMANHAP240S_11216501_A-b37/ReferencePanel_v2/MCMC1/chrMT_1kg_BDCHP-1X10-HUMANHAP240S_11216501_A-b37_imputed_MCMC1_haplogrep.vcf.gz"
#imp.info = "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/BDCHP-1X10-HUMANHAP240S_11216501_A-b37/ReferencePanel_v2/MCMC1/chrMT_1kg_BDCHP-1X10-HUMANHAP240S_11216501_A-b37_imputed_MCMC1_info"
imp_1kg.vcf <- read_tsv(imp.vcf, comment = '##', na = c(".", "", "NA"))
imp_1kg.info <- read_delim(imp.info, delim = " ")

### DATA MUNGING
## IMPUTED | WGS
#  obtain intersect of SNPs from imputed and wgs 
snp.intersect.imp <- wgs_1kg.vcf %>% 
  select(POS) %>% 
  semi_join(imp_1kg.vcf, by = 'POS') %>%
  mutate(POS = paste0('mt', POS))

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
imp_1kg <- as_tibble(sapply(imp_1kg, function(x) gsub('\\:.*', "", x)))
imp_1kg[imp_1kg == './.'] <- NA
imp_1kg[imp_1kg == '.'] <- NA
imp_1kg[imp_1kg == '0/0'] <- "0"
imp_1kg[imp_1kg == '0'] <- "0"
imp_1kg[imp_1kg == '1/1'] <- "1"
imp_1kg[imp_1kg == '1'] <- "1"
imp_1kg <- imp_1kg %>% 
  select(var_name, snp.intersect.imp$POS) %>% 
  mutate_at(vars(contains('mt')), as.integer)

## TYPED | WGS
# Obtain the intersect of SNPs from typed and wgs
snp.intersect.typ <- wgs_1kg.vcf %>% 
  select(POS) %>% 
  semi_join(typ_1kg.vcf, by = 'POS') %>%
  mutate(POS = paste0('mt', POS))

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
typ_1kg <- as_tibble(sapply(typ_1kg, function(x) gsub('\\:.*', "", x)))
typ_1kg[typ_1kg == './.'] <- NA
typ_1kg[typ_1kg == '.'] <- NA
typ_1kg[typ_1kg == '0/0'] <- "0"
typ_1kg[typ_1kg == '0'] <- "0"
typ_1kg[typ_1kg == '1/1'] <- "1"
typ_1kg[typ_1kg == '1'] <- "1"
typ_1kg <- typ_1kg %>% 
  select(var_name, snp.intersect.typ$POS) %>% 
  mutate_at(vars(contains('mt')), as.integer)

##  Munge WGS dataframe
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
wgs_1kg <- as_tibble(sapply(wgs_1kg, function(x) gsub('\\:.*', "", x)))
wgs_1kg[wgs_1kg == './.'] <- NA
wgs_1kg[wgs_1kg == '.'] <- NA
wgs_1kg[wgs_1kg == '0/0'] <- "0"
wgs_1kg[wgs_1kg == '0'] <- "0"
wgs_1kg[wgs_1kg == '1/1'] <- "1"
wgs_1kg[wgs_1kg == '1'] <- "1"
wgs_1kg.imp <- wgs_1kg %>% 
  select(var_name, snp.intersect.imp$POS) %>% 
  mutate_at(vars(contains('mt')), as.integer)
wgs_1kg.typ <- wgs_1kg %>% 
  select(var_name, snp.intersect.typ$POS) %>% 
  mutate_at(vars(contains('mt')), as.integer) 

## MUNGE INFO TABLE
imp_1kg.info <- mutate(imp_1kg.info, info_comb = ifelse(info_type0 == -1, info,info_type0 ))
imp_1kg.info <- mutate(imp_1kg.info, himc = ifelse(position %in% c(825, 1018, 1438, 1719, 1736, 2092, 3505, 3552, 3594, 4580, 4769, 4917, 4977, 5178, 5442, 6371, 7028, 8251, 8414, 8468, 8703, 9042, 9055, 9347, 9950, 10115, 10398, 10398, 10400, 10550, 11177, 11251, 11947, 12007, 12308, 12705, 13263, 13368, 13506, 13708, 13789, 14178, 14318, 14470, 14560, 14668, 14766, 15043, 15326, 15452, 15535, 16111, 16189, 16391), 'yes', 'no'))

### SNP Concordance Statisitics
# NOTE: SHOULD ONLY REALLY BE LOOKING AT TYPE 1 and 2 SNPS (in REF ONLY and SAMPLE ONLY)
# TYPE 0 SHOULD HAVE A MCC = 1
# TYPE 2 MAY BE PERFECT, SO LOOK AT THIS
# ONLY WANT TO LOOK AT MCC FOR IMPUTED SNPs

## calculate MCC between [imputed & wgs] and [typed & wgs] SNPs
mccr.geno.imp <- unlist(map2(wgs_1kg.imp[,2:ncol(wgs_1kg.imp)], imp_1kg[,2:ncol(imp_1kg)], function(a,b) mccr(a,b)))
mccr.geno.typ <- unlist(map2(wgs_1kg.typ[,2:ncol(wgs_1kg.typ)], typ_1kg[,2:ncol(typ_1kg)], function(a,b) mccr(a,b)))

##  Calculate concordance between imputed and typed SNPs
concodance.geno.imp <- unlist(map2(imp_1kg[,2:ncol(imp_1kg)], wgs_1kg.imp[,2:ncol(wgs_1kg.imp)], function(a, b) sum(a == b, na.rm = T)/length(a)))
concodance.geno.typ <- unlist(map2(typ_1kg[,2:ncol(typ_1kg)], wgs_1kg.typ[,2:ncol(wgs_1kg.typ)], function(a, b) sum(a == b, na.rm = T)/length(a)))

##  Calculate allele frequency of typed SNPs from WGS data
af.imp <- unlist(map(wgs_1kg.imp[,2:ncol(wgs_1kg.imp)], function(a) sum(a, na.rm = T)/length(a)))
af.typ <- unlist(map(wgs_1kg.typ[,2:ncol(wgs_1kg.typ)], function(a) sum(a, na.rm = T)/length(a)))

##  Dataframe for summary stats
summary.stats.imp <- tibble(mtSNP = colnames(wgs_1kg.imp[,2:ncol(wgs_1kg.imp)]),
                            pos = as.integer(gsub('mt', '', colnames(wgs_1kg.imp[,2:ncol(wgs_1kg.imp)]))),
                            af = af.imp,
                            mcc = mccr.geno.imp,
                            concodance = concodance.geno.imp)
summary.stats.typ <- tibble(mtSNP = colnames(wgs_1kg.typ[,2:ncol(wgs_1kg.typ)]),
                            pos = as.integer(gsub('mt', '', colnames(wgs_1kg.typ[,2:ncol(wgs_1kg.typ)]))),
                            af = af.typ,
                            mcc = mccr.geno.typ,
                            concodance = concodance.geno.typ)

##  merge on info.score file
summary.stats.imp <- left_join(summary.stats.imp, select(imp_1kg.info, c(-snp_id, -rs_id)), by = c('pos' = 'position')) 
#summary.stats.imp <- mutate(summary.stats.imp, info.cat = cut_width(summary.stats.imp$info, 0.25, boundary = 0))
#summary.stats.imp$info.cat = sub(",", "-", summary.stats.imp$info.cat)
if (length(unique(summary.stats.imp$info)) == 1) {
  tmp_val = as.numeric(summary.stats.imp$info[1])
  if (tmp_val > 0.75) {
    summary.stats.imp$info.cat = "(0.75-1]"
  } else if (tmp_val > 0.5 && tmp_val <= 0.75) {
    summary.stats.imp$info.cat = "(0.5-0.75]"
  } else if (tmp_val >= 0.25 && tmp_val <= 0.5) {
    summary.stats.imp$info.cat = "[0.25-0.5]"
  } else {
    summary.stats.imp$info.cat = "[0.0-0.25)"
  }
} else {
  summary.stats.imp <- mutate(summary.stats.imp, info.cat = cut_width(summary.stats.imp$info, 0.25, boundary = 0))
  summary.stats.imp$info.cat = sub(",", "-", summary.stats.imp$info.cat)
}
summary.stats.imp$info.cat = sub(",", "-", summary.stats.imp$info.cat)
##  basic summary stats
summary.stats.imp %>% count(af > 0.01)
summary.stats.imp %>% count(mcc > 0.4)
summary.stats.imp %>% count(concodance < 0.9)
summary.stats.imp %>% count(info > 0.3); summary.stats.imp %>% count(info > 0.5)

print(summary(summary.stats.imp))

##  merge on info.score file
summary.stats.typ <- left_join(summary.stats.typ, select(imp_1kg.info, c(-snp_id, -rs_id)), by = c('pos' = 'position'))
if (length(unique(summary.stats.typ$info)) == 1) {
  tmp_val = as.numeric(summary.stats.typ$info[1])
  if (tmp_val > 0.75) {
    summary.stats.typ$info.cat = "(0.75-1]"
  } else if (tmp_val > 0.5 && tmp_val <= 0.75) {
    summary.stats.typ$info.cat = "(0.5-0.75]"
  } else if (tmp_val >= 0.25 && tmp_val <= 0.5) {
    summary.stats.typ$info.cat = "[0.25-0.5]"
  } else {
    summary.stats.typ$info.cat = "[0.0-0.25)"
  }
} else {
  summary.stats.typ <- mutate(summary.stats.typ, info.cat = cut_width(summary.stats.typ$info, 0.25, boundary = 0))
  summary.stats.typ$info.cat = sub(",", "-", summary.stats.typ$info.cat)
}
summary.stats.typ$info.cat = sub(",", "-", summary.stats.typ$info.cat)
##  basic summary stats
#summary.stats.typ %>% count(af > 0.01)
#summary.stats.typ %>% count(mcc > 0.4)
#summary.stats.typ %>% count(concodance < 0.9)
#summary.stats.typ %>% count(info > 0.3); summary.stats.typ %>% count(info > 0.5)

##
write.csv(summary.stats.imp, out_file_imp, quote = F, row.names = F)
message(paste0("IMPUTED MCC STATS WRITTEN TO: ", out_file_imp))
write.csv(summary.stats.typ, out_file_typ, quote = F, row.names = F)
message(paste0("GENOTYPED MCC STATS WRITTEN TO: ", out_file_typ))
#summary.stats.imp
#summary.stats.typ

## END!
timer = T
if (timer == TRUE) {
  print("")
  print(printTime(proc.time() - start.time))
  print("")
}
