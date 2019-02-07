require(ggplot2)
require(tidyr)
require(emmeans)
require(dplyr)

#########################################################################
chip.table = read.table("~/GitCode/MitoImputePrep/scripts/INFORMATION_LISTS/b37_platforms.txt", header = F)
truth.table = read.table("~/GitCode/MitoImputePrep/metadata/HaploGrep_concordance/chrMT_1kg_diploid_haplogrep.txt", header = T)
container = "/Volumes/TimMcInerney/MitoImpute/data/HAPLOGROUPS/chips/" #BDCHP-1X10-HUMANHAP240S_11216501_A-b37/MCMC_Experiments/MCMC1"

truth.A = subset(truth.table, substr(truth.table$Haplogroup, 1, 1) == "L")
truth.A$Macrohaplogroup = substr(truth.A$Haplogroup, 1, 2)
truth.N = subset(truth.table, substr(truth.table$Haplogroup, 1, 1) != "L")
truth.N$Macrohaplogroup = substr(truth.N$Haplogroup, 1, 1)

truth.table = rbind(truth.A, truth.N)
truth.table = arrange(truth.table, truth.table$SampleID)

###############################################################################
## MAF 

exp.dir = "MCMC_Experiments"
exp.var = "MCMC1"
maf_df = chip.table
names(maf_df) = c("chip")

for (chip in 1:nrow(chip.table)) {
  tmp.file = paste0(container, chip.table$V1[chip], "/", exp.dir, "/", exp.var, "/", "chrMT_1kg_", chip.table$V1[chip], "_imputed_", exp.var, "_haplogrep.txt")
  if (file.exists(tmp.file)) {
    maf_df$imputed[chip] = T
  }
}

for (chip in 1:nrow(chip.table)) {
  tmp.file = paste0(container, chip.table$V1[chip], "/", exp.dir, "/", exp.var, "/", "chrMT_1kg_", chip.table$V1[chip], "_imputed_", exp.var, "_haplogrep.txt")
  if (file.exists(tmp.file)) {
    tmp.hg.table = read.table(tmp.file, header = T)
    
    ## CREATE COLUMN FOR MACRO HAPLOGROUPS
    # USE FIRST LETTER AND NUMBER FOR AFRICAN CLADES
    tmp.hg.table.A = subset(tmp.hg.table, substr(tmp.hg.table$Haplogroup, 1, 1) == "L")
    tmp.hg.table.A$Macrohaplogroup = substr(tmp.hg.table.A$Haplogroup, 1, 2)
    # USE FIRST LETTER FOR NON-AFRICAN CLADES
    tmp.hg.table.N = subset(tmp.hg.table, substr(tmp.hg.table$Haplogroup, 1, 1) != "L")
    tmp.hg.table.N$Macrohaplogroup = substr(tmp.hg.table.N$Haplogroup, 1, 1)
    # COMBINED AND REARRANGE
    tmp.hg.table = rbind(tmp.hg.table.A, tmp.hg.table.N)
    tmp.hg.table = arrange(tmp.hg.table, tmp.hg.table$SampleID)
    
    
  }
}












############ OLD

###############################################################################
## MAF 
f = "~/GitCode/MitoImputePrep/metadata/Concordance_tables/ConcordanceTables_MAF_Combined.csv"
x = read.csv(f, header = T)
x[,2] = factor(x[,2])
x2 = x

x = x[-(3:9)]
x = spread(x, names(x)[2], "Imputed.hg.Conc")
v = names(x)[2:ncol(x)]
#x$largest = NA
largest = c()
for (i in 1:nrow(x)) {
  if (apply(x[i,2:ncol(x)], 1, function(u) all(u %in% NA))) {
    #x$largest[i] = NA
    largest = c(largest, NA)
  } else {
    #x$largest[i] = v[which(unlist(x[i,2:ncol(x)]) == max(unlist(x[i,2:ncol(x)])))]
    largest = c(largest, paste(v[which(unlist(x[i,2:ncol(x)]) == max(unlist(x[i,2:ncol(x)]), na.rm = T))], collapse = ","))
  }
}
x$largest = largest
table(x$largest)


x2$diff = x2$Imputed.hg.Conc - x2$Typed.hg.Conc
x3 = x2
x2 = x2[-(3:10)]
x2 = spread(x2, names(x2)[2], "diff")
#x$largest = NA
largest = c()
for (i in 1:nrow(x2)) {
  if (apply(x2[i,2:ncol(x2)], 1, function(u) all(u %in% NA))) {
    #x$largest[i] = NA
    largest = c(largest, NA)
  } else {
    #x$largest[i] = v[which(unlist(x[i,2:ncol(x)]) == max(unlist(x[i,2:ncol(x)])))]
    largest = c(largest, paste(v[which(unlist(x2[i,2:ncol(x2)]) == max(unlist(x2[i,2:ncol(x2)]), na.rm = T))], collapse = ","))
  }
}
x2$largest = largest
table(x2$largest)

l = lm(diff ~ MAF, data = x3)
anova(l)
summary(l)
emmeans(l, pairwise ~ MAF)
summary(emmeans(l, pairwise ~ MAF))

boxplot(diff ~ MAF, data = x3)

l2 = lm(Imputed.hg.Conc ~ MAF, data = x3)
anova(l2)
summary(l2)
emmeans(l2, pairwise ~ MAF, type="response")
summary(emmeans(l2, pairwise ~ MAF))

boxplot(Imputed.hg.Conc ~ MAF, data = x3)

###############################################################################
## kHAP 
f = "~/GitCode/MitoImputePrep/metadata/Concordance_tables/ConcordanceTables_kHAP_Combined.csv"
x = read.csv(f, header = T)
x[,2] = factor(x[,2])
x2 = x

x = x[-(3:9)]
x = spread(x, names(x)[2], "Imputed.hg.Conc")
v = names(x)[2:ncol(x)]
#x$largest = NA
largest = c()
for (i in 1:nrow(x)) {
  if (apply(x[i,2:ncol(x)], 1, function(u) all(u %in% NA))) {
    #x$largest[i] = NA
    largest = c(largest, NA)
  } else {
    #x$largest[i] = v[which(unlist(x[i,2:ncol(x)]) == max(unlist(x[i,2:ncol(x)])))]
    largest = c(largest, paste(v[which(unlist(x[i,2:ncol(x)]) == max(unlist(x[i,2:ncol(x)]), na.rm = T))], collapse = ","))
  }
}
x$largest = largest
table(x$largest)


x2$diff = x2$Imputed.hg.Conc - x2$Typed.hg.Conc
x3 = x2
x2 = x2[-(3:10)]
x2 = spread(x2, names(x2)[2], "diff")
#x$largest = NA
largest = c()
for (i in 1:nrow(x2)) {
  if (apply(x2[i,2:ncol(x2)], 1, function(u) all(u %in% NA))) {
    #x$largest[i] = NA
    largest = c(largest, NA)
  } else {
    #x$largest[i] = v[which(unlist(x[i,2:ncol(x)]) == max(unlist(x[i,2:ncol(x)])))]
    largest = c(largest, paste(v[which(unlist(x2[i,2:ncol(x2)]) == max(unlist(x2[i,2:ncol(x2)]), na.rm = T))], collapse = ","))
  }
}
x2$largest = largest
table(x2$largest)

l = lm(diff ~ kHAP, data = x3)
anova(l)
summary(l)
emmeans(l, pairwise ~ kHAP)
summary(emmeans(l, pairwise ~ kHAP))

boxplot(diff ~ kHAP, data = x3)

l2 = lm(Imputed.hg.Conc ~ kHAP, data = x3)
anova(l2)
summary(l2)
emmeans(l2, pairwise ~ kHAP, type="response")
summary(emmeans(l2, pairwise ~ kHAP))

boxplot(Imputed.hg.Conc ~ kHAP, data = x3)

###############################################################################
## MCMC 
f = "~/GitCode/MitoImputePrep/metadata/Concordance_tables/ConcordanceTables_MCMC_Combined.csv"
x = read.csv(f, header = T)
x[,2] = factor(x[,2])
x2 = x

x = x[-(3:9)]
x = spread(x, names(x)[2], "Imputed.hg.Conc")
v = names(x)[2:ncol(x)]
#x$largest = NA
largest = c()
for (i in 1:nrow(x)) {
  if (apply(x[i,2:ncol(x)], 1, function(u) all(u %in% NA))) {
    #x$largest[i] = NA
    largest = c(largest, NA)
  } else {
    #x$largest[i] = v[which(unlist(x[i,2:ncol(x)]) == max(unlist(x[i,2:ncol(x)])))]
    largest = c(largest, paste(v[which(unlist(x[i,2:ncol(x)]) == max(unlist(x[i,2:ncol(x)]), na.rm = T))], collapse = ","))
  }
}
x$largest = largest
table(x$largest)


x2$diff = x2$Imputed.hg.Conc - x2$Typed.hg.Conc
x3 = x2
x2 = x2[-(3:10)]
x2 = spread(x2, names(x2)[2], "diff")
#x$largest = NA
largest = c()
for (i in 1:nrow(x2)) {
  if (apply(x2[i,2:ncol(x2)], 1, function(u) all(u %in% NA))) {
    #x$largest[i] = NA
    largest = c(largest, NA)
  } else {
    #x$largest[i] = v[which(unlist(x[i,2:ncol(x)]) == max(unlist(x[i,2:ncol(x)])))]
    largest = c(largest, paste(v[which(unlist(x2[i,2:ncol(x2)]) == max(unlist(x2[i,2:ncol(x2)]), na.rm = T))], collapse = ","))
  }
}
x2$largest = largest
table(x2$largest)

l = lm(diff ~ MCMC, data = x3)
anova(l)
summary(l)
emmeans(l, pairwise ~ MCMC)
summary(emmeans(l, pairwise ~ MCMC))

boxplot(diff ~ MCMC, data = x3)

l2 = lm(Imputed.hg.Conc ~ MCMC, data = x3)
anova(l2)
summary(l2)
emmeans(l2, pairwise ~ MCMC)
#emmeans(l2, pairwise ~ MCMC, type="response") # USE IF LOG SCALE IN lm FUNCTION
summary(emmeans(l2, pairwise ~ MCMC))

boxplot(Imputed.hg.Conc ~ MCMC, data = x3)
