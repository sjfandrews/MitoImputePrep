require(ggplot2)
require(tidyr)
require(emmeans)

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
