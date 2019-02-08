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
## MCMC 

exp.dir = "MCMC_Experiments"
exp.var = c("MCMC1", "MCMC5", "MCMC10", "MCMC20", "MCMC30")
maf_df = data.frame(matrix(ncol = 1, nrow = nrow(chip.table)))
names(maf_df) = c("chip")
maf_df$chip = chip.table$V1
maf_df$experiment = exp.dir

main_df = data.frame()

for (exp in 1:length(exp.var)) {
  out.file = paste0("/Volumes/TimMcInerney/MitoImpute/data/HAPLOGROUPS/combined/ConcordanceTables_", exp.var[exp], ".csv")
  maf_df$sub_experiment = exp.var[exp]
  for (chip in 1:nrow(chip.table)) {
    tmp.file = paste0(container, chip.table$V1[chip], "/", exp.dir, "/", exp.var[exp], "/", "chrMT_1kg_", chip.table$V1[chip], "_imputed_", exp.var[exp], "_haplogrep.txt")
    if (file.exists(tmp.file) == T) {
      maf_df$imputed[chip] = T
      chip.table$imputed[chip] = T
    } else {
      maf_df$imputed[chip] = F
      chip.table$imputed[chip] = F
    }
  }
  
  maf_df$typed_match[chip] = NA
  maf_df$typed_macro_match[chip] = NA
  maf_df$imputed_match[chip] = NA
  maf_df$imputed_macro_match[chip] = NA
  
  total = nrow(chip.table)
  total_imputed = nrow(subset(chip.table, chip.table$imputed == T))
  
  for (chip in 1:nrow(chip.table)) {
    message(paste0("WORKING ON CHIP FOR ", exp.var[exp], ":\t", chip, " / ", total))
    # TYPED FILE
    tmp1.file = paste0(container, chip.table$V1[chip], "/", exp.dir, "/", "chrMT_1kg_", chip.table$V1[chip], "_diploid.txt")
    if (file.exists(tmp1.file) == T) {
      tmp1.hg.table = read.table(tmp1.file, header = T)
      tmp1.hg.table$Range = NULL
      
      ## CREATE COLUMN FOR MACRO HAPLOGROUPS
      # USE FIRST LETTER AND NUMBER FOR AFRICAN CLADES
      tmp1.hg.table.A = subset(tmp1.hg.table, substr(tmp1.hg.table$Haplogroup, 1, 1) == "L")
      tmp1.hg.table.A$Macrohaplogroup = substr(tmp1.hg.table.A$Haplogroup, 1, 2)
      # USE FIRST LETTER FOR NON-AFRICAN CLADES
      tmp1.hg.table.N = subset(tmp1.hg.table, substr(tmp1.hg.table$Haplogroup, 1, 1) != "L")
      tmp1.hg.table.N$Macrohaplogroup = substr(tmp1.hg.table.N$Haplogroup, 1, 1)
      # COMBINED AND REARRANGE
      tmp1.hg.table = rbind(tmp1.hg.table.A, tmp1.hg.table.N)
      tmp1.hg.table = arrange(tmp1.hg.table, tmp1.hg.table$SampleID)
      
      tmp1.hg.table$typed_match = as.character(truth.table$Haplogroup) == as.character(tmp1.hg.table$Haplogroup)
      tmp1.hg.table$typed_macro_match = as.character(truth.table$Macrohaplogroup) == as.character(tmp1.hg.table$Macrohaplogroup)
      
      maf_df$typed_match[chip] = nrow(subset(tmp1.hg.table, tmp1.hg.table$typed_match == T)) / total_imputed
      maf_df$typed_macro_match[chip] = nrow(subset(tmp1.hg.table, tmp1.hg.table$typed_macro_match == T)) / total_imputed
    }
    
    # IMPUTED FILE
    tmp2.file = paste0(container, chip.table$V1[chip], "/", exp.dir, "/", exp.var[exp], "/", "chrMT_1kg_", chip.table$V1[chip], "_imputed_", exp.var[exp], "_haplogrep.txt")
    if (file.exists(tmp2.file) == T) {
      tmp2.hg.table = read.table(tmp2.file, header = T)
      tmp2.hg.table$Range = NULL
      
      ## CREATE COLUMN FOR MACRO HAPLOGROUPS
      # USE FIRST LETTER AND NUMBER FOR AFRICAN CLADES
      tmp2.hg.table.A = subset(tmp2.hg.table, substr(tmp2.hg.table$Haplogroup, 1, 1) == "L")
      tmp2.hg.table.A$Macrohaplogroup = substr(tmp2.hg.table.A$Haplogroup, 1, 2)
      # USE FIRST LETTER FOR NON-AFRICAN CLADES
      tmp2.hg.table.N = subset(tmp2.hg.table, substr(tmp2.hg.table$Haplogroup, 1, 1) != "L")
      tmp2.hg.table.N$Macrohaplogroup = substr(tmp2.hg.table.N$Haplogroup, 1, 1)
      # COMBINED AND REARRANGE
      tmp2.hg.table = rbind(tmp2.hg.table.A, tmp2.hg.table.N)
      tmp2.hg.table = arrange(tmp2.hg.table, tmp2.hg.table$SampleID)
      
      tmp2.hg.table$imputed_match = as.character(truth.table$Haplogroup) == as.character(tmp2.hg.table$Haplogroup)
      tmp2.hg.table$imputed_macro_match = as.character(truth.table$Macrohaplogroup) == as.character(tmp2.hg.table$Macrohaplogroup)
      
      maf_df$imputed_match[chip] = nrow(subset(tmp2.hg.table, tmp2.hg.table$imputed_match == T)) / total_imputed
      maf_df$imputed_macro_match[chip] = nrow(subset(tmp2.hg.table, tmp2.hg.table$imputed_macro_match == T)) / total_imputed
    }
  }
  write.csv(maf_df, out.file, row.names = F, quote = F)
  message(paste0("WROTE ", out.file, " TO DISK"))
  main_df = rbind(main_df, maf_df)
}
write.csv(main_df, paste0("/Volumes/TimMcInerney/MitoImpute/data/HAPLOGROUPS/combined/ConcordanceTables_", exp.dir,"_COMBINED.csv"), row.names = F, quote = F)

main_df$sub_experiment = factor(main_df$sub_experiment, levels = exp.var)

mcmc_box = ggplot(main_df, aes(x = sub_experiment, y = imputed_macro_match)) +
  geom_violin(fill = "#feb600", na.rm = T) +
  geom_boxplot(width = rel(0.25), notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428") +
  theme_bw() +
  labs(x = "Number of reference haplotypes used",
       y = "% concordance with resequenced dataset")
ggsave(filename = paste0("/Volumes/TimMcInerney/MitoImpute/data/HAPLOGROUPS/plots/ConcordanceTables_", exp.dir,".png"), plot = mcmc_box, units = "mm", width = 297, height = 210, dpi = 300)

## KHAP 

exp.dir = "kHAP_Experiments"
exp.var = c("kHAP100", "kHAP250", "kHAP500", "kHAP1000", "kHAP2500", "kHAP5000", "kHAP10000", "kHAP20000", "kHAP30000")
maf_df = data.frame(matrix(ncol = 1, nrow = nrow(chip.table)))
names(maf_df) = c("chip")
maf_df$chip = chip.table$V1
maf_df$experiment = exp.dir

main_df = data.frame()

for (exp in 1:length(exp.var)) {
  out.file = paste0("/Volumes/TimMcInerney/MitoImpute/data/HAPLOGROUPS/combined/ConcordanceTables_", exp.var[exp], ".csv")
  maf_df$sub_experiment = exp.var[exp]
  for (chip in 1:nrow(chip.table)) {
    tmp.file = paste0(container, chip.table$V1[chip], "/", exp.dir, "/", exp.var[exp], "/", "chrMT_1kg_", chip.table$V1[chip], "_imputed_", exp.var[exp], "_haplogrep.txt")
    if (file.exists(tmp.file) == T) {
      maf_df$imputed[chip] = T
      chip.table$imputed[chip] = T
    } else {
      maf_df$imputed[chip] = F
      chip.table$imputed[chip] = F
    }
  }
  
  maf_df$typed_match[chip] = NA
  maf_df$typed_macro_match[chip] = NA
  maf_df$imputed_match[chip] = NA
  maf_df$imputed_macro_match[chip] = NA
  
  total = nrow(chip.table)
  total_imputed = nrow(subset(chip.table, chip.table$imputed == T))
  
  for (chip in 1:nrow(chip.table)) {
    message(paste0("WORKING ON CHIP FOR ", exp.var[exp], ":\t", chip, " / ", total))
    # TYPED FILE
    tmp1.file = paste0(container, chip.table$V1[chip], "/", exp.dir, "/", "chrMT_1kg_", chip.table$V1[chip], "_diploid.txt")
    if (file.exists(tmp1.file) == T) {
      tmp1.hg.table = read.table(tmp1.file, header = T)
      tmp1.hg.table$Range = NULL
      
      ## CREATE COLUMN FOR MACRO HAPLOGROUPS
      # USE FIRST LETTER AND NUMBER FOR AFRICAN CLADES
      tmp1.hg.table.A = subset(tmp1.hg.table, substr(tmp1.hg.table$Haplogroup, 1, 1) == "L")
      tmp1.hg.table.A$Macrohaplogroup = substr(tmp1.hg.table.A$Haplogroup, 1, 2)
      # USE FIRST LETTER FOR NON-AFRICAN CLADES
      tmp1.hg.table.N = subset(tmp1.hg.table, substr(tmp1.hg.table$Haplogroup, 1, 1) != "L")
      tmp1.hg.table.N$Macrohaplogroup = substr(tmp1.hg.table.N$Haplogroup, 1, 1)
      # COMBINED AND REARRANGE
      tmp1.hg.table = rbind(tmp1.hg.table.A, tmp1.hg.table.N)
      tmp1.hg.table = arrange(tmp1.hg.table, tmp1.hg.table$SampleID)
      
      tmp1.hg.table$typed_match = as.character(truth.table$Haplogroup) == as.character(tmp1.hg.table$Haplogroup)
      tmp1.hg.table$typed_macro_match = as.character(truth.table$Macrohaplogroup) == as.character(tmp1.hg.table$Macrohaplogroup)
      
      maf_df$typed_match[chip] = nrow(subset(tmp1.hg.table, tmp1.hg.table$typed_match == T)) / total_imputed
      maf_df$typed_macro_match[chip] = nrow(subset(tmp1.hg.table, tmp1.hg.table$typed_macro_match == T)) / total_imputed
    }
    
    # IMPUTED FILE
    tmp2.file = paste0(container, chip.table$V1[chip], "/", exp.dir, "/", exp.var[exp], "/", "chrMT_1kg_", chip.table$V1[chip], "_imputed_", exp.var[exp], "_haplogrep.txt")
    if (file.exists(tmp2.file) == T) {
      tmp2.hg.table = read.table(tmp2.file, header = T)
      tmp2.hg.table$Range = NULL
      
      ## CREATE COLUMN FOR MACRO HAPLOGROUPS
      # USE FIRST LETTER AND NUMBER FOR AFRICAN CLADES
      tmp2.hg.table.A = subset(tmp2.hg.table, substr(tmp2.hg.table$Haplogroup, 1, 1) == "L")
      tmp2.hg.table.A$Macrohaplogroup = substr(tmp2.hg.table.A$Haplogroup, 1, 2)
      # USE FIRST LETTER FOR NON-AFRICAN CLADES
      tmp2.hg.table.N = subset(tmp2.hg.table, substr(tmp2.hg.table$Haplogroup, 1, 1) != "L")
      tmp2.hg.table.N$Macrohaplogroup = substr(tmp2.hg.table.N$Haplogroup, 1, 1)
      # COMBINED AND REARRANGE
      tmp2.hg.table = rbind(tmp2.hg.table.A, tmp2.hg.table.N)
      tmp2.hg.table = arrange(tmp2.hg.table, tmp2.hg.table$SampleID)
      
      tmp2.hg.table$imputed_match = as.character(truth.table$Haplogroup) == as.character(tmp2.hg.table$Haplogroup)
      tmp2.hg.table$imputed_macro_match = as.character(truth.table$Macrohaplogroup) == as.character(tmp2.hg.table$Macrohaplogroup)
      
      maf_df$imputed_match[chip] = nrow(subset(tmp2.hg.table, tmp2.hg.table$imputed_match == T)) / total_imputed
      maf_df$imputed_macro_match[chip] = nrow(subset(tmp2.hg.table, tmp2.hg.table$imputed_macro_match == T)) / total_imputed
    }
  }
  write.csv(maf_df, out.file, row.names = F, quote = F)
  message(paste0("WROTE ", out.file, " TO DISK"))
  main_df = rbind(main_df, maf_df)
}
write.csv(main_df, paste0("/Volumes/TimMcInerney/MitoImpute/data/HAPLOGROUPS/combined/ConcordanceTables_", exp.dir,"_COMBINED.csv"), row.names = F, quote = F)

main_df$sub_experiment = factor(main_df$sub_experiment, levels = exp.var)

k_hap_box = ggplot(main_df, aes(x = sub_experiment, y = imputed_macro_match)) +
  geom_violin(fill = "#feb600", na.rm = T) +
  geom_boxplot(width = 0.125, notch = T, fill = "#ea4e3c", na.rm = T, outlier.colour = "#802428") +
  theme_bw() +
  labs(x = "Number of reference haplotypes used",
       y = "% concordance with resequenced dataset")
ggsave(filename = paste0("/Volumes/TimMcInerney/MitoImpute/data/HAPLOGROUPS/plots/ConcordanceTables_", exp.dir,".png"), plot = k_hap_box, units = "mm", width = 297, height = 210, dpi = 300)

## MAF 

ref.panel = c("ReferencePanel_v2", "ReferencePanel_v4", "ReferencePanel_v3")
exp.var = c("MAF1%", "MAF0.5%", "MAF0.1%")
maf_df = data.frame(matrix(ncol = 1, nrow = nrow(chip.table)))
names(maf_df) = c("chip")
maf_df$chip = chip.table$V1
maf_df$experiment = exp.dir

main_df = data.frame()

for (exp in 1:length(exp.var)) {
  out.file = paste0("/Volumes/TimMcInerney/MitoImpute/data/HAPLOGROUPS/combined/ConcordanceTables_", exp.var[exp], ".csv")
  maf_df$sub_experiment = exp.var[exp]
  for (chip in 1:nrow(chip.table)) {
    tmp.file = paste0(container, chip.table$V1[chip], "/", ref.panel[exp], "/", "chrMT_1kg_", chip.table$V1[chip], "_imputed_", exp.var[exp], "_haplogrep.txt")
    if (file.exists(tmp.file) == T) {
      maf_df$imputed[chip] = T
      chip.table$imputed[chip] = T
    } else {
      maf_df$imputed[chip] = F
      chip.table$imputed[chip] = F
    }
  }
  
  maf_df$typed_match[chip] = NA
  maf_df$typed_macro_match[chip] = NA
  maf_df$imputed_match[chip] = NA
  maf_df$imputed_macro_match[chip] = NA
  
  total = nrow(chip.table)
  total_imputed = nrow(subset(chip.table, chip.table$imputed == T))
  
  for (chip in 1:nrow(chip.table)) {
    message(paste0("WORKING ON CHIP FOR ", exp.var[exp], ":\t", chip, " / ", total))
    # TYPED FILE
    tmp1.file = paste0(container, chip.table$V1[chip], "/", exp.dir, "/", "chrMT_1kg_", chip.table$V1[chip], "_diploid.txt")
    if (file.exists(tmp1.file) == T) {
      tmp1.hg.table = read.table(tmp1.file, header = T)
      tmp1.hg.table$Range = NULL
      
      ## CREATE COLUMN FOR MACRO HAPLOGROUPS
      # USE FIRST LETTER AND NUMBER FOR AFRICAN CLADES
      tmp1.hg.table.A = subset(tmp1.hg.table, substr(tmp1.hg.table$Haplogroup, 1, 1) == "L")
      tmp1.hg.table.A$Macrohaplogroup = substr(tmp1.hg.table.A$Haplogroup, 1, 2)
      # USE FIRST LETTER FOR NON-AFRICAN CLADES
      tmp1.hg.table.N = subset(tmp1.hg.table, substr(tmp1.hg.table$Haplogroup, 1, 1) != "L")
      tmp1.hg.table.N$Macrohaplogroup = substr(tmp1.hg.table.N$Haplogroup, 1, 1)
      # COMBINED AND REARRANGE
      tmp1.hg.table = rbind(tmp1.hg.table.A, tmp1.hg.table.N)
      tmp1.hg.table = arrange(tmp1.hg.table, tmp1.hg.table$SampleID)
      
      tmp1.hg.table$typed_match = as.character(truth.table$Haplogroup) == as.character(tmp1.hg.table$Haplogroup)
      tmp1.hg.table$typed_macro_match = as.character(truth.table$Macrohaplogroup) == as.character(tmp1.hg.table$Macrohaplogroup)
      
      maf_df$typed_match[chip] = nrow(subset(tmp1.hg.table, tmp1.hg.table$typed_match == T)) / total_imputed
      maf_df$typed_macro_match[chip] = nrow(subset(tmp1.hg.table, tmp1.hg.table$typed_macro_match == T)) / total_imputed
    }
    
    # IMPUTED FILE
    tmp2.file = paste0(container, chip.table$V1[chip], "/", exp.dir, "/", exp.var[exp], "/", "chrMT_1kg_", chip.table$V1[chip], "_imputed_", exp.var[exp], "_haplogrep.txt")
    if (file.exists(tmp2.file) == T) {
      tmp2.hg.table = read.table(tmp2.file, header = T)
      tmp2.hg.table$Range = NULL
      
      ## CREATE COLUMN FOR MACRO HAPLOGROUPS
      # USE FIRST LETTER AND NUMBER FOR AFRICAN CLADES
      tmp2.hg.table.A = subset(tmp2.hg.table, substr(tmp2.hg.table$Haplogroup, 1, 1) == "L")
      tmp2.hg.table.A$Macrohaplogroup = substr(tmp2.hg.table.A$Haplogroup, 1, 2)
      # USE FIRST LETTER FOR NON-AFRICAN CLADES
      tmp2.hg.table.N = subset(tmp2.hg.table, substr(tmp2.hg.table$Haplogroup, 1, 1) != "L")
      tmp2.hg.table.N$Macrohaplogroup = substr(tmp2.hg.table.N$Haplogroup, 1, 1)
      # COMBINED AND REARRANGE
      tmp2.hg.table = rbind(tmp2.hg.table.A, tmp2.hg.table.N)
      tmp2.hg.table = arrange(tmp2.hg.table, tmp2.hg.table$SampleID)
      
      tmp2.hg.table$imputed_match = as.character(truth.table$Haplogroup) == as.character(tmp2.hg.table$Haplogroup)
      tmp2.hg.table$imputed_macro_match = as.character(truth.table$Macrohaplogroup) == as.character(tmp2.hg.table$Macrohaplogroup)
      
      maf_df$imputed_match[chip] = nrow(subset(tmp2.hg.table, tmp2.hg.table$imputed_match == T)) / total_imputed
      maf_df$imputed_macro_match[chip] = nrow(subset(tmp2.hg.table, tmp2.hg.table$imputed_macro_match == T)) / total_imputed
    }
  }
  print(head(maf_df))
  write.csv(maf_df, out.file, row.names = F, quote = F)
  message(paste0("WROTE ", out.file, " TO DISK"))
  main_df = rbind(main_df, maf_df)
}
write.csv(main_df, paste0("/Volumes/TimMcInerney/MitoImpute/data/HAPLOGROUPS/combined/ConcordanceTables_", exp.dir,"_COMBINED.csv"), row.names = F, quote = F)



ggplot(main_df, aes(x = sub_experiment, y = imputed_macro_match)) +
  geom_boxplot()






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
