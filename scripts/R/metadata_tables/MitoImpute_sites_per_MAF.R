require(dplyr)
'%!in%' = function(x,y)!('%in%'(x,y))

MAF1pc = read.csv("~/GitCode/MitoImputePrep/metadata/MitoImpite_sites_per_MAF/MAF_1pc.csv", header = F)
MAF0.5pc = read.csv("~/GitCode/MitoImputePrep/metadata/MitoImpite_sites_per_MAF/MAF_0-5pc.csv", header = F)
MAF0.1pc = read.csv("~/GitCode/MitoImputePrep/metadata/MitoImpite_sites_per_MAF/MAF_0-1pc.csv", header = F)

sites = data.frame(MAF0.1pc$V1)
names(sites) = c("sites")
distinct_sites = distinct(sites)
multiallelic = sites$sites[duplicated(sites$sites) == T]

out_df = data.frame(matrix(NA, ncol = 1, nrow = nrow(distinct_sites)))
names(out_df) = c("site")
out_df$site = distinct_sites$sites
out_df$MAF_1pc = out_df$site %in% MAF1pc$V1
out_df$MAF_0.5pc = out_df$site %in% MAF0.5pc$V1
out_df$MAF_0.1pc = out_df$site %in% MAF0.1pc$V1
out_df$multiallelic = out_df$site %in% multiallelic

write.csv(out_df, "~/GitCode/MitoImputePrep/metadata/MitoImpite_sites_per_MAF/sites_per_MAF.csv", row.names = F, quote = F)
