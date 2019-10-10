HVRI  = c(1:577)
HVRII = c(16023:16569)

length(c(HVRI, HVRII))
rep("MT", length(c(HVRI, HVRII)))

df = data.frame(CHR = rep("MT", length(c(HVRI, HVRII))), SITE = c(HVRI, HVRII))

out_file = "/g/data1a/te53/MitoImpute/data/STRANDS/DLOOP_1-577_16023-16569/DLOOP_1-577_16023-16569_MT_snps.txt"

write.table(x = df,
            file = out_file,
            quote = F,
            row.names = F,
            col.names = F,
            sep = "\t")