x = entrez_search("nuccore", term="EF184582.1")
y = entrez_fetch("nuccore", id = x$ids, rettype = "native")
z = str_split(y, "\n")
z = lapply(z, str_trim)
z = unlist(z)

m = which("subtype country,"==z)
z[m[1]+1] %>%
  str_remove(., 'name \"') %>%
  str_remove(., '\"')

countries = c()
seq_list = read.table("~/Desktop/seq_country_test.txt", header = F)
names(seq_list) = c("seqID")
seq_list$Country = NA

for (i in 1:nrow(seq_list)) {
  if (i %% 100 == 0) {
    message(paste0(i , "  /  ", nrow(seq_list))) 
  }
  tmp_seq = as.character(seq_list$seqID[i])
  x = entrez_search("nuccore", term = tmp_seq)
  y = entrez_fetch("nuccore", id = x$ids, rettype = "native")
  z = str_split(y, "\n")
  z = lapply(z, str_trim)
  z = unlist(z)
  
  m = which("subtype country,"==z)
  tmp_country = z[m[1]+1] %>%
    str_remove(., 'name \"') %>%
    str_remove(., '\"')
  
  seq_list$Country[i] = tmp_country
}

write.csv(seq_list, "~/GitCode/MitoImputePrep/metadata/seq_country_list.csv", row.names = F, quote = F)
