require(rentrez)
require(tidyverse)

outFile = "~/GitCode/MitoImputePrep/metadata/seq_country_list.csv"

#x = entrez_search("nuccore", term="EF184582.1")
#y = entrez_fetch("nuccore", id = x$ids, rettype = "native")
#z = str_split(y, "\n")
#z = lapply(z, str_trim)
#z = unlist(z)
#
#m = which("subtype country,"==z)
#z[m[1]+1] %>%
#  str_remove(., 'name \"') %>%
#  str_remove(., '\"')
#
#countries = c()
seq_list = read.table("~/GitCode/MitoImputePrep/metadata/seq_country_test.txt", header = F)
names(seq_list) = c("seqID")
seq_list$Country = NA
seq_list$Checked = F

if (file.exists(outFile) == T) {
  seq_list = read.csv(outFile, header = T)
}

for (i in 1:nrow(seq_list)) { #nrow(seq_list)
  if (i %% 10 == 0) {
    message(paste0(i , "  /  ", nrow(seq_list))) 
  }
  
  if (seq_list$Checked[i] != T) {
    
    tmp_seq = as.character(seq_list$seqID[i])
    x = entrez_search("nuccore", term = tmp_seq)
    if (length(x$ids) > 0) {
      y = entrez_fetch("nuccore", id = x$ids, rettype = "native")
      z = str_split(y, "\n")
      z = lapply(z, str_trim)
      z = unlist(z)
      
      m = which("subtype country,"==z)
      tmp_country = z[m[1]+1] %>%
        str_remove(., 'name \"') %>%
        str_remove(., '\"') %>%
        str_replace(., pattern = ",", replacement = ";")
      
      seq_list$Country[i] = tmp_country
    }
    seq_list$Checked[i] = T
  }
  
  if (i %% 10 == 0) {
    write.csv(seq_list, outFile, row.names = F, quote = F)
  }
  
}

message(paste0("LAST i  =  ", i))

for (j in 1:nrow(seq_list)) {
  if (j %% 500 == 0) {
    message(paste0(j , "  /  ", nrow(seq_list))) 
  }
  if (seq_list$Checked[j] == T) {
    seq_list$Country[j] = str_replace(seq_list$Country[j], pattern = ",", replacement = ";")
  }
}

write.csv(seq_list, outFile, row.names = F, quote = F)

