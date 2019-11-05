library(readxl)
library(tidyverse)

supp_tables = read_xlsx("~/GitCode/MitoImputePrep/supplementary_information/MitoImpute_SupplementaryTables_abridged.xlsx", sheet = 1, skip = 1)
countries   = read_csv("~/GitCode/MitoImputePrep/metadata/seq_country_list.csv")
haplogroups = read_csv("/Volumes/TimMcInerney/MitoImpute/data/VCF/McInerney_Master_Alignment_July18_2018_ambigANDgap2missing_tagFilled_noInvar_HaploGrep.csv")

new_t1 = supp_tables

for (i in 1:nrow(supp_tables)) {
  cat(i, "  /   ", nrow(supp_tables), "\r")
  new_t1$Country[i] = countries$Country[match(supp_tables$Sequence[i], countries$seqID)]
  new_t1$Haplogroup[i] = haplogroups$Macrohaplogroup[match(supp_tables$Sequence[i], haplogroups$SampleID)]
}

new_t1
new_t2 = new_t1 %>% separate(Country, c("Country", "Region"), sep = ": ")
new_t2
write.csv(new_t2, "~/GitCode/MitoImputePrep/metadata/seq_country_list_44299.csv", row.names = F, quote = F)

new_t1 %>% filter(!is.na(Country))

new_t2 = new_t1[1743:1750,]
new_t2

new_t2 %>% separate(Country, c("Cnty", "Region"))
