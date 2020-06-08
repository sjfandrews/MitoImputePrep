library(tidyverse)
library(countrycode)

f = "~/GitCode/MitoImputePrep/metadata/seq_country_list.csv"

seq_list = read_csv(f, col_names = T)

seq_list$Region = gsub(".*:", "", seq_list$Country)
seq_list$Region = str_trim(seq_list$Region)
seq_list$Region[seq_list$Country == seq_list$Region] = NA
seq_list$Country_short = gsub("\\:.*","", seq_list$Country)
seq_list$Continent = countrycode(sourcevar = seq_list$Country_short, origin = "country.name", destination = "continent")

region_table    = as.data.frame(table(seq_list$Country))
country_table   = as.data.frame(table(seq_list$Country_short))
continent_table = as.data.frame(table(seq_list$Continent))

as_tibble(region_table)
as_tibble(country_table)
as_tibble(continent_table)


table(is.na(seq_list$Region))

write.csv(region_table, "~/GitCode/MitoImputePrep/metadata/country_information/region_table.csv", quote = F, row.names = F)
write.csv(country_table, "~/GitCode/MitoImputePrep/metadata/country_information/country_table.csv", quote = F, row.names = F)
write.csv(continent_table, "~/GitCode/MitoImputePrep/metadata/country_information/continent_table.csv", quote = F, row.names = F)
