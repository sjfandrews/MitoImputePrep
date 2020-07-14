library(tidyverse)

hg_file_1000g  = "~/GitCode/MitoImputePrep/metadata/HaploGrep_concordance/chrMT_1kg_diploid_haplogrep.txt"
out_file_1000g = sub(x = hg_file, pattern = ".txt", "_macro.csv")

hgs_1000g = read_tsv(hg_file_1000g)

hgs_1000g_mut = hgs %>%
  mutate(MacroHaplogroup = if_else(str_detect(Haplogroup, "^L"),
                                  substr(Haplogroup, start = 1, stop = 2),
                                  substr(Haplogroup, start = 1, stop = 1))
         ) %>%
  write_csv(path = out_file_1000g)
hgs_1000g_mut

hg_file_refpan  = "~/GitCode/MitoImputePrep/DerivedData/MasterAlignment/McInerney_Master_Alignment_July18_2018_ambigANDgap2missing_HaploGrep.txt"
out_file_refpan = sub(x = hg_file_refpan, pattern = ".txt", "_macro.csv")

hg_refpan = read_tsv(hg_file_refpan)

hgs_refpan_mut = hg_refpan %>%
  mutate(MacroHaplogroup = if_else(str_detect(Haplogroup, "^L"),
                                   substr(Haplogroup, start = 1, stop = 2),
                                   substr(Haplogroup, start = 1, stop = 1))
  ) %>%
  write_csv(path = out_file_refpan)
hgs_refpan_mut
