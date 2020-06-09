library(tidyverse)

umap_files <- list.files(path = "data/test/hyperparameters/", pattern = "umpap_*", full.names = TRUE)

umap_tabs <- umap_files %>%
  map(., read_tsv) %>%
  bind_rows() %>%
  group_by(umap_tabs, pcs, n_components, min_dist, n_neighbor) %>%
  nest() %>%
  arrange(pcs, min_dist, n_neighbor)

write_rds(umap_tabs, "data/test/umap_tune.rds.gz", compress = "gz")
