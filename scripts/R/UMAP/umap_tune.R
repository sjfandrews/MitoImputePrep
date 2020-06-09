

library(FactoMineR)
library(factoextra)
library(uwot)
library(tidyverse)
library(magrittr)
library(glue)

sessionInfo()

infile = snakemake@input[["infile"]]
outfile = snakemake@output[["outfile"]]
pcs = snakemake@params[["pcs"]]
n_components = as.numeric(snakemake@params[["n_components"]])
min_dist = as.numeric(snakemake@params[["min_dist"]])
n_neighbor = as.numeric(snakemake@params[["n_neighbor"]])

message(glue("Infile: {infile} \n Outfile: {outfile} \n Running: npc = {pcs}, n_components = {n_components}, min_dist = {min_dist}, n_neighbor = {n_neighbor} \n"))

ped <- read_tsv(infile)

message("\n PCA away \n")

set.seed(333)
res.pca <- ped %>% 
  select(starts_with("mt")) %>% 
  PCA(., scale.unit = F, ncp = 10, graph = FALSE) 

name_pcs <- as.character(glue("Dim.{npc}", npc = 1:pcs))
tab.pca <- res.pca %>% 
  get_pca_ind(res.pca) %>% 
  use_series(coord) %>% 
  as_tibble() %>% 
  select(., name_pcs)

message("\n uwot mate \n")

set.seed(333)
res.umap <- umap(tab.pca, min_dist = min_dist, spread = 1, n_neighbor = n_neighbor, n_components = n_components, 
             init = "agspectral", metric = "euclidean", verbose = FALSE) 

message("\n u-munged-wot \n")

out <- res.umap %>% 
  as.tibble() %>%
  rename(UMAP1 = V1, UMAP2 = V2) %>% 
#  bind_cols(res.pca) %>%
  bind_cols(select(ped, Individual, haplogroup, macro, supclu, Continent)) %>%
  mutate(pcs = pcs, 
         n_components = n_components, 
         min_dist = min_dist, 
         n_neighbor = n_neighbor) %>% 
  select(pcs, n_components, min_dist, n_neighbor, Individual, haplogroup, macro, supclu, Continent, UMAP1, UMAP2)

write_tsv(out, outfile)
