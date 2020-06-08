library(umapr)
library("FactoMineR")
library("factoextra")
library(plotly)
library(HiMC)

generate_snp_data_fixed <- function (map_file, ped_file) 
{
  map <- read.csv(map_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  header_row <- c("Family", "Individual", "Father", "Mother", 
                  "Sex", "Phenotype")
  snps = map[, 4]
  new_header = c(header_row, snps)
  ped <- read.csv(ped_file, sep = " ", header = FALSE, stringsAsFactors = FALSE, colClasses = 'character')
  range1 = seq(1, 6, by = 1)
  snp_data = data.frame(seq(1, nrow(ped), by = 1))
  for (i in range1) {
    snp_data[, i] = ped[, i]
  }
  range2 = seq(7, ncol(ped), by = 2)
  for (i in range2) {
    index = ((i - 7)/2) + 1
    snp_data[, index + 6] = paste(ped[, i], ped[, i + 1])
  }
  names(snp_data) <- new_header
  return(snp_data)
}

## Haplogroup assignments 
dat <- generate_snp_data_fixed("/Users/sheaandrews/GitCode/MitoImpute/ReferencePanel/ReferencePanel.map",
                               "/Users/sheaandrews/GitCode/MitoImpute/ReferencePanel/ReferencePanel.ped")

data(nodes)
classifications <- HiMC::getClassifications(dat)
classifications <- as_tibble(classifications[-nrow(classifications),]) #remove Mitochondria_Information row 

ped <- dat
colnames(ped) <- c("Family", "Individual", "Father", "Mother", 
                   "Sex", "Phenotype", paste0('mt', colnames(ped[, 7:ncol(ped)])))
ped <- as_tibble(ped) %>%
  na_if(., "0 0") %>% 
  mutate_at(vars(starts_with("mt")), as_factor) %>% 
  mutate_at(vars(starts_with("mt")), as.numeric) %>% 
  left_join(classifications) %>% 
  mutate(macro = ifelse(str_detect(haplogroup, "^L"),
                        substr(haplogroup, start = 1, stop = 2),
                        substr(haplogroup, start = 1, stop = 1))) %>% 
  select(Family, Individual, Father, Mother, Sex, Phenotype, haplogroup, macro, everything())


## PCA 
res.pca <- ped %>% 
  select(starts_with("mt")) %>% 
  PCA(., scale.unit = T, ncp = 10)

test <- get_pca_ind(res.pca)
ind.pca <- test$coord %>% as.tibble() %>%
  bind_cols(select(ped, Individual, haplogroup, macro)) 

ggplot(ind.pca, aes(x = Dim.1, y = Dim.2, colour = macro)) + geom_point() + theme_bw()

ind.pca %>% 
  select(starts_with("Dim")) %>%
  write_tsv("/Users/sheaandrews/Documents/gt-dimred/data/pcs.txt")

embedding <- umap(as.matrix(ind.pca[, 1:10]), min_dist = 0.5, n_neighbor = 30, n_components = 3)
embedding2 <- embedding %>% as.tibble() %>%
  bind_cols(select(ped, Individual, haplogroup, macro)) 

ggplot(embedding2, aes(x = UMAP1, y = UMAP2, color = macro)) + geom_point() + theme_bw()

plot_ly(embedding2, x = ~UMAP2, y = ~UMAP1, z = ~UMAP3, color = ~macro, type = 'scatter3d', mode = 'markers', marker = list(size = 3), 
        hoverinfo = 'text', text = ~paste('</br> ID: ', Individual, 
                                          '</br> Haplogroup: ', macro,
                                          '</br> PC1: ', round(UMAP1, 2),
                                          '</br> PC2: ', round(UMAP2, 2),
                                          '</br> PC3: ', round(UMAP3, 2)))

custom.config = umap.defaults
custom.config$n_neighbors = 30
custom.config$min_dist = 0.5

pca.umap <- ind.pca %>% 
  select(starts_with("Dim")) %>% 
  umap(., custom.config) 

pca.umap$layout %>% 
  as_tibble() %>%
  bind_cols(select(ped, X1)) %>%
  left_join(classifications, by = c("X1" = "Individual")) %>%
  ggplot(., aes(x = V1, y = V2, color = haplogroup)) + geom_point()

pca_umap <- read_table2('/Users/sheaandrews/Documents/gt-dimred/data/pcs_UMAP_PC10_NC2_NN15_MD0.5_euclidean_20200602_182358.txt', col_names = F) %>% 
  bind_cols(select(ped, Individual, haplogroup, macro)) 

ggplot(pca_umap, aes(x = X1, y = X2, color = macro)) + geom_point() + theme_bw()

plot_ly(ind.pca, x = ~Dim.2, y = ~Dim.1, z = ~Dim.3, color = ~haplogroup, type = 'scatter3d', mode = 'markers', marker = list(size = 3), 
        hoverinfo = 'text', text = ~paste('</br> ID: ', X1, 
                                          '</br> Haplogroup: ', haplogroup,
                                          '</br> PC1: ', round(Dim.1, 2),
                                          '</br> PC2: ', round(Dim.2, 2),
                                          '</br> PC3: ', round(Dim.3, 2)))



## other UMAP 
ped.umap <- ped %>% 
  select(starts_with("mt")) %>% 
  umap::umap()

ped.umap$layout %>% 
  as_tibble() %>%
  bind_cols(select(ped, X1)) %>%
  left_join(classifications, by = c("X1" = "Individual")) %>%
  ggplot(., aes(x = V1, y = V2, color = haplogroup)) + geom_point()

