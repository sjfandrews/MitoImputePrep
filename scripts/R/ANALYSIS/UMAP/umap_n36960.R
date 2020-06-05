library(umapr)
library(FactoMineR)
library(factoextra)
library(plotly)
library(HiMC)
library(tidyverse)
library(countrycode)

# A4 = 297 x 210
# PowerPoint widescreen = 252 x 116
# A4 = 297 x 210
# PowerPoint widescreen = 252 x 116
plot_multiplier = 1.5
wd = 250 * plot_multiplier
ht = 250 * plot_multiplier
DPI = 600

plot_dir = "~/Desktop/SANDBOX/MitoImpute/UMAP/plots/"

generate_snp_data_fixed <- function (map_file, ped_file) {
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

in_map_file = "~//Desktop/SANDBOX/MitoImpute/UMAP/data/ReferencePanel_v1_highQual_MAF0.01_filtered.map"
in_ped_file = "~//Desktop/SANDBOX/MitoImpute/UMAP/data/ReferencePanel_v1_highQual_MAF0.01_filtered.ped"

country_info_file = "~/GitCode/MitoImputePrep/metadata/seq_country_list.csv"
country_info = read_csv(country_info_file)

country_info$Region = gsub(".*:", "", country_info$Country)
country_info$Region = str_trim(country_info$Region)
country_info$Region[country_info$Country == country_info$Region] = NA
country_info$Country_short = gsub("\\:.*","", country_info$Country)
country_info$Continent = countrycode(sourcevar = country_info$Country_short, origin = "country.name", destination = "continent")

dat <- generate_snp_data_fixed(map_file = in_map_file, ped_file = in_ped_file)

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

ped$geographic_data = country_info$Country
ped$Continent       = country_info$Continent
ped$Country         = country_info$Country_short
ped$Region          = country_info$Region

ped = ped %>%
  select(Family, Individual, Father, Mother, Sex, Phenotype, haplogroup, macro, geographic_data, Continent, Country, Region, everything())

#write.table(ped, "~/GitCode/MitoImputePrep/DerivedData/ReferencePanel_v1_0.01/ReferencePanel_v1_highQual_MAF0.01_filtered_countryInfo.ped", row.names = F, quote = F, sep = "\t")

## PCA 
res.pca <- ped %>% 
  select(starts_with("mt")) %>% 
  PCA(., scale.unit = T, ncp = 10)

test <- get_pca_ind(res.pca)
ind.pca <- test$coord %>% 
  as_tibble() %>%
  bind_cols(select(ped, Individual, haplogroup, macro)) 

ggplot(ind.pca, aes(x = Dim.1, y = Dim.2, colour = macro)) +
  geom_point() +
  labs(x = "PCA1",
       y = "PCA2",
       colour = "Haplogroup") +
  guides(colour = guide_legend(title.position="top",
                               title.hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 0, size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.title.y = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.0)),
        legend.title = element_text(size = rel(1.25)),
        legend.position= "top",
        plot.title = element_blank(),
        strip.text.x = element_text(size = rel(5/5)),
        strip.text.y = element_text(size = rel(7/5)),
        panel.grid.major.x = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.grid.minor.x = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.grid.major.y = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.grid.minor.y = element_line(colour = "black", linetype = 2, size = rel(1/2)),
        panel.background = element_rect(fill = "transparent", colour = "black"),
        strip.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))
ggsave(paste0(plot_dir, "pca.png"), width = wd, height = ht, units = "mm", dpi = DPI)


ind.pca %>% 
  select(starts_with("Dim")) %>%
  write_tsv("~/Desktop/SANDBOX/MitoImpute/UMAP/pcs.txt")

embedding <- umap(as.matrix(ind.pca[, 1:10]), min_dist = 0.5, n_neighbor = 30, n_components = 3)
embedding2 <- embedding %>% as.tibble() %>%
  bind_cols(select(ped, Individual, haplogroup, macro)) 

ggplot(embedding2, aes(x = UMAP1, y = UMAP2, color = macro)) + geom_point() + theme_bw()

















# END!