library(tidyverse)
library(readr)

source("Data-cleaning/cover-1m2-cleaning.R")

# Read in the species data
JR_allcover <- JRm2cover %>%
  tbl_df() %>%
  mutate(species = as.character(species),
         species = ifelse(species == "TRAL" | species == "TRTR", "TRSP", species)) %>%
  group_by(year, species, treatment, trtrep, plot, uniqueID) %>%
  summarize(cover = sum(cover)) %>%
  tbl_df() %>%
  filter(species != "BARE", species != "ROCK") %>%
  mutate(plotyr = paste(uniqueID, year, spe = "_")) %>%
  spread(species, cover, fill = 0)


abundance <- as.matrix(vegdat[1:120, 11:25])
row.names(abundance) <- vegdat$plotyr

# traits <- as.matrix(traitdat[1:15, 2:11])
traits <- as.matrix(traitdat[1:15, c(2:ncol(traitdat))])
row.names(traits) <- traitdat$species

myrda <- rda(traits, scale = TRUE)
biplot(myrda)
cor(traits)