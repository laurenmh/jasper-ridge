##########################################
### Jasper Ridge Cover Data Formatting ###
##########################################

library(tidyverse)
library(readxl)


# Cover data by quadrat ---------------------------------------------------

## Cover data 
JRcover <- read_excel("~/Dropbox/California Data/percent cover database/percent cover vegetation Data 1983-2018.xlsx") %>%
  mutate(quadID = tolower(QuadratCode),
         species = SpeciesCode,
         year = SamplingYear,
         cover = PercentCover) %>%
  mutate(treatment=substr(quadID, 1, 1),
         trtrep=substr(quadID, 2,2),
         subplot=as.numeric(substr(quadID, 3, 4))) %>%
  select(-c(QuadratCode:PercentCover))

# write_csv(JRcover, "JR_cover.csv")


## Gopher data 
JRgopher1 <- read_csv("~/Dropbox/California Data/Gopher mound data/JR_gopher_1983_2015.csv") %>%
  select(-X1, -freq)

JRgopher2 <- read_excel("~/Dropbox/California Data/Gopher mound data/Gopher data 05-18.xlsx", sheet = 2) %>%
  mutate(quadID = tolower(QuadratCode),
         disturb = Numbersubquadratsdisturbed,
         year = SamplingYear) %>%
  select(-c(QuadratCode:Notes)) %>%
  filter(year > 2015 )

JRgopher <- rbind(JRgopher1, JRgopher2) %>%
  mutate(treatment=substr(quadID, 1, 1),
         trtrep=substr(quadID, 2,2),
         subplot=as.numeric(substr(quadID, 3, 4)))


# write_csv(JRgopher, "JR_gopher-disturbance.csv")


## Soil depth 
JRsoil <- read_csv("~/Dropbox/California Data/percent cover database/JR_soil_raw.csv") %>%
  mutate(quadID = tolower(QuadratCode)) %>%
  group_by(quadID) %>%
  summarize(depth = mean(SoilDepth)) %>%
  mutate(treatment=substr(quadID, 1, 1),
         trtrep=substr(quadID, 2,2),
         subplot=as.numeric(substr(quadID, 3, 4)))


write_csv(JRsoil, "JR_soil-depth.csv")



# By year -----------------------------------------------------------

JR_byyear <- JRcover %>%
  group_by(species, year) %>%
  summarize(cover = mean(cover))

# write_csv(JR_byyear, "JR_cover-by-year.csv")



# By plot -----------------------------------------------------------

JR_byplot <- JRcover %>%
  group_by(species, quadID) %>%
  summarize(cover = mean(cover))

# write_csv(JR_byplot, "JR_cover-by-plot.csv")

JRgopher_byplot <- JRgopher %>%
  group_by(quadID) %>%
  summarize(disturb = mean(disturb))

# write_csv(JRgopher, "JR_gopher-disturbance-by-plot.csv")

rm(JRgopher1, JRgopher2, JR_byplot, JR_byyear, JRgopher_byplot)
