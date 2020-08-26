source("Data-cleaning/times-since-gopher.R")
source("prism-cleaning.R")

library(tsvr)
library(cowplot)
library(grid)

threshold <- 1
upperthreshold <- 5

# plstcheck <- tog %>%
#   # only keep a record if the cover is above a threshold value in the initial year
#   mutate(keepdat = ifelse(rowdif == 0 & cover > threshold, 1, 0)) %>%
#   filter(keepdat == 1) %>%
#   
#   # subset the focal species
#   filter(species == "PLER" | species == "SIJU") %>%
#   
#   # only keep a plot if both species are above the threshold
#   group_by(quadID, treatment, trtrep, subplot, repnum2) %>%
#   mutate(nospp = length(unique((species)))) %>%
#   filter(nospp == 2) %>%
#   
#   # restrict to ts with at least 8 recrods
#   filter(ncount > 9) %>%
#   
#   # identify the starting year for each record
#   group_by(species, quadID, treatment, trtrep, subplot, repnum2) %>%
#   summarize(year = min(year))
# 
# plstcheck2 <- plstcheck %>%
#   group_by(species, year) %>%
#   summarize(repplots = n()) %>%
#   tbl_df() %>%
#   spread(species, repplots)
# 
# plstyears <- plstcheck2 %>%
#   filter(repplots > 9 & repplots < 40)
# 
# myyears <- plstyears$year

## in replicate years with at least 10 data points

tog2 <- tog %>%
  spread(species, cover)

tog2b <- left_join(tog2, annual)

tog3 <- tog %>%
  tbl_df() %>%
  select(quadID, species, year, cover) %>%
  mutate(species = paste(species, "prev", sep = "")) %>%
  mutate(year = year + 1) %>%
  spread(species, cover)



brpl <- left_join(tog2b, tog3) %>% 
  filter(BRMOprev > 0, BRMOprev < 3, PLERprev < 5)

ggplot(brpl, aes(x = ppt, y = BRMO/BRMOprev)) + geom_point() 


plbr <- left_join(tog2b, tog3) %>% 
  filter(PLERprev > 0, PLERprev < 3, BRMOprev > 30)

ggplot(plbr, aes(x = ppt, y = PLER/PLERprev)) + geom_point() # + scale_y_log10()


plbr <- left_join(tog2b, tog3) %>% 
  filter(LACAprev > 0, LACAprev < 3, PLERprev > 30, BRMOprev < 20)

ggplot(plbr, aes(x = ppt, y = LACA/LACAprev)) + geom_point() # + scale_y_log10()



plst <- tog %>%
  # only keep a record if the cover is above a threshold value in the initial year
  mutate(keepdat = ifelse(rowdif == 0 & cover > threshold & cover < upperthreshold, 1, 0)) %>%
  group_by(repnum2, species, quadID) %>%
  mutate(maxkeep = max(keepdat)) %>%
  filter(maxkeep > 0) %>%
  
  # subset the focal species
  filter(species == "PLER" | species == "SIJU") %>%
  
  # only keep a plot if both species are above the threshold
  group_by(quadID, treatment, trtrep, subplot, repnum2) %>%
  mutate(nospp = length(unique((species)))) %>%
  filter(nospp == 2) %>%
  
  # restrict to ts with at least 8 recrods
  filter(ncount > 9) %>%
  
  #identify the first year
  group_by(species, quadID, treatment, trtrep, subplot, repnum2) %>%
  mutate(minyear = min(year))



plst$uniqueID <- paste(plst$quadID, plst$minyear)

ggplot(subset(plst, rowdif < 10), aes(x=rowdif, y=cover, color = species, group = interaction(species,quadID))) + 
  geom_line() + facet_wrap(~interaction(minyear,quadID))

