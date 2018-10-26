#################################################
### Jasper Ridge Cover at 1m2 Data Formatting ###
#################################################

library(tidyverse)
library(readxl)

# Source the quadrat-level data
source("Data-cleaning/cover-cleaning.R")

# Key for 1m2 plots -------------------------------------------------------

# Separate the unique .5 m x .5 m plot name into its components
# treatment: c (control), g (gopher), r (rabbit)
# trtrep: the replicate set of treatments, 1:3
# subplot: the subplots within a treatment and trtrep, 1:24 (1:12 is the first block, 13:24 the second block)


# create a key to aggregate 0.5 m x 0.5 m subplots into 1 m x 1 m plots
top_block1 <- cbind(c(1,2,7,8), rep("block1_top", 4))
top_block2 <- cbind(c(13, 14,19, 20), rep("block2_top", 4))

middle_block1 <- cbind(c(3,4,9,10), rep("block1_middle", 4))
middle_block2 <- cbind(c(15, 16, 21, 22), rep("block2_middle", 4))

bottom_block1 <- cbind(c(5, 6, 11, 12), rep("block1_bottom", 4))
bottom_block2 <- cbind(c(17,18, 23, 24 ), rep("block2_bottom", 4))

key <-as.data.frame(rbind(top_block1, top_block2, 
                          middle_block1, middle_block2, 
                          bottom_block1, bottom_block2))
names(key)=c("subplot", "plot")

key <- key %>%
  tbl_df() %>%
  mutate(subplot=as.numeric(as.character(subplot)), 
         plot=as.character(plot))

rm(top_block1, top_block2, middle_block1, middle_block2, bottom_block1, bottom_block2)


# Cover data --------------------------------------------------------------

# merge JRcover and key; aggregate at the 1 m x 1 m plot scale
JRm2cover <- merge(JRcover, key) %>%
  tbl_df() %>%
  group_by(year, species, treatment, trtrep, plot) %>%
  summarize(cover=mean(cover)) %>%
  mutate(uniqueID=paste(treatment, trtrep, plot, sep="_")) %>%
  #removing the middle block for independence
  filter(plot!="block1_middle", plot!="block2_middle")

# write_csv(JRm2cover, "JR_cover_1mplot.csv")



# Soil depth data ---------------------------------------------------------

# Read in and merge JRsoil and key; aggregate at the 1m x 1m plot scale

# Merge JRcover and key; aggregate at the 1 m x 1 m plot scale
JRm2soil <- merge(JRsoil, key) %>%
  tbl_df() %>%
  group_by(treatment, trtrep, subplot, plot) %>%
  summarize(soil=mean(depth)) %>%
  mutate(uniqueID=paste(treatment, trtrep, plot, sep="_")) %>%
  tbl_df() %>%
  group_by(treatment, trtrep, plot) %>%
  summarize(soilDepth=mean(soil), 
            minSoilDepth=min(soil), 
            maxSoilDepth=max(soil)) %>%
  #removing the middle block for independence
  filter(plot != "block1_middle", plot != "block2_middle") 

## visualize relationships between soil categories
# ggplot(JRm2soil, aes(soilDepth, minSoilDepth)) + geom_point()
# ggplot(JRm2soil, aes(soilDepth, maxSoilDepth)) + geom_point()
# ggplot(JRm2soil, aes(maxSoilDepth, minSoilDepth)) + geom_point()


# write_csv(JRm2soil, "JR_soil_1mplot.csv")


# Gopher data -------------------------------------------------------------


# Merge JRgopher and key; aggregate at the 1 m x 1m plot scale

JRm2gopher <- merge(JRgopher, key) %>%
  tbl_df() %>%
  group_by(year, treatment, trtrep, plot) %>%
  summarize(disturb=mean(disturb)) %>%
  mutate(uniqueID=paste(treatment, trtrep, plot, sep="_")) %>%
  #removing the middle block for independence
  filter(plot!="block1_middle", plot!="block2_middle")

# write_csv(JRm2gopher, "JR_gopher_1mplot.csv")



# Zero check --------------------------------------------------------------


#look at zeros in gophers
ggplot(JRm2gopher, aes(x=disturb)) + geom_bar() + facet_wrap(~treatment)

#look at zeros in plants
JRcheck <- JR_aggcover %>%
  mutate(pres=ifelse(cover>0, "Pres", "Absent")) %>%
  mutate(count=1) %>%
  group_by(species, pres) %>%
  summarize(count=sum(count)) %>%
  #ggplot(JRcheck, aes(x=pres, fill=pres, y=count)) + geom_bar(stat="identity") + facet_wrap(~species)
  group_by(species) %>%
  mutate(totcount=sum(count), perpresent=count/totcount) %>%
  filter(pres=="Pres") %>%
  tbl_df() %>%
  mutate(perpresent=as.numeric(as.character(perpresent))) %>%
  arrange(-perpresent)

