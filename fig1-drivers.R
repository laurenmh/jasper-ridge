library(tidyverse)
library(corrplot)
library(lme4)
library(lmerTest)
library(nlme)
library(broom)
library(cowplot)

## Set ggplot2 theme
theme_set(theme_bw())
theme_update( panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              strip.background = element_blank(),
              text = element_text(size = 20),
              strip.text= element_text(size = 17), 
              axis.text = element_text(size = 14))

source("Data-cleaning/cover-1m2-cleaning.R")

JR_soil <- JRm2soil

JR_gopher <- JRm2gopher

source("prism-cleaning.R")

ggplot(JR_rain, aes(x=year, y=ppt)) + geom_line() + 
  geom_point(shape = 17) +
  geom_line(data = JR_rain, aes(x=year, y=fall_ppt), lty = "dotted") +
  geom_point(data = JR_rain, aes(x=year, y=fall_ppt), pch = 1) +
  geom_line(data = JR_rain, aes(x=year, y=winter_ppt), lty = "dashed") + 
  geom_point(data = JR_rain, aes(x=year, y=winter_ppt), pch = 15) +
  labs(x= "Sampling year", y = "Precipitation (mm)") + 
  scale_y_continuous(breaks=seq(0,1200,200)) + 
  scale_x_continuous(breaks=seq(1980,2020,5)) 

cor.test(JR_rain$ppt, JR_rain$winter_ppt)       
cor.test(JR_rain$winter_ppt, JR_rain$fall_ppt)       

ggplot(JR_rain, aes(x=ppt, y=winter_ppt)) + geom_point()
ggplot(JR_rain, aes(x=ppt, y=fall_ppt)) + geom_point()
ggplot(JR_rain, aes(x=winter_ppt, y=fall_ppt)) + geom_point()

ggplot(JR_gopher_mean_time, aes(x=year, y=meandisturb, lty = treatment)) + geom_line() 


ggplot(JR_gopher, aes(x=year, y=meandisturb, lty = treatment)) + geom_line() 
