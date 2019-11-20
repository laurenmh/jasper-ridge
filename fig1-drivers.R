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


cor.test(JR_rain$ppt, JR_rain$winter_ppt)       
cor.test(JR_rain$winter_ppt, JR_rain$fall_ppt)       

ggplot(JR_rain, aes(x=ppt, y=winter_ppt)) + geom_point()
ggplot(JR_rain, aes(x=ppt, y=fall_ppt)) + geom_point()
ggplot(JR_rain, aes(x=winter_ppt, y=fall_ppt)) + geom_point()

l <- lm(fall_ppt~year, data = JR_rain)
summary(l)

ggplot(JR_rain, aes(x=year, y=fall_ppt)) + geom_point() + geom_smooth(method = lm)


a <- ggplot(JR_rain, aes(x=year, y=ppt)) + geom_line() + 
  geom_point(shape = 17) +
  geom_line(data = JR_rain, aes(x=year, y=fall_ppt), lty = "dotted") +
  geom_point(data = JR_rain, aes(x=year, y=fall_ppt), pch = 1) +
  geom_line(data = JR_rain, aes(x=year, y=winter_ppt), lty = "dashed") + 
  geom_point(data = JR_rain, aes(x=year, y=winter_ppt), pch = 15) +
  labs(x= "Sampling year", y = "Precipitation (mm)") + 
  scale_y_continuous(breaks=seq(0,1200,200)) + 
  scale_x_continuous(breaks=seq(1980,2020,5)) 

JR_gopher_mean_time <- left_join(JR_gopher, JR_rain) %>%
  group_by( year, precip,  year, treatment) %>%
  summarize(meandisturb = mean(disturb), sedisturb = sd(disturb)/sqrt(length(disturb)))


ggplot(JR_gopher_mean_time, aes(x=year, y=meandisturb, lty = treatment)) + geom_line() 


JR_gopheroccur <- JRgopher %>%
  mutate(isdisturb = ifelse(disturb > 0, 1, 0)) %>%
  group_by(year, treatment, trtrep) %>%
  summarize(totdisturb = sum(isdisturb)) %>%
  mutate(treatment2 = as.factor(treatment),
         treatment2 = recode(treatment2, g = "Gopher excl", c = "Control", r = "Rabbit excl"),
         treatment2 = fct_relevel(treatment2, "Gopher excl","Control", "Rabbit excl"))


JR_gopheroccur2 <- JRgopher %>%
  mutate(isdisturb = ifelse(disturb > 0, 1, 0)) %>%
  group_by(quadID, treatment, trtrep) %>%
  summarize(totdisturb = sum(isdisturb)) %>%
  mutate(treatment2 = as.factor(treatment),
         treatment2 = recode(treatment2, g = "Gopher excl", c = "Control", r = "Rabbit excl"),
         treatment2 = fct_relevel(treatment2, "Gopher excl","Control", "Rabbit excl"))

ggplot(JR_gopheroccur, aes(x=year, y=totdisturb)) + geom_bar(stat = "identity") + facet_grid(treatment2 ~ trtrep)
ggplot(JR_gopheroccur, aes(x=year, y=totdisturb)) + geom_bar(stat = "identity") + facet_grid( trtrep ~ treatment2)

ggplot(JR_gopheroccur2, aes(x=totdisturb)) + geom_histogram() + facet_grid(treatment2 ~ trtrep)
ggplot(JR_gopheroccur2, aes(x=totdisturb))  + geom_histogram() + facet_grid( trtrep ~ treatment2)

ggplot(JR_gopher, aes(x=year, y=meandisturb, lty = treatment)) + geom_line() 

JRdrivers <- left_join(JR_gopheroccur2, JRsoil)
ggplot(JRdrivers, aes(x=depth, y=totdisturb)) + geom_point() + geom_smooth(method = lm) + facet_wrap(~treatment2)
ggplot(JRdrivers, aes(x=depth, y=totdisturb, color =treatment2)) + geom_point() + geom_smooth(method = lm) 

l <- lm(totdisturb ~ depth*treatment, data = JRdrivers)
summary(l)


ggplot(JR_gopher_mean, aes(x=soilDepth, y=meandisturb)) + geom_point() + geom_smooth(method = lm)  + facet_wrap(~treatment2)
l <- lm(meandisturb ~ soilDepth*treatment, data = JR_gopher_mean)
summary(l) 

