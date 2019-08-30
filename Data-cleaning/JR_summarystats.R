library(tidyr)
library(dplyr)
library(ggplot2)
library(corrplot)

## Set ggplot2 theme
theme_set(theme_bw())
theme_update( panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              strip.background = element_blank(),
              text = element_text(size = 20),
              strip.text= element_text(size = 17), 
              axis.text = element_text(size = 14))

source("Data-cleaning/cover-1m2-cleaning.R")

# Read in the species data
JR_aggcover <- JRm2cover %>%
  tbl_df() %>%
  mutate(species = as.character(species),
         species = ifelse(species == "TRAL" | species == "TRTR", "TRSP", species)) %>%
  group_by(year, species, treatment, trtrep, plot, uniqueID) %>%
  summarize(cover = sum(cover)) %>%
  tbl_df() %>%
  filter(species != "BARE", species != "ROCK") 

JR_meanall <- JR_aggcover %>%
  group_by(species) %>%
  summarize(meancover = mean(cover), maxcover = max(cover))

JR_soil <- JRm2soil

JR_gopher <- JRm2gopher

JR_rain <- read_csv("~/Dropbox/California Data/jrg_prism.csv") %>%
  select(-X1) %>%
  mutate(precip = ppt)


JR_spp <- read_csv("~/Dropbox/California Data/JR_species names.csv") %>%
  mutate(species = SpeciesCode)

JR_gopher_mean <- left_join(JR_gopher, JR_soil) %>%
  group_by(soilDepth, maxSoilDepth, minSoilDepth, plot, treatment, trtrep) %>%
  summarize(meandisturb = mean(disturb), sedisturb = sd(disturb)/sqrt(length(disturb)))


JR_tog0 <- left_join(JR_aggcover, JR_soil)
JR_tog1 <- left_join(JR_tog0, JR_rain)
JR_tog2 <- left_join(JR_tog1, JR_spp)
JR_tog3 <- left_join(JR_tog2, JR_gopher) 
JR_tog <- left_join(JR_tog3, JR_tog_prev)

JR_togmean <- left_join(JR_tog, JR_mean)

l <- lme(cover ~ treatment + soilDepth + disturb + ppt + prevprecip, random = ~1|trtrep/plot/uniqueID/year, 
         data = subset(JR_tog, species == "PLER"), na.action = na.omit)
summary(l)

library(lme4)
library(lmerTest)

l <- lmer(cover ~ treatment + soilDepth + disturb + ppt + prevprecip + (1|trtrep/uniqueID) + (1|year), 
          data = subset(JR_tog, species == "PLER"), na.action = na.omit)

l <- lmer(cover ~ soilDepth + disturb + ppt + (1|trtrep/uniqueID) + (1|year), 
          data = subset(JR_tog, species == "PLER"), na.action = na.omit)

summary(l)

dat <- JR_togmean %>%
  filter(species == "STPU") %>%
  filter(!is.na(cover), !is.na(precip), !is.na(prevprecip))

l <- lmer(cover ~  disturb + soilDepth+ precip + prevprecip  + meandisturb + 
            (1|trtrep/uniqueID) + 
            (1|year), 
          data = dat, na.action = na.omit)

summary(l)
step_res <- step(l)
final <- get_model(step_res)
summary(final)

anova(final)

summary(l)



a <-ggplot(JR_gopher_mean_time, aes(x=year, y=meandisturb/4*100)) +
  geom_point(aes(color = treatment, shape = treatment, 
                 group = treatment), size = 4) + geom_line(aes(color = treatment, 
                                                               group = treatment),
                                                           lwd = 1.1) + 
#  geom_errorbar(aes(ymin = meandisturb/4*100 - sedisturb/4*100,
#                    ymax = meandisturb/4*100 + sedisturb/4*100)) + 
  scale_color_manual(values = c("grey45",  "grey65",  "grey10")) + 
  theme(legend.position = "none") + labs(x="Year", y = "Percent disturbance") + 
  geom_smooth()

b <- ggplot(JR_rain, aes(x=year, y=precip)) + geom_line() + geom_smooth()

plot_grid(a,b, ncol = 1)

ggplot(JR_gopher_mean, aes(x=soilDepth, y=meandisturb/4*100)) +
  geom_smooth(se=F, method = "lm", color = "grey") + geom_point() +
  labs(x = "Soil depth (cm)", y = "Percent disturbance")


gophrain <- left_join(JR_gopher_mean_time, JR_tog_prev %>% select(year, prevprecip) %>%unique())


ggplot(gophrain, aes(x = precip, meandisturb, color = treatment)) + geom_point() +
  labs(x = "Rainfall (mm)", y = "Percent disturbance") + 
  geom_smooth(se = F, method = lm) + facet_wrap(~treatment)

ggplot(gophrain, aes(x = prevprecip, meandisturb, color = treatment)) + geom_point() +
  labs(x = "Rainfall (mm)", y = "Percent disturbance") + 
  geom_smooth(se = F, method = lm) + facet_wrap(~treatment)

l <- lm(meandisturb~ precip + prevprecip + treatment, data = gophrain)
summary(l)
lm1 <- lm(meandisturb~ precip + treatment + prevprecip, data = gophrain)
summary(lm1)
slm1 <- step(lm1)
summary(slm1)

JR_gopher2 <- left_join(JR_gopher, JR_rain)
JR_gopher3 <- left_join(JR_gopher2, JR_tog_prev %>% select(year, prevprecip))
JR_gopher4 <- left_join(JR_gopher3, JR_soil)

l <- lmer(disturb ~  soilDepth+ precip + prevprecip  +
            (1|trtrep/uniqueID) + 
            (1|year), 
          data = dat, na.action = na.omit)

summary(l)
step_res <- step(l)
final <- get_model(step_res)
summary(final)



JR_mean0 <- JR_tog %>%
  group_by(Genus, Species, LifeForm, species, soilDepth, maxSoilDepth, minSoilDepth, plot, treatment, trtrep) %>%
  summarize(meancover = mean(cover), secover = sd(cover)/sqrt(length(cover)),
            maxcover = max(cover)) 
  
JR_mean <- left_join(JR_mean0, JR_gopher_mean)

JR_tog_prev <- JR_tog %>%
  select(year, species, uniqueID, cover, precip) %>%
  mutate(year = year + 1) %>%
  mutate(prevcover = cover, prevprecip = precip) %>%
  select(-cover, - precip)

JR_mean_time_0 <- left_join(JR_tog, JR_tog_prev) %>%
  mutate(coverchange = (cover - prevcover)/prevcover) %>%
  group_by(Genus, Species, LifeForm, species, year, precip,  year, treatment, prevprecip) %>%
  summarize(meancover = mean(cover), secover = sd(cover)/sqrt(length(cover)),
            meancoverchange = mean(coverchange), secoverchange = sd(coverchange)/sqrt(length(coverchange)),
            meanprevcover = mean(prevcover))
  

JR_gopher_mean_time <- left_join(JR_gopher, JR_rain) %>%
  group_by( year, precip,  year, treatment) %>%
  summarize(meandisturb = mean(disturb), sedisturb = sd(disturb)/sqrt(length(disturb)))


JR_mean_time <- left_join(JR_mean_time_0, JR_gopher_mean_time)

ggplot(JR_mean_time, aes(x=meandisturb, y = meancover, color = treatment)) + geom_point() + facet_wrap(~species, scale = "free")



JR_x1 <- left_join(JR_tog, JR_tog_prev) %>%
  filter(species == "LOSU" | species == "ASGA" | species == "TRSP" | species == "BRHO") 

ggplot(subset(JR_tog, species == "BRMO"), aes(x=year, y =cover, color = soilDepth, group = uniqueID)) +
         geom_line(lwd=1.5) + theme_bw() + facet_grid(trtrep~treatment)

ggplot(subset(JR_tog, species == "BRMO"), aes(x=year, y =cover,  group = uniqueID)) +
  geom_line(lwd=1) + theme_bw() + facet_grid(trtrep~treatment) +
  geom_point(aes(color = disturb), size = 3) 


### CORRELATIONS ####
JR_wide <- JR_mean %>%
  group_by(species) %>%
  mutate(meanmeancov = mean(meancover)) %>%
  filter(meanmeancov > 0.47) %>%
    select(-secover, -maxcover, -meandisturb, -sedisturb, -meanmeancov) %>%
  # filter(species == "ASGA", species == "BRHO" | species == "BRSP" | species == )
  spread(species, meancover, fill = 0)

cormat<-cor(JR_wide[,7:ncol(JR_wide)])
rowSums(cormat)


##VISUALIZE THE CORREGRAM
corrplot(cormat,method="color",type="upper", tl.col="black", tl.cex = .8, diag=F)


cormatvis <- cormat %>%
  as.data.frame() %>%
  mutate(comp = row.names(cormat)) %>%
  gather(species, correlation, AGHE:VUMI) %>%
  filter(correlation != 1) %>%
  group_by(species) %>%
  mutate(meancor = mean(correlation)) 

cormatvis$species1 <- factor(cormatvis$species, levels = unique(cormatvis$species)[order(cormatvis$meancor)])

ggplot(cormatvis, aes(x=species1, y=correlation, color = species)) + geom_boxplot() +
  theme_bw() + geom_jitter(aes(color = comp),
   width = .01, height = 0, size = 2) 


ggplot(cormatvis, aes(x=species1, y=correlation)) + geom_jitter(aes(color = comp),
                                                                width = .01, height = 0, size = 2) +
  theme_bw() #+ geom_hline(yintercept = 0, color = "blue")



### CORRELATIONS OVER TIME ####
JR_wide_time <- JR_mean_time %>%
  group_by(species, year) %>%
  summarize(meancover = mean(meancover)) %>%
  tbl_df() %>%
  group_by(species) %>%
  mutate(meanmeancov = mean(meancover)) %>%
  filter(meanmeancov > 0.47) %>%
  select(species, year,  meancover) %>%
  # filter(species == "ASGA", species == "BRHO" | species == "BRSP" | species == )
  spread(species, meancover, fill = 0)

cormat_time<-cor(JR_wide_time[,2:ncol(JR_wide_time)])
rowSums(cormat_time)


##VISUALIZE THE CORREGRAM
corrplot(cormat_time,method="color",type="upper", tl.col="black", tl.cex = .8, diag=F)

cormatvis_time <- cormat_time %>%
  as.data.frame() %>%
  mutate(comp = row.names(cormat)) %>%
  gather(species, correlation, AGHE:VUMI) %>%
  filter(correlation != 1) %>%
  group_by(species) %>%
  mutate(meancor = mean(correlation)) 

cormatvis_time$species1 <- factor(cormatvis_time$species, levels = unique(cormatvis_time$species)[order(cormatvis_time$meancor)])

ggplot(cormatvis_time, aes(x=species1, y=correlation, color = species)) + geom_boxplot() +
  theme_bw() + geom_jitter(aes(color = comp),
                           width = .01, height = 0, size = 2) 


ggplot(cormatvis, aes(x=species1, y=correlation)) + geom_jitter(aes(color = comp),
                                                                width = .01, height = 0, size = 2) +
  theme_bw() #+ geom_hline(yintercept = 0, color = "blue")


#### FOCAL SPECIES ##


JR_mean_focal <- JR_mean %>%
  tbl_df() %>%
  filter(species%in%c("CAMU", "HELU", "LACA", "BRMO", "LOMU", "EVSP", "PLER", "BRSP", 
                      "VUMY", "LOSU", "SIJU", "STPU", "LAPL", "ORDE", "VUMI",
                       "CHPO", "BR__", "TIER", "MIDO")) %>%
  mutate(Genus = fct_relevel(Genus, c("Plantago", "Lasthenia", "Microseris", 
                                      "Vulpia","Bromus",  "Lolium", "Sitanion", "Stipa",
                                    "Calycadenia","Hemizonia",  "Lotus", "Orthocarpus",  "Layia",
                                    "Chlorogalum", "Brodiaea", "Evax", "Tillaea"))) %>%
  tbl_df() %>%
  mutate(dummy = "")

ggplot(JR_mean_focal, aes(x=treatment, y=meancover)) + facet_wrap(~Genus, scales = "free") + 
  geom_boxplot()






a <- ggplot(JR_mean_time_focal, aes(x=year, y=meancover, shape = treatment, lty = treatment)) + geom_line() + 
  geom_point() + facet_grid(Genus~dummy, scale = "free", switch = "y") + theme(legend.position = "none") + 
  labs(x="Year", y= "Average cover m2") + theme(strip.placement = "outside") 


c <- ggplot(JR_mean_focal, aes(x=meandisturb/4*100, y=meancover)) + facet_grid(Genus~dummy, scales = "free") + 
  geom_point() + 
  # geom_smooth(data = subset(JR_mean_focal, species == "BRMO"), aes(meandisturb, meancover),
  #             method = lm, se = FALSE) + 
  geom_smooth(data = subset(JR_mean_focal, species == "CAMU"), aes(meandisturb/4*100, meancover),
              method = lm, se = FALSE) + 
  geom_smooth(data = subset(JR_mean_focal, species == "LACA"), aes(meandisturb/4*100, meancover),
              method = lm, se = FALSE) + 
  geom_smooth(data = subset(JR_mean_focal, species == "LAPL"), aes(meandisturb/4*100, meancover),
              method = lm, se = FALSE) + 
  geom_smooth(data = subset(JR_mean_focal, species == "LOSU"), aes(meandisturb/4*100, meancover),
              method = lm, se = FALSE) + 
  geom_smooth(data = subset(JR_mean_focal, species == "PLER"), aes(meandisturb/4*100, meancover),
              method = lm, se = FALSE) + 
  geom_smooth(data = subset(JR_mean_focal, species == "BR__"), aes(meandisturb/4*100, meancover),
              method = lm, se = FALSE)+ 
  geom_smooth(data = subset(JR_mean_focal, species == "ORDE"), aes(meandisturb/4*100, meancover),
              method = lm, se = FALSE) + theme(strip.background = element_blank(), strip.text = element_blank()) + 
  labs(x="Percent disturbance over time", y = "")

l <- lm(meancover ~ meandisturb  +soilDepth, data = subset(JR_mean_focal, species == "TIER"))
summary(l)



b <- ggplot(JR_mean_focal, aes(x=soilDepth, y=meancover)) + facet_grid(Genus~dummy, scales = "free") + 
  geom_point() + 
#  geom_smooth(data = subset(JR_mean_focal, species == "BRMO"), aes(soilDepth, meancover),
 #             method = lm, se = FALSE) + 
  geom_smooth(data = subset(JR_mean_focal, species == "CAMU"), aes(soilDepth, meancover),
              method = lm, se = FALSE) + 
  geom_smooth(data = subset(JR_mean_focal, species == "LOSU"), aes(soilDepth, meancover),
              method = lm, se = FALSE) + 
  geom_smooth(data = subset(JR_mean_focal, species == "PLER"), aes(soilDepth, meancover),
              method = lm, se = FALSE)+ 
  geom_smooth(data = subset(JR_mean_focal, species == "EVSP"), aes(soilDepth, meancover),
              method = lm, se = FALSE)+ 
  geom_smooth(data = subset(JR_mean_focal, species == "VUMI"), aes(soilDepth, meancover),
              method = lm, se = FALSE)+ 
  geom_smooth(data = subset(JR_mean_focal, species == "CHPO"), aes(soilDepth, meancover),
              method = lm, se = FALSE)+ 
  geom_smooth(data = subset(JR_mean_focal, species == "STPU"), aes(soilDepth, meancover),
              method = lm, se = FALSE)+ 
  geom_smooth(data = subset(JR_mean_focal, species == "TIER"), aes(soilDepth, meancover),
              method = lm, se = FALSE)+ 
  geom_smooth(data = subset(JR_mean_focal, species == "HELU"), aes(soilDepth, meancover),
              method = lm, se = FALSE)+ 
  geom_smooth(data = subset(JR_mean_focal, species == "SIJU"), aes(soilDepth, meancover),
              method = lm, se = FALSE) + theme(strip.background = element_blank(), strip.text = element_blank()) + 
  labs(x="Soil depth (cm)", y = "")


l <- lm(meancover ~  meandisturb + soilDepth, data = subset(JR_mean_focal, species == "ORDE"))
summary(l)



JR_mean_time_focal <- JR_mean_time %>%
  filter(species%in%c("CAMU", "HELU", "LACA", "BRMO", "LOMU", "EVSP", "PLER", "BRSP", 
                      "VUMY", "LOSU", "SIJU", "STPU", "LAPL", "ORDE", "VUMI",
                      "CHPO", "BR__", "TIER", "MIDO"))  %>%
  tbl_df() %>%
  mutate(Genus = ifelse(Genus == "Orthocarpus", "Castilleja", Genus),
         Genus = ifelse(Genus == "Sitanion", "Elymus", Genus),
         Genus = ifelse(Genus == "Tillaea", "Crassula", Genus),
         Genus = ifelse(Genus == "Lotus", "Acmispon", Genus)) %>%
  mutate(Genus = fct_relevel(Genus, c("Plantago", "Lasthenia", "Microseris", 
                                      "Vulpia", "Bromus",  "Lolium", "Elymus", "Stipa",
                                      "Calycadenia","Hemizonia",  "Acmispon", "Castilleja",  "Layia",
                                      "Chlorogalum", "Brodiaea", "Evax", "Crassula"))) %>%
  tbl_df() %>%
  mutate(dummy = "")

JR_mean_time_focal2 <- JR_mean_time %>%
  filter(species%in%c("CAMU", "HELU", "LACA", "BRMO", "LOMU", "EVSP", "PLER", "BRSP", 
                      "VUMY", "LOSU", "SIJU", "STPU", "LAPL", "ORDE", "VUMI",
                      "CHPO", "BR__", "TIER", "MIDO")) %>%
  group_by(Genus, Species, LifeForm, year, species, precip) %>%
  summarize(meancover = mean(meancover)) %>%
  tbl_df() %>%
  mutate(Genus = fct_relevel(Genus, c("Plantago", "Lasthenia", "Microseris", 
                                      "Vulpia", "Bromus", "Lolium", "Sitanion", "Stipa",
                                      "Calycadenia", "Hemizonia",  "Lotus", "Orthocarpus",  "Layia",
                                      "Chlorogalum", "Brodiaea", "Evax", "Tillaea"))) %>%
  tbl_df() %>%
  mutate(dummy = "")




ggplot(JR_mean_time_focal, aes(x=precip, y=meancover)) + facet_wrap(~species, scales = "free") + 
  geom_point() + geom_smooth(se = F)

ggplot(JR_mean_time_focal, aes(x=year, y=meancover)) + geom_line(aes(shape = treatment, lty = treatment)) + 
  geom_point(aes(shape = treatment, lty = treatment)) + facet_wrap(~Genus, scale = "free") + geom_smooth(se=F, method = "lm")

ggplot(JR_mean_time_focal2, aes(x=precip, y=meancover, color)) + facet_grid(Genus~dummy, scales = "free") + 
  geom_point() + 
  # geom_smooth(data = subset(JR_mean_time_focal, species == "BRMO"), aes(precip, meancover),
  #             method = lm, se = FALSE) + 
  # geom_smooth(data = subset(JR_mean_time_focal2, species == "CAMU"), aes(precip, meancover),
  #             method = lm, se = FALSE) + 
  # geom_smooth(data = subset(JR_mean_time_focal2, species == "EVSP"), aes(precip, meancover),
  #             method = lm, se = FALSE) + 
  geom_smooth(data = subset(JR_mean_time_focal2, species == "LACA"), aes(precip, meancover),
              method = lm, se = FALSE) + 
  # geom_smooth(data = subset(JR_mean_time_focal2, species == "LAPL"), aes(precip, meancover),
  #             method = lm, se = FALSE) + 
  geom_smooth(data = subset(JR_mean_time_focal2, species == "SIJU"), aes(precip, meancover),
              method = lm, se = FALSE) + 
  geom_smooth(data = subset(JR_mean_time_focal2, species == "TIER"), aes(precip, meancover),
              method = lm, se = FALSE) #+
  # geom_smooth(data = subset(JR_mean_time_focal2, species == "VUMI"), aes(precip, meancover),
  #                                                   method = lm, se = FALSE) +
  # geom_smooth(data = subset(JR_mean_time_focal2, species == "CHPO"), aes(precip, meancover),
  #             method = lm, se = FALSE) + theme(strip.background = element_blank(), strip.text = element_blank()) + 
labs(x="Soil depth (cm)", y = "")


l <- lm(meancover ~ precip, data = subset(JR_mean_time_focal2, species == "LAPL"))
summary(l)




#### THE BIG PANEL FIG

## Set ggplot2 theme
theme_set(theme_bw())
theme_update( panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              strip.background = element_blank(),
              text = element_text(size = 20),
              strip.text= element_text(size = 16), 
              axis.text = element_text(size = 14))

a <- ggplot(JR_mean_time_focal, aes(x=year, y=meancover, lty = treatment, shape = treatment, color = treatment)) + geom_line() + 
  geom_point() + facet_grid(Genus~dummy, scale = "free", switch = "y") +
  labs(x="Year", y = expression(paste("Percent cover (m"^"2",")"))) + theme(strip.placement = "outside", strip.text = element_text(face = "italic")) +
  scale_color_manual(values = c("grey45",  "grey65",  "grey10")) + theme(legend.position = "none") 

l <- lm(meancover ~ year, data = subset(JR_mean_time_focal, species == "VUMI"))
summary(l)

b <- ggplot(JR_mean_focal, aes(x=soilDepth, y=meancover)) + facet_grid(Genus~dummy, scales = "free") + 
  #  geom_smooth(data = subset(JR_mean_focal, species == "BRMO"), aes(soilDepth, meancover),
  #             method = lm, se = FALSE) + 
  geom_smooth(data = subset(JR_mean_focal, species == "CAMU"), color = "lightgrey", lwd = .75, aes(soilDepth, meancover),
              method = lm, se = FALSE) + 
  geom_smooth(data = subset(JR_mean_focal, species == "LOSU"), color = "lightgrey", lwd = .75, aes(soilDepth, meancover),
              method = lm, se = FALSE) + 
  geom_smooth(data = subset(JR_mean_focal, species == "PLER"), color = "lightgrey", lwd = .75, aes(soilDepth, meancover),
              method = lm, se = FALSE)+ 
  geom_smooth(data = subset(JR_mean_focal, species == "EVSP"), color = "lightgrey", lwd = .75, aes(soilDepth, meancover),
              method = lm, se = FALSE)+ 
  geom_smooth(data = subset(JR_mean_focal, species == "VUMI"), color = "lightgrey", lwd = .75, aes(soilDepth, meancover),
              method = lm, se = FALSE)+ 
  geom_smooth(data = subset(JR_mean_focal, species == "CHPO"), color = "lightgrey", lwd = .75, aes(soilDepth, meancover),
              method = lm, se = FALSE)+ 
  geom_smooth(data = subset(JR_mean_focal, species == "STPU"), color = "lightgrey", lwd = .75, aes(soilDepth, meancover),
              method = lm, se = FALSE)+ 
  geom_smooth(data = subset(JR_mean_focal, species == "TIER"), color = "lightgrey", lwd = .75, aes(soilDepth, meancover),
              method = lm, se = FALSE)+ 
  geom_smooth(data = subset(JR_mean_focal, species == "HELU"), color = "lightgrey", lwd = .75, aes(soilDepth, meancover),
              method = lm, se = FALSE)+ 
  geom_smooth(data = subset(JR_mean_focal, species == "SIJU"), color = "lightgrey", lwd = .75, aes(soilDepth, meancover),
              method = lm, se = FALSE) + theme(strip.background = element_blank(), strip.text = element_blank()) + 
  geom_point() + 
  labs(x="Soil depth (cm)", y = "")


c <- ggplot(JR_mean_focal, aes(x=meandisturb/4*100, y=meancover)) + facet_grid(Genus~dummy, scales = "free") + 
 
  # geom_smooth(data = subset(JR_mean_focal, species == "BRMO"), aes(meandisturb, meancover),
  #             method = lm, se = FALSE) + 
  # geom_smooth(data = subset(JR_mean_focal, species == "CAMU"), color = "lightgrey", lwd = .75, aes(meandisturb/4*100, meancover),
  #             method = lm, se = FALSE) + 
  geom_smooth(data = subset(JR_mean_focal, species == "LACA"), color = "lightgrey", lwd = .75, aes(meandisturb/4*100, meancover),
              method = lm, se = FALSE) + 
  geom_smooth(data = subset(JR_mean_focal, species == "LAPL"), color = "lightgrey", lwd = .75, aes(meandisturb/4*100, meancover),
              method = lm, se = FALSE) + 
  # geom_smooth(data = subset(JR_mean_focal, species == "LOSU"), color = "lightgrey", lwd = .75, aes(meandisturb/4*100, meancover),
  #             method = lm, se = FALSE) + 
  # geom_smooth(data = subset(JR_mean_focal, species == "PLER"), color = "lightgrey", lwd = .75, aes(meandisturb/4*100, meancover),
  #             method = lm, se = FALSE) + 
  geom_smooth(data = subset(JR_mean_focal, species == "BR__"), color = "lightgrey", lwd = .75, aes(meandisturb/4*100, meancover),
              method = lm, se = FALSE)+ 
  geom_smooth(data = subset(JR_mean_focal, species == "ORDE"), color = "lightgrey", lwd = .75, aes(meandisturb/4*100, meancover),
              method = lm, se = FALSE) + theme(strip.background = element_blank(), strip.text = element_blank()) + 
  geom_point() + 
  labs(x="Percent disturbance", y = "")

l <- lm(meancover ~ meandisturb  +soilDepth, data = subset(JR_mean_focal, species == "BR__"))
summary(l)

d <- ggplot(JR_mean_time_focal2, aes(x=precip, y=meancover, color)) + facet_grid(Genus~dummy, scales = "free") + 
  # geom_smooth(data = subset(JR_mean_time_focal, species == "BRMO"), aes(precip, meancover),
  #             method = lm, se = FALSE) + 
  # geom_smooth(data = subset(JR_mean_time_focal2, species == "CAMU"), aes(precip, meancover),
  #             method = lm, se = FALSE) + 
  # geom_smooth(data = subset(JR_mean_time_focal2, species == "EVSP"), aes(precip, meancover),
  #             method = lm, se = FALSE) + 
  geom_smooth(data = subset(JR_mean_time_focal2, species == "LACA"), color = "lightgrey", lwd = .75, aes(precip, meancover),
              method = lm, se = FALSE) + 
  # geom_smooth(data = subset(JR_mean_time_focal2, species == "LAPL"), aes(precip, meancover),
  #             method = lm, se = FALSE) + 
  geom_smooth(data = subset(JR_mean_time_focal2, species == "SIJU"), color = "lightgrey", lwd = .75, aes(precip, meancover),
              method = lm, se = FALSE) + 
  geom_smooth(data = subset(JR_mean_time_focal2, species == "TIER"), color = "lightgrey", lwd = .75, aes(precip, meancover),
              method = lm, se = FALSE) +
  geom_point() + 
  labs(x="Rainfall (mm)", y = "")  + theme(strip.background = element_blank(), strip.text = element_blank()) 
# geom_smooth(data = subset(JR_mean_time_focal2, species == "VUMI"), aes(precip, meancover),
#                                                   method = lm, se = FALSE) +
# geom_smooth(data = subset(JR_mean_time_focal2, species == "CHPO"), aes(precip, meancover),
#             method = lm, se = FALSE) + theme(strip.background = element_blank(), strip.text = element_blank()) + 


#pdf("mastergraph.pdf", width = 16, height = 24)
plot_grid(a,b,c,d, nrow = 1, axis = "tb", align = "hv")
#dev.off()



JR_lateannual <- JR_tog %>%
  filter(cover > 0, (species == "CAMU" | species == "HELU")) %>%
  spread(species, cover, fill = 0)

cor.test(JR_lateannual$CAMU, JR_lateannual$HELU)

ggplot(JR_lateannual, 
       aes(x=CAMU, y = HELU)) + geom_point() + 
  facet_wrap(~treatment) + geom_smooth(method = "lm")


JR_N <- JR_tog %>%
  filter(cover > 0, (species == "LOSU" | species == "ASGA")) %>%
  spread(species, cover, fill = 0)


ggplot(JR_N, 
       aes(x=ASGA, y = LOSU)) + geom_point() + 
  facet_wrap(~treatment) + geom_smooth(method = "lm")

cor.test(JR_N$ASGA, JR_N$LOSU)



### GOPHER WITH SOIL DEPT ###
ggplot(JR_gopher_mean, aes(x=soilDepth, y=meandisturb)) + geom_point(aes(color =as.factor(trtrep))) + geom_smooth(se=F, method = "lm") +
  geom_errorbar(aes(ymin=meandisturb -sedisturb, ymax= meandisturb + sedisturb, color = as.factor(trtrep))) + 
  facet_wrap(~treatment)


ggplot(JR_gopher_mean, aes(x=soilDepth, y=meandisturb)) + geom_point(aes(color =as.factor(treatment))) + 
  geom_smooth(se=F, method = "lm") +
  geom_errorbar(aes(ymin=meandisturb -sedisturb, ymax= meandisturb + sedisturb, color = as.factor(treatment)))  



ggplot(JR_gopher_mean, aes(x=soilDepth, y=meandisturb)) + geom_point(aes(color =as.factor(treatment))) + 
  geom_smooth(aes(color = as.factor(treatment)), se=F, method = "lm") +
  geom_errorbar(aes(ymin=meandisturb -sedisturb, ymax= meandisturb + sedisturb, color = as.factor(treatment)))  


ggplot(JR_gopher_mean_time, aes(x=precip, y=meandisturb, color = as.factor(treatment))) + geom_point(size = 2) +
  geom_errorbar(aes(ymin=meandisturb -sedisturb, ymax= meandisturb + sedisturb, color = as.factor(treatment)))  +
  theme_bw() + facet_wrap(~treatment) + geom_smooth(se=F, method = "lm", formula = y ~ poly(x,3))


ggplot(JR_gopher_mean_time, aes(x=year, y=meandisturb, color = as.factor(treatment))) + geom_line(lwd = 2) + 
  geom_errorbar(aes(ymin=meandisturb -sedisturb, ymax= meandisturb + sedisturb, color = as.factor(treatment)))  +
  theme_bw() + facet_wrap(~treatment)



###### DOMINANT FORBS
ggplot(subset(JR_mean, species == "PLER"), 
       aes(x=soilDepth, y = meancover)) + geom_point(aes(color = as.factor(trtrep)), size =2) + 
  geom_smooth(se=F, method = "lm") + geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(trtrep))) + 
  facet_wrap(~treatment)

ggplot(subset(JR_mean_time, species == "PLER"), 
       aes(x=precip2, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_smooth(aes(color = treatment), se=F, method = "lm") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)


ggplot(subset(JR_mean_time, species == "PLER"), 
       aes(x=meandisturb, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_smooth(aes(color = treatment), se=F, method = "lm") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)


ggplot(subset(JR_mean_time, species == "PLER"), 
       aes(x=year, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_line(aes(color = treatment), lwd = 1) + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)


ggplot(subset(JR_tog, species == "MIDO" ), 
       aes(x=year, y = cover, group = uniqueID)) + geom_point(aes( color = disturb), size =4) + 
  geom_line( lwd = .1) + 
 theme_bw() +
  facet_wrap(trtrep~treatment)


ggplot(subset(JR_mean, species == "LACA"), 
       aes(x=soilDepth, y = meancover)) + geom_point(aes(color = as.factor(trtrep)), size =2) + 
  geom_smooth(se=F, method = "lm") + geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(trtrep))) + 
  facet_wrap(~treatment)

ggplot(subset(JR_mean_time, species == "LACA"), 
       aes(x=precip2, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_smooth(aes(color = treatment), se=F, method = "lm") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)

ggplot(subset(JR_mean_time, species == "LACA"), 
       aes(x=meandisturb, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_smooth(aes(color = treatment), se=F, method = "lm") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)

ggplot(subset(JR_mean_time, species == "LACA"), 
       aes(x=year, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_line(aes(color = treatment), lwd = 1) + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)




ggplot(subset(JR_mean, species == "MIDO"), 
       aes(x=soilDepth, y = meancover)) + geom_point(aes(color = as.factor(trtrep)), size =2) + 
  geom_smooth(se=F, method = "lm") + geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(trtrep))) + 
  facet_wrap(~treatment)

ggplot(subset(JR_mean_time, species == "MIDO"), 
       aes(x=precip2, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_smooth(aes(color = treatment), se=F, method = "lm") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)



ggplot(subset(JR_mean_time, species == "MIDO"), 
       aes(x=meandisturb, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_smooth(aes(color = treatment), se=F, method = "lm") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)

ggplot(subset(JR_mean_time, species == "MIDO"), 
       aes(x=year, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_line(aes(color = treatment), lwd = 1) + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)





###### ANNUAL GRASS
ggplot(subset(JR_mean, species == "BRHO"), 
       aes(x=soilDepth, y = meancover)) + geom_point(aes(color = as.factor(trtrep)), size =2) + 
  geom_smooth(se=F, method = "lm") + geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(trtrep))) + 
  facet_wrap(~treatment)

ggplot(subset(JR_mean_time, species == "BRHO"), 
       aes(x=precip2, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_smooth(aes(color = treatment), se=F, method = "lm") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)


 ggplot(subset(JR_mean_time, species == "BRHO" & meancover > 0), 
       aes(x=meandisturb, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_smooth(aes(color = treatment), se=F, method = "lm") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)
 
 ggplot(subset(JR_mean_time, species == "SIJU" & meanprevcover > 1), 
        aes(x=meandisturb, y = meancover)) + geom_point( size =4, aes(color = meanprevcover)) + 
   geom_smooth( se=F, method = "lm", formula = y ~ poly(x,3), color = "black") + 
   geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = meanprevcover)) + theme_bw() 
 
 ggplot(subset(JR_mean_time, species == "BRHO"), 
       aes(x=prevprecip2, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_smooth(aes(color = treatment), se=F, method = "lm") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)


ggplot(subset(JR_mean_time, species == "BRHO"), 
       aes(x=year, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_line(aes(color = treatment), lwd = 1) + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() #+
  facet_wrap(~treatment)


ggplot(subset(JR_mean, species == "VUMI"), 
       aes(x=soilDepth, y = meancover)) + geom_point(aes(color = as.factor(trtrep)), size =2) + 
  geom_smooth(se=F, method = "lm") + geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(trtrep))) + 
  facet_wrap(~treatment)

ggplot(subset(JR_mean_time, species == "VUMI"), 
       aes(x=precip2, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_smooth(aes(color = treatment), se=F, method = "lm") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)


ggplot(subset(JR_mean_time, species == "VUMI" & meancover > 0), 
       aes(x=meandisturb, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_smooth(aes(color = treatment), se=F, method = "lm") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)


ggplot(subset(JR_mean_time, species == "VUMI" & meancover > 0), 
       aes(x=meandisturb, y = meancover)) + geom_point( size =2) + 
  geom_smooth( se=F, method = "lm", formula = y ~ poly(x,3)) + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover)) + theme_bw() 

ggplot(subset(JR_mean_time, species == "VUMI"), 
       aes(x=year, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_line(aes(color = treatment), lwd = 1) + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() #+
facet_wrap(~treatment)

###### N-FIXER
ggplot(subset(JR_mean, species == "LOSU"), 
       aes(x=soilDepth, y = meancover)) + geom_point(aes(color = as.factor(trtrep)), size =2) + 
  geom_smooth(se=F, method = "lm") + geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(trtrep))) + 
  facet_wrap(~treatment)

ggplot(subset(JR_mean_time, species == "LOSU"), 
       aes(x=precip2, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_smooth(aes(color = treatment), se=F, method = "lm") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)


ggplot(subset(JR_mean_time, species == "LOSU" & meancover > 0), 
       aes(x=meandisturb, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_smooth(aes(color = treatment), se=F, method = "lm") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)


ggplot(subset(JR_mean_time, species == "LOSU" & meancover > 0), 
       aes(x=meandisturb, y = meancover)) + geom_point( size =2) + 
  geom_smooth( se=F, method = "lm", formula = y ~ poly(x,3)) + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover)) + theme_bw() 

ggplot(subset(JR_mean_time, species == "LOSU"), 
       aes(x=year, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_line(aes(color = treatment), lwd = 1) + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw()


ggplot(subset(JR_mean, species == "ASGA"), 
       aes(x=soilDepth, y = meancover)) + geom_point(aes(color = as.factor(trtrep)), size =2) + 
  geom_smooth(se=F, method = "lm") + geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(trtrep))) + 
  facet_wrap(~treatment)

ggplot(subset(JR_mean_time, species == "ASGA"), 
       aes(x=precip2, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_smooth(aes(color = treatment), se=F, method = "lm") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)



ggplot(subset(JR_mean_time, species == "ASGA" & meancover > 0), 
       aes(x=meandisturb, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_smooth(aes(color = treatment), se=F, method = "lm") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)


ggplot(subset(JR_mean_time, species == "ASGA" & meancover > 0), 
       aes(x=meandisturb, y = meancover)) + geom_point( size =4, aes(color = precip2)) + 
  geom_smooth( se=F, method = "lm", formula = y ~ poly(x,3)) + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = precip2)) + theme_bw() 


ggplot(subset(JR_mean_time, species == "ASGA"), 
       aes(x=year, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_line(aes(color = treatment), lwd = 1) + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw()


ggplot(subset(JR_mean, species == "TRSP"), 
       aes(x=soilDepth, y = meancover)) + geom_point(aes(color = as.factor(trtrep)), size =2) + 
  geom_smooth(se=F, method = "lm") + geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(trtrep))) + 
  facet_wrap(~treatment)

ggplot(subset(JR_mean_time, species == "TRSP"), 
       aes(x=precip2, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_smooth(aes(color = treatment), se=F, method = "lm") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)



ggplot(subset(JR_mean_time, species == "TRSP" ), 
       aes(x=meandisturb, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_smooth(aes(color = treatment), se=F, method = "lm") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)


ggplot(subset(JR_mean_time, species == "TRSP" ), 
       aes(x=meandisturb, y = meancover)) + geom_point( size =4, aes(color = precip2)) + 
  geom_smooth( se=F, method = "lm", formula = y ~ poly(x,3)) + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = precip2)) + theme_bw() 

ggplot(subset(JR_mean_time, species == "TRSP"), 
       aes(x=year, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_line(aes(color = treatment), lwd = 1) + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) +
 theme_bw()
# 
# a <- ggplot(subset(JR_mean_time, species == "TRSP"), 
#        aes(x=year, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
#   geom_line(aes(color = treatment), lwd = 1) + 
#   geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) +
#   facet_grid(~treatment) + theme_bw()
# 
# b <- ggplot(subset(JR_mean_time, species == "BRHO"), 
#             aes(x=year, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
#   geom_line(aes(color = treatment), lwd = 1) + 
#   geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + 
#   theme_bw() + facet_grid(~treatment)
# 
# grid.arrange(a,b)

## PERENNIALS 
ggplot(subset(JR_mean, species == "SIJU"), 
       aes(x=soilDepth, y = meancover)) + geom_point(aes(color = as.factor(trtrep)), size =2) + 
  geom_smooth(se=F, method = "lm") + geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(trtrep))) + 
  facet_wrap(~treatment)

ggplot(subset(JR_mean_time, species == "SIJU"), 
       aes(x=precip2, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_smooth(aes(color = treatment), se=F, method = "lm") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_grid(~treatment)

ggplot(subset(JR_mean_time, species == "SIJU" & meancover > 0), 
       aes(x=meandisturb, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_smooth(aes(color = treatment), se=F, method = "lm") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)

ggplot(subset(JR_mean, species == "SIJU"), 
       aes(x=meandisturb, y = meancover)) + geom_point(aes(color = as.factor(trtrep)), size =2) + 
  geom_smooth(se=F, method = "lm") + geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(trtrep))) + 
  facet_wrap(~treatment)


ggplot(subset(JR_mean_time, species == "SIJU" & meancover > 0), 
       aes(x=meandisturb, y = meancover)) + geom_point( size =4, aes(color = precip2)) + 
  geom_smooth( se=F, method = "lm", formula = y ~ poly(x,3), color = "black") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = precip2)) + theme_bw() 


ggplot(subset(JR_mean_time, species == "SIJU"), 
       aes(x=year, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_line(aes(color = treatment), lwd = 1) + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw()

ggplot(subset(JR_mean, species == "STPU"), 
       aes(x=soilDepth, y = meancover)) + geom_point(aes(color = as.factor(trtrep)), size =2) + 
  geom_smooth(se=F, method = "lm") + geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(trtrep))) + 
  facet_wrap(~treatment) 


ggplot(subset(JR_mean, species == "STPU"), 
       aes(x=meandisturb, y = meancover)) + geom_point(aes(color = as.factor(trtrep)), size =2) + 
  geom_smooth(se=F, method = "lm") + geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(trtrep))) + 
  facet_wrap(~treatment)


ggplot(subset(JR_mean_time, species == "STPU"), 
       aes(x=precip2, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_smooth(aes(color = treatment), se=F, method = "lm") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)



ggplot(subset(JR_mean_time, species == "STPU" & meancover > 0), 
       aes(x=meandisturb, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_smooth(aes(color = treatment), se=F, method = "lm") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)


ggplot(subset(JR_mean_time, species == "STPU" & meancover > 0), 
       aes(x=meandisturb, y = meancover)) + geom_point( size =4, aes(color = precip2)) + 
  geom_smooth( se=F, method = "lm", formula = y ~ poly(x,2), color = "black") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = precip2)) + theme_bw() 


ggplot(subset(JR_mean_time, species == "STPU"), 
       aes(x=year, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_line(aes(color = treatment), lwd = 1) + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw()




### ORTHOCARPUS
ggplot(subset(JR_mean_time, species == "ORDE"), 
       aes(x=precip2, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_smooth(aes(color = treatment), se=F, method = "lm") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)


ggplot(subset(JR_mean, species == "ORDE"), 
       aes(x=soilDepth, y = meancover)) + geom_point(aes(color = as.factor(trtrep)), size =2) + 
  geom_smooth(se=F, method = "lm") + geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(trtrep))) + 
  facet_wrap(~treatment)



ggplot(subset(JR_mean_time, species == "ORDE" & meancover > 0), 
       aes(x=meandisturb, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_smooth(aes(color = treatment), se=F, method = "lm") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)


ggplot(subset(JR_mean_time, species == "ORDE" & meancover > 0), 
       aes(x=meandisturb, y = meancover)) + geom_point( size =4, aes(color = precip2)) + 
  geom_smooth( se=F, method = "lm", formula = y ~ poly(x,2), color = "black") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = precip2)) + theme_bw() 


ggplot(subset(JR_mean_time, species == "ORDE"), 
       aes(x=year, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_line(aes(color = treatment), lwd = 1) + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw()




### Tillea
ggplot(subset(JR_mean_time, species == "TIER"), 
       aes(x=precip2, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_smooth(aes(color = treatment), se=F, method = "lm") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)


ggplot(subset(JR_mean, species == "TIER"), 
       aes(x=soilDepth, y = meancover)) + geom_point(aes(color = as.factor(trtrep)), size =2) + 
  geom_smooth(se=F, method = "lm") + geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(trtrep))) + 
  facet_wrap(~treatment)



ggplot(subset(JR_mean_time, species == "TIER" & meancover > 0), 
       aes(x=meandisturb, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_smooth(aes(color = treatment), se=F, method = "lm") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)


ggplot(subset(JR_mean_time, species == "TIER" & meancover > 0), 
       aes(x=meandisturb, y = meancover)) + geom_point( size =4, aes(color = precip2)) + 
  geom_smooth( se=F, method = "lm", formula = y ~ poly(x,3), color = "black") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = precip2)) + theme_bw() 



ggplot(subset(JR_mean_time, species == "TIER"), 
       aes(x=year, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_line(aes(color = treatment), lwd = 1) + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw()



### Brodiaea
ggplot(subset(JR_mean_time, species == "BRSP"), 
       aes(x=precip2, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_smooth(aes(color = treatment), se=F, method = "lm") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)


ggplot(subset(JR_mean, species == "BRSP"), 
       aes(x=soilDepth, y = meancover)) + geom_point(aes(color = as.factor(trtrep)), size =2) + 
  geom_smooth(se=F, method = "lm") + geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(trtrep))) + 
  facet_wrap(~treatment)


ggplot(subset(JR_mean, species == "BRSP"), 
       aes(x=meandisturb, y = meancover)) + geom_point(aes(color = as.factor(trtrep)), size =2) + 
  geom_smooth(se=F, method = "lm") + geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(trtrep))) + 
  facet_wrap(~treatment)


ggplot(subset(JR_mean_time, species == "BRSP" & meancover > 0), 
       aes(x=meandisturb, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_smooth(aes(color = treatment), se=F, method = "lm") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)


ggplot(subset(JR_mean_time, species == "BRSP" & meancover > 0), 
       aes(x=meandisturb, y = meancover)) + geom_point( size =4, aes(color = precip2)) + 
  geom_smooth( se=F, method = "lm", formula = y ~ poly(x,3), color = "black") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = precip2)) + theme_bw() 



ggplot(subset(JR_mean_time, species == "BRSP"), 
       aes(x=year, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_line(aes(color = treatment), lwd = 1) + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw()

### LATE_SEASON ANNUALS
## BOTH PUNKED OUT WITH THE LATEST DROUGHT
### Hemizonia
# likes rain
ggplot(subset(JR_mean_time, species == "HELU"), 
       aes(x=precip2, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_smooth(aes(color = treatment), se=F, method = "lm") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)


ggplot(subset(JR_mean, species == "HELU"), 
       aes(x=soilDepth, y = meancover)) + geom_point(aes(color = as.factor(trtrep)), size =2) + 
  geom_smooth(se=F, method = "lm") + geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(trtrep))) + 
  facet_wrap(~treatment)


ggplot(subset(JR_mean_time, species == "HELU" & meancover > 0), 
       aes(x=meandisturb, y = meancover)) + geom_point( size =4, aes(color = precip2)) + 
  geom_smooth( se=F, method = "lm", formula = y ~ poly(x,3), color = "black") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = precip2)) + theme_bw() 

# was never in the controls... 
ggplot(subset(JR_mean_time, species == "HELU"), 
       aes(x=year, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_line(aes(color = treatment), lwd = 1) + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw()


### Calycadenia
# likes rain
ggplot(subset(JR_mean_time, species == "CAMU"), 
       aes(x=precip2, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_smooth(aes(color = treatment), se=F, method = "lm") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)

# not big on soil depth
ggplot(subset(JR_mean, species == "CAMU"), 
       aes(x=soilDepth, y = meancover)) + geom_point(aes(color = as.factor(trtrep)), size =2) + 
  geom_smooth( se=F, method = "lm") +
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(trtrep))) + 
  facet_wrap(~treatment)

ggplot(subset(JR_mean, species == "CAMU"), 
       aes(x=soilDepth, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_smooth( se=F, method = "lm") +
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) 


ggplot(subset(JR_mean_time, species == "CAMU" ), 
       aes(x=meandisturb, y = meancover)) + geom_point( size =4, aes(color = precip2)) + 
  geom_smooth( se=F, method = "lm", formula = y ~ poly(x,3), color = "black") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = precip2)) + theme_bw() 


ggplot(subset(JR_mean_time, species == "CAMU"), 
       aes(x=year, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_line(aes(color = treatment), lwd = 1) + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw()

### CHPO
# likes rain
ggplot(subset(JR_mean_time, species == "CHPO"), 
       aes(x=precip2, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_smooth(aes(color = treatment), se=F, method = "lm") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)

# not big on soil depth
ggplot(subset(JR_mean, species == "CHPO"), 
       aes(x=soilDepth, y = meancover)) + geom_point(aes(color = as.factor(trtrep)), size =2) + 
  geom_smooth(se=F, method = "lm") + geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(trtrep))) + 
  facet_wrap(~treatment)


ggplot(subset(JR_mean_time, species == "CHPO" ), 
       aes(x=meandisturb, y = meancover)) + geom_point( size =4, aes(color = precip2)) + 
  geom_smooth( se=F, method = "lm", formula = y ~ poly(x,3), color = "black") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = precip2)) + theme_bw() 

ggplot(subset(JR_mean_time, species == "CHPO"), 
       aes(x=year, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_line(aes(color = treatment), lwd = 1) + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw()


### AGHE
# likes rain
ggplot(subset(JR_mean_time, species == "AGHE"), 
       aes(x=precip2, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_smooth(aes(color = treatment), se=F, method = "lm") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)

# not big on soil depth
ggplot(subset(JR_mean, species == "AGHE"), 
       aes(x=soilDepth, y = meancover)) + geom_point(aes(color = as.factor(trtrep)), size =2) + 
  geom_smooth(se=F, method = "lm") + geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(trtrep))) + 
  facet_wrap(~treatment)


ggplot(subset(JR_mean_time, species == "AGHE" ), 
       aes(x=meandisturb, y = meancover)) + geom_point( size =4, aes(color = precip2)) + 
  geom_smooth( se=F, method = "lm", formula = y ~ poly(x,3), color = "black") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = precip2)) + theme_bw() 


ggplot(subset(JR_mean_time, species == "AGHE"), 
       aes(x=year, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_line(aes(color = treatment), lwd = 1) + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw()



### LAYIA

ggplot(subset(JR_mean_time, species == "LAPL"), 
       aes(x=precip2, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_smooth(aes(color = treatment), se=F, method = "lm") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)


ggplot(subset(JR_mean, species == "LAPL"), 
       aes(x=soilDepth, y = meancover)) + geom_point(aes(color = as.factor(trtrep)), size =2) + 
  geom_smooth(se=F, method = "lm") + geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(trtrep))) + 
  facet_wrap(~treatment)


ggplot(subset(JR_mean_time, species == "LAPL" ), 
       aes(x=meandisturb, y = meancover)) + geom_point( size =4, aes(color = precip2)) + 
  geom_smooth( se=F, method = "lm", formula = y ~ poly(x,3), color = "black") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = precip2)) + theme_bw() 


# kinda like brome...
b<- ggplot(subset(JR_mean_time, species == "LAPL"), 
       aes(x=year, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_line(aes(color = treatment), lwd = 1) + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw()


a <- ggplot(subset(JR_mean_time, species == "BRHO"), 
       aes(x=year, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_line(aes(color = treatment), lwd = 1) + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw()


grid.arrange(a,b)
ggplot(subset(JR_mean_time, species == "MIDO"), 
       aes(x=precip2, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_smooth(aes(color = treatment), se=F, method = "lm") + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw() +
  facet_wrap(~treatment)



## try convergent cross-mapping with LAPL and BRHO
JR_x <- JR_aggcover %>%
  filter(species == "LAPL"| species == "BRHO")

ggplot(JR_x, aes(x=year, y=cover, color = species, group = interaction(species, uniqueID))) + geom_line() +
  theme_bw() + facet_grid(trtrep~treatment)# + scale_y_log10()




JR_x <- JR_aggcover %>%
  filter(species == "LOSU" | species == "ASGA" | species == "TRSP" | species == "BRHO") %>%
  mutate(Nfixer = ifelse(species != "BRHO", "yes", "no")) %>%
  group_by(year, treatment, Nfixer, treatment, trtrep, uniqueID) %>%
  summarize(cover = sum(cover))

ggplot(JR_x, aes(x=year, y=cover, color = Nfixer, group = interaction(Nfixer, uniqueID))) + geom_line() +
  theme_bw() + facet_grid(trtrep~treatment) + scale_y_log10()


JR_x <- JR_mean_time %>%
  filter(species == "LOSU" | species == "ASGA" | species == "TRSP" | species == "BRHO") %>%
  mutate(Nfixer = ifelse(species != "BRHO", "yes", "no")) %>%
  group_by(year, treatment, Nfixer, precip, precip2) %>%
  summarize(meancover = sum(meancover))


ggplot(subset(JR_x), 
       aes(x=year, y = meancover + 1, color = Nfixer)) + geom_point( size =2) + 
  geom_line( lwd = 1) + 
  #geom_errorbar(aes(ymin=meancover -secover + 1, ymax= meancover + secover + 1)) + 
  theme_bw() + 
  facet_wrap(~treatment)# + scale_y_log10()

ggplot(subset(JR_x), 
       aes(x=year, y = meancover + 1, color = species)) + geom_point( size =2) + 
  geom_line( lwd = 1) + 
  #geom_errorbar(aes(ymin=meancover -secover + 1, ymax= meancover + secover + 1)) + 
  theme_bw() + 
  facet_wrap(~treatment) + scale_y_log10()


a <- ggplot(subset(JR_mean_time, species == "LOSU"), 
            aes(x=year, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_line(aes(color = treatment), lwd = 1) + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw()


b <- ggplot(subset(JR_mean_time, species == "ASGA"), 
            aes(x=year, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_line(aes(color = treatment), lwd = 1) + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw()


c <- ggplot(subset(JR_mean_time, species == "TRSP"), 
            aes(x=year, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_line(aes(color = treatment), lwd = 1) + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw()


d <- ggplot(subset(JR_mean_time, species == "BRHO"), 
            aes(x=year, y = meancover)) + geom_point(aes(color = as.factor(treatment)), size =2) + 
  geom_line(aes(color = treatment), lwd = 1) + 
  geom_errorbar(aes(ymin=meancover -secover, ymax= meancover + secover, color = as.factor(treatment))) + theme_bw()

grid.arrange(a,b,c,d, ncol=1)


JR_mean_time_focal <- JR_mean_time %>%
  filter(species%in%c("CAMU", "HELU", "LACA", "BRMO", "LOMU", "EVSP", "PLER", "BRSP", 
                      "VUMY", "LOSU", "SIJU", "STPU", "LAPL", "ORDE", "VUMI",
                      "CHPO", "BR__", "TIER", "MIDO", "ASGA"))  %>%
  tbl_df() %>%
  mutate(Genus = ifelse(Genus == "Orthocarpus", "Castilleja", Genus),
         Genus = ifelse(Genus == "Sitanion", "Elymus", Genus),
         Genus = ifelse(Genus == "Tillaea", "Crassula", Genus),
         Genus = ifelse(Genus == "Lotus", "Acmispon", Genus)) %>%
  mutate(Genus = fct_relevel(Genus, c("Plantago", "Lasthenia", "Microseris", 
                                      "Vulpia", "Bromus",  "Lolium", "Elymus", "Stipa",
                                      "Calycadenia","Hemizonia",  "Acmispon", "Castilleja",  "Layia",
                                      "Chlorogalum", "Brodiaea", "Evax", "Crassula", "Astragalus"))) %>%
  tbl_df() %>%
  mutate(dummy = "")

JR_mean_time_focal2 <- JR_mean_time %>%
  filter(species%in%c("CAMU", "HELU", "LACA", "BRMO", "LOMU", "EVSP", "PLER", "BRSP", 
                      "VUMY", "LOSU", "SIJU", "STPU", "LAPL", "ORDE", "VUMI",
                      "CHPO", "BR__", "TIER", "MIDO", "ASGA")) %>%
  group_by(Genus, Species, LifeForm, year, species, precip) %>%
  summarize(meancover = mean(meancover)) %>%
  tbl_df() %>%
  mutate(Genus = fct_relevel(Genus, c("Plantago", "Lasthenia", "Microseris", 
                                      "Vulpia", "Bromus", "Lolium", "Sitanion", "Stipa",
                                      "Calycadenia", "Hemizonia",  "Lotus", "Orthocarpus",  "Layia",
                                      "Chlorogalum", "Brodiaea", "Evax", "Tillaea", "Astragalus"))) %>%
  tbl_df() %>%
  mutate(dummy = "")
