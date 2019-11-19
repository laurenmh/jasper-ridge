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



JR_spp <- read_csv("~/Dropbox/California Data/JR_species names.csv") %>%
  mutate(species = SpeciesCode)

JR_gopher_mean <- left_join(JR_gopher, JR_soil) %>%
  group_by(soilDepth, maxSoilDepth, minSoilDepth, plot, treatment, trtrep) %>%
  summarize(meandisturb = mean(disturb), sedisturb = sd(disturb)/sqrt(length(disturb))) %>%
  mutate(treatment2 = as.factor(treatment),
         treatment2 = reorder(treatment2, values = c(g,c, r)))


treatment2 = recode(treatment2, g = "Gopher excl", c = "Control", r = "Rabbit excl"),
treatment2 = stats::reorder(treatment2, `Gopher excl`,Control, `Rabbit excl`))


JR_tog0 <- left_join(JR_aggcover, JR_soil)
JR_tog1 <- left_join(JR_tog0, JR_rain)

JR_tog_prev <- JR_tog1 %>%
  select(year, species, uniqueID, cover, precip) %>%
  mutate(year = year + 1) %>%
  mutate(prevcover = cover, prevprecip = precip) %>%
  select(-cover, - precip)

JR_tog2 <- left_join(JR_tog1, JR_spp)
JR_tog3 <- left_join(JR_tog2, JR_gopher) 
JR_tog <- left_join(JR_tog3, JR_tog_prev)

JR_mean0 <- JR_tog %>%
  group_by(Genus, Species, LifeForm, species, soilDepth, maxSoilDepth, minSoilDepth, plot, treatment, trtrep) %>%
  summarize(meancover = mean(cover), secover = sd(cover)/sqrt(length(cover)),
            maxcover = max(cover)) 

JR_mean <- left_join(JR_mean0, JR_gopher_mean)

JR_togmean <- left_join(JR_tog, JR_mean)




JR_mean0 <- JR_tog %>%
  group_by(Genus, Species, LifeForm, species, soilDepth, maxSoilDepth, minSoilDepth, plot, treatment, trtrep) %>%
  summarize(meancover = mean(cover), secover = sd(cover)/sqrt(length(cover)),
            maxcover = max(cover)) 

JR_mean <- left_join(JR_mean0, JR_gopher_mean) %>%
  mutate(treatment2 = as.factor(treatment),
         treatment2 = recode(treatment2, g = "Gopher excl", c = "Control", r = "Rabbit excl"),
         treatment2 = fct_relevel(treatment2, "Gopher excl","Control", "Rabbit excl"))


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


myspp <- c("CAMU", "HELU", "LACA", "BRMO", "LOMU", "EVSP", "PLER",  
           "LOSU", "SIJU", "STPU", "LAPL", "ORDE", "VUMI",
           "CHPO", "BR__", "TIER", "MIDO", "ASGA")
lmoutput <- data.frame(variable = as.character(), estimate = as.numeric(), stderr = as.numeric(), df = as.numeric(), tval = as.numeric(), 
                       pval = as.numeric(), species = as.numeric(), Genus = as.numeric(), Species = as.numeric())
for(i in 1:length(myspp)) {
  
  dat <- subset(JR_togmean, species == myspp[i]) %>%
    filter(!is.na(cover), !is.na(precip), !is.na(prevprecip))
  
  l <- lmer(cover ~  soilDepth+ precip + prevprecip  + meandisturb + fall_ppt + winter_ppt + 
              (1|trtrep/uniqueID) + 
              (1|year), 
            data = dat, na.action = na.omit)
  
  
  step_res <- step(l)
  final <- get_model(step_res)
  
  lmout <- tidy(summary(final)$coefficients)
  names(lmout) = c("variable", "estimate", "stderr", "df", "tval", "pval")
  lmout$species <- unique(dat$species)
  lmout$Genus <- unique(dat$Genus)
  lmout$Species <- unique(dat$Species)
  lmoutput <- rbind(lmoutput, lmout)
}




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
