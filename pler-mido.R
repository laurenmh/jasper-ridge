source("Data-cleaning/times-since-gopher.R")
library(tsvr)

threshold <- 3

plmi <- tog %>%
  # only keep a record if the cover is above a threshold value in the initial year
  mutate(keepdat = ifelse(rowdif == 0 & cover > threshold, 1, 0)) %>%
  group_by(repnum2, species, quadID) %>%
  mutate(maxkeep = max(keepdat)) %>%
  filter(maxkeep > 0) %>%
  
  # subset the focal species
  filter(species == "PLER" | species == "MIDO") %>%
  
  # only keep a plot if both species are above the threshold
  group_by(quadID, treatment, trtrep, subplot, repnum2) %>%
  mutate(nospp = length(unique((species)))) %>%
  filter(nospp == 2) %>%
  
  # restrict to ts with at least 8 recrods
  filter(ncount > 9) %>%
  
  #identify the first year
  group_by(species, quadID, treatment, trtrep, subplot, repnum2) %>%
  mutate(minyear = min(year))




ggplot(subset(plmi, rowdif < 11), aes(x=rowdif, y=cover, color = species, group = interaction(species,quadID))) + 
  geom_line() + facet_wrap(~interaction(minyear,quadID))



plmi2 <- plmi %>%
  group_by(rowdif, quadID, species, treatment) %>%
  summarize(meancover = mean(cover)) %>%
  tbl_df() %>%
  group_by(rowdif, species) %>%
  summarize(cover = mean(meancover), secover = calcSE(meancover)) 


ggplot(subset(plmi2, rowdif < 15), aes(x=rowdif, y=cover, color = species)) + 
  geom_vline(xintercept = 1, color = "lightgrey", lwd = 10) + 
  geom_hline(yintercept = 11.2, color = "darkgrey", lty = "dashed") +
  geom_hline(yintercept = 20.4, color = "darkgrey", lty = "dashed") +
  geom_line() +
  geom_point() + 
  geom_errorbar(aes(ymin = cover - secover, ymax = cover + secover), width = .2) + 
  theme_classic() + labs(y="Percent cover", x = "Time point", color ="Species") + theme(text = element_text(size =14)) 
ggsave("Plantago-Microsersis_gopher-recover-focalyears.pdf", width = 6, height = 5)

plmi3 <- plmi %>%
  group_by(rowdif, quadID, species, treatment, minyear) %>%
  summarize(meancover = mean(cover)) %>%
  tbl_df() %>%
  group_by(rowdif, species, minyear) %>%
  summarize(cover = mean(meancover), secover = calcSE(meancover)) 


ggplot(subset(plmi3, rowdif < 15), aes(x=rowdif, y=cover, color = species)) + 
  geom_vline(xintercept = 1, color = "lightgrey", lwd = 5) + 
  geom_line() +
  geom_point() + 
  geom_errorbar(aes(ymin = cover - secover, ymax = cover + secover), width = .2) + 
  theme_classic() + labs(y="Percent cover", x = "Time point", color ="Species") + theme(text = element_text(size =14)) + facet_wrap(~minyear)
ggsave("Plantago-Microsersis_gopher-recover-focalyears-faceted.pdf", width = 10, height = 7)



## calculate the VRs within each quadID
outnames<-c("quadID", "treatment", "classicVR", "longVR", "shortVR")
siteout<-as.data.frame(matrix(nrow=0, ncol=5))
names(siteout)<-outnames
quads <- unique(plmi$quadID)

for (i in 1:length(quads)){
  
  subber <- subset(plmi, quadID == quads[i] &  rowdif < 15) %>%
    tbl_df()
  
  subber2 <- subber %>%
    select(year, species, cover) %>%
    spread(species, cover, fill = 0)
  
  d <- t(as.matrix(subber2[,2:dim(subber2)[2]]))
  
  subdat <-subber %>%
    select(quadID, treatment) %>%
    unique()
  
  res0 <- vreq_classic(d)
  subdat$classicVR <- res0[[3]]
  
  res <-tsvreq_classic(d)
  
  aggresLong<-aggts(res,res$ts[res$ts>=4])
  aggresShort<-aggts(res,res$ts[res$ts<4])
  
  subdat$longVR <- aggresLong[[3]]
  subdat$shortVR <- aggresShort[[3]]
  
  siteout<-rbind(siteout, subdat)
  
}

# aggregate and summarize the ts values
siteout2 <- siteout %>%
  gather(metric, value, classicVR:shortVR) %>%
  group_by(metric) %>%
  summarize(meanval = mean(value), seval = sd(value)/sqrt(n()))
siteout2$metric2 <- c("Classic", "Long 
Time-Scale", "Short 
Time-Scale")
siteout2$facorder <- c(3,2,1)


# plot the average and se of the
#ggplot(siteout2, aes(x=metric, y = value, color = metric)) + geom_boxplot() + facet_wrap(~treatment)
ggplot(siteout2, aes(x=reorder(metric2, facorder), y = meanval, color = metric)) +   
  geom_hline(yintercept = 1, color = "grey", lty = "dashed") + 
  geom_point(size = 3) + 
  geom_errorbar(aes(ymin = meanval - seval, ymax = meanval + seval), width = .2)  +
  theme_bw() + theme(legend.position = "none", text = element_text(size = 14)) + 
  scale_color_manual(values = c( "turquoise4", "darkgreen",  "darkblue"))  + 
  labs(x = "", y="Variance Ratio")
ggsave("Plantago-Microseris_tsvr_focalyears-repnumunder15.pdf", width = 8, height = 5)  
