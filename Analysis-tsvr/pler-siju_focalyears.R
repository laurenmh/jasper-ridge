source("Data-cleaning/times-since-gopher.R")
library(tsvr)
library(cowplot)
library(grid)

threshold <- 3
upperthreshold <- 50

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



plst2 <- plst %>%
  group_by(rowdif, quadID, species, treatment, minyear, uniqueID) %>%
  summarize(meancover = mean(cover)) %>%
  tbl_df() %>%
  group_by(rowdif, species) %>%
  summarize(cover = mean(meancover), secover = calcSE(meancover)) 


a <- ggplot(subset(plst2, rowdif < 9), aes(x=rowdif, y=cover, color = species)) + 
  geom_vline(xintercept = 1, color = "lightgrey", lwd = 10) + 
 # geom_hline(yintercept = 9.88, color = "darkgrey", lty = "dashed") +
#  geom_hline(yintercept = 21.6, color = "darkgrey", lty = "dashed") +
  geom_line() +
  geom_point() + 
  geom_errorbar(aes(ymin = cover - secover, ymax = cover + secover), width = .2) + 
  theme_classic() + labs(y="Percent cover", x = "Time point", color ="Species") + theme(text = element_text(size =14))
#ggsave("Plantago-Sitanion_gopher-recover-focalyears.pdf", width = 6, height = 5)



## calculate the VRs within each quadID
outnames<-c("quadID", "treatment", "classicVR", "longVR", "shortVR")
siteout<-as.data.frame(matrix(nrow=0, ncol=5))
names(siteout)<-outnames
quads <- unique(plst$uniqueID)

for (i in 1:length(quads)){
  
  subber <- subset(plst, uniqueID == quads[i]) %>%
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
b <- ggplot(siteout2, aes(x=reorder(metric2, facorder), y = meanval, color = metric)) +   
  geom_hline(yintercept = 1, color = "grey", lty = "dashed") + 
  geom_point(size = 3) + 
  geom_errorbar(aes(ymin = meanval - seval, ymax = meanval + seval), width = .2) + ylim(0.75,1.25) +
  theme_classic() + theme(legend.position = "none", text = element_text(size = 14)) + 
  scale_color_manual(values = c( "lightseagreen", "darkgreen",  "blue"))  + 
  labs(x = "", y="Variance Ratio")
#ggsave("Plantago-Sitanion_tsvr_focalyears.pdf", width = 8, height = 5)  

#ggsave("Plantago-Microseris_tsvr_focalyears-repnumunder15.pdf", width = 8, height = 5)  



legd <- legendGrob(c("Short-term Driver VR", "Long-term Driver VR","Classic VR", expression(italic("Plantago")), expression(italic("Elymus"))), do.lines=TRUE,
                   gp=gpar(col = c("blue", "darkgreen", "lightseagreen", "darkorchid1", "orange3"), lwd = 2), nrow = 1)


c <- plot_grid(a + theme(legend.position = "none") + annotate("text", x=0, y =32, label="a)", size = 5) +
                 scale_color_manual(values = c("darkorchid1", "orange3")) + labs(x="Time", y = "Percent Cover") + 
                 theme(panel.border = element_rect(colour = "black", fill=NA, size=.75)),
               b  + annotate("text", x=.6, y =1.25, label="b)", size = 5) +
                 theme(panel.border = element_rect(colour = "black", fill=NA, size=.75)),
               align = c("hv"))

pdf("annual_perennial_empirical2.pdf", width = 10, height = 5)
plot_grid(c, legd, nrow = 2, rel_heights = c(9.5,.5))
dev.off()




### only in never disturbed plots

plst <- tog %>%
  # filter(species == "PLER" | species == "SIJU") %>%
  mutate(keepdat = ifelse(contyear == 0 & cover > 3, 1, 0)) %>%
  mutate(maxkeep = max(keepdat)) %>%
  filter(maxkeep > 0) %>%
  mutate(keepdat2 = ifelse(contyear == 0 & disturb == 0, 1, 0)) %>%
  mutate(maxkeep2 = max(keepdat2))  %>%
  filter(maxkeep2 > 0) %>%
  # mutate(keepdat = ifelse(rowdif == 0 & species == "PLER" & cover > 5, 1, 0),
  #        keepdat = ifelse(rowdif == 0 & species == "SIJU" & cover > 5, 1, keepdat)) %>%
  # mutate(keepdat2 = ifelse(timepoint == 0 & species == "PLER" & cover < 5, 1, 0),
  #        keepdat2 = ifelse(timepoint == 0 & species == "SIJU" & cover < 4, 1, keepdat2)) %>%
  group_by(repnum2, species, quadID) %>%
  mutate(maxkeep = max(keepdat)) %>%
  filter(maxkeep > 0) %>%
  tbl_df() %>%
  filter(ncount > 4) #%>%
# filter(keepdat2 > 0)

plst2 <- plst %>%
  group_by(rowdif, quadID, species, treatment) %>%
  summarize(meancover = mean(cover)) %>%
  tbl_df() %>%
  group_by(rowdif, species) %>%
  summarize(cover = mean(meancover), secover = calcSE(meancover)) 

ggplot(plst2, aes(x=rowdif, y=cover, color = species)) + geom_point()  +
  geom_line() + xlim(0,8) + facet_wrap(~species)

plst3 <- plst %>%
  filter(species == "PLER" | species == "SIJU") %>%
  group_by(rowdif, quadID, species, treatment) %>%
  summarize(meancover = mean(cover)) %>%
  tbl_df() %>%
  group_by(species, quadID) %>%
  # mutate(maxcover = ifelse(rowdif == 0, meancover, NA),
  #        maxcover = max(maxcover, na.rm=T),
  #        meancover = (meancover -maxcover)/maxcover) %>%
  tbl_df() %>%
  group_by(rowdif, species) %>%
  summarize(cover = mean(meancover), secover = calcSE(meancover)) 


ggplot(subset(plst3, rowdif < 9), aes(x=rowdif, y=cover, color = species)) + geom_point() +
  geom_hline(yintercept = 1, color = "grey", lty = "dashed") +
  geom_line() +
  geom_errorbar(aes(ymin = cover - secover, ymax = cover + secover), width = .2)


plst3 <- plst %>%
  filter(species == "PLER" | species == "SIJU") %>%
  group_by(rowdif, quadID, species, treatment) %>%
  summarize(meancover = mean(cover)) %>%
  tbl_df() %>%
  group_by(rowdif, species) %>%
  summarize(cover = mean(meancover), secover = calcSE(meancover)) 


ggplot(subset(plst3, rowdif < 9), aes(x=rowdif, y=cover, color = species)) + geom_point() + 
  geom_vline(xintercept = 1, color = "lightgrey", lwd = 10) + 
  geom_hline(yintercept = 9.47, color = "darkgrey", lty = "dashed") +
  geom_hline(yintercept = 23.0, color = "darkgrey", lty = "dashed") +
  geom_line() +
  geom_errorbar(aes(ymin = cover - secover, ymax = cover + secover), width = .2) + 
  theme_classic() + labs(y="Percent cover", x = "Time point", color ="Species") + theme(text = element_text(size =14))
ggsave("Plantago-Sitanion_gopher-recover.pdf", width = 6, height = 5)


## all data
plstall <- JRcover %>%
  filter(species == "PLER" | species == "SIJU",
         treatment != "r")

plstall2 <- plstall %>%
  group_by(quadID, species) %>%
  mutate(totalcover = sum(cover)) %>%
  tbl_df() %>%
  group_by(quadID) %>%
  mutate(mincover = min(totalcover)) %>%
  filter(mincover > 0)


outnames<-c("quadID", "treatment", "classicVR", "longVR", "shortVR")
siteout<-as.data.frame(matrix(nrow=0, ncol=5))
names(siteout)<-outnames
quads <- unique(plstall2$quadID)

for (i in 1:length(quads)){
  
  subber <- subset(plstall2, quadID == quads[i]) %>%
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

siteout2 <- siteout %>%
  gather(metric, value, classicVR:shortVR) %>%
  group_by(treatment, metric) %>%
  summarize(meanval = mean(value), seval = sd(value)/sqrt(n()))
siteout2$metric2 <- c("Classic", "Long 
                      Time-Scale", "Short 
                      Time-Scale", "Classic", "Long 
                      Time-Scale", "Short 
                      Time-Scale")
siteout2$facorder <- c(3,2,1, 3, 2,1)
siteout2$treatment2 <- c(rep("Control", 3), rep("Gopher Excluded", 3))


#ggplot(siteout2, aes(x=metric, y = value, color = metric)) + geom_boxplot() + facet_wrap(~treatment)
ggplot(siteout2, aes(x=reorder(metric2, facorder), y = meanval, color = metric)) +   
  geom_hline(yintercept = 1, color = "grey", lty = "dashed") + 
  geom_point(size = 3) + 
  facet_wrap(~treatment2) +
  geom_errorbar(aes(ymin = meanval - seval, ymax = meanval + seval), width = .2) + ylim(0.75,1.25) +
  theme_bw() + theme(legend.position = "none", text = element_text(size = 14)) + 
  scale_color_manual(values = c( "turquoise4", "darkgreen",  "darkblue"))  + 
  labs(x = "", y="Variance Ratio")
ggsave("Plantago-Sitanion_tsvr_alldat.pdf", width = 8, height = 5)  






demo <- subset(plst3, rowdif < 9) %>%
  select(rowdif, species, cover) %>%
  spread(species, cover)

d <- t(as.matrix(demo[,2:dim(demo)[2]]))


res0 <- vreq_classic(d)
classicVR <- res0[[3]]

res <-tsvreq_classic(d)

aggresLong<-aggts(res,res$ts[res$ts>=4])
aggresShort<-aggts(res,res$ts[res$ts<4])

longVR <- aggresLong[[3]]
shortVR <- aggresShort[[3]]

demoout <- data.frame(rbind(classicVR, longVR, shortVR))
names(demoout) = "value"
demoout$metric <- row.names(demoout)
demoout$metric2 <- c("Classic", "Long 
                     Time-Scale", "Short 
                     Time-Scale")
demoout$facorder <- c(3,2,1)

ggplot(demoout, aes(x=reorder(metric2, facorder), y=value, color = metric2)) + 
  geom_point(size = 3) + ylim(.5, 2) +
  theme_bw() + 
  theme(legend.position = "none", text = element_text(size = 14)) + labs(y="Variance Ratio", x = "") + 
  geom_hline(yintercept = 1, color = "grey", lty = "dashed") + 
  scale_color_manual(values = c( "turquoise4", "darkgreen",  "darkblue")) 
ggsave("Plantago-Sitanion_tsvr_demodat.pdf", width = 5, height = 5)  

## now for two annuals with different responses

plmi <- tog %>%
  filter(species == "PLER" | species == "MIDO") %>%
  mutate(keepdat = ifelse(rowdif == 0 & cover > 3, 1, 0)) %>%
  group_by(repnum2, species, quadID) %>%
  mutate(maxkeep = max(keepdat)) %>%
  filter(maxkeep > 0) %>%
  tbl_df() %>%
  filter(ncount > 4) 


plmi2 <- plmi %>%
  group_by(rowdif, quadID, species, treatment) %>%
  summarize(meancover = mean(cover)) %>%
  tbl_df() %>%
  group_by(rowdif, species) %>%
  summarize(cover = mean(meancover), secover = calcSE(meancover)) 

ggplot(plmi2, aes(x=rowdif, y=cover, color = species)) + geom_line()



ggplot(subset(plmi2, rowdif < 20), aes(x=rowdif, y=cover, color = species)) + geom_point() + 
  geom_vline(xintercept = 1, color = "lightgrey", lwd = 10) + 
  geom_hline(yintercept = 9.47, color = "darkgrey", lty = "dashed") +
  geom_hline(yintercept = 23.0, color = "darkgrey", lty = "dashed") +
  geom_line() +
  geom_errorbar(aes(ymin = cover - secover, ymax = cover + secover), width = .2) + 
  theme_classic() + labs(y="Percent cover", x = "Time point", color ="Species") + theme(text = element_text(size =14))
ggsave("Plantago-Microseris_gopher-recover.pdf", width = 6, height = 5)


demo <- subset(plmi2, rowdif < 9) %>%
  select(rowdif, species, cover) %>%
  spread(species, cover)

d <- t(as.matrix(demo[,2:dim(demo)[2]]))


res0 <- vreq_classic(d)
classicVR <- res0[[3]]

res <-tsvreq_classic(d)

aggresLong<-aggts(res,res$ts[res$ts>=4])
aggresShort<-aggts(res,res$ts[res$ts<4])

longVR <- aggresLong[[3]]
shortVR <- aggresShort[[3]]

demoout <- data.frame(rbind(classicVR, longVR, shortVR))
names(demoout) = "value"
demoout$metric <- row.names(demoout)
demoout$metric2 <- c("Classic", "Long 
                     Time-Scale", "Short 
                     Time-Scale")
demoout$facorder <- c(3,2,1)


ggplot(demoout, aes(x=reorder(metric2, facorder), y=value, color = metric2)) + 
  geom_point(size = 3) + ylim(0, 2) +
  theme_bw() + 
  theme(legend.position = "none", text = element_text(size = 14)) + labs(y="Variance Ratio", x = "") + 
  geom_hline(yintercept = 1, color = "grey", lty = "dashed") + 
  scale_color_manual(values = c( "turquoise4", "darkgreen",  "darkblue")) 


## now by treatment

plstall <- JRcover %>%
  filter(species == "MIDO" | species == "PLER",
         treatment != "r")

plstall2 <- plstall %>%
  group_by(quadID, species) %>%
  mutate(totalcover = sum(cover)) %>%
  tbl_df() %>%
  group_by(quadID) %>%
  mutate(mincover = min(totalcover)) %>%
  filter(mincover > 0)


outnames<-c("quadID", "treatment", "classicVR", "longVR", "shortVR")
siteout<-as.data.frame(matrix(nrow=0, ncol=5))
names(siteout)<-outnames
quads <- unique(plstall2$quadID)

for (i in 1:length(quads)){
  
  subber <- subset(plstall2, quadID == quads[i]) %>%
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



siteout2 <- siteout %>%
  gather(metric, value, classicVR:shortVR) %>%
  group_by(treatment, metric) %>%
  summarize(meanval = mean(value), seval = sd(value)/sqrt(n()))

siteout2$metric2 <- c("Classic", "Long 
                      Time-Scale", "Short 
                      Time-Scale", "Classic", "Long 
                      Time-Scale", "Short 
                      Time-Scale")
siteout2$facorder <- c(3,2,1, 3, 2,1)
siteout2$treatment2 <- c(rep("Control", 3), rep("Gopher Excluded", 3))


#ggplot(siteout2, aes(x=metric, y = value, color = metric)) + geom_boxplot() + facet_wrap(~treatment)
ggplot(siteout2, aes(x=reorder(metric2, facorder), y = meanval, color = metric)) +   
  geom_hline(yintercept = 1, color = "grey", lty = "dashed") + 
  geom_point(size = 3) + 
  facet_wrap(~treatment2) +
  geom_errorbar(aes(ymin = meanval - seval, ymax = meanval + seval), width = .2) + 
  theme_bw() + theme(legend.position = "none", text = element_text(size = 14)) + 
  scale_color_manual(values = c( "turquoise4", "darkgreen",  "darkblue"))  + 
  labs(x = "", y="Variance Ratio")

ggsave("Plantago-Microseris_tsvr_alldat.pdf", width = 8, height = 5)  


quadinfo <- tog %>%
  tbl_df() %>%
  select(quadID, yrsdist, treatment) %>%
  unique()

siteout3 <- left_join(siteout, quadinfo) %>%
  gather(metric, value, classicVR:shortVR) %>%
  mutate(freqdist = ifelse(yrsdist > 5, 1, 0))  %>%
  group_by(freqdist, metric) %>%
  summarize(meanval = mean(value), seval = sd(value)/sqrt(n()))


siteout3$metric2 <- c("Classic", "Long 
                      Time-Scale", "Short 
                      Time-Scale", "Classic", "Long 
                      Time-Scale", "Short 
                      Time-Scale")
siteout3$facorder <- c(3,2,1, 3, 2,1)

ggplot(siteout3, aes(x=reorder(metric2, facorder), y = meanval, color = metric)) +   
  geom_hline(yintercept = 1, color = "grey", lty = "dashed") + 
  geom_point(size = 3) + 
  facet_wrap(~freqdist) +
  geom_errorbar(aes(ymin = meanval - seval, ymax = meanval + seval), width = .2) + 
  theme_bw() + theme(legend.position = "none", text = element_text(size = 14)) + 
  scale_color_manual(values = c( "turquoise4", "darkgreen",  "darkblue"))  + 
  labs(x = "", y="Variance Ratio")
ggsave("Plantago-Microseris_tsvr_demodat.pdf", width = 5, height = 5)  


ggplot(quadinfo, aes(x=yrsdist)) + geom_histogram() + facet_wrap(~treatment)
