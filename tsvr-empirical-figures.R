source("Data-cleaning/times-since-gopher.R")
library(tsvr)
library(cowplot)
library(grid)

###########################
# Plantago and Microseris #
##########################

# set a threshold minimum abundance
threshold <- 3

# subset the data
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



# # plot the raw data
# ggplot(subset(plmi, rowdif < 11), aes(x=rowdif, y=cover, color = species, group = interaction(species,quadID))) + 
#   geom_line() + facet_wrap(~interaction(minyear,quadID))


# aggregate for the over time panel
plmi2 <- plmi %>%
  group_by(rowdif, quadID, species, treatment) %>%
  summarize(meancover = mean(cover)) %>%
  tbl_df() %>%
  group_by(rowdif, species) %>%
  summarize(cover = mean(meancover), secover = calcSE(meancover)) 

# make the over time panel
a <- ggplot(subset(plmi2, rowdif < 15), aes(x=rowdif, y=cover, color = species)) + 
  geom_vline(xintercept = 1, color = "lightgrey", lwd = 10) + 
  # geom_hline(yintercept = 11.2, color = "darkgrey", lty = "dashed") +
  #  geom_hline(yintercept = 20.4, color = "darkgrey", lty = "dashed") +
  geom_line() +
  geom_point() + 
  geom_errorbar(aes(ymin = cover - secover, ymax = cover + secover), width = .2) + 
  theme_classic() + labs(y="Percent cover", x = "Time point", color ="Species") + theme(text = element_text(size =14)) 

# # check what patterns look like within each replicate
# plmi3 <- plmi %>%
#   group_by(rowdif, quadID, species, treatment, minyear) %>%
#   summarize(meancover = mean(cover)) %>%
#   tbl_df() %>%
#   group_by(rowdif, species, minyear) %>%
#   summarize(cover = mean(meancover), secover = calcSE(meancover)) 
# 
# 
# ggplot(subset(plmi3, rowdif < 15), aes(x=rowdif, y=cover, color = species)) + 
#   geom_vline(xintercept = 1, color = "lightgrey", lwd = 5) + 
#   geom_line() +
#   geom_point() + 
#   geom_errorbar(aes(ymin = cover - secover, ymax = cover + secover), width = .2) + 
#   theme_classic() + labs(y="Percent cover", x = "Time point", color ="Species") + theme(text = element_text(size =14)) + facet_wrap(~minyear)
# #ggsave("Plantago-Microsersis_gopher-recover-focalyears-faceted.pdf", width = 10, height = 7)



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
Timescale", "Short 
Timescale")
siteout2$facorder <- c(3,2,1)


# plot the average and se of the
#ggplot(siteout2, aes(x=metric, y = value, color = metric)) + geom_boxplot() + facet_wrap(~treatment)
b <- ggplot(siteout2, aes(x=reorder(metric2, facorder), y = meanval, color = metric)) +   
  geom_hline(yintercept = 1, color = "grey", lty = "dashed") + 
  geom_point(size = 3) + 
  geom_errorbar(aes(ymin = meanval - seval, ymax = meanval + seval), width = .2)  +
  theme_classic() + theme(legend.position = "none", text = element_text(size = 14)) + 
  scale_color_manual(values = c( "lightseagreen", "darkgreen",  "blue"))  + 
  labs(x = "", y="Variance Ratio")


# # defunct legend creation, defunct plotting
# legd <- legendGrob(c("Short-term Driver VR", "Long-term Driver VR","Classic VR", expression(italic("Microseris")), expression(italic("Plantago"))), do.lines=TRUE,
#                    gp=gpar(col = c("blue", "darkgreen", "lightseagreen", "black", "darkgrey"), lwd = 2), nrow = 1)

# 
# c <- plot_grid(a + theme(legend.position = "none") + annotate("text", x=0, y =32, label="a)", size = 5) +
#                  scale_color_manual(values = c("black","darkgrey")) + labs(x="Time", y = "Percent Cover") + 
#                  theme(panel.border = element_rect(colour = "black", fill=NA, size=.75))+
#                  annotate("text", x = 12, y = 22, label = "Plantago erecta",fontface = 'italic', color = "darkgrey") +
#                  annotate("text", x = 7, y = 11, label = "Elymus glaucus",fontface = 'italic', color = "black"),
#                b  + annotate("text", x=.6, y =1.25, label="b)", size = 5) +
#                  theme(panel.border = element_rect(colour = "black", fill=NA, size=.75)),
#                align = c("hv"))
# 
# # old one-scenario outout
# # pdf("pler-mido-2panels.pdf", width = 10, height = 5)
# # plot_grid(c, legd, nrow = 2, rel_heights = c(9.5,.5))
# # dev.off()


#######################
# Plantago and Elymus #
#######################

# set thresholds for abundance values
threshold <- 3
upperthreshold <- 50

# subset the data
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



# plot raw data as a check
# ggplot(subset(plst, rowdif < 10), aes(x=rowdif, y=cover, color = species, group = interaction(species,quadID))) + 
#   geom_line() + facet_wrap(~interaction(minyear,quadID))

# aggregate by year for the first panel
plst2 <- plst %>%
  group_by(rowdif, quadID, species, treatment, minyear) %>%
  summarize(meancover = mean(cover)) %>%
  tbl_df() %>%
  group_by(rowdif, species) %>%
  summarize(cover = mean(meancover), secover = calcSE(meancover)) 


a2 <- ggplot(subset(plst2, rowdif < 9), aes(x=rowdif, y=cover, color = species)) + 
  geom_vline(xintercept = 1, color = "lightgrey", lwd = 10) + 
  # geom_hline(yintercept = 9.88, color = "darkgrey", lty = "dashed") +
  #  geom_hline(yintercept = 21.6, color = "darkgrey", lty = "dashed") +
  geom_line() +
  geom_point() + 
  geom_errorbar(aes(ymin = cover - secover, ymax = cover + secover), width = .2) + 
  theme_classic() + labs(y="Percent cover", x = "Time point", color ="Species") + theme(text = element_text(size =14)) 


## calculate the VRs within each quadID
outnames<-c("quadID", "treatment", "classicVR", "longVR", "shortVR")
siteout<-as.data.frame(matrix(nrow=0, ncol=5))
names(siteout)<-outnames
quads <- unique(plst$quadID)

for (i in 1:length(quads)){
  
  subber <- subset(plst, quadID == quads[i]) %>%
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
Timescale", "Short 
Timescale")
siteout2$facorder <- c(3,2,1)


# plot the average and se of the
#ggplot(siteout2, aes(x=metric, y = value, color = metric)) + geom_boxplot() + facet_wrap(~treatment)
b2 <- ggplot(siteout2, aes(x=reorder(metric2, facorder), y = meanval, color = metric)) +   
  geom_hline(yintercept = 1, color = "grey", lty = "dashed") + 
  geom_point(size = 3) + 
  geom_errorbar(aes(ymin = meanval - seval, ymax = meanval + seval), width = .2) + ylim(0.75,1.25) +
  theme_classic() + theme(legend.position = "none", text = element_text(size = 14)) + 
  scale_color_manual(values = c( "lightseagreen", "darkgreen",  "blue"))  + 
  labs(x = "", y="Variance Ratio")

# # obsolete code to add a legend
# legd <- legendGrob(c("Short-term Driver VR", "Long-term Driver VR","Classic VR", expression(italic("Plantago")), expression(italic("Elymus"))), do.lines=TRUE,
#                    gp=gpar(col = c("blue", "darkgreen", "lightseagreen", "darkorchid1", "orange3"), lwd = 2), nrow = 1)


  

# put it all together
c <- plot_grid(a + theme(legend.position = "none") + annotate("text", x=0, y =32, label="a)", size = 5) +
                 scale_color_manual(values = c("black","darkgrey")) + labs(x="Time", y = "Percent Cover") + 
                 theme(panel.border = element_rect(colour = "black", fill=NA, size=.75))+
                 annotate("text", x = 12, y = 22, label = "Plantago erecta",fontface = 'italic', color = "darkgrey") +
                 annotate("text", x = 7, y = 12, label = "Microseris douglasii",fontface = 'italic', color = "black") +
                 scale_x_continuous(breaks = c(0,2,4,6,8,10,12)) + 
                 labs(x=""),
               b  + annotate("text", x=.6, y =1.25, label="b)", size = 5) +
                 theme(panel.border = element_rect(colour = "black", fill=NA, size=.75)) + labs (y=""),
               align = c("hv"))

c2 <- plot_grid(a2 + theme(legend.position = "none") + annotate("text", x=0, y =32, label="c)", size = 5) +
                 scale_color_manual(values = c("darkorchid1", "orange3")) + labs(x="Time", y = "Percent Cover") + 
                 theme(panel.border = element_rect(colour = "black", fill=NA, size=.75)) +
                 annotate("text", x = 6, y = 31, label = "Plantago erecta",fontface = 'italic', color = "darkorchid1") +
                 annotate("text", x = 6, y = 10, label = "Elymus glaucus",fontface = 'italic', color = "orange3"),
               b2  + annotate("text", x=.6, y =1.25, label="d)", size = 5) +
                 theme(panel.border = element_rect(colour = "black", fill=NA, size=.75)),
               align = c("hv"))

pdf("jr-empirical.pdf", width = 10, height = 8)
plot_grid(c, c2, nrow = 2)
 dev.off()
