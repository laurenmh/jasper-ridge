source("Data-cleaning/cover-cleaning.R")


#A function that calculates SE
calcSE<-function(x){
  x2<-na.omit(x)
  sd(x2)/sqrt(length(x2))
}

## A function that orders timepoints based on gopher disturbance
## Sets a timepoint with gophers as 0 and counts up from there to the next timepoint with gophers
## Describes each  successive iteration as "repnum"
gopher_time<-function(data1, site="quadID", year="year", disturbno = 0){
  outnames<-c("quadID", "year", "contyear", "timpoint", "repnum")
  siteout<-as.data.frame(matrix(nrow=0, ncol=5))
  names(siteout)<-outnames
  t=0
  contyr=-1
  repnum=0
  yrs=unique(data1$year)
  for (i in 1:length(yrs)){
    subyear<-data1%>%
      filter(year==yrs[i])
    t=ifelse(subyear$disturb > disturbno, 0, t+1)
    contyr=contyr+1
    repnum=ifelse(subyear$disturb > disturbno, repnum+1, repnum)
    subyear$contyear<-contyr
    subyear$timepoint<-t
    subyear$repnum<-repnum
    siteout<-rbind(siteout, subyear)
  }
  return(siteout)
}

## Apply gopher_time to whenever the site was bombed out (all quads disturbed)
X <- split(JRgopher, JRgopher$quadID)
out <- lapply(X, FUN=gopher_time, "quadID", "year", 3)
gophbigtime <- cbind(do.call("rbind", out)) %>%
  tbl_df() %>%
  select(quadID, disturb, year, contyear, timepoint, repnum) %>%
  mutate(timepoint = ifelse(timepoint > contyear, NA, timepoint)) %>%
  group_by(quadID) %>%
  mutate(rowdif =lead(timepoint),
         repnum2 = lead(repnum)) %>%
  tbl_df() %>%
  group_by(quadID, repnum2) %>%
  mutate(ncount = n())
 

tog <- left_join(JRcover, gophbigtime) 

plst <- tog %>%
  # filter(species == "PLER" | species == "SIJU") %>%
   mutate(keepdat = ifelse(rowdif == 0 & cover > 3, 1, 0)) %>%
  
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

plstall <- tog %>%
  filter(species == "PLER" | species == "SIJU") 
