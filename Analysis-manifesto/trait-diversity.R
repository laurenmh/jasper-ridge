library(tidyverse)
library(readr)
library(vegan)
library(FD)
library(cowplot)

## Set ggplot2 theme
theme_set(theme_bw())
theme_update( panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
              strip.background = element_blank(),
              text = element_text(size = 16),
              strip.text= element_text(size = 17), 
              axis.text = element_text(size = 14))

source("Data-cleaning/cover-1m2-cleaning.R")


trait.dat0 <- read_csv("~/Dropbox/California Data/trait_all.csv") 
names(trait.dat0)
names(trait.dat0) = c("species", "code", "status", "habit", "longevity", "Nfixer", "rootDepth",
                      "SRL", "LMA", "leafN", "A", "TcorA", "PS2", "WUE", "height", "germDays", "RS",
                      "seed", "NUE", "NUEadj") 

trait.dat <- trait.dat0 %>%
  mutate(code = ifelse(code == "brho", "brmo", code),
         code = ifelse(code == "dica", "br__", code),
         code = ifelse(code == "crco", "tier", code),
         code = ifelse(code == "hesp", "evsp", code),
         code = ifelse(code == "heco", "helu", code),
         code = ifelse(code == "lowr", "losu", code),
         code = ifelse(code == "napu", "stpu", code),
         code = ifelse(code == "elmu", "siju", code),
         code = ifelse(code == "cade", "orde", code),
         code = ifelse(code == "pose", "posc", code))

traitcodes <- unique(trait.dat$code)
covercodes <- unique(tolower(JRm2cover$species))

tog <- intersect(traitcodes, covercodes)
sort(tog)
sort(covercodes)
sort(traitcodes)

trait.dat2 <- trait.dat %>%
  filter(species%in%c( "Layia platyglossa","Calycadenia multiglandulosa",
                       "Microseris douglasii", "Castilleja densiflora",
                       "Hemizonia congesta",
                       "Lotus wrangelianus" ,  "Trifolium albopurpureum" , "Plantago erecta",
                       "Astragalus gambelianus", "Lasthenia californica", "Hesperevax sparsiflora" ,
                       "Crassula connata", "Elymus multisedus", "Nassella pulchra", "Lolium multiflorum", "Bromus hordaceous", "Chlorogalum pomeridianum var. pomeridianum","Dichelostemma capitatum", "Triteleia laxa",  "Vulpia microstachys var pauciflora")) %>%
  filter(species != "Crassula connata") %>% 
  mutate(RS = ifelse(species == "Triteleia laxa", 3.68009626, RS)) %>%
  mutate(SRL = ifelse(species == "Dichelostemma capitatum", 0.02259358, RS)) %>%
  mutate(NUE = as.numeric(as.character(NUE)))%>%
  select(-status, -habit, -longevity, -Nfixer, -A, -germDays, -seed) %>%
  filter(code != "tral") %>%
  arrange(code)

traitinfo <- sort(unique(trait.dat2$code))
tog2 <- intersect(traitinfo, covercodes)

# Read in the species data
JR_allcover <- JRm2cover %>%
  tbl_df() %>%
  mutate(species = as.character(species),
         species = ifelse(species == "TRAL" | species == "TRTR", "TRSP", species),
         species = tolower(species)) %>%
  filter(species%in%c(tog2)) %>%
  group_by(year, species, treatment, trtrep, plot, uniqueID) %>%
  summarize(cover = sum(cover)) %>%
  tbl_df() %>%
  filter(species != "BARE", species != "ROCK") %>%
  mutate(plotyr = paste(uniqueID, year, sep = "_")) %>%
  spread(species, cover, fill = 0)


abundance <- as.matrix(JR_allcover[1:1332, 7:ncol(JR_allcover)])
row.names(abundance) <- JR_allcover$plotyr


traitdat <- trait.dat2 %>%
  filter(code%in%c(tog2))

traits <- as.matrix(traitdat[1:17, c(3:ncol(traitdat))])
row.names(traits) <- traitdat$code

sort(colnames(abundance))
sort(rownames(traits))

myrda <- rda(traits, scale = TRUE)
biplot(myrda)
cor(traits)


# sppout<-as.data.frame(scores(myrda, choices=c(1,2), display=c("species")))
# sppout$type<-"species"
# sppout$name<-rownames(sppout)
# siteout<-as.data.frame(scores(myrda, choices=c(1,2), display=c("sites")))
# siteout$code<-rownames(siteout)


# extract values

# extract values
siteout <- as.data.frame(scores(myrda, choices=c(1,2), display=c("sites")))
siteout$code<-rownames(siteout)
siteout$name <- siteout$code

enviroout<-as.data.frame(scores(myrda, choices=c(1,2), display=c("species")))
enviroout$type<-"traits"
enviroout$name<-rownames(enviroout)



PCtraits <- as.data.frame(siteout[,1:2])


PCresults<-dbFD(PCtraits, abundance, corr="cailliez")

PCresults <- data.frame(PCresults) 
PCresults$plotyr <- rownames(PCresults)
PCresults <- PCresults %>%
  tbl_df() %>%
  separate(plotyr, c("treatment", "rep", "block", "position", "year")) 

PCresultagg <- PCresults %>%
  group_by(year, treatment) %>%
  summarize_all(funs(mean)) %>%
  tbl_df() %>%
  mutate(year = as.numeric(year))

A <- ggplot(PCresultagg, aes(x=year, y=CWM.PC1, color = treatment)) + geom_line() +
  labs(x = "Year", y = "CWM of PC1")

B <- ggplot(PCresultagg, aes(x=year, y=CWM.PC2, color = treatment)) + geom_line()  +
  labs(x = "Year", y = "CWM of PC2")

plot_grid(A + labs(x="") + theme(axis.text.x = element_blank()), B, ncol = 1)

tog <- right_join(trait.dat, siteout) %>%
  mutate(func = paste(status, habit, longevity, sep = "_")) %>%
  mutate(func = tolower(func)) %>%
  mutate(species = ifelse(species == "Vulpia microstachys var pauciflora", "Vulpia microstachys", species),
         species = ifelse(species == "Dichelostemma capitatum", "Brodiaea sp.", species),
         species = ifelse(species == "Nassella pulchra", "Stipa pulchra", species),
         species = ifelse(species == "Chlorogalum pomeridianum var. pomeridianum",
                          "Chlorogalum pomeridianum", species))


#pdf("TraitPCA-justforbs.pdf", width = 9, height = 7.5)
#pdf("TraitPCA_noLegumes.pdf", width = 9, height = 7.5)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")

C <- ggplot(tog, aes(x=PC1, y=PC2))+ 
  geom_hline(aes(yintercept=0), color="grey") + 
  geom_vline(aes(xintercept=0), color="grey") +
  geom_text(aes(label = code, color = func), size = 4) +
  # scale_color_manual(values = c("grey20", "grey70")) +
  geom_segment(data = enviroout,
               aes(x = 0, xend =  PC1,
                   y = 0, yend =  PC2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "black") + #grid is required for arrow to work.
  geom_text(data = enviroout,
            aes(x=  PC1*1.2, y =  PC2*1.2, #we add 10% to the text to push it slightly out from arrows
                label = name), #otherwise you could use hjust and vjust. I prefer this option
            size = 5,
            hjust = 0.5, 
            color="black") + 
  theme_bw() +theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                    text=element_text(size = 16))+ 
  xlab(paste("Axis 1 (",sprintf("%.1f",myrda$CA$eig["PC1"]/myrda$tot.chi*100,3),"%)",sep="")) +
  ylab(paste("Axis 2 (",sprintf("%.1f",myrda$CA$eig["PC2"]/myrda$tot.chi*100,3),"%)",sep="")) +
  scale_color_manual(values = cbPalette) + theme(legend.position = "none")

#pdf("trait-pca-patterns-code.pdf", width = 10, height = 5)
plot_grid(C + annotate("text", x = -1.7, y = 1.6, label = "(A)", size = 5),
          plot_grid(A + labs(x="") + theme(axis.text.x = element_blank()) +
                      annotate("text", x= 1984, y=.35, label = "(B)", size = 5), 
                    B + annotate("text", x= 1984, y=.3, label = "(C)", size = 5), ncol = 1))
#dev.off()
siteout

JR_allcover2 <- JRm2cover %>%
  tbl_df() %>%
  mutate(species = as.character(species),
         species = ifelse(species == "TRAL" | species == "TRTR", "TRSP", species),
         species = tolower(species)) %>%
  filter(species%in%c(tog2)) %>%
  group_by(year, species, treatment, trtrep, plot, uniqueID) %>%
  summarize(cover = sum(cover)) %>%
  tbl_df() %>%
  filter(species != "BARE", species != "ROCK") %>%
  mutate(code = tolower(species))

covtogPC2 <- right_join(JR_allcover2, siteout) %>%
  mutate(PC1cat = ifelse(PC1 > 0, "high", "low"),
         PC2cat = ifelse(PC2 > 0, "high", "low")) %>%
  group_by(year, PC2cat, treatment, trtrep, plot, uniqueID) %>%
  summarize(cover = sum(cover)) %>%
  group_by(year, PC2cat) %>%
  summarize(meancover = mean(cover))


covtogPC1 <- right_join(JR_allcover2, siteout) %>%
  mutate(PC1cat = ifelse(PC1 > 0, "high", "low"),
         PC2cat = ifelse(PC2 > 0, "high", "low")) %>%
  group_by(year, PC1cat, treatment, trtrep, plot, uniqueID) %>%
  summarize(cover = sum(cover)) %>%
  group_by(year, PC1cat) %>%
  summarize(meancover = mean(cover))


library(codyn)
ggplot(covtogPC2, aes(x=year, y=meancover, color = PC2cat)) + geom_line()

variance_ratio(subset(covtogPC2, year < 2007), "year", "PC2cat", "meancover", 500)
variance_ratio(subset(covtogPC2, year >= 2007), "year", "PC2cat", "meancover", 500)


variance_ratio(subset(covtogPC1, year < 2007), "year", "PC1cat", "meancover", 500)
variance_ratio(subset(covtogPC1, year >= 2007), "year", "PC1cat", "meancover", 500)


variance_ratio(subset(JR_allcover2), "year", "code", "cover", 10, "uniqueID")

covtogPC1all <- right_join(JR_allcover2, siteout) %>%
  mutate(PC1cat = ifelse(PC1 > 0, "high", "low"),
         PC2cat = ifelse(PC2 > 0, "high", "low")) %>%
  group_by(year, PC1cat, treatment, trtrep, plot, uniqueID) %>%
  summarize(cover = sum(cover))

variance_ratio(covtogPC1all, "year", "PC1cat", "cover", 10, "uniqueID")
variance_ratio(subset(covtogPC1all, year <= 2002), "year", "PC1cat", "cover", 10, "uniqueID")
variance_ratio(subset(covtogPC1all, year > 2002), "year", "PC1cat", "cover", 10, "uniqueID")


covtogPC2all <- right_join(JR_allcover2, siteout) %>%
  mutate(PC1cat = ifelse(PC1 > 0, "high", "low"),
         PC2cat = ifelse(PC2 > 0, "high", "low")) %>%
  group_by(year, PC2cat, treatment, trtrep, plot, uniqueID) %>%
  summarize(cover = sum(cover))


variance_ratio(covtogPC2all, "year", "PC2cat", "cover", 10, "uniqueID")
variance_ratio(subset(covtogPC2all, year <= 2002), "year", "PC2cat", "cover", 10, "uniqueID")
variance_ratio(subset(covtogPC2all, year > 2002), "year", "PC2cat", "cover", 10, "uniqueID")

ggplot(covtogPC2all, aes(x=year, y=cover, color = PC2cat)) + geom_point() + facet_wrap(~PC2cat)


JR_rain <- read_csv("~/Dropbox/California Data/jrg_prism.csv") %>%
  select(-X1) %>%
  mutate(precip = ppt)

covPC1rain <- left_join(covtogPC1all, JR_rain) %>%
  mutate(plotrep = paste(treatment, trtrep, sep ="_"))

ggplot(covPC1rain, aes(x=growing_season_ppt, y=cover)) + geom_point() + facet_wrap(~PC1cat) + 
  geom_smooth(method = "lm")

covPC2rain <- left_join(covtogPC2all, JR_rain) %>%
  mutate(plotrep = paste(treatment, trtrep, sep ="_"))

ggplot(covPC2rain, aes(x=growing_season_ppt, y=cover)) + geom_point() + facet_wrap(~PC2cat) + 
  geom_smooth(method = "lm")


covPC1rain <- left_join(covtogPC1all, JR_rain)

ggplot(covPC1rain, aes(x=precip, y=cover)) + geom_point() + facet_wrap(~PC1cat) + 
  geom_smooth(method = "lm")


ggplot(covPC2rain, aes(x=growing_season_ppt, y=cover, color = treatment)) + geom_point() + facet_wrap(~PC2cat) + 
  geom_smooth(method = "lm") 

library(nlme)


a<-lme(cover ~ growing_season_ppt*treatment, random=~1|plotrep/uniqueID, 
       data = subset(covPC2rain, PC2cat == "high"), na.action = na.omit)
summary(a)



a<-lme(cover ~ growing_season_ppt + treatment, random=~1|plotrep/uniqueID, 
       data = subset(covPC1rain, PC1cat == "high"), na.action = na.omit)
summary(a)


# 
# 
# 
# mytrts <- names(siteout)[-3]
# results_onetrt <- as.data.frame(cbind(matrix(nrow = 0, ncol=10)))
# names(results_onetrt) = c("nbsp",      "sing.sp",   "FRic",      "qual.FRic", "FEve",    "FDis",  "RaoQ",     "Trait",     "plotyr",    "trait")
# for (i in 1:length(mytrts)){
#   PCtraitssub <- as.data.frame(PCtraits[,i])
#   names(PCtraitssub)="Trait"
#   row.names(PCtraitssub) <- row.names(PCtraits)
#   results_onetrt0<-dbFD(PCtraitssub, abundance, corr="cailliez")
#   results_onetrt0<-data.frame(results_onetrt0)
#   results_onetrt0$plotyr = row.names(results_onetrt0)
#   results_onetrt0$trait <- mytrts[i]
#   results_onetrt <- rbind(results_onetrt, results_onetrt0)
# }
