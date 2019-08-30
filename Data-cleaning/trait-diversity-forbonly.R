library(tidyverse)
library(readr)
library(vegan)
library(FD)

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
                       "Astragalus gambelianus", "Lasthenia californica", "Hesperevax sparsiflora")) %>%
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

traits <- as.matrix(traitdat[1:10, c(3:ncol(traitdat))])
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
  group_by(year) %>%
  summarize_all(funs(mean)) %>%
  mutate(year = as.numeric(year))

ggplot(PCresultagg, aes(x=year, y=CWM.PC1)) + geom_line() + geom_smooth(se=F, method = "lm")
ggplot(PCresultagg, aes(x=year, y=CWM.PC2)) + geom_line() + geom_smooth(se=F, method = "lm")



tog <- right_join(trait.dat, siteout) %>%
  mutate(func = paste(status, habit, sep = "_")) %>%
  mutate(func = tolower(func))


#pdf("TraitPCA-justforbs.pdf", width = 9, height = 7.5)
#pdf("TraitPCA_noLegumes.pdf", width = 9, height = 7.5)

ggplot(tog, aes(x=PC1, y=PC2))+ 
  geom_hline(aes(yintercept=0), color="grey") + 
  geom_vline(aes(xintercept=0), color="grey") +
  geom_text(aes(label = name, color = func), size = 4) +
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
                    text=element_text(size = 20))+ 
  xlab(paste("Axis 1 (",sprintf("%.1f",myrda$CA$eig["PC1"]/myrda$tot.chi*100,3),"%)",sep="")) +
  ylab(paste("Axis 2 (",sprintf("%.1f",myrda$CA$eig["PC2"]/myrda$tot.chi*100,3),"%)",sep="")) 
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
