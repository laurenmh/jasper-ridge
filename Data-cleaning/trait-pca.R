library(tidyverse)
library(readr)
library(vegan)

trait.dat <- read_csv("~/Dropbox/California Data/trait_all.csv") 
names(trait.dat)
names(trait.dat) = c("species", "code", "status", "habit", "longevity", "Nfixer", "rootDepth",
               "SRL", "LMA", "leafN", "A", "TcorA", "PS2", "WUE", "height", "germDays", "RS",
               "seed", "NUE", "NUEadj")

trait.dat2 <- trait.dat %>%
  filter(species%in%c("Elymus multisedus", "Layia platyglossa","Calycadenia multiglandulosa",
                      "Nassella pulchra", "Microseris douglasii", "Castilleja densiflora",
                      "Vulpia microstachys var pauciflora", "Hemizonia congesta",
                      "Lotus wrangelianus" ,  "Trifolium albopurpureum" , "Dichelostemma capitatum",
                      "Triteleia laxa","Lolium multiflorum", "Plantago erecta","Chlorogalum pomeridianum var. pomeridianum",
                      "Bromus hordaceous", "Astragalus gambelianus", "Lasthenia californica", "Hesperevax sparsiflora" ,
                      "Crassula connata")) %>%
  filter(species != "Crassula connata") %>% 
  mutate(RS = ifelse(species == "Triteleia laxa", 3.68009626, RS)) %>%
  mutate(SRL = ifelse(species == "Dichelostemma capitatum", 0.02259358, RS)) %>%
  mutate(NUE = as.numeric(as.character(NUE)))%>%
  select(-code, -status, -habit, -longevity, -Nfixer, -A, -TcorA, -germDays, -seed) 
 
# filter(!species%in%c("Lepidium nitidum", "Dodecatheon hendersonii", "Anagalis arvensis", "Hordeum marinum",
#                       "Hordeum brachyantherum ssp. Californicum"))

  # "Micropus californicus", "Poa secunda", "Eschscholzia californica" , "Calandrinia ciliata", "Gilia clivorum"   
  
# matrix for PCA
traits <- as.matrix(trait.dat2[,c(2:ncol(trait.dat2))])
row.names(traits) <- trait.dat2$species

# run PCA
myrda <- rda(traits, scale = TRUE)

# extract values
siteout <- as.data.frame(scores(myrda, choices=c(1,2), display=c("sites")))
siteout$species<-rownames(siteout)
siteout$name <- siteout$species

enviroout<-as.data.frame(scores(myrda, choices=c(1,2), display=c("species")))
enviroout$type<-"traits"
enviroout$name<-rownames(enviroout)

# merge PC axes with trait data
tog <- left_join(trait.dat, siteout) %>%
  mutate(func = paste(status, habit, sep = "_"))



#pdf("TraitPCA.pdf", width = 9, height = 7.5)
#pdf("TraitPCA_noLegumes.pdf", width = 9, height = 7.5)

ggplot(tog, aes(x=PC1, y=PC2))+ 
  geom_hline(aes(yintercept=0), color="grey") + 
  geom_vline(aes(xintercept=0), color="grey") +
  geom_text(aes(label = name, color = func), size = 5) +
  # scale_color_manual(values = c("grey20", "grey70")) +
  geom_segment(data = enviroout,
               aes(x = 0, xend =  PC1,
                   y = 0, yend =  PC2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "black") + #grid is required for arrow to work.
  geom_text(data = enviroout,
            aes(x=  PC1*1.2, y =  PC2*1.2, #we add 10% to the text to push it slightly out from arrows
                label = name), #otherwise you could use hjust and vjust. I prefer this option
            size = 6,
            hjust = 0.5, 
            color="black") + 
  theme_bw() +theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                    text=element_text(size = 20))+ 
  xlab(paste("Axis 1 (",sprintf("%.1f",myrda$CA$eig["PC1"]/myrda$tot.chi*100,3),"%)",sep="")) +
  ylab(paste("Axis 2 (",sprintf("%.1f",myrda$CA$eig["PC2"]/myrda$tot.chi*100,3),"%)",sep="")) 

#dev.off()

