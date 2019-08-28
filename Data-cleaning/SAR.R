source("Data-cleaning/cover-cleaning.R")
source("Data-cleaning/cover-1m2-cleaning.R")

library(vegan)

JRcoverall <- JRcover 

JRcover <- JRcover %>%
  filter(species !="ROCK", species != "BARE") %>%
  tbl_df()

JRcover$trtreptog <- paste(JRcover$treatment, JRcover$trtrep, sep = "")
JRcover$trtrepyear <- paste(JRcover$treatment, JRcover$trtrep, JRcover$year, sep = "_")


SAR <- function(dat, trtrepyear, quadID, treatment, species, cover, year, trtreptog){
subber <- dat %>%
  arrange(quadID) %>%
  select(quadID, species, cover) %>%
  spread(species, cover, fill = 0)
subquads <- subber$quadID
subber$quadID <-NULL 
sp1 <- specaccum(subber)
plotarea <- seq(.25, 6, by=.25)
#plot(log(plotarea), log(sp1$richness), add = TRUE, col = "red")

intercept <- coef(lm(log(sp1$richness)~log(plotarea)))[1]
slope <- coef(lm(log(sp1$richness)~log(plotarea)))[2]
year <- unique(dat$year)
trtrep <- unique(dat$trtreptog)
treatment <- unique(dat$treatment)
out <- cbind(trtrep, treatment, year, intercept, slope)
names(out) = c("plot", "treatment", "year", "intercept", "slope")
outout <- data.frame(out)
return(outout)
}


X <- split(JRcover, JRcover["trtrepyear"])

out <- lapply(X, SAR)
reps <- unique(JRcover["trtrepyear"])
output <- cbind(reps, do.call("rbind", out)) %>%
  tbl_df() %>%
  mutate(intercept = as.numeric(as.character(intercept)),
         slope = as.numeric(as.character(slope)),
         year = as.numeric(as.character(year)))



ggplot(output, aes(x=treatment, y=intercept)) + geom_boxplot()
ggplot(output, aes(x=treatment, y=slope)) + geom_boxplot()


ggplot(output, aes(x=as.factor(year), y=intercept)) + geom_boxplot()
ggplot(output, aes(x=as.factor(year), y=slope)) + geom_boxplot()
ggplot(output, aes(x=as.factor(year), y=slope/intercept)) + geom_boxplot()

JRgopher2 <- JRgopher %>%
  group_by(year, treatment, trtrep) %>%
  summarize(disturb = mean(disturb)) %>%
  mutate(trtrep = paste(treatment, trtrep, sep = ""))

JRsoil2 <- JRsoil %>%
  group_by(treatment, trtrep) %>%
  summarize(depth = mean(depth)) %>%
  mutate(trtrep = paste(treatment, trtrep, sep = ""))

JRrain <- read_csv("~/Dropbox/California Data/jrg_prism.csv") %>%
  select(-X1)
  
  
  
output2 <- left_join(output, JRgopher2)
output3 <- left_join(output2, JRsoil2)
output4 <- left_join(output3, JRrain)


ggplot(output4, aes(x=disturb, y=slope)) + geom_point() + geom_smooth(se=F, method = 'lm')
ggplot(output4, aes(x=disturb, y=intercept)) + geom_point() + geom_smooth(se=F, method = 'lm')


ggplot(output4, aes(x=depth, y=slope)) + geom_point() + geom_smooth(se=F, method = 'lm')
ggplot(output4, aes(x=depth, y=exp(intercept))) + geom_point(aes(color = treatment)) + geom_smooth(se=F, method = 'lm')


ggplot(output4, aes(x=ppt, y=slope)) + geom_point() + geom_smooth(se=F, method = 'lm')
ggplot(output4, aes(x=ppt, y=intercept)) + geom_point(aes(color = treatment)) + geom_smooth(se=F, method = 'lm')


l <- lm(intercept~disturb, data = output4)
summary(l)

l <- lm(slope~disturb, data = output4)
summary(l)

l <- lm(intercept~ppt, data = output4)
summary(l)

l <- lm(slope~ppt, data = output4)
summary(l)


l <- lm(intercept~depth, data = output4)
summary(l)

l <- lm(slope~depth, data = output4)
summary(l)

l <- lm(intercept ~ depth + disturb + ppt, data = output4 )
summary(l)

l <- lm(slope ~ depth + disturb + ppt, data = output4 )
summary(l)


### rare species in aggregate
JRcoverbyrare <- JRcover %>%
  group_by(species) %>%
  mutate(totalmean = mean(cover)) %>%
  mutate(cat = ifelse(totalmean > 3, "dominant", "rare")) %>%
  group_by(quadID, year, treatment, trtrep, subplot, trtrepyear, cat) %>%
  summarize(cover = sum(cover)) %>%
  group_by(year, treatment, trtrep, cat) %>%
  summarize(meancover = mean(cover))

ggplot(JRcoverbyrare, aes(x=year, y=meancover, color = cat, group = interaction(cat, treatment, trtrep))) + geom_point()+ geom_line()  +facet_wrap(trtrep~treatment)

JRcoverbyrare_rain <- left_join(JRcoverbyrare, JRrain)
ggplot(JRcoverbyrare_rain, aes(x=ppt, y=meancover, color = cat)) + geom_point()  +
  facet_grid(cat~treatment, scale = "free") + geom_smooth(method = "lm")



### now output the actual curves

SARcurve <- function(dat, trtrepyear, quadID, treatment, species, cover, year, trtreptog){
  subber <- dat %>%
    arrange(quadID) %>%
    select(quadID, species, cover) %>%
    spread(species, cover, fill = 0)
  subquads <- subber$quadID
  subber$quadID <-NULL 
  sp1 <- specaccum(subber)
  plotarea <- seq(.25, 6, by=.25)
  #plot(log(plotarea), log(sp1$richness), add = TRUE, col = "red")
  richness <- sp1$richness
  sd <- sp1$sd
  
  out <- as.data.frame(cbind(plotarea, richness, sd))
  out$year <- unique(dat$year)
  out$trtrep <- unique(dat$trtreptog)
  out$treatment <- unique(dat$treatment)
  #out$trtrepyear <- unique(dat$trtrepyear)
  return(out)
}


X <- split(JRcover, JRcover["trtrepyear"])

out <- lapply(X, SARcurve)
out2 <- unlist(out)
reps <- unique(JRcover["trtrepyear"])
outputcurve <- cbind(reps, do.call("rbind", out)) %>%
  tbl_df()


outputcurve2 <- left_join(outputcurve, JRgopher2)
outputcurve3 <- left_join(outputcurve2, JRsoil2)
outputcurve4 <- left_join(outputcurve3, JRrain)

curvedepth <- outputcurve4 %>%
  mutate(depthquantile = ntile(depth, 3)) %>%
  group_by(depthquantile, plotarea) %>%
  summarize(rich = mean(richness), serich = sd(richness)/sqrt(length(richness)))

ggplot(curvedepth, aes(x=plotarea, y=rich, color = depthquantile, group = as.factor(depthquantile))) + geom_point() + geom_line() + 
  geom_errorbar(aes(ymin = rich - serich, ymax = rich + serich))

ggplot(curvedepth, aes(x=log(plotarea), y=log(rich), color = depthquantile, group = as.factor(depthquantile))) + geom_point() + geom_line() +
  geom_errorbar(aes(ymin = log(rich - serich), ymax = log(rich + serich)))

curvedist<- outputcurve4 %>%
  mutate(disturbquantile = ntile(depth, 3)) %>%
  group_by(disturbquantile, plotarea) %>%
  summarize(rich = mean(richness), serich = sd(richness)/sqrt(length(richness)))

ggplot(curvedist, aes(x=plotarea, y=rich, color = disturbquantile, group = as.factor(disturbquantile))) + geom_point() + geom_line() + 
  geom_errorbar(aes(ymin = rich - serich, ymax = rich + serich))

ggplot(curvedist, aes(x=log(plotarea), y=log(rich), color = disturbquantile, group = as.factor(disturbquantile))) + geom_point() + geom_line() +
  geom_errorbar(aes(ymin = log(rich - serich), ymax = log(rich + serich)))


curverain <- outputcurve4 %>%
  filter(!is.na(growing_season_ppt)) %>%
  mutate(rainquantile = ntile(growing_season_ppt, 3)) %>%
  group_by(rainquantile, plotarea) %>%
  summarize(rich = mean(richness), serich = sd(richness)/sqrt(length(richness))) 

ggplot(curverain, aes(x=(plotarea), y=(rich), color = rainquantile, group = as.factor(rainquantile))) + geom_point() + geom_line() +
  geom_errorbar(aes(ymin = (rich - serich), ymax = (rich + serich)))


ggplot(curverain, aes(x=log(plotarea), y=log(rich), color = rainquantile, group = as.factor(rainquantile))) + geom_point() + geom_line() +
  geom_errorbar(aes(ymin = log(rich - serich), ymax = log(rich + serich)))




curvetime <- outputcurve4 %>%
  mutate(yearquantile = ntile(year, 3)) %>%
  group_by(yearquantile, plotarea) %>%
  summarize(rich = mean(richness), serich = sd(richness)/sqrt(length(richness))) 

ggplot(curvetime, aes(x=(plotarea), y=(rich), color = yearquantile, group = as.factor(yearquantile))) + geom_point() + geom_line() +
  geom_errorbar(aes(ymin = (rich - serich), ymax = (rich + serich)))


ggplot(curvetime, aes(x=log(plotarea), y=log(rich), color = yearquantile, group = as.factor(yearquantile))) + geom_point() + geom_line() +
  geom_errorbar(aes(ymin = log(rich - serich), ymax = log(rich + serich)))
