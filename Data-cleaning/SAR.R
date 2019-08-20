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

  
output2 <- left_join(output, JRgopher2)
ggplot(output2, aes(x=disturb, y=slope)) + geom_point() + geom_smooth(se=F, method = 'lm')
ggplot(output2, aes(x=disturb, y=intercept)) + geom_point() + geom_smooth(se=F, method = 'lm')

output3 <- left_join(output2, JRsoil2)
ggplot(output3, aes(x=depth, y=slope)) + geom_point() + geom_smooth(se=F, method = 'lm')
ggplot(output3, aes(x=depth, y=intercept)) + geom_point(aes(color = treatment)) + geom_smooth(se=F, method = 'lm')
