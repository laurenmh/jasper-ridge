library(tidyverse)
library(lubridate)

## Read in the jrn_prism data
jrg_prism <- read_csv("~/Dropbox/California Data/PRISM_ppt_tmin_tmean_tmax_tdmean_vpdmin_vpdmax_provisional_4km_19810801_20190401_37.4040_-122.2228.csv",
                      skip = 10) %>%
  tbl_df() 
names(jrg_prism) = c("Date", "pptin", "tminF", "tmeanF", "tmaxF", "tdmean", "vpdmin", "vpdmax")
# Create year and month columns to calculate growing season values for temp and ppt
jrg_prism$year<-as.numeric(format(as.Date(jrg_prism$Date),"%Y"))
jrg_prism$month<-as.numeric(format(as.Date(jrg_prism$Date),"%m"))
jrg_prism$ppt <- jrg_prism$pptin*25.4
jrg_prism$tmin <- (jrg_prism$tminF -32)*(5/9)
jrg_prism$tmean <- (jrg_prism$tmeanF -32)*(5/9)
jrg_prism$tmax <- (jrg_prism$tmaxF -32)*(5/9)

seasonal_clim <- jrg_prism %>%
  filter(month%in%c(9,10,11,12,1,2,3,4)) %>%
  mutate(year = ifelse(month > 8, year + 1, year)) %>%
  mutate(season = ifelse(month%in%c(9,10,11), "fall", "winter"),
         season = ifelse(month%in%c(3,4), "spring", season)) %>%
  group_by(season, year) %>%
  summarize(ppt = sum(ppt),
            mintmin = min(tmin),
            meantmin = mean(tmin),
            tmean = mean(tmean),
            maxtmax = max(tmax),
            meantmax = mean(tmax)) %>%
  gather(variable, value, mintmin:meantmax, ppt) %>%
  tbl_df() %>%
  mutate(season_var = paste(season, variable, sep = "_")) %>%
  select(-season, -variable) %>%
  spread(season_var, value)

growing_season <- jrg_prism %>%
  filter(month%in%c(9,10,11,12,1,2,3,4)) %>%
  mutate(year = ifelse(month > 8, year + 1, year)) %>%
  mutate(season = ifelse(month%in%c(9,10,11), "fall", "winter"),
         season = ifelse(month%in%c(3,4), "spring", season)) %>%
  group_by(year) %>%
  summarize(ppt = sum(ppt),
            mintmin = min(tmin),
            meantmin = mean(tmin),
            tmean = mean(tmean),
            maxtmax = max(tmax),
            meantmax = mean(tmax))

annual <- jrg_prism %>%
  mutate(year = ifelse(month > 6, year + 1, year)) %>%
  group_by(year) %>%
  summarize(ppt = sum(ppt))

JR_rain <- left_join(growing_season, seasonal_clim) %>%
  mutate(precip = ppt)

write_csv(JR_rain, "JR_rain.csv")


## play around with date of first rain, etc
daily_season <- jrg_prism %>%
  filter(month%in%c(9,10,11,12,1,2,3,4)) %>%
  mutate(year = ifelse(month > 8, year + 1, year)) %>%
  mutate(season = ifelse(month%in%c(9,10,11), "fall", "winter"),
         season = ifelse(month%in%c(3,4), "spring", season)) %>%
  arrange(Date) %>%
  group_by(year) %>%
  mutate(raintally = cumsum(ppt)) %>%
  mutate(daymonth = yday(Date)) %>%
  mutate(daymonth = ifelse(month < 8, daymonth + 365, daymonth),
         daymonth = daymonth - 244)

firstrain <- daily_season %>%
  filter(raintally > 0) %>%
  group_by(year) %>%
  mutate(firstrain = min(Date)) %>%
  mutate(tokeep = ifelse(raintally == ppt, 1, 0)) %>%
  filter(tokeep == 1)

ggplot(firstrain, aes(x=year, y=daymonth)) + geom_point() + geom_line()
ggplot(firstrain, aes(x=daymonth, y = ppt)) + geom_point()

# rain accumulation curves
ggplot(subset(daily_season), aes(x=daymonth, y=raintally, color = as.factor(year))) + geom_line()
ggplot(subset(daily_season), aes(x=daymonth, y=raintally)) + geom_line() + facet_wrap(~year)
