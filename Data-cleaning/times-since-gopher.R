source("Data-cleaning/cover-cleaning.R")
library(tsvr)

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
  mutate(ncount = n()) %>%
  mutate(yrsdist = max(repnum))
 

tog <- left_join(JRcover, gophbigtime) 

