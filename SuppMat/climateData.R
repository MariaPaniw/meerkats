
### Paniw et al. XXXXX

### R script to construct rainfall and temperature deviations from raw data 

# load necessary packages


library(lubridate)
library(dplyr)
library(stringr)
library(ggplot2)
library(plyr)
library(scales)

# Function to calculate number of days in a month

numberOfDays <- function(date) {
  m <- format(date, format="%m")
  
  while (format(date, format="%m") == m) {
    date <- date + 1
  }
  
  return(as.integer(format(date - 1, format="%d")))
}

## RAINFALL
## current conditions 1997-2006 used to fit vital-rate GAMs (see meerkat IPMs)

rainfall=read.csv("/Users/maria/Dropbox/Meerkats/SuppMat/rainfall_GPCP.csv") # average daily rainfall

rainfall$date=as.Date(rainfall$date,format="%m.%d.%y") 
rainfall$month=month(rainfall$date)
rainfall$year=year(rainfall$date)

rain.sub=rainfall[rainfall$lon==21.25&rainfall$lat==(-26.25),] # retain rainfall values for the appropriate coordinates
rain.sub=rain.sub[rain.sub$date>"1996-11-01"&rain.sub$date<"2017-02-01",] 

rain.sub$rainSUMlag=NA # empty vector to save lagged (1.5 months) summed rainfall

for(i in 2:nrow(rain.sub)){
  
  rain.sub$rainSUMlag[i]=30.4*rain.sub$rain[i-1]+15*rain.sub$rain[i]
  
}

rain.sub$rainSUMlag=as.numeric(scale(rain.sub$rainSUMlag)) # scale to get standardized deviation 

mean.rain=aggregate(rainSUMlag~month,data=rain.sub,mean,na.rm=T) # monthly mean deviation

colnames(mean.rain)[2]="rain.mu"

rain.sub$rainMean=left_join(rain.sub,mean.rain,by="month")$rain.mu
rain.sub$rainSD=rain.sub$rainSUMlag-rain.sub$rainMean # interannual monthly value substracted from montly value averaged over all years (1997-2016)

plot(rain.sub$date,rain.sub$rainSD,type="b")


## Rainfall (RCP scenarios 2017-2066) used to project IPMs

# note that rainfall rainfall simulation prior to 2006 are considered historic simulations,
# while simulation from 2006 on apply the Representative Concentration Pathways (RCPs), 
# which provide concentrations of atmospheric greenhouse gas (GHG) and the trajectory that is taken over time to reach those concentrations. 
# These RCPs are named according to the level of radiative forcing (enhanced greenhouse effect or warming) that they produce by the year 2100.
# For more detail, see https://gisclimatechange.ucar.edu/gis-data-ar5

rainfall.rcp=read.csv("/Users/maria/Dropbox/Meerkats/SuppMat/rainfall_RCP_scen.csv")
rainfall.rcp$date=as.Date(rainfall.rcp$date,format="%d.%m.%y") 

# Go through scenarios
scen=c("RCP2.6","RCP4.5","RCP6.0","RCP8.5")

proj.sd=NULL
for(i in 1:length(scen)){
 
  # get historic mean of monthly lagged value
 hist=rainfall.rcp[rainfall.rcp$date<"2017-01-15"&rainfall.rcp$scen==scen[i],] 
 hist$rainSUMlag=NA
 
 for(j in 2:nrow(hist)){
   
   hist$rainSUMlag[j]=30.4*(hist$rain[j-1]/numberOfDays(hist$date[j-1]))+15*(hist$rain[j]/numberOfDays(hist$date[j]))

 }
 
 hist$rainSUMlag=as.numeric(scale(hist$rainSUMlag))
 
 # Create new average rainfall: 
 mean.hist=aggregate(rainSUMlag~month,data=hist,mean,na.rm=T)
 
 colnames(mean.hist)[2]="rain.mu"
 
 # join average monthly values to main data frame
 hist$rainMean=left_join(hist,mean.hist,by="month")$rain.mu
 hist$rainSD=hist$rainSUMlag-hist$rainMean
 
 ### Do the same for future projection (2017 -2066)
 
 proj=rainfall.rcp[rainfall.rcp$date>="2016-12-14"&rainfall.rcp$scen==scen[i],] 
 proj$rainSUMlag=NA
 
 for(j in 2:nrow(proj)){
   
   proj$rainSUMlag[j]=30.4*(proj$rain[j-1]/numberOfDays(proj$date[j-1]))+15*(proj$rain[j]/numberOfDays(proj$date[j]))
   
 }
 
 proj$rainSUMlag=as.numeric(scale(proj$rainSUMlag))
 
 # join average monthly values from historic conditions
 proj$rainMean=left_join(proj,mean.hist,by="month")$rain.mu
 
 # lastly, calculate standardized deviations 
 
 proj$rainSD=proj$rainSUMlag-proj$rainMean
  
 proj.sd=rbind(proj.sd,hist[-which(is.na(hist$rainSD)),],proj[-which(is.na(proj$rainSD)),])
}

### To make deviations comparable between rainfall values used to estimate GAM and ones from the RCP future projections, 
# we simply calculate relative deviations as follows

current=proj.sd[proj.sd$date<"2017-01-15",]
future=proj.sd[proj.sd$date>"2017-01-14"&proj.sd$date<"2077-01-14",]

stoch.month.Rain=array(NA,c(20,12,50,length(scen)))

for(x in 1:length(scen)){
  
  for(i in 1:12){
    
    sub=current[current$month==i&current$scen==scen[x],] 
    
    # for each current month (20 in total), find deviations from 2017-2066
    for(j in 1:nrow(sub)){
      
      prop=future[future$month==i&future$scen==scen[x],"rainSD"]-sub$rainSD[j]
      stoch.month.Rain[j,i,,x]=prop[-51] 
      
    }
    
  }
}

## plot 

dev=adply(stoch.month.Rain,c(1,2,3,4))
colnames(dev)=c("currentY","month","futureY","scen","rel.dev")
levels(dev$currentY)=1997:2016
levels(dev$scen)=scen
levels(dev$futureY)=2017:2066
dev$date=as.Date(paste(dev$futureY,dev$month,"15",sep="-"))

ggplot(data=dev,aes(date,rel.dev,col=currentY))+
  geom_line()+
  facet_wrap(~scen)+
  theme_bw()+
  xlab("Time") +
  ylab("Relative rainfall deviation")+
  theme(panel.grid = element_blank())+
  theme(axis.text.y = element_text(size=18,colour="black"))+
  theme(axis.title = element_text(size=25))+
  theme(axis.title.x = element_text(vjust=-0.4),
        axis.title.y = element_text(vjust=0.9),
        plot.title = element_text(size=22,lineheight=.8, face="bold"))+
  
  theme(plot.margin = unit(c(0.2,0.1,0.5,0.1), "cm"))+
  scale_x_date(date_breaks = "12 month", labels = date_format("%Y"), expand = c(0, 0))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=12,color="black"))

#### plot rainfall deviation used in population projections (adding relative deviation to observed deviation)
rain.sub$year=factor(rain.sub$year)
dev$month=as.numeric(dev$month)
dev$obs.dev=left_join(dev,rain.sub,by=c("currentY"="year","month"="month"))$rainSD
dev$future.dev=dev$obs.dev+dev$rel.dev

ggplot(data=dev,aes(date,future.dev,col=currentY))+
  geom_line()+
  facet_wrap(~scen)+
  theme_bw()+
  xlab("Time") +
  ylab("Projected rainfall deviation")+
  theme(panel.grid = element_blank())+
  theme(axis.text.y = element_text(size=18,colour="black"))+
  theme(axis.title = element_text(size=25))+
  theme(axis.title.x = element_text(vjust=-0.4),
        axis.title.y = element_text(vjust=0.9),
        plot.title = element_text(size=22,lineheight=.8, face="bold"))+
  
  theme(plot.margin = unit(c(0.2,0.1,0.5,0.1), "cm"))+
  scale_x_date(date_breaks = "12 month", labels = date_format("%Y"), expand = c(0, 0))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=12,color="black"))

## MAXIMUM TEMPERATURE

## current conditions 1997-2006 used to fit vital-rate GAMs (see meerkat IPMs)

temperature=read.csv("/Users/maria/Dropbox/Meerkats/SuppMat/tempNOAA_CPC.csv") # average daily temperature

temperature$date=as.Date(temperature$date) 
temperature$month=month(temperature$date)
temperature$year=year(temperature$date)

sub=temperature[temperature$date>="1996-12-01"&temperature$date<"2017-02-05",] 

start=as.Date("01/15/1997",format="%m/%d/%Y")
end=as.Date("01/15/2017",format="%m/%d/%Y")

dates=seq(start, end , by="month")

tempMAXlag=data.frame(lag=rep(NA,length(dates)),
                      date=rep(NA,length(dates)))# empty data frame to save lagged (1.5 months) mean maximum temperatures



for(i in 1:length(dates)){
  
  tempMAXlag$lag[i]=mean(sub$temp[sub$date<=dates[i]&sub$date>=(dates[i]-months(1)-days(14))]) # mean maximum temperature
  tempMAXlag$date[i]=as.character(dates[i])
}

tempMAXlag$lag=as.numeric(scale(tempMAXlag$lag)) # scale 
tempMAXlag$month=month(as.Date(tempMAXlag$date))
tempMAXlag$year=year(as.Date(tempMAXlag$date))
mean.max=aggregate(lag~month,data=tempMAXlag,mean,na.rm=T) # average monthly values

colnames(mean.max)[2]="max.mu"

tempMAXlag$maxMean=left_join(tempMAXlag,mean.max,by="month")$max.mu
tempMAXlag$maxSD=tempMAXlag$lag-tempMAXlag$maxMean # standardized deviation in maximun temperature values
tempMAXlag$date=as.Date(tempMAXlag$date)


plot(tempMAXlag$date,tempMAXlag$maxSD,type="b")

## temperature (RCP scenarios 2017-2066) used to project IPMs

# note that temperature temperature simulation prior to 2006 are considered historic simulations,
# while simulation from 2006 on apply the Representative Concentration Pathways (RCPs), 
# which provide concentrations of atmospheric greenhouse gas (GHG) and the trajectory that is taken over time to reach those concentrations. 
# These RCPs are named according to the level of radiative forcing (enhanced greenhouse effect or warming) that they produce by the year 2100.
# For more detail, see https://gisclimatechange.ucar.edu/gis-data-ar5

temperature.rcp=read.csv("/Users/maria/Dropbox/Meerkats/SuppMat/temperature_RCP_scen.csv")
temperature.rcp$date=as.Date(temperature.rcp$date,format="%d.%m.%y") 

# Go through scenarios
scen=c("RCP2.6","RCP4.5","RCP6.0","RCP8.5")

proj.sd=NULL
for(i in 1:length(scen)){
  
  # get historic mean of monthly lagged value
  hist=temperature.rcp[temperature.rcp$date<"2017-01-15"&temperature.rcp$scen==scen[i],] 
  hist$tempMAXlag=NA
  
  for(j in 2:nrow(hist)){
    
    hist$tempMAXlag[j]=mean(hist$temp[j-1], hist$temp[j])
    
  }
  
  hist$tempMAXlag=as.numeric(scale(hist$tempMAXlag))
  
  # Create new average temperature: 
  mean.hist=aggregate(tempMAXlag~month,data=hist,mean,na.rm=T)
  
  colnames(mean.hist)[2]="temp.mu"
  
  # join average monthly values to main data frame
  hist$tempMean=left_join(hist,mean.hist,by="month")$temp.mu
  hist$tempSD=hist$tempMAXlag-hist$tempMean
  
  ### Do the same for future projection (2017 -2066)
  
  proj=temperature.rcp[temperature.rcp$date>="2016-12-14"&temperature.rcp$scen==scen[i],] 
  proj$tempMAXlag=NA
  
  for(j in 2:nrow(proj)){
    
    proj$tempMAXlag[j]=mean(proj$temp[j-1], proj$temp[j])
    
  }
  
  proj$tempMAXlag=as.numeric(scale(proj$tempMAXlag))
  
  # join average monthly values from historic conditions
  proj$tempMean=left_join(proj,mean.hist,by="month")$temp.mu
  
  # lastly, calculate standardized deviations 
  
  proj$tempSD=proj$tempMAXlag-proj$tempMean
  
  proj.sd=rbind(proj.sd,hist[-which(is.na(hist$tempSD)),],proj[-which(is.na(proj$tempSD)),])
}

### To make deviations comparable between temperature values used to estimate GAM and ones from the RCP future projections, 
# we simply calculate relative deviations as follows

current=proj.sd[proj.sd$date<"2017-01-15",]
future=proj.sd[proj.sd$date>"2017-01-14"&proj.sd$date<"2077-01-14",]

stoch.month.temp=array(NA,c(20,12,50,length(scen)))

for(x in 1:length(scen)){
  
  for(i in 1:12){
    
    sub=current[current$month==i&current$scen==scen[x],] 
    
    # for each current month (20 in total), find deviations from 2017-2066
    for(j in 1:nrow(sub)){
      
      prop=future[future$month==i&future$scen==scen[x],"tempSD"]-sub$tempSD[j]
      stoch.month.temp[j,i,,x]=prop[-51] 
      
    }
    
  }
}

## plot 

dev=adply(stoch.month.temp,c(1,2,3,4))
colnames(dev)=c("currentY","month","futureY","scen","rel.dev")
levels(dev$currentY)=1997:2016
levels(dev$scen)=scen
levels(dev$futureY)=2017:2066
dev$date=as.Date(paste(dev$futureY,dev$month,"15",sep="-"))

ggplot(data=dev,aes(date,rel.dev,col=currentY))+
  geom_line()+
  facet_wrap(~scen)+
  theme_bw()+
  xlab("Time") +
  ylab("Relative temperature deviation")+
  theme(panel.grid = element_blank())+
  theme(axis.text.y = element_text(size=18,colour="black"))+
  theme(axis.title = element_text(size=25))+
  theme(axis.title.x = element_text(vjust=-0.4),
        axis.title.y = element_text(vjust=0.9),
        plot.title = element_text(size=22,lineheight=.8, face="bold"))+
  
  theme(plot.margin = unit(c(0.2,0.1,0.5,0.1), "cm"))+
  scale_x_date(date_breaks = "12 month", labels = date_format("%Y"), expand = c(0, 0))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=12,color="black"))

#### plot rainfall deviation used in population projections (adding relative deviation to observed deviation)
tempMAXlag$year=factor(tempMAXlag$year)
dev$month=as.numeric(dev$month)
dev$obs.dev=left_join(dev,tempMAXlag,by=c("currentY"="year","month"="month"))$maxSD
dev$future.dev=dev$obs.dev+dev$rel.dev

ggplot(data=dev,aes(date,future.dev,col=currentY))+
  geom_line()+
  facet_wrap(~scen)+
  theme_bw()+
  xlab("Time") +
  ylab("Projected temperature deviation")+
  theme(panel.grid = element_blank())+
  theme(axis.text.y = element_text(size=18,colour="black"))+
  theme(axis.title = element_text(size=25))+
  theme(axis.title.x = element_text(vjust=-0.4),
        axis.title.y = element_text(vjust=0.9),
        plot.title = element_text(size=22,lineheight=.8, face="bold"))+
  
  theme(plot.margin = unit(c(0.2,0.1,0.5,0.1), "cm"))+
  scale_x_date(date_breaks = "12 month", labels = date_format("%Y"), expand = c(0, 0))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=12,color="black"))

