
### Paniw et al. XXXXX

### R script to construct meerkat IPMs

# load necessary packages

library(lubridate)
library(boot)
library(plyr)
library(ggplot2)
library(Cairo)
library(scales)

# load Vital-rate GAMs 

file.names <- list.files(path = "/Users/maria/Dropbox/Meerkats/GAMs")

for(i in 1:length(file.names)){
  
  load(paste("/Users/maria/Dropbox/Meerkats/GAMs/",file.names[i],sep=""))
}


###################################################
# put together IPM

##### STEP 1: create functions
# x = mass (discretized, 100 bins)
# ageM = age in months
# rain/tempSD = rainfall/temperature standardized deviation (variation)
# density = population density
# pregCat = pregnancy category 
# massM = mass of mother

# Pup survival 
f.pup.surv=function(x,ageM,month,rainSD,tempSD,density){
  
  new.data=expand.grid(mass=x,ageM=ageM,month=month,rainSD=rainSD,
                       tempSD=tempSD,density=density)
  new.data$pred=predict(pup.surv,newdata = new.data,type="response")
  return(new.data)
}

# Pup growth 
f.pup.gr=function(x,y,ageM,month,rainSD,tempSD,density,year){
  
  new.data=expand.grid(mass=x,ageM=ageM,month=month,rainSD=rainSD,
                       tempSD=tempSD,density=density,year=year)
  mu=predict(pup.gr,newdata = new.data)
  var=2*exp(as.numeric(predict(pup.gr.var,newdata = new.data))) # note that we multiply by 2 here, as well as for the other relevant models, because model predictions consistently underestimated observed variance in size 
  temp1 <- sqrt(2*pi*var)
  temp2 <- ((y-mu)^2)/(2*var)
  
  return(exp(-temp2)/temp1)
}

# Juvenile survival 
f.juv.surv=function(x,density,rainSD,tempSD){
  
  new.data=expand.grid(mass=x,density=density,rainSD=rainSD,
                       tempSD=tempSD)
  new.data$pred=predict(juv.surv,newdata = new.data,type="response")
  return(new.data)
}

# Juvenile growth 
f.juv.gr=function(x,y,ageM,month,rainSD,tempSD,density,year){
  new.data=expand.grid(mass=x,ageM=ageM,month=month,rainSD=rainSD,
                       tempSD=tempSD,density=density,year=year)
  mu=predict(juv.gr,newdata = new.data)
  var=2*exp(as.numeric(predict(juv.gr.var,newdata = new.data)))
  temp1 <- sqrt(2*pi*var)
  temp2 <- ((y-mu)^2)/(2*var)
  
  return(exp(-temp2)/temp1)
}

# Subadult survival 
f.sub.surv=function(x,ageM,month,density,rainSD,tempSD){
  
  new.data=expand.grid(mass=x,ageM=ageM,month=month,rainSD=rainSD,
                       tempSD=tempSD,density=density)
  new.data$pred=predict(sub.surv,newdata = new.data,type="response")
  return(new.data)
}

# Subadult growth 
f.sub.gr=function(x,y,ageM,month,rainSD,tempSD,density,year){
  
  new.data=expand.grid(mass=x,ageM=ageM,month=month,rainSD=rainSD,
                       tempSD=tempSD,density=density,year=year)
  mu=predict(sub.gr,newdata = new.data)
  var=2*exp(as.numeric(predict(sub.gr.var,newdata = new.data)))
  temp1 <- sqrt(2*pi*var)
  temp2 <- ((y-mu)^2)/(2*var)
  
  return(exp(-temp2)/temp1)
}

# Adult helper survival 
f.help.surv=function(x,month,density,rainSD,tempSD,year){
  
  new.data=expand.grid(mass=x,month=month,rainSD=rainSD,
                       tempSD=tempSD,density=density,year=year)
  new.data$pred=predict(help.surv,newdata = new.data,type="response")
  return(new.data)
}

# Adult helper growth
f.help.gr=function(x,y,pregCat,month,rainSD,tempSD,density,year){
  
  new.data=expand.grid(mass=x,pregCat=pregCat,month=month,rainSD=rainSD,
                       tempSD=tempSD,density=density,year=year)
  mu=predict(help.gr,newdata = new.data)
  var=2*exp(as.numeric(predict(help.gr.var,newdata = new.data)))
  temp1 <- sqrt(2*pi*var)
  temp2 <- ((y-mu)^2)/(2*var)
  return(exp(-temp2)/temp1)
}

# Emigration 
f.help.emig=function(x,month,density,rainSD,tempSD,year){
  
  new.data=expand.grid(mass=x,month=month,rainSD=rainSD,
                       tempSD=tempSD,density=density,year=year)
  new.data$pred=predict(help.emig,newdata = new.data,type="response")
  return(new.data)
}

# Helper transition to dominant  
f.toDom=function(x,month,density,rainSD){
  
  new.data=expand.grid(mass=x,month=month,density=density,rainSD=rainSD)
  new.data$pred=predict(toDom,newdata = new.data,type="response")
  return(new.data)
}

# Non-pregnant helper transition to 1-month pregnant helper
f.help.NPH=function(x,month,density,rainSD,tempSD,year){
  
  new.data=expand.grid(mass=x,month=month,rainSD=rainSD,
                       tempSD=tempSD,density=density,year=year)
  new.data$pred=predict(help.transNPH,newdata = new.data,type="response")
  return(new.data)
}

# Non-pregnant helper transition to 1-month pregnant dominant
f.help.NPD=function(x,month){
  
  new.data=expand.grid(mass=x,month=month)
  new.data$pred=predict(help.transNPD,newdata = new.data,type="response")
  return(new.data)
}

# First month pregnant helper transition to 2-month pregnant helper or aborts
f.help.FH=function(x,month,density,rainSD,tempSD,year){
  
  new.data=expand.grid(mass=x,month=month,rainSD=rainSD,
                       tempSD=tempSD,density=density,year=year)
  new.data$pred=predict(help.transFH,newdata = new.data,type="response")
  return(new.data)
}

# If aborted at first month of pregnancy, becoming 1-month pregnant immediately or remaining non-pregnant
f.help.FH2=function(x,month,density,rainSD,tempSD){
  
  new.data=expand.grid(mass=x,month=month,rainSD=rainSD,
                       tempSD=tempSD,density=density)
  new.data$pred=predict(help.transFH2,newdata = new.data,type="response")
  return(new.data)
}

# Second month pregnant helper transition to  helper with litter or aborts
f.help.SH=function(x,month,density,year){
  
  new.data=expand.grid(mass=x,month=month,density=density,year=year)
  new.data$pred=predict(help.transSH,newdata = new.data,type="response")
  return(new.data)
}

# If aborted at 2 months of pregnancy, becoming 1-month pregnant immediately or remaining non-pregnant
f.help.SH2=function(x,month,density,rainSD,tempSD,year){
  
  new.data=expand.grid(mass=x,month=month,rainSD=rainSD,
                       tempSD=tempSD,density=density,year=year)
  new.data$pred=predict(help.transSH2,newdata = new.data,type="response")
  return(new.data)
}

# Helper with weaning litter becoming 1-month pregnant 
f.help.BH=function(x,month){
  
  new.data=expand.grid(mass=x,month=month)
  new.data$pred=predict(help.transBH,newdata = new.data,type="response")
  return(new.data)
}

# Number of pups born to helper mother 
f.help.pups=function(x,month,density,rainSD){
  
  new.data=expand.grid(mass=x,month=month,rainSD=rainSD,
                       density=density)
  new.data$pred=predict(help.pups,newdata = new.data,type="response")
  return(new.data)
}

# Mass of pups born to helper mother 
f.help.off.mass=function(x,y,month,rainSD,tempSD,density,year){
  
  new.data=expand.grid(massM=x,month=month,rainSD=rainSD,
                       tempSD=tempSD,density=density,year=year)
  mu=predict(help.off.mass,newdata = new.data)
  
  var=2*exp(as.numeric(predict(help.off.mass.var,newdata = new.data)))
  temp1 <- sqrt(2*pi*var)
  temp2 <- ((y-mu)^2)/(2*var)
  
  return(exp(-temp2)/temp1)
}

# Dominant survival 
f.dom.surv=function(x,month,density,rainSD,tempSD){
  
  new.data=expand.grid(mass=x,month=month,rainSD=rainSD,
                       tempSD=tempSD,density=density)
  new.data$pred=predict(dom.surv,newdata = new.data,type="response")
  return(new.data)
}

# Dominant growth  
f.dom.gr=function(x,y,pregCat,month,rainSD,tempSD,density,year){
  
  new.data=expand.grid(mass=x,pregCat=pregCat,month=month,rainSD=rainSD,
                       tempSD=tempSD,density=density,year=year)
  mu=predict(dom.gr,newdata = new.data)
  var=2*exp(as.numeric(predict(dom.gr.var,newdata = new.data)))
  temp1 <- sqrt(2*pi*var)
  temp2 <- ((y-mu)^2)/(2*var)
  return(exp(-temp2)/temp1)
  
}

# Non-pregnant dominant transition to 1-month pregnant
f.dom.NP=function(x,month,density,rainSD,tempSD,year){
  
  new.data=expand.grid(mass=x,month=month,rainSD=rainSD,
                       tempSD=tempSD,density=density,year=year)
  new.data$pred=predict(dom.transNP,newdata = new.data,type="response")
  return(new.data)
}

# First-month pregnant dominant transition to 2-month pregnant or aborts
f.dom.F=function(x,month,density,rainSD,tempSD){
  
  new.data=expand.grid(mass=x,month=month,rainSD=rainSD,
                       tempSD=tempSD,density=density)
  new.data$pred=predict(dom.transF,newdata = new.data,type="response")
  return(new.data)
}

# If aborted at first month of pregnancy, becoming 1-month pregnant immediately or remaining non-pregnant
f.dom.F2=function(x,month,density,rainSD,tempSD,year){
  
  new.data=expand.grid(mass=x,month=month,rainSD=rainSD,tempSD=tempSD,density=density,year=year)
  new.data$pred=predict(dom.transF2,newdata = new.data,type="response")
  return(new.data)
}

# Second-month pregnant dominant giving birth to pups or aborts
f.dom.Sec=function(x,month,density,rainSD,tempSD,year){
  
  new.data=expand.grid(mass=x,month=month,rainSD=rainSD,
                       tempSD=tempSD,density=density,year=year)
  new.data$pred=predict(dom.transSec,newdata = new.data,type="response")
  return(new.data)
}

# If aborted at second month of pregnancy, becoming 1-month pregnant immediately or remaining non-pregnant
f.dom.Sec2=function(x,tempSD){
  
  new.data=expand.grid(mass=x,tempSD=tempSD)
  new.data$pred=predict(dom.transSec2,newdata = new.data,type="response")
  return(new.data)
}

# Dominant with litter remaining non-pregnant or becoming pregnant 
f.dom.B=function(x,month,density,rainSD,tempSD,year){
  
  new.data=expand.grid(mass=x,month=month,rainSD=rainSD,
                       tempSD=tempSD,density=density,year=year)
  new.data$pred=predict(dom.transB,newdata = new.data,type="response")
  return(new.data)
}

# If becoming prgnant when with litter, 1-month pregnant vs 2-month pregnant 
f.dom.B2=function(x,month,density,rainSD,tempSD,year=year){
  
  new.data=expand.grid(mass=x,month=month,rainSD=rainSD,
                       tempSD=tempSD,density=density,year=year)
  new.data$pred=predict(dom.transB2,newdata = new.data,type="response")
  return(new.data)
}

# Number of pups recruited per dominant female
f.dom.pups=function(x,density){
  
  new.data=expand.grid(mass=x,density=density)
  
  new.data$pred=predict(dom.pups,newdata = new.data,type="response")
  return(new.data)
}

# Mass of the pups 
f.dom.off.mass=function(x,y,month,rainSD,tempSD,density,year){
  
  new.data=expand.grid(massM=x,month=month,rainSD=rainSD,
                       tempSD=tempSD,density=density,year=year)
  mu=predict(dom.off.mass,newdata = new.data)
  
  var=2*exp(as.numeric(predict(dom.off.mass.var,newdata = new.data)))
  temp1 <- sqrt(2*pi*var)
  temp2 <- ((y-mu)^2)/(2*var)
  # vr=residuals(dom.off.mass)^2
  # sd=sqrt(mean(vr))
  # temp1 <- sqrt(2*pi)*sd
  # temp2 <- ((y-mu)^2)/(2*sd^2)
  return(exp(-temp2)/temp1)
  
}

###############################
#### STEP 2 FILL IN IPM KERNELS

###################################################################

# Define IPM parameters (for mid-point integration)

minMass=3.89 # minimum observed mass
maxMass=6.88 # maximum observed mass

n.bins = 100; n.stage = 20; # number of bins and total number of life-cycle stages considered 
b <- minMass+c(0:n.bins)*(maxMass-minMass)/n.bins 
z <- 0.5*(b[1:n.bins]+b[2:(n.bins+1)]) # bin midpoint 
h <- (maxMass - minMass)/n.bins # bin width

# Intial female population vector (Here Dec 1996)
load("/Users/maria/Dropbox/Meerkats/SuppMat/N0.rda")

dens0=sum(N0)/23.8034 # initial density - total number of females/population range (Dec 1996)

# Current observed (1997-2016) temperature and rainfall variation

load("/Users/maria/Dropbox/Meerkats/SuppMat/clim.obs.rda")

### Start simulations
av_month=1:12
av_year=as.character(seq(1997,2016))
ageM.pup=c("1","2","3")
ageM.juv=c("4","5","6")
ageM.sub=c("7","8","9","10","11","12")

# Empty arrays to hold simulation results in 
mass.dist.obs=array(0,c(length(N0),length(av_year),length(av_month))) # distribution of x-mass individuals per stage
PR.obs=array(0,c(length(av_year),length(av_month))) # simulated population range
densSim.obs=array(0,c(length(av_year),length(av_month))) # simulated density
dominants.obs=array(0,c(length(av_year),length(av_month))) # simulated number of dominant females

# Begin building IPMs 
dens1=dens0
N1=N0

for(j in 1:length(av_year)){ # loop over years - 1997-2016 in sequence 
  for(i in 1:length(av_month)){ # loop over months - 1-12 in sequence
    IPM=array(0,c(n.bins*n.stage,n.bins*n.stage)) 
   
    tempSD=clim.obs$tempSD[clim.obs$year==av_year[j]&clim.obs$month==av_month[i]] # temperature standardized deviation for each month and year
    rainSD=clim.obs$rainSD[clim.obs$year==av_year[j]&clim.obs$month==av_month[i]] # rainfall standardized deviation for each month and year
    density=dens1
    
    ### Age 1 Pups
    
    G=h*t(outer(z,z,f.pup.gr,ageM.pup[1],av_month[i],rainSD,tempSD,density,av_year[j]))
    
    # control for eviction
    G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
    G[is.na(G)]=0
    S=diag(f.pup.surv(z,ageM.pup[1],av_month[i],rainSD,tempSD,density)$pred)
    
    IPM[(n.bins+1):(2*n.bins),1:n.bins]=G%*%S
    
    ### Age 2 Pups
    
    G=h*t(outer(z,z,f.pup.gr,ageM.pup[2],av_month[i],rainSD,tempSD,density,av_year[j]))
    # control for eviction
    G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
    G[is.na(G)]=0
    S=diag(f.pup.surv(z,ageM.pup[2],av_month[i],rainSD,tempSD,density)$pred)
    
    IPM[(2*n.bins+1):(3*n.bins),(n.bins+1):(2*n.bins)]=G%*%S
    
    ### Age 3 Pups
    
    G=h*t(outer(z,z,f.pup.gr,ageM.pup[3],av_month[i],rainSD,tempSD,density,av_year[j]))
    # control for eviction
    G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
    G[is.na(G)]=0
    S=diag(f.pup.surv(z,ageM.pup[3],av_month[i],rainSD,tempSD,density)$pred)
    
    IPM[(3*n.bins+1):(4*n.bins),(2*n.bins+1):(3*n.bins)]=G%*%S
    
    ### Age 4 Juvenile
    
    G=h*t(outer(z,z,f.juv.gr,ageM.juv[1],av_month[i],rainSD,tempSD,density,av_year[j]))
    # control for eviction
    G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
    G[is.na(G)]=0
    S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
    
    IPM[(4*n.bins+1):(5*n.bins),(3*n.bins+1):(4*n.bins)]=G%*%S
    
    ### Age 5 Juvenile
   
    G=h*t(outer(z,z,f.juv.gr,ageM.juv[2],av_month[i],rainSD,tempSD,density,av_year[j]))
    # control for eviction
    G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
    G[is.na(G)]=0
    S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
    
    IPM[(5*n.bins+1):(6*n.bins),(4*n.bins+1):(5*n.bins)]=G%*%S
    
    ### Age 6 Juvenile
    
    G=h*t(outer(z,z,f.juv.gr,ageM.juv[3],av_month[i],rainSD,tempSD,density,av_year[j]))
    # control for eviction
    G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
    G[is.na(G)]=0
    S=diag(f.juv.surv(z,density,rainSD,tempSD)$pred)
    
    IPM[(6*n.bins+1):(7*n.bins),(5*n.bins+1):(6*n.bins)]=G%*%S
    
    ### Age 7 Subadult
   
    G=h*t(outer(z,z,f.sub.gr,ageM.sub[1],av_month[i],rainSD,tempSD,density,av_year[j]))
    # control for eviction
    G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
    G[is.na(G)]=0
    S=diag(f.sub.surv(z,ageM.sub[1],av_month[i],density,rainSD,tempSD)$pred)
    
    IPM[(7*n.bins+1):(8*n.bins),(6*n.bins+1):(7*n.bins)]=G%*%S
    
    ### Age 8 Subadult
    
    G=h*t(outer(z,z,f.sub.gr,ageM.sub[2],av_month[i],rainSD,tempSD,density,av_year[j]))
    # control for eviction
    G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
    G[is.na(G)]=0
    S=diag(f.sub.surv(z,ageM.sub[2],av_month[i],density,rainSD,tempSD)$pred)
    
    IPM[(8*n.bins+1):(9*n.bins),(7*n.bins+1):(8*n.bins)]=G%*%S
    
    ### Age 9 Subadult
    
    G=h*t(outer(z,z,f.sub.gr,ageM.sub[3],av_month[i],rainSD,tempSD,density,av_year[j]))
    # control for eviction
    G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
    G[is.na(G)]=0
    S=diag(f.sub.surv(z,ageM.sub[3],av_month[i],density,rainSD,tempSD)$pred)
    
    IPM[(9*n.bins+1):(10*n.bins),(8*n.bins+1):(9*n.bins)]=G%*%S
    
    ### Age 10 Subadult
    
    G=h*t(outer(z,z,f.sub.gr,ageM.sub[4],av_month[i],rainSD,tempSD,density,av_year[j]))
    # control for eviction
    G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
    G[is.na(G)]=0
    S=diag(f.sub.surv(z,ageM.sub[4],av_month[i],density,rainSD,tempSD)$pred)
    
    IPM[(10*n.bins+1):(11*n.bins),(9*n.bins+1):(10*n.bins)]=G%*%S
    
    ### Age 11 Subadult
    
    G=h*t(outer(z,z,f.sub.gr,ageM.sub[5],av_month[i],rainSD,tempSD,density,av_year[j]))
    # control for eviction
    G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
    G[is.na(G)]=0
    S=diag(f.sub.surv(z,ageM.sub[5],av_month[i],density,rainSD,tempSD)$pred)
    
    IPM[(11*n.bins+1):(12*n.bins),(10*n.bins+1):(11*n.bins)]=G%*%S
    
    ### Age 12 Subadult
    
    
    G=h*t(outer(z,z,f.sub.gr,ageM.sub[6],av_month[i],rainSD,tempSD,density,av_year[j]))
    # control for eviction
    G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
    G[is.na(G)]=0
    S=diag(f.sub.surv(z,ageM.sub[6],av_month[i],density,rainSD,tempSD)$pred)
    
    IPM[(12*n.bins+1):(13*n.bins),(11*n.bins+1):(12*n.bins)]=G%*%S
    
    ### HELPER 
    
    S=diag(f.help.surv(z,av_month[i],density,rainSD,tempSD,av_year[j])$pred)
    
    E=diag(f.help.emig(z,av_month[i],density,rainSD,tempSD,av_year[j])$pred)
    
    Ph_d=diag(f.toDom(z,av_month[i],density,rainSD)$pred)
    
    ### not pregnant -> staying not pregnant helper
    pregCat="np"
    G=h*t(outer(z,z,f.help.gr,pregCat,av_month[i],rainSD,tempSD,density,av_year[j]))
    G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
    G[is.na(G)]=0
    
    Thh.np.p1=diag(f.help.NPH(z,av_month[i],density,rainSD,tempSD,av_year[j])$pred)
    
    IPM[(12*n.bins+1):(13*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*(diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.np.p1))
    
    ### not pregnant -> pregnant helper 1 month 
    
    IPM[(13*n.bins+1):(14*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*(diag(1,n.bins)-Ph_d)*Thh.np.p1)
    
    ### not pregnant -> not pregnant dominant
    Thd.p1=diag(f.help.NPD(z,av_month[i])$pred)
    
    IPM[(16*n.bins+1):(17*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*Ph_d*(diag(1,n.bins)-Thd.p1))
    
    ### not pregnant -> pregnant dominant 1 month
    
    IPM[(17*n.bins+1):(18*n.bins),(12*n.bins+1):(13*n.bins)]=G%*%(S*(diag(1,n.bins)-E)*Ph_d*Thd.p1)
    
    ### pregnant helper 1 month -> pregnant helper 2 month 
    
    pregCat="first"
    G=h*t(outer(z,z,f.help.gr,pregCat,av_month[i],rainSD,tempSD,density,av_year[j]))
    G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
    G[is.na(G)]=0
    
    Thh.p1.p2=diag(f.help.FH(z,av_month[i],density,rainSD,tempSD,av_year[j])$pred)
    
    IPM[(14*n.bins+1):(15*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*Thh.p1.p2)
    
    ### pregnant helper 1 month -> pregnant helper 1 month
    Thh.p1.p1=diag(f.help.FH2(z,av_month[i],density,rainSD,tempSD)$pred)
    
    IPM[(13*n.bins+1):(14*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p1.p2)*Thh.p1.p1)
    
    ### pregnant helper 1 month -> non-pregnant helper 
    IPM[(12*n.bins+1):(13*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p1.p2)*(diag(1,n.bins)-Thh.p1.p1))
    
    ### pregnant helper 1 month -> pregnant dominant 2 months 
    IPM[(18*n.bins+1):(19*n.bins),(13*n.bins+1):(14*n.bins)]=G%*%Ph_d
    
    ### pregnant helper 2 month -> birth
    pregCat="second"
    G=h*t(outer(z,z,f.help.gr,pregCat,av_month[i],rainSD,tempSD,density,av_year[j]))
    G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
    G[is.na(G)]=0
    
    Thh.p2.b=diag(f.help.SH(z,av_month[i],density,av_year[j])$pred)
    
    IPM[(15*n.bins+1):(16*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*Thh.p2.b)
    
    ### pregnant helper 2 month -> abortion and back to pregnant helper 1 month
    Thh.p2.p1=diag(f.help.SH2(z,av_month[i],density,rainSD,tempSD,av_year[j])$pred)
    
    IPM[(13*n.bins+1):(14*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p2.b)*Thh.p2.p1)
    
    ### pregnant helper 2 month -> abortion and back to non-pregnant helper
    
    IPM[(12*n.bins+1):(13*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*(diag(1,n.bins)-Thh.p2.b)*(diag(1,n.bins)-Thh.p2.p1))
    
    ### pregnant helper 2 month -> birth as dominant
    
    IPM[(19*n.bins+1):(20*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%(Ph_d*diag(5/16,n.bins))
    
    ### pregnant helper 2 month -> non-pregnant dominant
    
    IPM[(16*n.bins+1):(17*n.bins),(14*n.bins+1):(15*n.bins)]=G%*%(Ph_d*diag(11/16,n.bins))
    
    ### helper with litter -> non-pregnant helper
    pregCat="birth"
    G=h*t(outer(z,z,f.help.gr,pregCat,av_month[i],rainSD,tempSD,density,av_year[j]))
    G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
    G[is.na(G)]=0
    
    Thh.b.p1=as.numeric(f.help.BH(z,av_month[i])$pred)
    
    IPM[(12*n.bins+1):(13*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*diag(1-Thh.b.p1,n.bins))
    
    ### helper with litter -> pregnant helper
    
    IPM[(13*n.bins+1):(14*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%((diag(1,n.bins)-Ph_d)*diag(Thh.b.p1,n.bins))
    
    ### helper with litter -> non-pregnant dominant
    
    IPM[(16*n.bins+1):(17*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%(Ph_d*0.5)
    
    ### helper with litter -> pregnant dominant
    
    IPM[(17*n.bins+1):(18*n.bins),(15*n.bins+1):(16*n.bins)]=G%*%(Ph_d*0.5)
    
    
    ### DOMINANT
    S=diag(f.dom.surv(z,av_month[i],density,rainSD,tempSD)$pred)
    
    ### non-pregnant -> non-pregnant
    pregCat="np"
    G=h*t(outer(z,z,f.dom.gr,pregCat,av_month[i],rainSD,tempSD,density,av_year[j]))
    G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
    G[is.na(G)]=0
    
    Td.np.p1=diag(f.dom.NP(z,av_month[i],density,rainSD,tempSD,av_year[j])$pred)
    
    IPM[(16*n.bins+1):(17*n.bins),(16*n.bins+1):(17*n.bins)]=G%*%(S*(diag(1,n.bins)-Td.np.p1))
    
    ### non-pregnant -> pregnant
    
    IPM[(17*n.bins+1):(18*n.bins),(16*n.bins+1):(17*n.bins)]=G%*%(S*Td.np.p1)
    
    ### pregnant month 1  -> pregnant month 2
    pregCat="first"
    G=h*t(outer(z,z,f.dom.gr,pregCat,av_month[i],rainSD,tempSD,density,av_year[j]))
    G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
    G[is.na(G)]=0
    
    Td.p1.p2=diag(f.dom.F(z,av_month[i],density,rainSD,tempSD)$pred)
    
    IPM[(18*n.bins+1):(19*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%Td.p1.p2
    
    ### pregnant month 1  -> pregnant month 1
    if(av_year[j]=="1997"|av_year[j]=="1998"){
      
      Td.p1.p1=diag(f.dom.F2(z,av_month[i],density,rainSD,tempSD,"1999")$pred)
    }else{ Td.p1.p1=diag(f.dom.F2(z,av_month[i],density,rainSD,tempSD,av_year[j])$pred) }
    
    
    IPM[(17*n.bins+1):(18*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%((diag(1,n.bins)-Td.p1.p2)*Td.p1.p1)
    
    ### pregnant month 1  -> non-pregnant
    
    IPM[(16*n.bins+1):(17*n.bins),(17*n.bins+1):(18*n.bins)]=G%*%((diag(1,n.bins)-Td.p1.p2)*(diag(1,n.bins)-Td.p1.p1))
    
    ### pregnant month 2 -> birth
    pregCat="second"
    G=h*t(outer(z,z,f.dom.gr,pregCat,av_month[i],rainSD,tempSD,density,av_year[j]))
    G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
    G[is.na(G)]=0
    
    Td.p2.b=diag(f.dom.Sec(z,av_month[i],density,rainSD,tempSD,av_year[j])$pred)
    
    IPM[(19*n.bins+1):(20*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%Td.p2.b
    
    ### pregnant month 2 -> pregnant 1 month
    Td.p2.p1=diag(f.dom.Sec2(z,tempSD)$pred)
    
    IPM[(17*n.bins+1):(18*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%((diag(1,n.bins)-Td.p2.b)*Td.p2.p1)
    
    ### pregnant month 2 -> non-pregnant
    
    IPM[(16*n.bins+1):(17*n.bins),(18*n.bins+1):(19*n.bins)]=G%*%((diag(1,n.bins)-Td.p2.b)*(diag(1,n.bins)-Td.p2.p1))
    
    ### with litter -> non-pregnant
    
    pregCat="birth"
    G=h*t(outer(z,z,f.dom.gr,pregCat,av_month[i],rainSD,tempSD,density,av_year[j]))
    G=G/matrix(as.vector(apply(G,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
    G[is.na(G)]=0
    
    Td.b.np=diag(f.dom.B(z,av_month[i],density,rainSD,tempSD,av_year[j])$pred)
    
    IPM[(16*n.bins+1):(17*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%Td.b.np
    
    ### with litter -> pregnant 1 month
    Td.b.p1=diag(f.dom.B2(z,av_month[i],density,rainSD,tempSD,av_year[j])$pred)
    
    IPM[(17*n.bins+1):(18*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%((diag(1,n.bins)-Td.b.np)*Td.b.p1)
    
    ### with litter -> pregnant 2 month
    
    IPM[(18*n.bins+1):(19*n.bins),(19*n.bins+1):(20*n.bins)]=G%*%((diag(1,n.bins)-Td.b.np)*(diag(1,n.bins)-Td.b.p1))
    
    ###### RECRUITMENT
    
    # from helper
    R=diag(as.numeric(f.help.pups(z,av_month[i],density,rainSD)$pred),n.bins)
    D=h*t(outer(z,z,f.help.off.mass,av_month[i],rainSD,tempSD,density,av_year[j]))
    D=D/matrix(as.vector(apply(D,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
    D[is.na(D)]=0
    
    IPM[1:n.bins,(15*n.bins+1):(16*n.bins)]=D%*%R
    
    # from dominant
    R=diag(f.dom.pups(z,density)$pred)
    D=h*t(outer(z,z,f.dom.off.mass,av_month[i],rainSD,tempSD,density,av_year[j]))
    D=D/matrix(as.vector(apply(D,2,sum)),nrow=n.bins,ncol=n.bins,byrow=TRUE)
    D[is.na(D)]=0
    
    IPM[1:n.bins,(19*n.bins+1):(20*n.bins)]=D%*%R
    
    N1 <- IPM%*%as.numeric(N1) # multiply IPM by vector 
    
    t.sub=as.Date(paste(as.numeric(av_year[j]),av_month[i],15,sep="-"))%m+% months(1)
    if(t.sub=="2017-01-15") t.sub<-as.Date("2016-12-15")
    doms=sum(N1[(length(N1)-n.bins*4+1):length(N1)])
    poprange=predict(m.pop.range$gam, newdata = data.frame(popSize=sum(N1),domNum=doms,
                                                   year=as.character(year(t.sub)),
                                                   month=month(t.sub)), se.fit = F)
    
    dens1=2*(sum(N1)/poprange) # multiply by two to add males 
    
    mass.dist.obs[,j,i]= N1
    PR.obs[j,i]=poprange
    densSim.obs[j,i]=dens1
    dominants.obs[j,i]=doms
  }
  
}

dyn.obs=list(mass.dist.obs=mass.dist.obs,PR.obs=PR.obs,densSim.obs=densSim.obs,dominants.obs=dominants.obs) 


######### PLOTS

### NOTE: one can either use the results after running the code above (but it will take ca. 30 minutes),
## or just load simulations that we included in the supplmentary material

load("/Users/maria/Dropbox/Meerkats/SuppMat/dyn.obs")
load("/Users/maria/Dropbox/Meerkats/SuppMat/dens.obs.rda")
load("/Users/maria/Dropbox/Meerkats/SuppMat/dom.obs.rda")
load("/Users/maria/Dropbox/Meerkats/SuppMat/pop.range.obs.rda")

### PLOT DENSITY
dens.sim=dyn.obs$densSim.obs

dens.sim=adply(dens.sim,c(1,2))
colnames(dens.sim)=c("year","month","density")
levels(dens.sim$year)=1997:2016
dens.sim$date=as.Date(paste(dens.sim$year,dens.sim$month,"15",sep="-"))
dens.sim=dens.sim[order(dens.sim$date),]
dens.sim$cat="Projected"
dens.obs$cat="Observed"

cor.test(dens.sim$density,dens.obs$density)# significant correlation 

all=rbind(dens.obs,dens.sim[,c(4,3,5)])

ggplot(data=all,aes(date,density,col=cat,group=cat))+
  geom_line(size=1.5)+
  scale_color_manual(name="",values=c("black","darkgreen"))+
  theme_bw()+
  theme(legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=18),
        legend.key.size = unit(2, "lines"),
        legend.background =element_rect(fill = "transparent"),
        legend.position = c(0.1,.9))+
  guides(color=guide_legend(ncol=1))+
  
  xlab("Time") +
  ylab("Density")+
  theme(panel.grid = element_blank())+
  theme(axis.text.y = element_text(size=18,colour="black"))+
  theme(axis.title = element_text(size=25))+
  theme(axis.title.x = element_text(vjust=-0.4),
        axis.title.y = element_text(vjust=0.9),
        plot.title = element_text(size=22,lineheight=.8, face="bold"))+
  
  theme(plot.margin = unit(c(0.2,0.1,0.5,0.1), "cm"))+
  scale_x_date(date_breaks = "12 month", labels = date_format("%Y"), expand = c(0, 0))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=16,color="black"))


### PLOT PROPORTION OF DOMINANTS

tot.females=apply(dyn.obs$mass.dist.obs,c(2,3),sum)
dom.sim=dyn.obs$dominants.obs/tot.females


dom.sim=adply(dom.sim,c(1,2))
colnames(dom.sim)=c("year","month","dom")
levels(dom.sim$year)=1997:2016
dom.sim$date=as.Date(paste(dom.sim$year,dom.sim$month,"15",sep="-"))
dom.sim=dom.sim[order(dom.sim$date),]
dom.sim$cat="Projected"
dom.obs$cat="Observed"

cor.test(dom.sim$dom,dom.obs$dom)# significant correlation 

all=rbind(dom.obs,dom.sim[,c(4,3,5)])

ggplot(data=all,aes(date,dom,col=cat,group=cat))+
  geom_line(size=1.5)+
  scale_color_manual(name="",values=c("black","darkgreen"))+
  theme_bw()+
  theme(legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=18),
        legend.key.size = unit(2, "lines"),
        legend.background =element_rect(fill = "transparent"),
        legend.position = c(0.1,.9))+
  guides(color=guide_legend(ncol=1))+
  
  xlab("Time") +
  ylab("Proportion ominants")+
  theme(panel.grid = element_blank())+
  theme(axis.text.y = element_text(size=18,colour="black"))+
  theme(axis.title = element_text(size=25))+
  theme(axis.title.x = element_text(vjust=-0.4),
        axis.title.y = element_text(vjust=0.9),
        plot.title = element_text(size=22,lineheight=.8, face="bold"))+
  
  theme(plot.margin = unit(c(0.2,0.1,0.5,0.1), "cm"))+
  scale_x_date(date_breaks = "12 month", labels = date_format("%Y"), expand = c(0, 0))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size=16,color="black"))




