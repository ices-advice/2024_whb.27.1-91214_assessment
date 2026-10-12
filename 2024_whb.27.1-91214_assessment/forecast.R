## Run analysis, write model results

## Before:
## After:

library(icesTAF)
library(stockassessment)

mkdir("forecast")

# Load model run for forecast
load(file="model/model_forecast.RData")

# TAC year and years to average across
assessmentYear      <- fit.forecast$data$years[(length(fit.forecast$data$years)-1)]
TACyear             <- assessmentYear + 1
#ave.years vector of years to average for weights, maturity, M and such
ave.years           <- c((assessmentYear-2):assessmentYear)   # for mean weights etc. (copy of preliminary data year is not used, used since 2019)

# Set years for recruitment to average across - depending on taking full time series or high recruitment regime only
# rec.years           <-1981:2020  # years to calculate GM recruitment (full time series minus 1)
rec.years           <-1996:2023  # years to calculate GM recruitment (high recruitment regime) 

# Set reference points
Fpa                 <-0.32
FMSY                <-0.32
Flim                <-0.88
Blim                <-1500000
Bpa                 <-2250000
MSYBtrigger         <-2250000

# Set latest catch advice and latest catch
last.advice         <-1529754  # Published 29 September 2023
last.catch          <-tail(catchtable(fit.forecast,obs.show=TRUE),1)[1,4] 

# Extract latest F  
Fsq                 <-fbartable(fit.forecast)[as.character(TACyear-1),1]  #F status quo = average F from last year - used for one of the catch options

#if (!extended) stop('Forecast must be done using the "extended dataset" to get the right (process corrected) initial stock numbers')

# Get recruitment table and extract latest recruitment value
rec                 <-rectable(fit.forecast)  # recruitment
rec
rec.last            <-tail(rec,1)[1,1]

# Estimate the geometric mean recruitment based on the selected years
GM                  <-exp(mean(log(rectable(fit.forecast)[as.character(rec.years),1])))  # all years in rec.year
GM

# Set reference year for recruitment
recRefYear          <-TACyear-1
refRecruit          <-rectable(fit.forecast)[as.character(recRefYear),1]

# use SAM estimate for age 1 in terminal year (year of the preliminary data)
my.nscale           <- c(GM/refRecruit,1,GM/refRecruit)  # scaling of N: age 1 first year, age 2 first year and age 1 second+ year

# change the SAM estimate for age 1 in the terminal year (year of the preliminary data), in this case use GM for all years
#my.nscale<- c(GM/refRecruit,GM/refRecruit,GM/refRecruit)  # scaling of N: age 1 first year, age 2 first year and age 1 second+ year

my.nscale


source("forecast_deterministic.R")

#just testing
tst                 <-forecast_deterministic(fit.forecast, 
                                             ave.years =ave.years,
                                             fscale    =c(NA,NA,NA),
                                             fval      =c(0.32,0.32,0.32), 
                                             nscale    =my.nscale,
                                             label     =paste0("test"),
                                             rec.years =TACyear-1)

tab                 <-attr(tst, "tab")
tab
ssb.start.TAC.year  <-tab[1,'ssb:median']

# Set up optimizing function 
SAM_optimize<-function(x,speci) {
  FF        <-x
  my.fval   <-c(FF,FF)
  set.seed(12345)
  a         <-forecast_deterministic(fit.forecast, 
                                     ave.years  =ave.years, 
                                     fscale     =rep(NA,2),
                                     fval       =my.fval, 
                                     nscale     =my.nscale,
                                     label      =paste0("Status quo F in first year followed by F",FF),
                                     rec.years  =TACyear-1)
  b         <-attr(a, "tab")
  out       <-b[speci$y,speci$target]
  target    <-(speci$value-out)^2
  return(target)
}

# Set-up scenario table
scen<- list(
  list(lab="Long-term management strategy (F=FMSY)", target='F',value=FMSY,y=0),
  list(lab="MSY approach: FMSY",                 target='F',value=FMSY,y=0),
  list(lab="F = 0",                              target='F',value=1E-6,y=0),
  list(lab="Fpa",                                target='F',value=Fpa,y=0),
  list(lab="Flim",                               target='F',value=Flim,y=0),
  list(lab=paste0("SSB (",TACyear+1,") = Blim"), target='ssb:median',value=Blim,y=2),
  list(lab=paste0("SSB (",TACyear+1," = Bpa"),   target='ssb:median',value=Bpa,y=2),
  list(lab=paste0("SSB (",TACyear+1,") = MSY Btrigger"), target='ssb:median',value=MSYBtrigger,y=2),
  list(lab=paste0("F = F (",TACyear-1,")"), target='F',value=Fsq,y=0),
  list(lab=paste0("SSB (",TACyear+1,") = SSB (",TACyear,")"), target='ssb:median',value=ssb.start.TAC.year,y=2),
  list(lab=paste0("Catch (",TACyear,") = Catch (",TACyear-1,")"), target='catch:median',value=last.catch,y=1),
  list(lab=paste0("Catch (",TACyear,") = Catch (",TACyear-1,") -20 %"),target='catch:median',value=last.catch*0.8,y=1),
  list(lab=paste0("Catch (",TACyear,") = Catch (",TACyear-1,")  +25%"),target='catch:median',value=last.catch*1.25,y=1),
  list(lab=paste0("Catch (",TACyear,") = Advice (",TACyear-1,") -20 %"),target='catch:median',value=last.advice*0.8,y=1),
  list(lab=paste0("Catch (",TACyear,") = Advice (",TACyear-1,") +25%"),target='catch:median',value=last.advice*1.25,y=1)
)

# Get catch options for all the scenarios from the optimizing function
b <-lapply(scen,function(x){ if (x$target=='F') #EB: what is done here?
  return(list(minimum=x$value,objective=0)) 
  else optimize(SAM_optimize,c(0,7),speci=x)})

# should be small ( less than around 1000)
unlist(lapply(b,function(x) x$objective)) #EB: what are we checking here?

targetF             <-unlist(lapply(b,function(x) x$minimum))
labels              <-unlist(lapply(scen,function(x) x$lab))

names(targetF)      <-labels
t(t(targetF))

# fsq, SSB=Bpa, SSB= Blim, SSB=last year's SSB catch(y+1)=catch(y)
add.targetF         <-c(0.05,0.10,seq(from=0.15,to=0.35,by=.01),0.45,0.50)
names(add.targetF)  <-paste('F =',add.targetF)
targetF             <-c(targetF,add.targetF)
round(t(t(targetF)),3)

my.fscale           <-c(NA,NA)   # F status quo the first year 

FC                  <-list()
set.seed(12345)

for (FF in targetF) { # EB: what is done here?
  my.fval            <-c(FF,FF)
  #print(my.fval)
  set.seed(12345)
  FC[[length(FC)+1]] <- forecast_deterministic(fit.forecast, 
                                               ave.years    =ave.years, 
                                               fscale       =rep(NA,2),
                                               fval         =my.fval, 
                                               nscale       =my.nscale,
                                               label        =paste0("Status quo F in first year followed by F",FF),
                                               rec.years    =TACyear-1)
  #cat(FF,'\n')
}

save(FC,gem,Fpa,FMSY,Flim,Blim,Bpa,MSYBtrigger,last.advice,last.catch,Fsq,rec.last,ssb.start.TAC.year,scen,b,targetF,labels,file="forecast/forecast.RData")
