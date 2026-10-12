## Prepare plots and tables for report

## Before:
## After:

library(icesTAF)
library(stockassessment)
library(tidyverse)

mkdir("report")

# Load model runs and forecast results
load("model/model_base.RData")
load("model/model_forecast.RData")
load("forecast/forecast.RData")

# Load data used for forecast
load("data/data_forecast.RData")

# Load reporting functions
source("common.R")
source("MV_functions.R")

# Set assessment year and TAC year
assessmentYear      <- fit.base$data$years[length(fit.base$data$years)]
TACyear             <- assessmentYear + 1


############################################################################################################-
# Base model run ----

# Set base run as 'main' fit
fit <- fit.base

# Run script that makes first large bunch of plots and tables for the report
source("plotscript.R")

# Table parameter estimates
CommonThings(fit.base)
tp<-tableParameters(fit.base)
print(round(tp,2),na.print = "")
xtab(tp, caption=paste('Table 11. Parameter estimates','.',sep=''), cornername='Parameter \ Year', 
     file=file.path("report",paste(stamp,'_tab11a.html',sep='')), dec=rep(2,ncol(ftab)))

# Table parameter estimates last 5 years based on retrospective analysis
a<-cbind(
  tableParameters(RETRO[[4]],doCV=FALSE),
  tableParameters(RETRO[[3]],doCV=FALSE),
  tableParameters(RETRO[[2]],doCV=FALSE),
  tableParameters(RETRO[[1]],doCV=FALSE),
  tableParameters(fit,doCV=FALSE))

colnames(a)<-as.character((assessmentYear-4):assessmentYear)
xtab(a, caption=paste('Table 11. Parameter estimates, retrospective','.',sep=''), cornername='Parameter \ Year', 
     file=file.path("report",paste(stamp,'_tab11.html',sep='')), dec=rep(2,ncol(ftab)))



############################################################################################################-
# Forecast results and catch advice ----

# Just to check: first catch advice according to long-term management plan
FC[[1]]

# Check N at age
FbarAges<-c(3,7)
test<-FC[[1]]
n1<-exp(test[[1]]$sim) # Stock at the beginning of the TAC year
n2<-exp(test[[2]]$sim) # Stock at the beginning at next year
n<-rbind(n1[1,],n2[1,]) #what are the different columns
N<-n[,1:10]
cat(paste0("Stock N at the beginning of the TAC year(",TACyear,") and the following year(",TACyear+1,")\n:"));round(N)

# Check F at age and exploitation pattern
FF<-cbind(n[,11:19],n[,19])
round(FF,3)
expl<-FF[2,]
expl<-expl/mean(expl[FbarAges[1]:FbarAges[2]])
expl

# Get mean stock weight, proportion mature and natural mortality used in the forecast
mw<-attr(test,"fit")[['data']]$stockMeanWeight
mw[as.character(TACyear),]
pm<-attr(test,"fit")[['data']]$propMat
nm<-attr(test,"fit")[['data']]$natMor

# Table with input to short-term forecast
tab<-t(rbind(mw[as.character(TACyear-1),],gem$sw ,dat.forecast$propMat[as.character(TACyear-1),],dat.forecast$natMor[as.character(TACyear-1),],expl,N[1,]))
rownames(tab)<-paste('Age',attr(test,"fit")[['conf']]$minAge:attr(test,"fit")[['conf']]$maxAge)
colnames(tab)<-c(paste0('Mean weight in the stock and catch (kg) in ',as.character(TACyear-1)),paste0('Mean weight in the stock and catch (kg) in ', as.character(TACyear),'+'),'Proportion mature','Natural mortality','Exploitation pattern',paste0('Stock number(',TACyear,') (thousands)'))

stamp<-gsub('-[[:digit:]]{4}$','',gsub(':','.',gsub(' ','-',gsub('^[[:alpha:]]{3} ','',date()))))

xtab(tab, caption=paste('Table 20. Input to short term forecast','.',sep=''), cornername='Age', 
     file=file.path('report',paste(stamp,'_tab20.html',sep='')), dec=c(3,3,2,2,3,0))

xtab(tab, caption=paste('UNROUNDED Table 21. Input to short term forecast','.',sep=''), cornername='Age', 
     file=file.path('report',paste(stamp,'_tab21.html',sep='')), dec=c(5,5,5,5,5,2))

# Get SSB and catch in intermediate year
ft<-attr(FC[[1]],"tab")
refSSB<-ft[as.character(TACyear),'ssb:median']
refCatch<-catchtable(fit.forecast,obs.show = TRUE)[as.character(TACyear-1),'sop.catch']

# Table with intermediate year assumptions
ft<-rbind(Fsq, 
          refSSB,
          rec.last,
          ft[as.character(TACyear),'rec:median'],
          ft[as.character(TACyear+1),'rec:median'],
          refCatch)

rownames(ft)<-c(paste0('F ages ',FbarAges[1],'-',FbarAges[2],' (',TACyear-1,')'),
                paste0('SSB (',TACyear,')'),
                paste0('R age 1 (',TACyear-1,')'),
                paste0('R age 1 (',TACyear,')'),     
                paste0('R age 1 (',TACyear+1,')'),
                paste0('Total catch (',TACyear-1,')'))  
ft<-cbind(ft,rep(NA,dim(ft)[[1]]))
colnames(ft)<-c('Value','Notes')

stamp<-gsub('-[[:digit:]]{4}$','',gsub(':','.',gsub(' ','-',gsub('^[[:alpha:]]{3} ','',date()))))

xtab(ft, caption=paste('Table 19. Intermediate year input to short term forecast','.',sep=''), cornername='Values', 
     file=file.path('report',paste(stamp,'_tab19.html',sep='')), dec=c(3,0))

# Calculate change in SSB, catch adn advice
SSByear<-TACyear+1
if(exists("FC")){  
  res<-lapply(FC, function(f){
    tf<-attr(f,"tab");
    c(tf[as.character(TACyear),'catch:median'],tf[as.character(TACyear),'fbar:median'],tf[as.character(SSByear),'ssb:median'])
  })
  
  res<-cbind(matrix(unlist(res),ncol=3,byrow=TRUE),0,0,0)
  res[,4]<-(res[,3]-refSSB)/refSSB*100   
  res[,5]<-(res[,1]-refCatch)/refCatch*100 
  res[,6]<-(res[,1]-last.advice)/last.advice*100 
  head(res)
}

# Table with catch options
rownames(res)<-names(targetF)
colnames(res)<-c(paste0('Catch(',TACyear,')'), paste0('F(',TACyear,')'),paste0('SSB(',SSByear,')'),'% SSB change*', '% Catch change**','% Advice change***')
head(res,14)
stamp<-gsub('-[[:digit:]]{4}$','',gsub(':','.',gsub(' ','-',gsub('^[[:alpha:]]{3} ','',date()))))
xtab(res, caption=paste('Table 22.  short term forecast','.',sep=''), cornername='Basis', 
     file=file.path('report',paste(stamp,'_tab22.html',sep='')), dec=c(0,3,0,1,1,1))
#res

# Check the catch proportions by age for the next year for the first scenario (long-term management strategy)
catchAge <-  FC[[1]][[1]][["catchByAge"]][,1]
catchProp <- catchAge / sum(catchAge)
ctprop <- matrix(round(catchAge,0))
ctprop <- cbind(ctprop, round(catchProp*100,1))
colnames(ctprop) <- c("Catch at age (tonnes)","% of catch")
rownames <- c(fit.forecast$data$minAgePerFleet[1]:fit.forecast$data$maxAgePerFleet[1])
ctprop

stamp<-gsub('-[[:digit:]]{4}$','',gsub(':','.',gsub(' ','-',gsub('^[[:alpha:]]{3} ','',date()))))
xtab(ctprop, caption=paste('Table 23.  short term forecast catch proportions by age for TAC year','.',sep=''), cornername='Age', 
     file=file.path('report',paste(stamp,'_tab23.html',sep='')), dec=c(0,0,1))
