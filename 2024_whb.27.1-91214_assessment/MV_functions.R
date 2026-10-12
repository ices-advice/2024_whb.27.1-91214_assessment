# various variables, used later on
CommonThings<-function(fit) {
  minAge<<-fit$conf$minAge
  maxAge<<-fit$conf$maxAge
  nAges<<-maxAge-minAge+1
  ages<<-paste('Age',minAge:maxAge)
  aIndex<<-(minAge:maxAge)-minAge+1
  nFleet<<-sum(fit$data$fleetTypes==2)
  surveys<-fit$data$surveys
  surveyNames<<-names(surveys)
}

logMigEst<-function(fit) {
  scale<-exp(fit$pl$logMig)
  upper<-exp(fit$pl$logMig+2*fit$plsd$logMig)
  lower<-exp(fit$pl$logMig-2*fit$plsd$logMig)
  scale<-matrix(c(scale,lower,upper),nrow=3,byrow=T)
  colnames(scale)<-paste('Age',fit$data$minAgePerFleet[1]:fit$data$maxAgePerFleet[1])
  rownames(scale)<-c('estimate','lower','upper')
  return(scale)
}
 
catchScaling<-function(fit,byYear=TRUE) {
  scale<-exp(fit$pl$logScale)
  upper<-exp(fit$pl$logScale+2*fit$plsd$logScale)
  lower<-exp(fit$pl$logScale-2*fit$plsd$logScale)
  scale<-matrix(c(scale,lower,upper),nrow=3,byrow=T)
  if(byYear) colnames(scale)<-fit$conf$keyScaledYears
  rownames(scale)<-c('estimate','lower','upper')
  return(scale)
}

Ppar<-function(fit) {
  a<-fit$rep$P; colnames(a)<-as.character(fit$data$years);
  rownames(a)<-ages
  return(a)
}

PrintProcessNoise<-function(fits) {
  cat('\nprocess noise N,\n')
  for (i in (1:length(fits))) cat(names(fits)[i],'\t',round(exp(fits[[i]]$pl$logSdLogN),2),'\n')
}
PrintlogLik<-function(fits) { 
 cat('\n,\n log Likelihood\n');for (i in (1:length(fits))) { cat(names(fits)[i],'\t',paste0(round(logLik(fits[[i]]),1),'\t(df=',length(fits[[i]]$sdrep$par.fixed),')\n')) }
}

PrintAIC<-function(fits) {
  cat('\nAIC,\n')
  for (i in (1:length(fits))) cat(names(fits)[i],'\t',round(AIC(fits[[i]]),2),'\n')
}


AICc<-function(fit) {
  n<-nobs(fit)
  k<-dim(partable(fit))[[1]]
  return(AIC(fit)+(2*k*(k+1))/(n-k-1))
}


catchScalingPlot<-function(fit,byYear=TRUE,meanF=FALSE) {
  CommonThings(fit)
  scale<-exp(fit$pl$logScale)
  upper<-exp(fit$pl$logScale+2*fit$plsd$logScale)
  lower<-exp(fit$pl$logScale-2*fit$plsd$logScale)
  scale<-matrix(c(scale,lower,upper),nrow=3,byrow=T)
  if (byYear) colnames(scale)<-fit$conf$keyScaledYears
  rownames(scale)<-c('estimate','lower','upper')
  round(scale,2)
  
  
  # F table
  tabF<-exp(fit$pl$logF)
  #rownames(tabF)<-makeLabels(fit$conf$keyLogFsta[1,],fill=NULL)
  lastF<-dim(tabF)[[1]]
  for (a in ((dim(tabF)[[1]]+1):nAges)) tabF<-rbind(tabF,tabF[lastF,])
  rownames(tabF)<-ages
  colnames(tabF)<-as.character(fit$data$years)
  round(tabF,3)
  
  sy<-as.character(fit$conf$keyScaledYears)
  s<-exp(fit$pl$logScale)
  s
  
  SC<-fit$conf$keyParScaledYA
  for (i in (1:length(s))) {
    SC[(i-1)==SC]<-s[i] 
  }  
  SC[SC<=0]<- 1
  rownames(SC)<-as.character(fit$conf$keyScaledYears)
  SC
  
  scalingFactor<-t(tabF)
  scalingFactor[]<-1
  scalingFactor[sy,]<-SC[sy, ]
  scalingFactor<-t(scalingFactor)
  
  trueF<-tabF/scalingFactor
  unacounted<-tabF-trueF
  
  
  X11(width=8, height=6); par<-par(mfrow = c(2, 2))
  maxF<-max(tabF)
  matplot(y=t(tabF),x=as.numeric(colnames(tabF)),type='b',main='F including contributions from catch scaling',xlab='',ylab='Mortality rate',ylim=c(0,maxF),lwd=2)
  abline(v=min(fit$conf$keyScaledYears),lty=4,col=1);   abline(v=max(fit$conf$keyScaledYears),lty=4,col=1)
  
  matplot(y=t(unacounted),x=as.numeric(colnames(scalingFactor)),type='b',main='Unacounted mortality',xlab='',ylab='Mortality rate',ylim=c(0,maxF),lwd=2)
  abline(v=min(fit$conf$keyScaledYears),lty=4,col=1);   abline(v=max(fit$conf$keyScaledYears),lty=4,col=1)
  
  matplot(y=t(trueF),x=as.numeric(colnames(scalingFactor)),type='b',main='F excluding contributions from catch scaling',xlab='',ylab='Mortality rate',ylim=c(0,maxF),lwd=2)
  abline(v=min(fit$conf$keyScaledYears),lty=4,col=1);   abline(v=max(fit$conf$keyScaledYears),lty=4,col=1)
  
  matplot(y=t(scalingFactor),x=as.numeric(colnames(scalingFactor)),type='b',main='Scaling factor',xlab='',ylab='Scaling factor',lwd=2)
  abline(v=min(fit$conf$keyScaledYears),lty=4,col=1);   abline(v=max(fit$conf$keyScaledYears),lty=4,col=1)
  
  
  FtoM<-t(unacounted)
  return(FtoM)
}



makeLabels<-function(k,pre=NULL,fill=' ---') {
  # test  k<-fit$conf$keyVarObs[5,];k
  kIndex<-aIndex[k>=0]
  kUnique<-!duplicated(k[kIndex])
  kNames<-ages[kUnique]
  ageUnique<-ages[k>=0][kUnique]
  fl<-matrix(' ',ncol=2,nrow=length(ageUnique))
  fl[,1]<- ageUnique
  i<-0;j<-0
  for (a in kUnique) {j<-j+1; if (a) i<-i+1 else  fl[i,2]<-paste0('-',as.character(j))}
  lab<-paste0(fl[,1],fl[,2])
  if (!is.null(pre)) lab[1]<-paste(pre,lab[1])
  if (length(lab)>1 & !is.null(fill)) for ( i in (2:length(lab))) lab[i]<-paste(fill,lab[i])
  return(lab)
}



#observation variance
obsVar<-function(fit,roundV=2) {
  CommonThings(fit)
  #print(fnames)
  pp<-t(fit$conf$keyVarObs)
  colnames(pp)<-c('Catch',surveyNames)
  rownames(pp)<-ages
  pp[pp>=0]<- exp(fit$pl$logSdLogObs)[pp+1]
  pp<-t(pp)
  pp[pp== -1]<-NA
  return(round(pp,roundV))
}



doPlot<-function(fit,byYear=TRUE,plotP=FALSE){
  
  plotfit<-function(Type='SSB',inp) {
    for (i in (1:length(inp))) {
      x<-inp[[i]]
      if (Type=='SSB')  { Y<-x[,'SSB'];       maxY<-maxVal["SSB"] }
      if (Type=='R')    { Y<-x[,'R(age 1)'];  maxY<-maxVal['R(age 1)'] }
      if (Type=='Fbar') { Y<-x[,'Fbar(3-5)']; maxY<-maxVal['Fbar(3-5)'] }
      years<-as.numeric(rownames(x))
      if (i==1) plot(x=years,y=Y,col=cols[i],pch=as.character(i),lwd=2,xlab='Year', ylab=Type,main=Type,type='b',ylim=c(0,maxY)) else points(x=years,y=Y,col=cols[i],lwd=2,pch=as.character(i),type='b')  }
  }
  
  plotscale<-function(inp) {
    mult=1
    for (i in (1:length(inp))) {
      x<-inp[[i]]
      Y<-exp(x$pl$logScale)
      if (length(x$conf$keyScaledYears)==length(Y)) {
        if (i==1) plot(x=x$conf$keyScaledYears,y=Y,ylim=c(0,maxVal),col=cols[i],pch=as.character(i),lwd=2,main='Catch scaling',xlab='Scaling year',ylab='Catch scaling', type='b') else points(x=x$conf$keyScaledYears,y=Y,col=cols[i],lwd=2,pch=as.character(i),type='b') 
      }
    } 
    abline(h=1,lty=3,col=1,lwd=2)
  } 
  
  plotPval<-function(inp) {
    for (i in (1:length(inp))) {
      xx<-inp[[i]]
      if (i==1) plot(x=xx$data$years,y=xx$rep$P[1,],type='b',ylim=c(0,1),xlab='Year',ylab='proportion',main='Proportion Kattegat cod at age 1',pch=as.character(i),lwd=2,col=cols[i]) 
      if (i>1)  points(x=xx$data$years,y=xx$rep$P[1,],type='b',pch=as.character(i),lwd=2,col=cols[i])
       } 
  }  
  
 
  
  fs<-lapply(fit,summary)
  fss<-do.call('rbind',fs)
  maxVal<-apply(fss,2,max)

  cols<-rainbow(length(fs))
  cols<-palette()
 # X11(width=7,height=5)
#  par(mfrow = c(2, 2))
#  par(mar=c(4, 4, 4, 2) + 0.1)  #  c(bottom, left, top, right)
 
  
  X11(width=10,height=8)
  des<-matrix(c(0,1,0,0.10,
                0,1,0.10,1.0 ), nrow=2, ncol=4, byrow=TRUE)
  split.screen(des)     # split display into two screens (screen 1-2)
  split.screen(c(2,2), screen = 2) # now split the bottom half into 2x2 (screen 3-6)
  
  screen(3);  plotfit(Type='SSB',inp=fs)
  screen(4);  plotfit(Type='R',inp=fs)
  screen(5);  plotfit(Type='Fbar',inp=fs)
  if (byYear) {
    maxVal<-max(unlist(lapply(fit,function(x)max(exp(x$pl$logScale)))))
    screen(6);  plotscale(inp=fit)
  } 
  if (plotP) {
    screen(6);  plotPval(inp=fit)
  }
   screen(1); 
   par(mar=c(0, 0, 0, 0))   #  c(bottom, left, top, right)
   if (length(fs)>6) yl<-3 else yl<-2
   plot(1,1, type="n", yaxt="n", xaxt='n',xlab='',ylab="",ylim=c(-0.2,yl),xlim=c(0,3),main=NULL,frame.plot=F)
   for (i in (1:length(fit))) {
     #cat('\nx=',(i-1) %%  3,' y=',(i*2-1) %/% 6,'\n')
     text(x=(i-1) %%  3,y=(i*2-1) %/% 6,labels=names(fit)[i],adj=c(0,0),offset=0,col=cols[i],cex=1.5)
   }
   box()
  close.screen(all = TRUE)    # exit split-screen mode

}

ZwithMigration<-function(fit,adjustZ=TRUE) {
  
  tabN<-exp(fit$pl$logN)
  lastN<-dim(tabN)[[1]]
  
  tabF<-exp(fit$pl$logF)
  lastF<-dim(tabF)[[1]]
  for (a in ((dim(tabF)[[1]]+1):lastN)) tabF<-rbind(tabF,tabF[lastF,])

  colnames(tabF)<-as.character(fit$data$years)
  rownames(tabF)<-paste('Age',fit$conf$minAge:fit$conf$maxAge) 
    
  #round(tabF,3)

  M<-t(fit$data$natMor)
  P<-fit$rep$P;  
  MM<-t(fit$data$migMor)
  mm<-(1-P)*MM
  Z<-tabF+M
  if (adjustZ) {
   if (fit$data$Nmodel==1 | fit$data$Nmodel==3) Z<- -log( exp(-(tabF+M)) *(P+(1-P)*exp(-MM)) )
   if (fit$data$Nmodel==2 |fit$data$Nmodel==4) Z<- -log( exp(-(tabF+M)) *(P+(1-P)*MM )) 
  }
 return(Z)
}


ProcessNoise<-function(fit,doPlot=0,main1='Joint sample residuals log(N)',main2='Joint sample residuals log(F)',newGraph=TRUE) {
  resp <- procres(fit)
  if (doPlot>0) {
      resNF<-data.frame(year=resp$year,age=resp$age,fleet=resp$fleet,res=resp$residual)
      resN<-subset(resNF,fleet==doPlot)
      residual<-tapply(resN$res,list(resN$year,resN$age),sum)
      if (newGraph) X11(h=8,w=8)
      matplot(x=rownames(residual),y=residual,type='b',lwd=2, xlab='Year',main=ifelse(doPlot==1,main1,main2))
      abline(h=0,lty=3)
 
  }
  return(resp)   
}

### process noise
processNoiseZ<-function(fit,main='Negative values show that more fish are removed\n than calculated from M and F',all=FALSE,old=TRUE) {
  a<-fit$rep$resNZ
  rownames(a)<-paste(fit$conf$minAge:fit$conf$maxAge) 
  colnames(a)<-head(fit$data$years,-1)
  a<-t(a)
  if (all) {X11(w=10,h=8);matplot(y=a,x=head(fit$data$years,-1),xlab='year',ylab='process Z',main=main,lwd=2,type='b',pch=colnames(a));abline(h=0,lty=3,lwd=2)}
  if (old) {X11(w=10,h=8);matplot(y=a[,-1],x=head(fit$data$years,-1),xlab='year',ylab='process Z',main=main,lwd=2,type='b',pch=colnames(a[,-1]));abline(h=0,lty=3,lwd=2)}
  return(t(a))  
}
#processNoiseZ(fit)


tableParameters<-function(fit,doCV=TRUE,fnames=names(surveys)) {
  
  tabP<-matrix(
    c(NA,exp(fit$pl$logSdLogFsta),NA,exp(fit$pl$logSdLogN),NA,exp(fit$pl$logSdLogObs),NA,exp(fit$pl$logFpa)*1000,NA,exp(fit$pl$itrans_rho)/(1+exp(fit$pl$itrans_rho)))
               ,ncol=1)
  if (doCV) {
    tabP2<-matrix(      c (NA,fit$plsd$logSdLogFsta, NA,fit$plsd$logSdLogN, NA,fit$plsd$logSdLogObs, NA,fit$plsd$logFpa,NA,fit$plsd$itrans_rho) ,ncol=1)
    tabP<-cbind(tabP,tabP2)
    colnames(tabP)<-c('value','CV')
  } else colnames(tabP)<-c('value')
  
  makeLabs<-function(k,pre=NULL,fill=' ---') {
    lab<-makeLabels(k,fill=fill,pre=pre)
    llab<-length(lab)
    rnames[rowNo:(rowNo+llab-1)]<<-lab
    rowNo<<-rowNo+llab+1
  }
  
  rnames<-rep('--',nrow(tabP))
  rnames[is.na(tabP[,1])]<-c('Random walk variance','Process error','Observation variance','Survey catchability','Rho')
  rowNo<-2
  
  makeLabs(fit$conf$keyVarF[1,],pre='-F')
  makeLabs(fit$conf$keyVarLogN,pre='-log(N)')
  makeLabs(fit$conf$keyVarObs[1,],pre='-Catch')
  nFl<-dim(fit$conf$keyVarObs)[[1]]-1
  for (i in (1:nFl)) { rowNo<-rowNo-1;makeLabs(fit$conf$keyVarObs[i+1,],paste0('-',pre=fnames[i]))}
  rowNo<-rowNo+1
  for (i in (1:nFl)) { rowNo<-rowNo-1;makeLabs(fit$conf$keyLogFpar[i+1,],paste0('-',pre=fnames[i]))}
  rownames(tabP)<-rnames
  return(tabP)
}



# Catchability
Catchability<-function(fit,addUncertanty=FALSE) {
  pp<-fit$conf$keyLogFpar
  pp<-t(pp[-1,])
  colnames(pp)<-names(surveys)
  rownames(pp)<-ages
  pp.cv<-pp
  pp[pp>=0]<-(exp(fit$pl$logFpar)*1000)[pp+1]
  pp<-t(pp); 
  pp[pp== -1]<-NA
  round(pp,1)
  

  pp.cv[pp.cv>=0]<-fit$plsd$logFpar[pp.cv+1]
  pp.cv<-t(pp.cv)
  pp.cv[pp.cv== -1]<-NA
  round(pp.cv,2)
 
  cat('Catchability\n')
  print(round(pp,2))
  if (addUncertanty) {
    cat('Catchability CV\n')
    print(round(pp.cv,2))
    cat('\nCatchability lower\n')
    print(round(pp-(pp.cv*pp*2),2))
    cat('\nCatchability upper\n')
    print(round(pp+(pp.cv*pp*2),2))
  }
  matplot(t(pp),type='b')
  pp<-pp/rowMeans(pp,na.rm=TRUE)
  matplot(t(pp),type='b', ylab='Standardized catchability',xlab='Age',lwd=2)
  legend("topright",legend=rownames(pp),pch=as.character(1:dim(pp)[[1]]))
  ### 
}


FleetAssess<-function(incl) {
  incl<-unlist(incl)
  incl<-incl[!is.na(incl)]
  excl<-setdiff(1:(nFleet+1),incl)
  a<-runwithout(fit,fleet=excl,silent = TRUE)
  return(list(incl=incl,excl=excl,fit=a))
}


doPjuvenile<-function(fit,trans=TRUE) {
  Pjuv<-fit$pl$Pjuv; 
  Pjuvsd<-fit$plsd$Pjuv;
  PP<-matrix(0,ncol=length(Pjuv),nrow=3)
  colnames(PP)<-as.character(head(fit$data$years,-1));rownames(PP)<-c('estimate','lower','upper')
  PP[1,]<-Pjuv
  PP[2,]<-Pjuv-1.95*Pjuvsd
  PP[3,]<-Pjuv+1.95*Pjuvsd
  if (trans) PP<-1.0/(1.0+exp(-PP))
  print(round(PP,2))
}

plotPjuv<-function(fits) {
  a<-lapply(fits,function(x)return(1/(1+exp(-x$pl$Pjuv))))
  runs<-names(fits)
  b<-t(matrix(unlist(a),byrow=TRUE,nrow=length(a)))
  
  X11(width=9,height=7)
  cols<-palette()
  des<-matrix(c(0,1,0,0.10,
                0,1,0.10,1.0 ), nrow=2, ncol=4, byrow=TRUE)
  split.screen(des)     # split display into two screens (screen 1-2)
  screen(2)
  matplot(b,x=head(fits[[1]]$data$years,-1),lwd=2,type='b',xlab='Year',ylab='Proportion Kattegat cod')
  screen(1); 
  par(mar=c(0, 0, 0, 0))   #  c(bottom, left, top, right)
  if (length(runs)>6) yl<-3 else yl<-2
  plot(1,1, type="n", yaxt="n", xaxt='n',xlab='',ylab="",ylim=c(-0.2,yl),xlim=c(0,3),main=NULL,frame.plot=F)
  for (i in (1:length(runs))) {
    #cat('\nx=',(i-1) %%  3,' y=',(i*2-1) %/% 6,'\n')
    text(x=(i-1) %%  3,y=(i*2-1) %/% 6,labels=runs[i],adj=c(0,0),offset=0,col=cols[i],cex=1.5)
  }
  box()
  close.screen(all = TRUE)    # exit split-screen mode
  
} 

doPlotF<-function(fit,plotP=FALSE){
  
  plotZ<-function(age='Age 1',inp) {
    for (i in (1:length(inp))) {
      y<-inp[[i]][age,]
      years<-as.numeric(names(y))
      if (i==1) {
        plot(x=years,y=y,col=cols[i],pch=as.character(i),lwd=2,xlab='Year', ylab='F at age',main=a,type='b',ylim=c(0,maxVal[age])) 
      } else points(x=years,y=y,col=cols[i],lwd=2,pch=as.character(i),type='b')  
      
    }
  }
  
  
  fs<-lapply(lapply(fit,faytable),t)
  maxVal<-apply(as.data.frame(fs),1,max)
  
  cols<-rainbow(length(fs))
  cols<-palette()
  
  X11(width=8,height=10)
  
  des<-matrix(c(0,1,0,0.10,
                0,1,0.10,1.0 ), nrow=2, ncol=4, byrow=TRUE)
  split.screen(des)     # split display into two screens (screen 1-2)
  split.screen(c(3,2), screen = 2) # now split the bottom half into 2x2 (screen 3-6)
  
  sc<-2
  for (a in names(maxVal)) {
    sc<-sc+1; screen(sc)
    par(mar= c(3, 4, 4, 2) + 0.1) #(bottom, left, top, right) 
    plotZ(age=a,inp=fs)
  }
  
  screen(1); 
  par(mar=c(0, 0, 0, 0))   #  c(bottom, left, top, right)
  if (length(fs)>6) yl<-3 else yl<-2
  plot(1,1, type="n", yaxt="n", xaxt='n',xlab='',ylab="",ylim=c(-0.2,yl),xlim=c(0,3),main=NULL,frame.plot=F)
  for (i in (1:length(fit))) {
    #cat('\nx=',(i-1) %%  3,' y=',(i*2-1) %/% 6,'\n')
    text(x=(i-1) %%  3,y=(i*2-1) %/% 6,labels=names(fit)[i],adj=c(0,0),offset=0,col=cols[i],cex=1.5)
  }
  box()
  close.screen(all = TRUE)    # exit split-screen mode
  
}


doPlotZ<-function(fit,plotP=FALSE,adjustZ=TRUE){
  
  plotZ<-function(age='Age 1',inp) {
    for (i in (1:length(inp))) {
      y<-inp[[i]][age,]
      years<-as.numeric(names(y))
      if (i==1) {
         plot(x=years,y=y,col=cols[i],pch=as.character(i),lwd=2,xlab='Year', ylab='Z',main=a,type='b',ylim=c(0,maxVal[age])) 
      } else points(x=years,y=y,col=cols[i],lwd=2,pch=as.character(i),type='b')  
      
      }
  }
 
  
  fs<-lapply(fit,ZwithMigration,adjustZ=adjustZ)
  maxVal<-apply(as.data.frame(fs),1,max)

  cols<-rainbow(length(fs))
  cols<-palette()
 
  X11(width=8,height=10)
  
  des<-matrix(c(0,1,0,0.10,
                0,1,0.10,1.0 ), nrow=2, ncol=4, byrow=TRUE)
  split.screen(des)     # split display into two screens (screen 1-2)
  split.screen(c(3,2), screen = 2) # now split the bottom half into 2x2 (screen 3-6)
  
  sc<-2
  for (a in names(maxVal)) {
    sc<-sc+1; screen(sc)
    par(mar= c(3, 4, 4, 2) + 0.1) #(bottom, left, top, right) 
    plotZ(age=a,inp=fs)
  }
 
  screen(1); 
  par(mar=c(0, 0, 0, 0))   #  c(bottom, left, top, right)
  if (length(fs)>6) yl<-3 else yl<-2
  plot(1,1, type="n", yaxt="n", xaxt='n',xlab='',ylab="",ylim=c(-0.2,yl),xlim=c(0,3),main=NULL,frame.plot=F)
  for (i in (1:length(fit))) {
    #cat('\nx=',(i-1) %%  3,' y=',(i*2-1) %/% 6,'\n')
    text(x=(i-1) %%  3,y=(i*2-1) %/% 6,labels=names(fit)[i],adj=c(0,0),offset=0,col=cols[i],cex=1.5)
  }
  box()
  close.screen(all = TRUE)    # exit split-screen mode
  
}


extractCatchN<-function(fit){
  cobs<-NULL
  for (i in (1:(dim(fit$data$idx1)[[2]]-1)))  cobs<-c(cobs,(fit$data$idx1[1,i]:fit$data$idx2[1,i]))
  catage<-matrix(exp(fit$data$logobs[cobs+1]),byrow=TRUE,ncol=nage)
  rownames(catage)<-head(fit$data$years,-1)
  colnames(catage)<-colnames((fit$data$propMat))
  return(catage)
}

