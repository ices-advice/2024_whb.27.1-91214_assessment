Extended<-FALSE  # MV addition, to specify if (dummy) survey data exists for the year after with most recent (preliminary) catch data)

basefit<-NULL
if(file.exists("boot/data/previous_assessment/model.RData")){
  local({load("boot/data/previous_assessment/model.RData"); basefit<<-fit})
  if(abs(logLik(basefit)-logLik(fit.base))<1.0e-9)basefit<<-NULL
}


plotcounter<-1
tit.list<-list()

setcap<-function(title="", caption=""){   
  tit.list[length(tit.list)+1]<<-paste("# Title",plotcounter,"\n")
  tit.list[length(tit.list)+1]<<-paste(title,"\n")
  tit.list[length(tit.list)+1]<<-paste("# Caption",plotcounter,"\n")
  tit.list[length(tit.list)+1]<<-paste(caption,"\n")
  plotcounter<<-plotcounter+1 
}


############################## plots ##############################
plots<-function(){
  par(cex.lab=1, cex.axis=1, mar=c(5,5,1,1))
  
  if(exists("fit")){
    if (is.null(basefit)) fits<-fit else fits <- c(fit,basefit)
    
    ssbplot(fits, addCI=TRUE,drop=0)
    stampit()
   
    setcap("Spawning stock biomass", "Spawning stock biomass. 
           Estimates from the current run and point wise 95% confidence 
           intervals are shown by black line and shaded area.")
    
    fbarplot(fits, addCI=TRUE,drop=Extended)
    stampit()
    
    recplot(fits, addCI=TRUE, las=0,drop=Extended)
    stampit()
    
    catchplot(fits, addCI=TRUE)
    stampit()
    
    par(mfrow=c(2,2))
    catchplot(fits, addCI=TRUE)
    recplot(fits, las=0, drop=Extended)
    fbarplot(fits, las=0, drop=Extended)
    ssbplot(fits, las=0, drop=0)
    
    stampit()
    par(mfrow=c(1,1))
     par(mfrow=c(3,1))
    par(mar=c(c(3, 4, 2.0, 2) + 0.1)) #c(bottom, left, top, right) 
    recplot(fits, las=0, drop=Extended)
    fbarplot(fits, las=0, drop=Extended)
    ssbplot(fits, las=0, drop=0)
    stampit()
    par(mfrow=c(1,1))
    par(cex.lab=1, cex.axis=1, mar=c(5,5,1,1))
    
    
    plot(fits)
    stampit()
    
    plot(ypr(fit, aveYears = 20))
    stampit()
    
    obscorrplot(fit)
    stampit()
    
    
    for(f in 1:fit$data$noFleets){
      fitplot(fit, fleets=f)
      stampit()
    }
    
    
    FF<-faytable(fit)
    if (Extended) FF<-head(FF,-1)
    
    matplot(rownames(FF),FF,lty=1:ncol(FF),col=1:ncol(FF), type='l', xlab='Year', ylab='F',lwd=4)
    legend('topleft', col=1:ncol(FF), lty=1:ncol(FF), legend=colnames(FF), bty='n', lwd=4)
    stampit()
    
    # MV exploitation pattern  ( F divided by Fbar)
    Fbar<-fbartable(fit)[,1]
    if (Extended) Fbar<-head(Fbar,-1)
    
    FF1<-FF/ rep(Fbar,times=dim(FF)[[2]])
    matplot(rownames(FF1),FF1,lty=1:ncol(FF),col=1:ncol(FF), type='l', xlab='Year', lwd=4,ylab='F scaled to Fbar')
    legend('topleft', col=1:ncol(FF), lty=1:ncol(FF), legend=colnames(FF), bty='n', lwd=4,ncol=2)
    stampit()
    
    
    # MV exploitation pattern  ( F divided by average F
    FF2<-FF/ rep(apply(FF,1,mean),times=dim(FF)[[2]])
    matplot(rownames(FF2),FF2,lty=1:ncol(FF),col=1:ncol(FF), type='l', xlab='Year', lwd=4,ylab='F scaled to average F')
    legend('topleft', col=1:ncol(FF), lty=1:ncol(FF), legend=colnames(FF), bty='n', lwd=4,ncol=2)
    stampit()
    
    op2<-par(mar=c(5, 4, 2, 2) + 0.1, mfcol=c(2,2))
    matplot(rownames(FF),FF,lty=1:ncol(FF),col=1:ncol(FF), type='l', xlab='Year', ylab='F',lwd=4)
    legend('topleft', col=1:ncol(FF), lty=1:ncol(FF), legend=colnames(FF), bty='n', lwd=4)
    matplot(rownames(FF1),FF1,lty=1:ncol(FF),col=1:ncol(FF), type='l', xlab='Year', lwd=4,ylab='F scaled to Fbar')
    legend('topleft', col=1:ncol(FF), lty=1:ncol(FF), legend=colnames(FF), bty='n', lwd=4,ncol=2)
    matplot(rownames(FF2),FF2,lty=1:ncol(FF),col=1:ncol(FF), type='l', xlab='Year', lwd=4,ylab='F scaled to average F')
    legend('topleft', col=1:ncol(FF), lty=1:ncol(FF), legend=colnames(FF), bty='n', lwd=4,ncol=2)
    stampit()
    
    par(mfrow=c(1,1))
    par(cex.lab=1, cex.axis=1, mar=c(5,5,1,1))
    
    # Selectivity of the Fishery
    
    FF<-faytable(fit)
    if (Extended) FF<-head(FF,-1)
    
    sel <-FF/ rep(Fbar,times=dim(FF)[[2]])
    sel <- t(sel)
    age.sel <- as.numeric(dimnames(sel)[[1]])
    
    par(mfrow=c(3,3), mai=c(0.4,0.5,0.2,0.2))
    for(i in seq(1,dim(sel)[2],5)){
      plot(age.sel, sel[,i], type="l", xlab="", ylab="", lwd=1.5, ylim=c(0,max(sel)))
      if (i+1<dim(sel)[2])try(lines(age.sel, sel[,i+1], col="red", lwd=1.5))
      if (i+2<dim(sel)[2]) try(lines(age.sel, sel[,i+2], col="blue", lwd=1.5))
      if (i+3<dim(sel)[2]) try(lines(age.sel, sel[,i+3], col="green", lwd=1.5))
      if (i+4<dim(sel)[2]) try(lines(age.sel, sel[,i+4], col="blue3", lwd=1.5))
      
      legend("bottomright",paste(c(colnames(sel)[i],colnames(sel)[i+1],colnames(sel)[i+2],colnames(sel)[i+3],colnames(sel)[i+4])), lty=rep(1,5), col=c("black","red","blue","green","blue3"),bty="n")
    }
    mtext("Age", 1, outer=T, line=1)
    mtext("F/Fbar", 2, outer=T, line=1)
    stampit(); 
    
    
    par(mfrow=c(1,1), mar=c(4,5,2,2))
    mn.sel <- apply(sel,1,mean)
    sd.sel <- apply(sel,1,quantile,probs=c(0.025,0.975))
    #age.sel <- minAge:(nrow(sel)-1+minAge)
    plot(age.sel, mn.sel, ylim=c(0,max(sd.sel)), type="l", xlab='Age', ylab="F/Fbar", lwd=2)
    for(i in 1:length(age.sel)){
      lines(rep(age.sel[i],2), sd.sel[,i])
    }
    stampit();setcap("Selection pattern", "Average (and standard deviation of) selection pattern all years")
    
    
    par(cex.lab=1, cex.axis=1, mar=c(5,5,1,1),mfrow=c(1,1)) 
    
    
    
    #par(op)
    
    
    #Q<-fit$pl$logFpar
    #Qsd<-fit$plsd$logFpar
    #key<-fit$conf$keyLogFpar
    #fun<-function(x)if(x<0){NA}else{Q[x+1]}
    #FF<-Vectorize(fun)
    #ages<-fit$conf$minAge:fit$conf$maxAge
    #matplot(ages, exp(t(matrix(FF(key), nrow=5))), type="l", lwd=5, lty="solid", xlab="Ages", ylab="Q")
    #legend("topright", lwd=5, col=2:5, legend=attr(fit$data, "fleetNames")[2:5])
    #stampit()
    
  }  
  
  if(exists("RES")){  
    plot(RES)
    par(mfrow=c(1,1))
    stampit()
    
  }
  
  
  if(exists("RESP")){  
    plot(RESP)
    par(mfrow=c(1,1))
    stampit()
  } 
  
  
  if(exists("LO")){  
    ssbplot(LO)
    stampit()
    
    fbarplot(LO)
    stampit()
    
    recplot(LO)
    stampit()
    
    catchplot(LO)
    stampit()
    
    par(mfrow=c(2,2))
    recplot(LO)
    fbarplot(LO)
    ssbplot(LO)
    catchplot(LO)
    stampit()
    par(mfrow=c(1,1))
    
  } 
  
  if(exists("RETRO")){ 
  
    sdState<-function(fit, y=max(fit$data$years)-1:0){
      idx <- names(fit$sdrep$value) == "logR"
      sdLogR<-fit$sdrep$sd[idx][fit$data$years%in%y]
      idx <- names(fit$sdrep$value) == "logssb"
      sdLogSSB<-fit$sdrep$sd[idx][fit$data$years%in%y]
      idx <- names(fit$sdrep$value) == "logfbar"
      sdLogF<-fit$sdrep$sd[idx][fit$data$years%in%y]
      ret<-cbind(sdLogR, sdLogSSB, sdLogF)
      rownames(ret)<-y
      colnames(ret)<-c("sd(log(R))", "sd(log(SSB))", "sd(log(Fbar))")
      return(ret)
    }
    
    sdtab<<-sdState(fit,y=max(fit$data$years)-20:0)
    
    
    mho_Rho<-function() {
      fs<-summary(fit)[,c(1,4,7)]
      n<-length(RETRO)
      rho<-matrix(0,nrow=n,ncol=3)
      colnames(rho)<-colnames(fs)
      rownames(rho)<-head(tail(rownames(fs),n+1),n)
      for (i in (1:n)) {
        rs<-summary(RETRO[[i]])[,c(1,4,7)] 
        y<-tail(rownames(rs),1)
        for (j in (1:3)) rho[y,j]<-(rs[y,j]-fs[y,j])/fs[y,j]
      }
       return(rho)

    }
     if (!Extended) {
      rho<-mho_Rho()
      rho.mean<-colMeans(rho)
      rho.out<<-rbind(rho,rho.mean)
    }
     
    
    # SSB plot 
  
    ssbplot(RETRO, las=0, drop=0)
    if (!Extended) legend("topright", legend=paste("Mohn's rho =",round(rho.mean[2],3)), bty="n")
    stampit()
    
    
    fbarplot(RETRO, las=0, drop=0)
    if (!Extended) legend("topright", legend=paste("Mohn's rho =",round(rho.mean[3],3)), bty="n")
    stampit()
    
    recplot(RETRO, las=0, drop=0)
    if (!Extended) legend("topright", legend=paste("Mohn's rho =",round(rho.mean[1],3)), bty="n")
    
    stampit()
    
    catchplot(RETRO)
    stampit()
    
    par(mfrow=c(2,2))
    recplot(RETRO, las=0, drop=0)
    if (!Extended) legend("topright", legend=paste("Mohn's rho =",round(rho.mean[1],3)), bty="n")
    
    fbarplot(RETRO, las=0, drop=0)
    if (!Extended) legend("topright", legend=paste("Mohn's rho =",round(rho.mean[3],3)), bty="n")
    
    ssbplot(RETRO, las=0, drop=0)
    if (!Extended) legend("topright", legend=paste("Mohn's rho =",round(rho.mean[2],3)), bty="n")
    
    catchplot(RETRO)
    stampit()
    par(mfrow=c(1,1))
  } 
  
  if(exists("FC")){  
    #MV  lapply(FC, function(f){plot(f); title(attr(f,"label"), outer=TRUE, line=-1); stampit()})
  }  
  if(exists("RETRO")){ 
    matplot(x=rownames(sdtab),y=sdtab,xlab='Year', ylab="CV",pch='RSF',type='b')
    stampit()
  }
}


setwd('report')
file.remove(dir(pattern='png$'))
stamp<-gsub('-[[:digit:]]{4}$','',gsub(':','.',gsub(' ','-',gsub('^[[:alpha:]]{3} ','',date()))))
png(filename = paste(stamp,"_%03d.png", sep=''), width = 480, height = 480,
    units = "px", pointsize = 10, bg = "white")
plots()    
dev.off()

writeLines(unlist(tit.list),'titles.cfg') 

png(filename = paste("big_",stamp,"_%03d.png", sep=''), width = 1200, height = 1200, 
    units = "px", pointsize = 20, bg = "white")
plots()    
dev.off()

#pdf(onefile=FALSE, width = 8, height = 8)
#  plots()    
#dev.off()

file.remove(dir(pattern='html$'))

tsb<-tsbtable(fit)
colnames(tsb)<-c("TSB","Low", "High")
tab.summary <- cbind(summary(fit), tsb)
xtab(tab.summary, caption=paste('Table 1. Estimated recruitment, spawning stock biomass (SSB), 
                                and average fishing mortality','.',sep=''), cornername='Year', 
     file=paste(stamp,'_tab1.html',sep=''), dec=c(0,0,0,0,0,0,3,3,3,0,0,0))


ftab <- faytable(fit)
xtab(ftab, caption=paste('Table 2. Estimated fishing mortality at age','.',sep=''), cornername='Year \ Age', 
     file=paste(stamp,'_tab2.html',sep=''), dec=rep(3,ncol(ftab)))

ntab <- ntable(fit)
xtab(ntab, caption=paste('Table 3. Estimated stock numbers at age','.',sep=''), cornername='Year \ Age', 
     file=paste(stamp,'_tab3.html',sep=''), dec=rep(0,ncol(ntab)))

ptab <- partable(fit)
xtab(ptab, caption=paste('Table 4. Table of model parameters','.',sep=''), cornername='Parameter name', 
     file=paste(stamp,'_tab4.html',sep=''), dec=rep(3,ncol(ptab)))


cwtab<-fit$data$catchMeanWeight

cwtab<-cwtab[,,1] # inserted by MV 29 August 2024
xtab(x=cwtab, caption=paste('Table 5. Mean weigh at age in the catch (kg)','.',sep=''), cornername='Year \ Age', 
     file=paste(stamp,'_tab5.html',sep=''), dec=rep(3,ncol(cwtab)))


cntab<-t(matrix(exp(fit$data$logobs[fit$data$aux[,'fleet']==1]),nrow=10))
dimnames(cntab)<-dimnames(cwtab)
xtab(cntab, caption=paste('Table 6. Catch at age','.',sep=''), cornername='Year \ Age', 
     file=paste(stamp,'_tab6.html',sep=''), dec=rep(0,ncol(cntab)))




yieldtab<-catchtable(fit,obs.show=TRUE)
xtab(yieldtab, caption=paste('Table 7. Estimated and observed (sum of product of C*CW) catch','.',sep=''), cornername='Year', 
     file=paste(stamp,'_tab7.html',sep=''), dec=rep(0,ncol(yieldtab)))



mtab <- modeltable(c(Current=fit, base=basefit))
mdec <- c(2,0,2,6)[1:ncol(mtab)]
xtab(mtab, caption=paste('Table 8. Model fitting','.',sep=''), cornername='Model', 
     file=paste(stamp,'_tab8.html',sep=''), dec=mdec)
if (!Extended) xtab(rho.out, caption=paste("Table 10. Mohn'r rho",'.',sep=''), cornername='Last data year', 
     file=paste(stamp,'_tab10.html',sep=''), dec=rep(3,ncol(rho.out)))

sdState<-function(fit, y=max(fit$data$years)-1:0){
  idx <- names(fit$sdrep$value) == "logR"
  sdLogR<-fit$sdrep$sd[idx][fit$data$years%in%y]
  idx <- names(fit$sdrep$value) == "logssb"
  sdLogSSB<-fit$sdrep$sd[idx][fit$data$years%in%y]
  idx <- names(fit$sdrep$value) == "logfbar"
  sdLogF<-fit$sdrep$sd[idx][fit$data$years%in%y]
  ret<-cbind(sdLogR, sdLogSSB, sdLogF)
  rownames(ret)<-y
  colnames(ret)<-c("sd(log(R))", "sd(log(SSB))", "sd(log(Fbar))")
  return(ret)
}

sdtab<-sdState(fit,y=max(fit$data$years)-20:0)

xtab(sdtab, caption=paste('Table 12. Model CV','.',sep=''), cornername=' ', 
     file=paste(stamp,'_tab12.html',sep=''), dec=rep(2,3))

# matplot(x=rownames(sdtab),y=sdtab,xlab='Year', ylab="CV",pch='RSF',type='b')



if(exists("FC") & FALSE){  
  ii<-0
  lapply(FC, function(f){
    ii<<-ii+1;
    tf<-attr(f,"tab");
    dec<-c(3,3,3,rep(0,ncol(tf)-3));
    xtab(tf, caption=paste0('Forecast table ',ii,'. ', attr(f,"label"),'.'), 
         cornername='Year', file=paste(stamp,'_tabX',ii,'.html',sep=''), dec=dec);      
  })
}  

setwd("..") 

