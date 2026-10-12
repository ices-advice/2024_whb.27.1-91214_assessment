## Preprocess data, write TAF data tables

## Before:
## After:

library(icesTAF)
library(stockassessment)

mkdir("data")

# Add data to boot/data folder
taf.boot(software = FALSE)

# Function to open sAM data files from boot/data folder
read.ices.taf <- function(...) {
  read.ices(taf.data.path(...))
}

# Load data
cn<-read.ices.taf("cn.dat")
cw<-read.ices.taf("cw.dat")
dw<-read.ices.taf("dw.dat")
lf<-read.ices.taf("lf.dat")
lw<-read.ices.taf("lw.dat")
mo<-read.ices.taf("mo.dat")
nm<-read.ices.taf("nm.dat")
pf<-read.ices.taf("pf.dat")
pm<-read.ices.taf("pm.dat")
sw<-read.ices.taf("sw.dat")
surveys<-read.ices.taf("survey.dat")

# Set up data for base run and forecast run

# 1. Base run includes only data up to the intermediate year
# Remove the last year from data files containing this additional year
mo.base <- mo[-nrow(mo),]
nm.base <- nm[-nrow(nm),]
pf.base <- pf[-nrow(pf),]
pm.base <- pm[-nrow(pm),]
sw.base <- sw[-nrow(sw),]
surveys.base <- surveys
surveys.base$IBWSS <- surveys.base$IBWSS[-nrow(surveys.base$IBWSS),]
attributes(surveys.base$IBWSS)$time <- attributes(surveys$IBWSS)$time #add again attribute
attributes(surveys.base$IBWSS)$twofirst <- attributes(surveys$IBWSS)$twofirst #add again attribute

# Setup as SAM object
dat.base<-setup.sam.data(surveys=surveys.base,
                         residual.fleet=cn, 
                         prop.mature=mo.base, 
                         stock.mean.weight=sw.base, 
                         catch.mean.weight=cw, 
                         dis.mean.weight=dw, 
                         land.mean.weight=lw,
                         prop.f=pf.base, 
                         prop.m=pm.base, 
                         natural.mortality=nm.base, 
                         land.frac=lf)

# 2. Forecast run includes additional year for some of the files, so they're ready as they are
# Setup as SAM object
dat.forecast<-setup.sam.data(surveys=surveys,
                    residual.fleet=cn, 
                    prop.mature=mo, 
                    stock.mean.weight=sw, 
                    catch.mean.weight=cw, 
                    dis.mean.weight=dw, 
                    land.mean.weight=lw,
                    prop.f=pf, 
                    prop.m=pm, 
                    natural.mortality=nm, 
                    land.frac=lf)

# Save data
save(dat.base, file="data/data_base.RData")
save(dat.forecast, file="data/data_forecast.RData")
