## Run analysis, write model results

## Before:
## After:

library(icesTAF)
library(stockassessment)

mkdir("model")

# Load model input data
load("data/data_base.RData")
load("data/data_forecast.RData")

# Load configuration file
conf <- loadConf(dat.base, "boot/data/model.cfg")

# Set up default parameter values
par.base <- defpar(dat.base,conf) 
par.forecast <- defpar(dat.forecast,conf)

# Run base model and forecast model
fit.base<-sam.fit(dat.base,conf,par.base)
fit.forecast<-sam.fit(dat.forecast,conf,par.forecast)

# Residuals for base run
RES              <-residuals(fit.base)

# Process noise for base run
RESP             <-procres(fit.base)

# Retrospective runs for base run
RETRO             <-stockassessment::retro(fit.base,year=5)

# Save model output for base run and forecast run
save(fit.base,RES,RESP,RETRO,file="model/model_base.RData")
save(fit.forecast,file="model/model_forecast.RData")
