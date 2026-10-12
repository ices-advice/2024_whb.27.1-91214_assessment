## Extract results of interest, write TAF output tables

## Before:
## After:

library(icesTAF)
library(stockassessment)

mkdir("output")

# Load base model
load("model/model_base.RData")

pdf("output/diagnostics.pdf", 6, 6)
ssbplot(fit.base)
fbarplot(fit.base)
recplot(fit.base)
catchplot(fit.base)
fitplot(fit.base, fleet=1)
fitplot(fit.base, fleet=2)
plot(RESP)
plot(RES)

ssbplot(RETRO,main="RETRO SSB")
fbarplot(RETRO,partial=F,main="RETRO Fbar")
recplot(RETRO,main = "RETRO Recruitment" , las=0, drop=1)
catchplot( RETRO, main = "RETRO Catch")

dev.off()
