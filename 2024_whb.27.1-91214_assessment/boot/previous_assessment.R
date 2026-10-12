## Download model run from previous assessment year for comparison later on in reporting script

library(icesTAF)

sam_assessment <- "BW-2023"

sam_dir <-
  paste0(
    "https://stockassessment.org/datadisk/stockassessment/userdirs/user3/",
    sam_assessment,
    "/run/"
  )

files <- "model.RData"

for (file in files) {
  download(paste0(sam_dir, file))
}