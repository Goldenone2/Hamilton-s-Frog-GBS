# set working directory
setwd("~/uoo04306/frogs_gbs/TrioML")

# Clear workspace
rm(list = ls())

# Load Packages
library("vcfR")
library("ggplot2")
library("related")
library("adegenet")
library("dplyr")

message("Packages_Loaded")

# Run TrioML
BoatBay_results <- coancestry("BoatBay.txt", trioml = 1, trioml.num.reference = 20)
#RDS means we can save the R object for plotting later :) 
saveRDS(BoatBay_results, file = "BoatBay_results.rds")

message("Boat Bay done")

TePakeka_results <- coancestry("TePakeka.txt", trioml = 1, trioml.num.reference = 20)
saveRDS(TePakeka_results, file = "TePakeka_results.rds")

message("Te Pakeka done")

Motuara_results <- coancestry("Motuara.txt", trioml = 1, trioml.num.reference = 27)
saveRDS(Motuara_results, file = "Motuara_results.rds")

message("Motuara done")

Takapourewa_results <- coancestry("Takapourewa.txt", trioml = 1, trioml.num.reference = 15)
saveRDS(Takapourewa_results, file = "Takapourewa_results.rds")

message("Takapourewa done")


warnings()
