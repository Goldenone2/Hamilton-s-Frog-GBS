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
# these are great to check in the terminal my R code is working

# Run TrioML
BoatBay_results <- coancestry("BoatBay.txt", dyadml =1, wang =1 )
#RDS means we can save the R object for plotting later :) 
saveRDS(BoatBay_results, file = "BoatBay_results.rds")

message("Boat Bay done")

TePakeka_results <- coancestry("TePakeka.txt", dyadml =1, wang =1)
saveRDS(TePakeka_results, file = "TePakeka_results.rds")

message("Te Pakeka done")

Motuara_results <- coancestry("Motuara.txt", dyadml =1, wang =1)
saveRDS(Motuara_results, file = "Motuara_results.rds")

message("Motuara done")

Takapourewa_results <- coancestry("Takapourewa.txt", dyadml =1, wang =1)
saveRDS(Takapourewa_results, file = "Takapourewa_results.rds")

message("Takapourewa done")


warnings()
