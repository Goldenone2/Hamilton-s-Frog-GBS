# set working directory
setwd("~/uoo04306/frogs_gbs/TrioML")

# Clear workspace
rm(list = ls())

# Load Packages
library("vcfR")
library("ggplot2")
library("related")
library("adegenet")

# Run TrioML
results <- coancestry("HamGeno.txt", trioml = 1, trioml.num.reference = 82)
#RDS means we can save the R object for plotting later :) 
saveRDS(results, file = "TrioML_results.rds")