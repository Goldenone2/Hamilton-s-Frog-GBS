Swab.Tissue.Jitter
================
Hadley Muller
2025-04-07

I have completed SNP calling for the Hamilton’s Frog data, which
contains both tissues and, buccal swabs.

Because this is the first time anyone has utilised buccal swabbing on
*Leiopelma*, I’d like to produce box (or jitter) plots comparing
coverage depth and total sites.

### Data Import and Setup

``` r
# Clear workspace
rm(list = ls())

# Load required Packages
library(ggplot2)
library(dplyr)
library(patchwork)

# Load data
Dep <- read.table("Depth.idepth", header = TRUE)

#Add sample information to the dataframe
Data <- Dep %>%
  mutate(Tissue = case_when(
    grepl("^B", INDV) ~ "Toe Clip",   
    grepl("^T", INDV) ~ "Toe Clip",    
    grepl("^M_", INDV) ~ "Toe Clip",
    grepl("^MT", INDV) ~ "Buccal Swab",   
    grepl("^G", INDV) ~ "Toe Clip"
    ))
```
#Simple t-test
```r
#Subset data
Dtoeclip <- (Data[Data$Tissue == "Toe Clip", "MEAN_DEPTH"])
Dswab <- (Data[Data$Tissue == "Buccal Swab", "MEAN_DEPTH"])

Ntoeclip <- (Data[Data$Tissue == "Toe Clip", "N_SITES"])
Nswab <- (Data[Data$Tissue == "Buccal Swab", "N_SITES"])

#Welch's t-test
t.test(Dtoeclip, Dswab, var.equal=FALSE)

t.test(Ntoeclip, Nswab, var.equal=FALSE)
```

#Create Jitterplot

``` r
#Total Sites
plot1 <- ggplot(Data, aes(x = Tissue, y = N_SITES)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, colour = "cornflowerblue", alpha = 0.3) +
  theme_light() +
  labs(y= "Total sites", x = "Sample material")

#Coverage Depth
plot2 <- ggplot(Data, aes(x = Tissue, y = MEAN_DEPTH)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, colour = "darkorange", alpha = 0.3) +
  theme_light() +
  labs(y= "Mean Depth", x = "Sample material")

(plots = wrap_plots(plot1,plot2)) +
  plot_layout(axis_titles = "collect")
```

![](Swab.Tissue.Jitter_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->
