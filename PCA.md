PCA
================
Hadley Muller
2025-04-16

# Analysis Overview

I will investigate the natural divergence or structure between
Hamilton’s frog populations (esp. Takapourewa & Te Pakeka) using PCA, a
dimensionalality reduction analysis that takes our high-dimensional data
and reduces it to a few principal components.

### Data Import and Setup

``` r
# Clear workspace
rm(list = ls())

#Load required packages
library(ggplot2)
library(ggrepel)
```

    ## Warning: package 'ggrepel' was built under R version 4.4.3

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
# Load eigenvectors and eigenvalues
eigenvec <- read.delim("PCA_Ham.eigenvec", header = FALSE, sep = " ")
eigenval <- read.delim("PCA_Ham.eigenval", header = FALSE)
```

I’ll tidy this data set, PLINK’s output aren’t informative and I need to
add population information.

``` r
#remove extra ID column
eigenvec <- eigenvec[,-1]

#set names
names(eigenvec)[1] <- "INDV"
names(eigenvec)[2:ncol(eigenvec)] <- paste0("PC", 1:(ncol(eigenvec)-1))

#Provide population information
eigenvec <- eigenvec %>%
  mutate(Pop = case_when(
    grepl("^B", INDV) ~ "Boat Bay",   
    grepl("^T", INDV) ~ "Stephens Island",    
    grepl("^M_", INDV) ~ "Maud Island",
    grepl("^MT", INDV) ~ "Motuara",   
    grepl("^G", INDV) ~ "Maud Island"
    ))
```

\#Visualise our Data

\###EigenValues

``` r
#Calculate the percentage variance explained
PVE <- data.frame(PC = 1:20, pve = eigenval$V1/sum(eigenval$V1)*100)

ggplot(PVE, aes(x = PC, y = pve)) + 
  geom_bar(stat="identity") +
  theme_light()
```

![](PCA_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
ggplot(eigenvec, aes(x = PC1, y = PC2, colour = Pop)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(
    y = paste0("Principal component 2 (",round(PVE$pve[2], 2)," %)"), 
    x = paste0("Principal component 1 (",round(PVE$pve[1], 2)," %)")) + 
  theme_light()
```

![](PCA_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->
