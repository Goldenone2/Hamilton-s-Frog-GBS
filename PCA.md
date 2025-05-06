# Principal Component Analysis
I will investigate the natural divergence or structure between Hamilton’s frog populations (esp. Takapourewa & Te Pakeka) using PCA, a dimensionalality reduction analysis that takes our high-dimensional data and reduces it to a few principal components.

### Set up
First, we need to convert our .vcf in the relevant plink formats. Plink was made for human genomic data and thus expects: chromosome location information, and pedigree information: --allow-extra-chr allows me to proceed with non-standard chromosome names and numbers,  --double-id allows me to duplicate the ID of my samples for both "family" and "individual" ID

Notes from the [PLINK Documentation](https://www.cog-genomics.org/plink/1.9/input) that if you're dealing with a draft assembly with lots of contigs, rather than actual autosomes—the standard PLINK build can handle that if you name your contigs 'contig1', 'contig2', etc. and use the --allow-extra-chr flag!
```sh
# mkdir PCA
cd PCA
awk '{if($0 !~ /^#/) print "contig"$0; else print $0}' ../HamFrogR08maxsnps1DP5.recode.vcf > PLINKvcf_with_contig.vcf
```
```sh
plink --vcf PLINKvcf_with_contig.vcf --make-bed --out PLINK_Ham --allow-extra-chr --double-id
```
### Peform PCA
```sh
plink --bfile PLINK_Ham --pca --out PCA_Ham --allow-extra-chr --double-id
```
## Peform PCA on a subset: translocation data only
To see whether there is more intricate structure between the Muad Island (natural source populaiton) and, Boat Bay or Motuara (translocated populations) I will run PCA on this small subset. I will make a .txt file with the individuals I want to --keep
```sh
cd /home/mulha552/uoo04306/frogs_gbs/PCA/PCA_Focused
module load VCFtools
vcftools --vcf ../PLINKvcf_with_contig.vcf --keep Maud.txt --recode --recode --out MaudIslandOnly
```

```sh
plink --vcf MaudIslandOnly.recode.vcf plink --make-bed --out PLINK_Maud --allow-extra-chr --double-id
plink --bfile PLINK_Maud --pca --out PCA_Maud --allow-extra-chr --double-id
```

# Data Visualisation in R
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
