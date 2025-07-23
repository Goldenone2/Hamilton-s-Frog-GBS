# Phylogenetics 
I will do a phylogenetic analysis using IQTREE; to understand the natural divergence between Maud Island and Stephen's Island; population level dynamics between founder and source populations will not have an effect here.
```sh
mdkir tree
cd tree
```
First , I need convert the SNPs in my .vcf file to a phylogenetic input format for IQtree: phylip, using [vcf2phylip](https://github.com/edgardomortiz/vcf2phylip)
https://github.com/edgardomortiz/vcf2phylip

```sh
git clone https://github.com/edgardomortiz/vcf2phylip.git
vcf2phylip/vcf2phylip.py --input ../HamFrogR08maxsnps1DP5.recode.vcf
```
I have no idea what subsitution model is best for my data, so I'll try a model selection method outlined [here](http://www.iqtree.org/doc/Tutorial). MFP is a specially model specifier *modelfinderplus* and looks for the model that gives the lowest bayesian information criterion.

```sh
module load IQ-TREE
iqtree2 -s HamFrogR08maxsnps1DP5.recode.min4.phy -m MFP -st DNA
```
Let's have a look at the results:
```
Best-fit model: TIM+F+I+R2 chosen according to BIC
```
Now, I can run my final tree, with the TIM+F+I+R2 model and 1000 bootstraps. Our parameters specify: -TIM trnasitional model, allows different rates for transitions and transversion, +F use tghe emprirical base frequencies rather than assuming equal (frogs have higher GC proportions than other species), +I account for sites which don't evolve at all and +R2 use the FreeRate model with 2 discrete categories (modelling how some sites evolve faster than others). 

```sh
iqtree2 -nt 16 -s HamFrogR08maxsnps1DP5.recode.min4.phy -st DNA -m TIM+F+I+R2 -bb 1000 -pre Final.FrogTree
```

## Visualisations in R


Righto, I’ve now done my maximum-likelihood phylogenetic analysis using
IQTREE2. I will visualise this using the package ‘ggtree’ loosely
following [this
tutorial](https://arftrhmn.net/creating-a-publication-quality-phylogeny-using-ggtree/).

### Data Import and Setup

``` r
# Clear workspace
rm(list = ls())

#Load required packages

# install.packages("phangorn")
# install.packages("BiocManager")
# BiocManager::install("ggtree")
# BiocManager::install("treeio")

library(phangorn)
```

    ## Loading required package: ape

``` r
library(ggtree)
```

    ## ggtree v3.16.0 Learn more at https://yulab-smu.top/contribution-tree-data/
    ## 
    ## Please cite:
    ## 
    ## Guangchuang Yu, David Smith, Huachen Zhu, Yi Guan, Tommy Tsan-Yuk Lam.
    ## ggtree: an R package for visualization and annotation of phylogenetic
    ## trees with their covariates and other associated data. Methods in
    ## Ecology and Evolution. 2017, 8(1):28-36. doi:10.1111/2041-210X.12628

    ## 
    ## Attaching package: 'ggtree'

    ## The following object is masked from 'package:ape':
    ## 
    ##     rotate

``` r
library(treeio)
```

    ## treeio v1.32.0 Learn more at https://yulab-smu.top/contribution-tree-data/
    ## 
    ## Please cite:
    ## 
    ## LG Wang, TTY Lam, S Xu, Z Dai, L Zhou, T Feng, P Guo, CW Dunn, BR
    ## Jones, T Bradley, H Zhu, Y Guan, Y Jiang, G Yu. treeio: an R package
    ## for phylogenetic tree input and output with richly annotated and
    ## associated data. Molecular Biology and Evolution. 2020, 37(2):599-603.
    ## doi: 10.1093/molbev/msz240

``` r
library(ggnewscale)

# Load data
Tree <- read.iqtree('Final.FrogTree.treefile')
meta <- read.csv('meta.csv', fileEncoding = "UTF-8")

#Set factor levels
meta$Population <- factor(meta$Population, levels = c("Takapourewa", "Te Pākeka", "Boat Bay", "Motuara" ))

# Midpoint root the phylo *inside* the tree data object, keeping other metadata.
Tree@phylo <- midpoint(Tree@phylo)
```

## Visualisations

``` r
Tree1 <- ggtree(Tree) %<+% meta + 
  geom_tippoint(aes(colour = Population))

Tree1
```

![](Hamilton.Tree_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

### Bootstrap Support

IQTree2 comes with bootstrap support in both SH_aLRT and UFBoot values.
In order to annotate or filter our tree, we want to know which nodes are
supported by both values by adding a new column to our tree data.

``` r
# Add bootstrap support
Tree1$data$bootstrap <- '0'

Tree1$data[which(Tree1$data$SH_aLRT >= 70 & Tree1$data$UFboot >= 70),]$bootstrap <- '1'

Tree1$data$bootstrap <- as.factor(Tree1$data$bootstrap)

Tree1 <- Tree1 + new_scale_color() +
  geom_tree(aes(color=bootstrap)) +
  scale_colour_manual(name = 'Bootstrap', values = c('1' = 'black', '0' = 'grey'), guide = "none")

Tree1
```

![](Hamilton.Tree_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
summary(Tree1$data$SH_aLRT)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    ##   24.00   58.50   70.00   69.01   81.50  100.00      84

``` r
summary(Tree1$data$UFboot)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    ##   24.00   58.50   70.00   69.01   81.50  100.00      84
