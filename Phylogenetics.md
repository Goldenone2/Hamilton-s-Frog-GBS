# Phylogenetic Modelling
I undertook a simple phylogenetic analysis using IQTREE. Primarily, to visualise the natural divergence between Maud Island and Stephen's Island; population level dynamics between founder and source populations will not have an effect here.
```sh
mdkir tree
cd tree
```
First , I  convert the SNPs in my .vcf file to a phylogenetic input format for IQtree: phylip, using [vcf2phylip](https://github.com/edgardomortiz/vcf2phylip)

```sh
git clone https://github.com/edgardomortiz/vcf2phylip.git
vcf2phylip/vcf2phylip.py --input ../HamFrogR08maxsnps1DP5.recode.vcf
```
I had no idea what substitution model is best for my data, so I tried a model selection method outlined [here](http://www.iqtree.org/doc/Tutorial). MFP is a specially model specifier *modelfinderplus* and looks for the model that gives the lowest Bayesian information criterion.

```sh
module load IQ-TREE
iqtree2 -s HamFrogR08maxsnps1DP5.recode.min4.phy -m MFP -st DNA
```
Let's have a look at the results:
```
Best-fit model: TIM+F+I+R2 chosen according to BIC
```
I ran my final tree, with the TIM+F+I+R2 model and 1000 bootstraps. 
Our parameters specify: -TIM trnasitional model, allows different rates for transitions and transversion, +F use tghe emprirical base frequencies rather than assuming equal (frogs have higher GC proportions than other species), +I account for sites which don't evolve at all and +R2 use the FreeRate model with 2 discrete categories (modelling how some sites evolve faster than others). 

```sh
iqtree2 -nt 16 -s HamFrogR08maxsnps1DP5.recode.min4.phy -st DNA -m TIM+F+I+R2 -bb 1000 -pre Final.FrogTree
```

# Phylogenetic Tree
Using the package ‘ggtree,’ I loosely followed [this tutorial](https://arftrhmn.net/creating-a-publication-quality-phylogeny-using-ggtree/) to produce a midpoint rooted visualisation.

## Data Import and Setup
Install BiocManager related package if required.
``` r
# install.packages("BiocManager")
# BiocManager::install("ggtree")
# BiocManager::install("treeio")
```
``` r
# Clear workspace
rm(list = ls())

# Load Packages
library(ggtree)
```

    ## ggtree v3.16.0 Learn more at https://yulab-smu.top/contribution-tree-data/
    ## 
    ## Please cite:
    ## 
    ## Guangchuang Yu.  Data Integration, Manipulation and Visualization of
    ## Phylogenetic Trees (1st edition). Chapman and Hall/CRC. 2022,
    ## doi:10.1201/9781003279242, ISBN: 9781032233574

``` r
library(phangorn)
```

    ## Loading required package: ape

    ## 
    ## Attaching package: 'ape'

    ## The following object is masked from 'package:ggtree':
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
library(patchwork)
```

    ## Warning: package 'patchwork' was built under R version 4.5.3

``` r
library(ggplot2)
```

    ## Warning: package 'ggplot2' was built under R version 4.5.3

``` r
# Load data
Tree <- read.iqtree('Final.FrogTree.treefile')
meta <- read.csv('meta.csv', fileEncoding = "UTF-8")

# Set factor levels
meta$Population <- factor(meta$Population, levels = c("Takapourewa", "Te Pākeka", "Boat Bay", "Motuara" ))

# Midpoint root the the tree data object, preserving node data
Tree@phylo$node.label <- as.character(Tree@data$UFboot)

# Midpoint root
Tree@phylo <- midpoint(Tree@phylo)
```

## Visualisations
First a basic tree, coloured by population. Using *branch.length =‘none’* allows a cladogram without branch length scaling; useful to fit bootstrap values on the tree later.
``` r
Tree1 <- 
  ggtree(Tree, branch.length = "none") %<+% meta +
  geom_tippoint(aes(colour = Population)) +
  scale_colour_manual(values = c("cornflowerblue", "darkorange","#F763E0", "#44AA99")) +
  theme(legend.text = element_text(size=13)) +
  theme(legend.title = element_text(size=14))
```

    ## Warning: `aes_()` was deprecated in ggplot2 3.0.0.
    ## ℹ Please use tidy evaluation idioms with `aes()`
    ## ℹ The deprecated feature was likely used in the ggtree package.
    ##   Please report the issue at <https://github.com/YuLab-SMU/ggtree/issues>.
    ## This warning is displayed once per session.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning in fortify(data, ...): Arguments in `...` must be used.
    ## ✖ Problematic arguments:
    ## • as.Date = as.Date
    ## • yscale_mapping = yscale_mapping
    ## • hang = hang
    ## ℹ Did you misspell an argument name?

``` r
Tree1
```

![](Phylogenetics_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

### Bootstrap Support
IQTree2 provides node support with:

- SH_aLRT, a likelihood measure
- UFBoot, a more *traditional* bootstrap measure.

For the purpose of this simple tree, I  annotate with UFBoot. However, in other cases both measures can be used to evaluate node support. I need to moved the legend to the top left for patchwork to create an effective visualisation.

``` r
# Add bootstrap support
Tree_UF <- Tree1 + geom_text2(aes(subset = !isTip, label = label),size = 2.2,hjust = -0.3, fontface = "bold") +
  theme(
   legend.position = c(0.02, 0.95),
  legend.justification = c(0, 1) # anchor top left
  )

Tree_UF
```

![](Phylogenetics_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

### Export final plot
This phylogeny and the PCAs were run locally, and plotted together using patchwork. The final plots with red stars indicating nodes with >80% bootstrap support, were added in Photoshop. 

``` r
Plot_Manuscript <- 
  ((PCA_All | PCA_Maud)/ Tree_UF) +
  plot_layout(heights  = c(1, 2), 
    axis_titles = "collect") +
  plot_annotation(tag_levels = "A")

Plot_Manuscript

ggsave(filename="Figure2.pdf", plot = Plot_Manuscript, dpi = 300, width = 11, height = 13, device = cairo_pdf)
```
