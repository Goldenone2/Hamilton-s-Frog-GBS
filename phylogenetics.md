## phylogenetics 
I will do a phylogenetic analysis using IQTREE; to understand the natural populaiton structure between Maud Island and Stephen's Island; population level dynamics between founder and source populations will not have an effect here.
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
I have no idea what subsitution model is best for my data, so I'll try a model selection method outlined [here](http://www.iqtree.org/doc/Tutorial). MFP is a specially model specifier *modelfinderplus* and looks for the model that gives the lowest bayesian information criterion. I ran this as a quick SLURM job :) 

```sh
module load IQ-TREE
iqtree2 -s HamFrogR08maxsnps1DP5.recode.min4.phy -m MFP -st DNA
```
Then I run IQtree using my own conda environment, the GTR+G model with 1000 bootstraps.

```sh
iqtree2 -nt 16 -s HamFrogR08maxsnps1DP5.recode.min4.phy -st DNA -m GTR+G -bb 1000  -pre inferred
```
download FIGtree to visualise it.
