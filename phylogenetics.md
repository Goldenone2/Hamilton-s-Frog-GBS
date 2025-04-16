## Phylogenetics 
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

UFBoot does not converge after 1000 bootstraps, so I will trial with 3,000!

```sh
iqtree2 -nt 16 -s HamFrogR08maxsnps1DP5.recode.min4.phy -st DNA -m TIM+F+I+R2 -bb 3000  -pre Final.FrogTree
```
Download FIGtree to visualise it.
