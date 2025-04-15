# Principal Component Analysis
I will intially investigate populaiton structure using a PCA, a simple dimensionality reduction technique can help to understand my complex GBS dataset without modelling!

## Set up
First, we need to convert our .vcf in the relevant plink formats; note plink was/is made for human genomic data and thus expects chromosome location information. --allow-extra-chr allows me to proceed with non-standard chromosome names, and numbers. 
```Bash
plink --vcf ../HamFrogR08maxsnps1DP5.recode.vcf --make-bed --out plink_Ham --allow-extra-chr
```
