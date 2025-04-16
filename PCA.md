# Principal Component Analysis
I will intially investigate populaiton structure using a PCA, a simple dimensionality reduction technique can help to understand my complex GBS dataset without modelling!

## Set up
First, we need to convert our .vcf in the relevant plink formats. Plink was made for human genomic data and thus expects: chromosome location information, and pedigree information: --allow-extra-chr allows me to proceed with non-standard chromosome names and numbers,  --double-id allows me to duplicate the ID of my samples for both "family" and "individual" ID

Notes from the [PLINK Documentation]{https://www.cog-genomics.org/plink/1.9/input}: if you're dealing with a draft assembly with lots of contigs, rather than actual autosomes—the standard PLINK build can handle that if you name your contigs 'contig1', 'contig2', etc. and use the --allow-extra-chr flag!
```sh
# mkdir PCA
cd PCA
awk '{if($0 !~ /^#/) print "contig"$0; else print $0}' ../HamFrogR08maxsnps1DP5.recode.vcf > PLINKvcf_with_contig.vcf
```
```sh
plink --vcf ../HamFrogR08maxsnps1DP5.recode.vcf --make-bed --out PLINK_Ham --allow-extra-chr --double-id
```
