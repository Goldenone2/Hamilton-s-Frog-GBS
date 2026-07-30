# Filtering our VCF for TrioML
To reduce computational load and improve the reliability of this analysis, I filtered singletons. [bcftools](https://samtools.github.io/bcftools/bcftools.html#common_options) made my life much easier here. Vcftools doesn't have the filtering I required readily available.

```bash
cd TrioML
module load VCFtools
module load BCFtools
```
### singletons first
Vcftools --singletons created an output file listing all singletons and private doubletons.
```bash
vcftools --vcf HamFrogR08maxsnps1DP5.recode.vcf --singletons --out flagged
```
Cut takes the first two columns i.e. the chrom and pos (exclduing indv). Tail removes the header.
```bash
tail -n +2 flagged.singletons | cut -f1,2 > singletons_to_remove.txt
```
View tells bcftools that we want to subset/filter our .vcf file, -T is the list of targets with ^ telling us to exclude them, -O v says we want an uncompressed vcf file, -o is our file name.
```bash
bcftools view -T ^singletons_to_remove.txt HamFrogR08maxsnps1DP5.recode.vcf -O v -o HamFrog_nosingletons.vcf
bcftools view -H HamFrog_nosingletons.vcf | wc -l
```
