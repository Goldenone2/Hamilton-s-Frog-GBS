# Filtering our VCF for TrioML
TrioML is both computationally arduous. I'll do a little additonal filtering here to make our lives easier down the road by removing singletons.
[bcftools](https://samtools.github.io/bcftools/bcftools.html#common_options) can make our lives much easier here; where vcftools doesn't have the filtering we'de like.

````bash
cd TrioML
module load VCFtools
module load BCFtools
```

### singletons first
Vcftools --singeltons creates an output file listing all singletons and private doubletons.
```bash
vcftools --vcf HamFrogR08maxsnps1DP5.recode.vcf --singletons --out flagged
```
Cut takes the first two columns i.e. the chrom and pos (exclduing indv). Tail removes the header.
```bash
tail -n +2 flagged.singletons | cut -f1,2 > singletons_to_remove.txt
```
Bcftools, view tells it that we want to subset/filter our vcf file, -T is the list of targets with ^ telling us to exclude them, -O v says we want an uncompressed vcf file, -o is our file name.
```bash
bcftools view -T ^singletons_to_remove.txt HamFrogR08maxsnps1DP5.recode.vcf -O v -o HamFrog_nosingletons.vcf
bcftools view -H HamFrog_nosingletons.vcf | wc -l
```
I get a lot of warnings but seems like filtering worked (4724 SNPs left).
