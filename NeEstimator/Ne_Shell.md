# NeEstimator
Convert vcf to genpop fomrat for NeEstimator.
```bash
# mkdir Ne
cd Ne
```
First, we need to convert our .vcf file into a format servicable for NeEstimator, namely *genpop*.
```bash
module load Python
./vcf2genepop.pl vcf=../HamFrogR08maxsnps1DP5.recode.vcf pops="^MT_,^(M_|G_),^B_,^T_" > HamFrog.gen
```
I'll also produce a dataset based on 4,000 random SNPS. I'll breakdown this code:
- {} is a command group to be redirected together to >
- grep '^#' is find lines which start with # which in a vcf file is the metadata
- grep -v '^#' then extracts only the SNP data -v is saying take the inverse of this
- | is pipe, so the SNP data is piped to shuffle and then head which take the first 4,000 lines
```bash
{ grep '^#' ../HamFrogR08maxsnps1DP5.recode.vcf; grep -v '^#' ../HamFrogR08maxsnps1DP5.recode.vcf | shuf | head -n 4000; } > Random4kSNPsforNe.vcf
./vcf2genepop.pl vcf=Random4kSNPsforNe.vcf pops="^MT_,^(M_|G_),^B_,^T_" > HamFrog4kSNPs.gen
```
I've then created the output paramter files on my laptop which Ludo has run on his computer ...... the Ne2-1L file is made for a 32-bit operating system (it hasn't been upgraded since 2008) making in incompatible with NesI without extra effort. Using the Java GUI on a simple computer is much easier in this case.

