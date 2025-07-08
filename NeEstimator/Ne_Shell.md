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
