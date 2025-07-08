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
I've then creating my output paramter files using the GUI on my laptop using only LinKage Disequilibium, these info and option files can then be run from the command line at a directory containing "info" "option" and "Ne2-1L."
```bash
module load Java
# chmod +x Ne2-1L; commented out but marks the file as executable.
./Ne2-1L i:info o:option
```
