# Per-population filtering
We applied a heterozygosity filter of 0.65 to our dataset globally using 'populations.' Because we have uneven sample sizes between Takapourewa and Te Pakeka + associated translocations, it's possible to have missed some Ho > 0.65 sites within the Takapourewa population. Quickly, I will verify our results aren't a result of erroneous filtering. An additional analysis suggested after peer-review, and run on the University of Otago's Aoraki HPC environment.

I will split our original vcf file by population groups, and calculate the per-SNP observed heterozygosity expressed as a proportion (rather than number of individuals). Code based on [this BioStars Post](https://www.biostars.org/p/291147/#291167). 

## Set up
Load a vcf/bcftools conda environment:
```bash
source /projects/sciences/zoology/rawlence_lab/conda/etc/profile.d/conda.sh
conda activate bcftools_hadley 
```
Create separate .bcf files for Takapourewa, and Te Pakeka + associated translocations:
```bash
grep -v '^#' HamFrogR08maxsnps1DP5.recode.vcf | cut -f1 | sort -u > contigs.txt
awk '{print "##contig=<ID="$1">"}' contigs.txt > contig_header.txt

bcftools annotate -h contig_header.txt -Ov -o HamFrogR08maxsnps1DP5.fixed.vcf HamFrogR08maxsnps1DP5.recode.vcf

bcftools view -S Takapourewa.txt -Ob -o Takapourewa_Only.bcf HamFrogR08maxsnps1DP5.fixed.vcf
bcftools view -S ^Takapourewa.txt -Ob -o Te_Pakeka_Only.bcf HamFrogR08maxsnps1DP5.fixed.vcf
```
Calculate per SNP heterozygosity, switch input / output as needed:
```bash
paste <(bcftools view Te_Pakeka_Only.bcf \
    | awk -F"\t" 'BEGIN {print "CHR\tPOS\tID\tREF\tALT"} !/^#/ \
    {print $1"\t"$2"\t"$3"\t"$4"\t"$5}') \
  <(bcftools query -f '[\t%GT]\n' Te_Pakeka_Only.bcf | \
    awk 'BEGIN {print "nHet\tnCalled\tHo"} {
        nHet=0
        nCalled=0
        for (i=1; i<=NF; i++) {
            if ($i != "./.") {
                nCalled++
                if ($i ~ /^0\/1$/ || $i ~ /^1\/0$/)
                    nHet++}}
        if (nCalled > 0)
            print nHet "\t" nCalled "\t" nHet/nCalled
        else
            print "0\t0\tNA"
    }') > Te_Pakeka_Per_SNP_Het.txt
```
Remember total SNPs is 21250; count the number of SNPs > 0.65. 
```bash
awk 'NR > 1 && $8 > 0.65 {print $1 "\t" $2}' Takapourewa_Per_SNP_Het.txt > Takapourewa_highHet.txt
wc -l Takapourewa_highHet.txt
# Result: 256

awk 'NR > 1 && $8 > 0.65 {print $1 "\t" $2}' Te_Pakeka_Per_SNP_Het.txt > Te_Pakeka_highHet.txt
wc -l Te_Pakeka_highHet.txt
# Result: 26

cat Takapourewa_highHet.txt Te_Pakeka_highHet.txt | sort -k1,1 -k2,2n -u > all_highHet.txt
wc -l all_highHet.txt
#Result: 282
```

## Filtering
```bash
bcftools view -T ^all_highHet.txt -Ov -o HamFrogR08maxsnps1DP5.Under065.vcf HamFrogR08maxsnps1DP5.fixed.vcf
vcftools --vcf HamFrogR08maxsnps1DP5.Under065.vcf --het --out het_Under065
```

