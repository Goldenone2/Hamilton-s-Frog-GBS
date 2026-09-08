# Per-population filtering
We applied a heterozygosity filter of 0.65 to our dataset globally using 'populations.' Because we have uneven sample sizes between Takapourewa and Te Pakeka + associated translocations, it's possible to have missed some Ho > 0.65 sites within the Takapourewa population. Quickly, I will verify our results aren't a result of erroneous filtering. An additional analysis suggested after peer-review, and run on the University of Otago's Aoraki HPC environment.

I will split our original vcf file by population groups, and calculate the per-SNP observed heterozygosity expressed as a proportion (rather than number of individuals). Code based on [this BioStars Post](https://www.biostars.org/p/291147/#291167). 

## Set up
```bash
# Load vcftools and bcftools:
source /projects/sciences/zoology/rawlence_lab/conda/etc/profile.d/conda.sh
conda activate bcftools_hadley 
```
Create separate bcf files for Takapourewa, and Te Pakeka + associated translocations.
```bash
grep -v '^#' HamFrogR08maxsnps1DP5.recode.vcf | cut -f1 | sort -u > contigs.txt
awk '{print "##contig=<ID="$1">"}' contigs.txt > contig_header.txt

vcftools --vcf HamFrogR08maxsnps1DP5.recode.vcf --remove Takapourewa.txt --recode-bcf --contigs contig_header.txt --out Te_Pakeka_Only
vcftools --vcf HamFrogR08maxsnps1DP5.recode.vcf --keep Takapourewa.txt --recode-bcf --contigs contig_header.txt --out Takapourewa_Only
```

## Filtering
Change input file out output name as needed:
```bash
paste <(bcftools view Takapourewa_Only.recode.bcf |\
    awk -F"\t" 'BEGIN {print "CHR\tPOS\tID\tREF\tALT"} \
      !/^#/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5}') \
    \
  <(bcftools query -f '[\t%GT]\n' 1000Genomes.Norm.bcf |\
    awk 'BEGIN {print "nHet\tnCalled\tHo"} \
      {nHet=0
        nCalled=0
        for (i=1; i<=NF; i++) {
          if ($i != "./.") {
            nCalled++
            if ($i ~ /^0\/1$/ || $i ~ /^1\/0$/)
              nHet++
          }
        }
        print nHet "\t" nCalled "\t" nHet/nCalled
      }') > Takapourewa_Per_SNP_Het.txt
``
