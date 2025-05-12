# Pairwise Fst between populaitons 
First I'd like to caclulate pairwise Fst between populations, to do this using vcf tools I need to sequentially run --weri-fst-pop for each pair of my four populations. I must specify the populations by creating .txt files listing the relevant individuals. Later, I can visualise these result as a matrix in R. 

```bash
vcftools --vcf input_data.vcf --weir-fst-pop population_1.txt --weir-fst-pop population_2.txt --out pop1_vs_pop2

```
