# Pairwise F<sub>st</sub>
I'd like to calculate pairwise FF<sub>st</sub> across my four populations in this analysis. Although [VCFtools](https://vcftools.github.io/documentation.html) does allow you to do this, it is inefficient and I'd have to do so manually for each pair. Instead, I'll try to streamline this anlaysis in R with the package SNPRelate.

Let's first take our popmap.txt we've used for denovo_map.pl in Stacks, and create a new tab separted .txt file that *actually* contains the relevant 'pop' information

```bash
cd /home/mulha552/uoo04306/frogs_gbs
awk '{
  if ($1 ~ /^B/)   $2 = "Boat Bay";
  else if ($1 ~ /^T/)  $2 = "Stephens Island";
  else if ($1 ~ /^M_/) $2 = "Maud Island";
  else if ($1 ~ /^MT/) $2 = "Motuara";
  else if ($1 ~ /^G/)  $2 = "Maud Island";
  print $1, $2;
}' popmap.txt > meta.txt
```
