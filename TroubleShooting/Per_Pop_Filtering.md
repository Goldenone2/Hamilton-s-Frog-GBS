# Per-population filtering
We applied a heterozygosity filter of 0.65 to our dataset globally using 'populations.' Because we have uneven sample sizes between Takapourewa and Maud Island + associated translocations, it's possible to have missed some Ho > 0.65 sites within the Takapourewa population. Quickly, I will verify our results aren't a result of erroneous filtering. An additional analysis suggested after peer-review, and run on the University of Otago's Aoraki HPC environment.

```bash
# Load vcftools:
source /projects/sciences/zoology/rawlence_lab/conda/etc/profile.d/conda.sh
conda activate basic
```
