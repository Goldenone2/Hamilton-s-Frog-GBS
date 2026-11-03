### Admixture
Based on: https://speciationgenomics.github.io/ADMIXTURE/
Run on Aoraki Otago University's HPC for Supplementary Information based on reviewer suggestions :)

```bash
## Already contains necessary packages to have re-run all other analyses with additional filtering:
source /projects/sciences/zoology/rawlence_lab/conda/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1
conda activate bcftools_hadley

conda install -c conda-forge -c bioconda admixture

mkdir admixture
cd admixture
```

## File conversion
```bash
#!/bin/bash
#SBATCH --job-name=Filter_plink
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=02:00:00
#SBATCH --output=logs/Filter_plink_%j.log
#SBATCH --chdir=
##----------Set Up-------------##
set -euo pipefail
mkdir -p logs

source /projects/sciences/zoology/rawlence_lab/conda/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1
conda activate bcftools_hadley


##----------Code-------------##
awk '{if($0 !~ /^#/) print "contig"$0; else print $0}' ../HamFrogR08maxsnps1DP5.Under065.vcf > PLINKvcf_with_contig.vcf

echo "filtering singletons and invariant sites..."
vcftools --vcf PLINKvcf_with_contig.vcf \
	--mac 2 \
	--min-alleles 2 --max-alleles 2 \
	--recode --recode-INFO-all \
	--stdout | bgzip -c > filtered_renamed.vcf.gz

bcftools index -t filtered_renamed.vcf.gz

echo "converting to plink bed/bim/fam..."
plink --vcf filtered_renamed.vcf.gz --make-bed --out filtered_renamed --allow-extra-chr --double-id

echo "zeroing out chromosome codes in bim..."
awk '{$1="0";print $0}' filtered_renamed.bim > filtered_renamed.bim.tmp
mv filtered_renamed.bim.tmp filtered_renamed.bim

echo "done"
```

Note: the run of admixture below is an array that will run jobs for each value of K between 2 and 12 (see line 57 '''#SBATCH --array=2-12''')
```
#!/bin/bash
#SBATCH --job-name=Frog_admixture
#SBATCH --array=2-12
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=6G
#SBATCH --time=24:00:00
#SBATCH --output=logs/admixture_%a_%j.log
##----------Set Up-------------##
set -euo pipefail

source /projects/sciences/zoology/rawlence_lab/conda/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1
conda activate bcftools_hadley

FILE=filtered_renamed.bed

admixture --cv -j${SLURM_CPUS_PER_TASK} $FILE ${SLURM_ARRAY_TASK_ID} > log${SLURM_ARRAY_TASK_ID}.out
```
