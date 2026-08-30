### admixture

based on: https://speciationgenomics.github.io/ADMIXTURE/

```bash
module load stack/.2024-06-silent  gcc/12.2.0        plink/1.9-beta6.27
  conda activate admixtureconda install bioconda::admixture
mkdir ADMIXTURE
cd ADMIXTURE
```

File conversion
```bash
#!/bin/bash
#SBATCH --job-name=filter_plink
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=02:00:00
#SBATCH --output=logs/filter_plink_%j.log
set -euo pipefail

mkdir -p logs

source /cluster/project/gdc/shared/stack/GDCstack.sh
module load vcftools
module load bcftools
module load plink

VCF=freebayes_calling/Cyanistes_caeruleus_NOoutgroup_rm_call_minQ30_SNP_minDP3_maxmiss_06_azure_autosomes_filt_renamed_norm.vcf.gz


echo "filtering singletons and invariant sites..."
vcftools --gzvcf ${VCF} \
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

Note that the run of admixture below is an array that will run jobs for each value of K between 2 and 12 (see line 57 '''#SBATCH --array=2-12''')
'''
```
#!/bin/bash
#SBATCH --job-name=admixture
#SBATCH --array=2-12
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=6G
#SBATCH --time=24:00:00
#SBATCH --output=logs/admixture_%a_%j.log
source ~/miniconda3/etc/profile.d/conda.sh

conda activate admixture

FILE=filtered_renamed.bed

admixture --cv -j${SLURM_CPUS_PER_TASK} $FILE ${SLURM_ARRAY_TASK_ID} > log${SLURM_ARRAY_TASK_ID}.out
FILE=filtered_renamed.bed

admixture --cv -j${SLURM_CPUS_PER_TASK} $FILE.bed ${SLURM_ARRAY_TASK_ID} > log${SLURM_ARRAY_TASK_ID}.out
```
