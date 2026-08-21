# Single Nucleotide Polymorphism Calling for Hamilton's frog *Leiopelma hamiltoni*
SNP calling with [Stacks](https://catchenlab.life.illinois.edu/stacks/) for a [Genotyping_by_Sequencing](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0019379) dataset of *Leiopelma hamiltoni*

## Quality Control and Adaptor Trimming
### Testing: QC
```bash
# Create a subset of 250,000 from the raw reads: 
zcat NS0229_S1_L004_R1_001.fastq.gz | head -n 1000000 > testR1.fastq 
zcat NS0229_S1_L004_R2_001.fastq.gz | head -n 1000000 > testR2.fastq

# Run QC
module load FastQC 
fastqc *.fastq
```

### Testing: adapter trimming
We need to remove quite a few adapter's from our sequences, especially at high read positions. So, based on the FastQC, we are trimming off adapters and remove reads shorter than 125p with Cutadapt.

Stacks requires all reads to be the same length. Cutting too high the would loose too many reads after the adaptors are removed, but if we made it too low then for those, in our case ~90%, of reads with no adaptors will be shortened too much loosing information (Remember Ilumina reads are only 150bp long).

Parameters were: 
- -a specifies the adaptor sequence for forward reads (R1), and -A for reverse (R2) reads;
- -q specifies trimming low quality bases below a QC score of 25 from the 3' end;
- -o (-p) specifies the output file for forward reads (and for reverse reads);
- --minimum length specifies that trimmed reads must be 125 bp long or should be discarded;
- --length shortens all reads to the required 125 length.
```bash 
module load cutadapt 
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC  -q 25 -o trimmed_testR1.fastq  --minimum-length 125:125 --length 125  -p  trimmed_testR2.fastq testR1.fastq  testR2.fastq

fastqc trimmed*fastq
```
To check the adapters were gone, we look in the FastQC reports of these trimmed test data. Specifically, we looked for "bases proceeding adaptors" as these should be random, if a single nucleotide were to be overrepresented it may indicate the entire adaptor wasn't removed. For the Hamilton's frog data, we do see there is a problematic 'per base sequence content.' However, this was caused by the presence of a limited the number of barcodes in the sequences, ignored for now.

### CutAdapt (an example SLURM Jobs on NeSi)
NeSi or the New Zealand eScience Infrastructure is the national HPC platform. Including an exemplar on submitting jobs etc 😊
```bash
# Create a new job file with nano
nano Frog_trim.sl
```
Specify arguments with #SBATCH, and then our code:
```bash 
#!/bin/bash -e
#SBATCH --job-name=Frogtrim # Job name (shows up in the queue)
#SBATCH --time=48:00:00      # Walltime (HH:MM:SS), if our job finishes before this no worries, give ample time in case
#SBATCH --mem=8G          # Memory in G, minimum 1G per core.
#SBATCH --cpus-per-task=8 # Number of cores for our job

cd /home/mulha552/uoo04306/frogs_gbs/source_files
module load cutadapt
cutadapt -j 8 -a AGATCGGAAGAGC -A AGATCGGAAGAGC  -q 25  -o trimmed_NS0229_Hamiltons_S1_R1_001.fastq  --minimum-length 125:125  --length 125  -p trimmed_NS0229_Hamiltons_S1_R2_001.fastq  NS0229_S1_L004_R1_001.fastq.gz  NS0229_S1_L004_R2_001.fastq.gz
cd ..
```
Useful SLURM commands:
```bash
# In an interactive session, we can check code works (i.e. no syntax errors) before we submit:
bash Frog_trim.sl

# Use Ctrl+C to abort the local processing; and Submit the job:
sbatch Frogtrim.sl

# Check the queue:  
squeue -u mulha552 

# or cancel the job:
scancel Frogtrim

# Inspect .log file, with less, tail, or cat:
cat Frogtrim.log
```
Run trimming utilising the exemplar SLURM job, and check with a final FastQC.
```bash 
fastqc trimmed_NS0229_Hamiltons_S1_R1_001.fastq trimmed_NS0229_Hamiltons_S1_R2_001.fastq
```

## Demultiplexing
The file Barcode.txt contains all three plates sequenced at AgResearch on a Single NovaSeq lane with the unique PstI barcodes identifying  wells on each plate. I changed the names of the wells with Hamilton's frogs to give an informative name on population / individual ID.

Demultiplexing is run as a SLURM job. Stacks has great documentation, including useful information on [process_radtags](https://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php).
```bash 
cd /home/mulha552/uoo04306/frogs_gbs/source_files/raw
module load Stacks #2.61
process_radtags -P   -p ../raw/ -o ../samples/ -b ../barcodes.txt -e pstI -r -c  --inline-inline
```
Results:
 ```
3666129448  total sequences  
489768784 barcode not found drops (13.4%)
889709 low quality read drops (0.0%)
0 poly-G run drops (0.0%)
43575443 RAD cutsite not found drops (1.2%)
3131895512 retained reads (85.4%)
```

## Concatenate Reads
We want to produce a single file per sample, excluding the skinks, snails, tuatara etc. sequenced on this same sequencing run. Popmap.txt has a first column as the frogIDs taken from Barcode.txt, with the second column, in our case, just all filled with 'Pop.'See the [information](http://catchenlab.life.illinois.edu/stacks/manual/#popmap) available in the Stacks documentation.

### Python (an example SLURM Job)
Run concatenation as a SLURM job, but as a Python Environment.
```bash
# Create a new job file with nano
nano nano frogconcat.py
```
We use a python shebang for this file, and cannot specify requirements using #SBATCH as with bash or shell. Popmap.txt needs to be in the current working directory.
```python
#!/usr/bin/env python3

import os
os.chdir('/home/mulha552/uoo04306/frogs_gbs') #python for 'cd'
allfiles_to_concat = os.listdir("samples/")
alltokeep=[]
with open("popmap.txt") as f:
	for line in f:
		print(line)
		samplename=line.split("\t")[0]
		#print (samplename) # should be commented out in the slurm job, but I can use this as a check things are running.
		checkfiles=[filename for filename in allfiles_to_concat  if filename.startswith(samplename)]
		if len (checkfiles)!=4: # that was weirdly complicated because some sample name are contained in others different ways, but the vcheck above solve it uysing the rem file
			print(line)
			raise Exception
		else:
			os.system("zcat "+" ".join(["samples/"+checkfile for checkfile in checkfiles])+"|  gzip -c > samples_concat/"+samplename+".fq.gz" )

```
```bash
# In an interactive session, we can check code works (i.e. no syntax errors) before we submit:
module load Python
python3 frogconcat.py

# Submit the a job, with arguments specified in command line:
sbatch -A uoo04306 -t 48:00:00 -c 8 --mem 8G --job-name=frogconcat frogconcat.py
```

## Subsampling
```
# Check how many million reads we have per sample:
cd /home/mulha552/uoo04306/frogs_gbs//samples_concat
for fq in *.fq.gz; do
   sample=$(basename "$fq" .fq.gz)
    reads=$(zgrep -c "^@" "$fq")
     reads_mill=$(awk -v r="$reads" 'BEGIN {printf "%.2f", r/1000000}')
     echo -e "${sample}\t${reads}\t${reads_mill}M"
done > ../read_counts_all_samples.txt
```

Because our dataset has such a huge number of reads, we subsample to take the first 5 millions reads (i.e. 20 million lines). We need a workable amount of sequencing for our computing capacity. For any samples where we may have < 5million reads, we analyse to take all of the data. 

Run as a python SLURM job
```python
import os
os.chdir('/home/mulha552/uoo04306/frogs_gbs') #python for 'cd'
with open("popmap.txt") as f:
	for line in f:
		sample  =line.split("\t")[0]
		print(sample)
		os.system("zcat samples_concat/"+sample+".fq.gz |  head -n 20000000  | gzip -c > samples_subsampled/"+sample+".fq.gz ") 
```

## Parameters optimisation
Parameter optimisation is about deciding with which stacks parameters (-M, -m, -o, -T etc.) we continue to use when running on the entire dataset. You can read a relatively clear tutorial [here](http://catchenlab.life.illinois.edu/stacks/param_tut.php), which describes the parameters dictating stack and SNP calling de-Novo. We roughly followed optimisations methods described in this [paper](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12775).
```bash
# Testing on 30 sample subset:
shuf ../popmap.txt | head -n 30 > popmap_opti.txt
```
This loop automates the creation of directories, scripts, and submission of a SLURM job for processing the specified random 30 samples using the denovo_map.pl commando; each job is submitted with different parameters, based on the loop variable i:
```bash
for i in 2 3 4 5 6 7 8
do mkdir -p M$i; echo '#!/bin/sh' > runM$i.sh
echo cd /home/mulha552/uoo04306/frogs_gbs/source_files/para_opti >> runM$i.sh
echo "module load Stacks/2.61-gimkl-2022a" >> runM$i.sh
echo "denovo_map.pl --samples ../samples_subsampled/ --popmap popmap_opti.txt  -o M$i  -M $i -n $i -m 3 -T 8" >> runM$i.sh
sbatch -A uoo04306 -t 2-00:00:00 -J M$i -c 8 --mem=64G runM$i.sh
done
```
We take the SNP data from denovo_map.pl and compute a variety of useful statistics using [populations](http://catchenlab.life.illinois.edu/stacks/comp/populations.php). We are interest in the number of polymorphic loci at -R 0.8 i.e. loci covered in 80% of individuals. 
```bash
module load Stacks/2.61-gimkl-2022a

for i in 2 3 4 5 6 7 8 
do
populations -P M$i -R 0.8 -M popmap_opti.txt --vcf
done

# See how many loci are left by looking in the log of populations;
# cut -f 1, is looking only a column one which lists fragment ID, and uniq means we are only counting unique fragment IDs:

for i in 2 3 4 5 6 7 8
do
echo $i
cat M$i/populations.snps.vcf | grep -v '#' | cut -f 1 | uniq | wc -l
done
```
Hamilton's frog data have the most polymorphic loci with M= 2 and n = 2 ... these parameters will be run on the entire dataset!


## Final SNP Calling
```bash
#!/bin/bash
cd /home/mulha552/uoo04306/frogs_gbs
module load Stacks/2.61-gimkl-2022a
denovo_map.pl --samples samples_subsampled/ --popmap popmap.txt  -o M2_final  -M 2 -n 2 -m 3 -T 8
```
SLURM job, given the amount of data:
```bash
sbatch -A uoo04306 -t 5-00:00:00 -J M2Final -c 8 --mem=300G -p hugemem M2Final.sl
```

## Population Statistics
Run populations at a lower threshold (-R 0.5) allowing us to check for low quality samples.
```bash
# You could reduce the -R value here if you had a poor read depth dataset:
populations -P M2_final -M popmap.txt  --vcf --structure --plink --treemix --max-obs-het 0.65 -R 0.8  -O M2_Final

# View and sort individual frogs by missingness:
module load VCFtools 
vcftools --vcf M2_Final/populations.snps.vcf --missing-indv
sort -k 4n out.imiss

# Only MT_24_FAU has more than 60% missing data. Remove from popmap.txt, to create popmap_clean.txt:
grep -v "^MT_24_FAU\\s" popmap.txt > popmap_clean.txt
```
Run populations, but filtering positions found in less than 80%, keeping a maximum of one SNP per locus, and a maximum observed heterozygosity to 0.65. 

We keep one SNP per locus because our analysis assume that SNPs are independent, and those on the same loci may not. Also because it can reduce the impact of erroneous loci made up of repetitive regions, although not a primary reason for this filtering in our method. 

Maximum observed heterozygosity is set to 0.65, because we'd never reasonably expect heterozygosity to be > 0.5 in a natural populations. Were Hamilton's frog to have a duplication in its genome then these may assemble together with all individuals called at Heterozygotes at a mutation site.
```bash
populations -P M2_final/ -M popmap_clean.txt  --vcf --structure --plink --treemix --max-obs-het 0.65 -R 0.8  --write-single-snp -O M2_final
```
Results:
```
Removed 4068568 loci that did not pass sample/population constraints from 4090567 loci.
Kept 21999 loci, composed of 2605505 sites; 37125 of those sites were filtered, 21250 variant sites remained.
Mean genotyped sites per locus: 118.44bp (stderr 0.01).

Population summary statistics (more detail in populations.sumstats_summary.tsv):
  pop: 70.591 samples per locus; pi: 0.028234; all/variant/polymorphic sites: 2605505/21250/21250; private alleles: 0
Populations is done.
```
```bash
# Save with a meaningful name:
cp M2_final/populations.snps.vcf HamFrogR08maxsnps1.vcf
```
### Coverage
If we have a 'true' heterozygote, but at low coverage, say 2x, then we have a 50% chance of incorrectly calling a homozygote in our data. Potentially we could see a relationship between coverage depth and heterozygosity, at low depth. We want to ensure to filter SNPs so that we loose this relationship. 
```bash
for i in 2 3 4 5 6
do
vcftools --vcf ../HamFrogR08maxsnps1.vcf --minDP $i --recode --out filtered_Depth$i
vcftools --vcf filtered_Depth$i.recode.vcf --het --out het_Depth$i
vcftools --vcf filtered_Depth$i.recode.vcf --depth --out depth$i
echo "Complete for minDP = $i"
done
```
I ran [Troubleshooting](./TroubleShooting/) in R, and even with a minimum coverage depth of two there is no relationship at all. The results suggest we have a high coverage dataset, but we will still utilise data cut of at a minimum depth of five for analyses.
