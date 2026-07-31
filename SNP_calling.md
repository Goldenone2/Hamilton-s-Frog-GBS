# Single Nucleotide Polymorphism Calling 
Create a folder for the (raw data) source files my source files.
```bash 
mkdir source_files
```
Take a subset of 1'00000 lines i.e. 250'000 reads. Use this subset to trial code on the first run through, and check there are no issues with my sequences before submitting a full SLURM job. 
```bash 
cd source_files
zcat NS0229_S1_L004_R1_001.fastq.gz | head -n 1000000 > testR1.fastq 
zcat NS0229_S1_L004_R2_001.fastq.gz | head -n 1000000 > testR2.fastq
```
## Quality Control and Adaptor Trimming

### Testing: QC
```bash 
module load FastQC 
fastqc *.fastq
```
Looking at the generated testR1.html and testR2.html, there is quite a few adapter's in there, especially in high read positions.

### Testing: adapter trimming
Trim off adapters and remove reads shorter than 125p with Cutadapt; 125bp was chosen based on the our Fastq quality check. 

Stacks requires al reads to be the same length. Cutting too high the will loose too many reads after the adaptors are removed, but if we make it too low then for those, in our case ~90%, of reads with no adaptors will be shortened too much loosing information (Remember Ilumina reads are only 150bp long).

Parameters are: -a specifies the adaptor sequence for forward reads (R1), and -A for reverse (R2) reads, -q specifies trimming low quality bases below a QC score of 25 from the 3' end, -o (-p) specifies the output file for forward reads (and for reverse reads), minimum length specifies that reads must be 125 bp long or should be discarded, length shortens all reads to the required 125 length.
```bash 
cd source_files
module load cutadapt # 
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC  -q 25 -o trimmed_testR1.fastq  --minimum-length 125:125 --length 125  -p  trimmed_testR2.fastq testR1.fastq  testR2.fastq
cd ..
```
Check the adapters are gone, we need to make this decision based on the output from the Cutadapt testing. Specifically, looks for "bases proceeding adaptors" these should be random, if a single nucleotide were to be overrepresented it may indicate the entire adaptor wasn't removed.
```bash 
fastqc trimmed*fastq
```
Checking the .html files to see whether that worked. For the Hamilton's frog data, we do see there is a problematic 'per base sequence content.' However this is caused by the presence of a limited the number of barcodes in the sequences, ignore for now.

### SLURM Job
Run that on all reads. Note: as my first analysis on NeSi (New Zealand's national HPC), I will include some superfluous detail on submitting jobs etc. 😊

Create a new job file:
```bash 
nano Frogtrim.sl
```
Specify arguments with #SBATCH, and afterwards include code:
```bash 
#!/bin/bash -e
#SBATCH --job-name=Frogtrim # job name (shows up in the queue)
#SBATCH --time=48:00:00      # Walltime (HH:MM:SS), if our job finishes before this no worries, give ample time in case
#SBATCH --mem=8G          # Memory in G, minimum 1G per core.
#SBATCH --cpus-per-task=8 #number of cores for our job

cd /home/mulha552/uoo04306/frogs_gbs/source_files
module load cutadapt
cutadapt -j 8 -a AGATCGGAAGAGC -A AGATCGGAAGAGC  -q 25  -o trimmed_NS0229_Hamiltons_S1_R1_001.fastq  --minimum-length 125:125  --length 125  -p trimmed_NS0229_Hamiltons_S1_R2_001.fastq  NS0229_S1_L004_R1_001.fastq.gz  NS0229_S1_L004_R2_001.fastq.gz
cd ..
```
Check code works (i.e. no syntax errors) before we submit:
```bash 
sh Frogtrim.sl
```
Kill the local processing with Ctrl+C, and Submit the job:
```bash 
sbatch Frogtrim.sl
```
Check the queue: We can also cancel by ' scancel '<jobid>" '. Use squeue to get the ID. 
```bash 
squeue -u mulha552 # -u specifies just my user's jobs
```
or cancel the job:
```bash
scancel frogtrim
```
SLURM jobs produce a .log file, inspect the log file with less, cat or nano to verify code has run correctly.

Once trimmed, make a final Fastq check. 
```bash 
fastqc trimmed_NS0229_Hamiltons_S1_R1_001.fastq trimmed_NS0229_Hamiltons_S1_R2_001.fastq
```

## Demultiplexing
Barcode.txt contains all three plates sequenced at Ag Research on a Single NovaSeq lane with the unique PstI barcodes which identify the wells on each plate. I changed the names of the wells with frogs to give an informative name on population and other metadata etc.

Make directory:
```bash 
cd .. #assuming I'm in source files
mkdir raw
```
Copy files:
```bash 
cd raw
cp /home/mulha552/uoo04306/frogs_gbs/source_files/trimmed_NS0229_Hamiltons_S1_R1_001.fastq /home/mulha552/uoo04306/frogs_gbs/source_files/trimmed_NS0229_Hamiltons_S1_R2_001.fastq /home/mulha552/uoo04306/frogs_gbs/raw
cd ..
```
Run demultiplexing as a SLURM job. Useful info can be found [here](https://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php):
```bash 
cd /home/mulha552/uoo04306/frogs_gbs/source_files/raw
module load Stacks #2.61
process_radtags -P   -p ../raw/ -o ../samples/ -b ../barcodes.txt -e pstI -r -c  --inline-inline # NO -q often used for process-radtags gives me an error because of it, but no worries, cutadapatalready took care of this
 ```
Results:
 ```
3666129448 total sequences
 489768784 barcode not found drops (13.4%)
    889709 low quality read drops (0.0%)
         0 poly-G run drops (0.0%)
  43575443 RAD cutsite not found drops (1.2%)
3131895512 retained reads (85.4%)
```

## Concatenate Reads
To produce one file per sample, excluding the skinks, snails, tuatara etc. sequenced on this same sequencing run

Popmap.txt has a first column as the frogIDs taken from Barcode.txt, with the second column just all filled with 'Pop';[Info](http://catchenlab.life.illinois.edu/stacks/manual/#popmap).

Make directory:
```bash 
mkdir samples_concat
```
### Python SLURM
Run concatenation as a SLURM job. Because this is within a Python Environment, submission looks different. 

Create a new job file:
```bash 
nano frogconcat.py #create python script
```
We use a python shebang for this file (rather than a shell/bash one) and we cannot specify requirements using #SBATCH, this is now done later. 
- Copy the code below into frogconcat.py.
- popmap.txt needs to be in the current working directory.

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
Check code works (i.e. no syntax errors) before we submit:
```bash
module load Python #need to load a python environment//ipython console
python3 frogconcat.py
```
Submit a job. For python specify arguments here in the command line: -A specifies the account to 'charge' the job, -t specifies time, -c specifies the number of cpus per task, --mem specifies memory, --job-name specifies a job name
```bash
sbatch -A uoo04306 -t 48:00:00 -c 8 --mem 8G --job-name=frogconcat frogconcat.py
squeue -u mulha552 #check its working
```
Check of how many million reads we have per sample:
```
cd /home/mulha552/uoo04306/frogs_gbs//samples_concat
for fq in *.fq.gz; do
   sample=$(basename "$fq" .fq.gz)
    reads=$(zgrep -c "^@" "$fq")
     reads_mill=$(awk -v r="$reads" 'BEGIN {printf "%.2f", r/1000000}')
     echo -e "${sample}\t${reads}\t${reads_mill}M"
done > ../read_counts_all_samples.txt
```
## Subsampling
*Hamilton's frog have a large genome so I have special code*
We will subsample from the concatenated reads to take the first 5 million reads (code reads 20,000,000 lines because each read is 4 lines in a .fastq) so we are working with an actually workable amount of sequence. We still run the previous code, because for those samples where we may have <5million reads we need to take all of it.

 We're also running as a pyhton slurm job, as above...
```bash
mkdir samples_subsampled
```
```python
import os
os.chdir('/home/mulha552/uoo04306/frogs_gbs') #python for 'cd'
with open("popmap.txt") as f:
	for line in f:
		sample  =line.split("\t")[0]
		print(sample)
		os.system("zcat samples_concat/"+sample+".fq.gz |  head -n 20000000  | gzip -c > samples_subsampled/"+sample+".fq.gz ") 
```
Then, run denovo_map with samples_concat as input.

## Parameters optimisation
Ok, so, paramater optimisation is about deciding with stacks paramters (-M, -m, -o, -T etc) we will continue to use when running on the entire dataset. You can read a relatively clear tutorial [here](http://catchenlab.life.illinois.edu/stacks/param_tut.php), which describes how stacks without a reference genome works and, what these paramters dictate. Not included here are -o, output directory and, -t, number of threads.

At this stage, I'll make a popmap and exclude all samples with less than 100'000 reads (combining forward and reverse, i.e. 50k retained reads). First check number of reads (remember .fq files have 4 lines for each read), then create popmap.

```bash
cd samples_subsampled #assuming in frogs_gbs
output_file="read_counts.txt"
> "$output_file" # clears previous contents
for fq_file in *.fq.gz; do # Loop through all fq.gz files in the current directory
lines=$(zcat "$fq_file" | wc -l)  # retrieve number of lines in the file
reads=$((lines / 4)) # Calculate the number of reads
echo "$fq_file: $reads reads" >> "$output_file"  # echo the file name and number of reads to the output file
done

# read_counts.txt #optionally delete text file, I'm unsure whether this inerfere's with future script...
```

Great, it looks like we don't need to remove any frogs from our analysis!

## Parameter optimisation

Run on a random 30 samples.
```bash
#mkdir para_opti
cd para_opti
shuf ../popmap.txt | head -n 30 > popmap_opti.txt # not done if no optimisation
```
This code is a "for loop." This loop automates the creation of directories, scripts, and slurm job for processing the specified random 30 samples using the denovo_map.pl tool. Each job is submitted with different parameters (based on the loop variable i), and each job will run with different configurations of parameters (see where variable i is included in the code). 

Note from future Hadley, Ludo has used echo which doesn't execute the code rather repeats it in the terminal; Ludo's told it to 'echo' in the slurm (.sh) files. This technique is different from something you could just do in the terminal....

I believe this method of paramater optimisation, roughly, follow this [paper](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12775).

```bash
for i in 2 3 4 5 6 7 8
do mkdir -p M$i; echo '#!/bin/sh' > runM$i.sh
echo cd /home/mulha552/uoo04306/frogs_gbs/source_files/para_opti
echo "module load Stacks/2.61-gimkl-2022a" >> runM$i.sh
echo "denovo_map.pl --samples ../samples_subsampled/ --popmap popmap_opti.txt  -o M$i  -M $i -n $i -m 3 -T 8" >> runM$i.sh
sbatch -A uoo04306 -t 2-00:00:00 -J M$i -c 8 --mem=64G runM$i.sh
done
```
```bash
squeue -u mulha552
```
The *populations* command takes the the SNP data from *denovo_map.pl* and, can compute a variety of useful population genetics statistics, See [here](http://catchenlab.life.illinois.edu/stacks/comp/populations.php) for more information. 
Focus on the number of loci at -R 0.8 (loci covered in 80% of inds):
```bash
module load Stacks/2.61-gimkl-2022a
```
```bash
for i in 2 3 4 5 6 7 8 
do
populations -P M$i -R 0.8 -M popmap_opti.txt --vcf
done
```
See how many loci are left by looking in the log of populations:

cut -f 1, is looking only a column one which lists fragment ID, and uniq means we are only counting unique fragment IDs....
```bash
for i in 2 3 4 5 6 7 8
do
echo $i
cat M$i/populations.snps.vcf | grep -v '#' | cut -f 1 | uniq | wc -l
done
```
We have the most polymoprhic loci with M=n of 2, so that is now what we'll run for the enitre study!!


## Run on the whole dataset
I will run this as SLURM job, we have a lot of data...
```bash
mkdir M2_final
```
```bash
#!/bin/bash
cd /home/mulha552/uoo04306/frogs_gbs
module load Stacks/2.61-gimkl-2022a
denovo_map.pl --samples samples_subsampled/ --popmap popmap.txt  -o M2_final  -M 2 -n 2 -m 3 -T 8
```
```bash
sbatch -A uoo04306 -t 5-00:00:00 -J M2Final -c 8 --mem=300G -p hugemem M2Final.sl
```
First lets run populations with a lower threshold (-R 0.5) aallowing us to check for low quality samples. 
```bash
populations -P M2_final -M popmap.txt  --vcf --structure --plink --treemix --max-obs-het 0.65 -R 0.8  -O M2_Final
```
You can reduce the -R value here if you had a poor read depth dataset. *For your info sort is a shell command -k is saying to sort in order by column 4, so it shows everything sorted by the amount of missing data*
```bash
module load VCFtools 
vcftools --vcf M2_Final/populations.snps.vcf --missing-indv
sort -k 4n out.imiss
```
The following sample has more than 60% missing data: MT_24_FAU
```
Let's remove 'FAU' from popmap.txt and create popmap_clean.txt by using -v which prints all the lines that *do not* match waht I've specified. ^ ensures our match is at the beginning of a line and \\s account for any whitespace; Grep is a tool used to search and manipulate text within files.
```bash
grep -v "^MT_24_FAU\\s" popmap.txt > popmap_clean.txt
```
Now I'll re-run populations, but filtering positions found in less than 80% and keeping a maximum of one SNP per locus. I keep one SNP per locus because ou anlaysis assume that SNPs are independant, and thsoe on the same loci are not (also because it can reduce the impact of erroneous loci made up of repetitive regions....but not a reason for my methods). I've also set the max-obs-het to 0.65, becasue we'd never expect het to be >0.5. If we have a duplication in the Hamilton's frog genome then these may assemble together with all individuals called at Heterozygotes at a mutation site, max-obvs-het allows us to filter this out. 
```bash
populations -P M2_final/ -M popmap_clean.txt  --vcf --structure --plink --treemix --max-obs-het 0.65 -R 0.8  --write-single-snp -O M2_final
```
My results:
```
Removed 4068568 loci that did not pass sample/population constraints from 4090567 loci.
Kept 21999 loci, composed of 2605505 sites; 37125 of those sites were filtered, 21250 variant sites remained.
Mean genotyped sites per locus: 118.44bp (stderr 0.01).

Population summary statistics (more detail in populations.sumstats_summary.tsv):
  pop: 70.591 samples per locus; pi: 0.028234; all/variant/polymorphic sites: 2605505/21250/21250; private alleles: 0
Populations is done.
```
Only SNPs found in 80% of individuals are kept, so 21,250 variable sites! 

We also have indv missigness (using code above) <0.23 which is great give we subsampled, reaffirms our choice.

Save it with a meaningful name:

```bash
 cp M2_final/populations.snps.vcf HamFrogR08maxsnps1.vcf
```
## Coverage
I must filter my dataset by coverage; remeber, if we have a 'true' heterozygote, but at low coverage (say two), then we have a 50% chance of incorrectly calling a homozygote in our data.... there will be an inherent relationship between coverage depth and heterozygosity (at low depth). I want to filter my SNPs so that we loose this relationship. 

I've madea bsic a 'for' loop. I filter for minimum depth using --minDP from values of 2 to 6, --recode generate a new .vcf files for each filter. Then I calculte heterozygosity  using --het, and depth using --depth.
```bash
module load VCFtools
mkdir filtered_coverage
```
```bash
for i in 2 3 4 5 6
do
vcftools --vcf ../HamFrogR08maxsnps1.vcf --minDP $i --recode --out filtered_Depth$i
vcftools --vcf filtered_Depth$i.recode.vcf --het --out het_Depth$i
vcftools --vcf filtered_Depth$i.recode.vcf --depth --out depth$i
echo "Complete for minDP = $i"
done
```
Right, so, if we have a lack at the Coverage Plots we can see that even with a minimum coverage depth of two there is no relationship at all! Overall, the results usggest we have a high coverage dataset. I'll still make the choice to cut of at a miniumum depth of five; although we don't see any relationship this is still a good choice and, given our dataset we will not use much data!

### Re-run populations
For some analyses, specifcally for Stairway Plot, we need to filter using *populaitons* not *VCFTools.* When we filter with VCFTools we may remove the variable individuals for some of the called SNPs making them monomorphic, which can be problematic. Addtionally, stairway plot 2 requires accurate information of the number of *sites* both monomoprhic and polymorphics; it's easiest to get this infromation directly from population's output. 

**in the final published version of this code this should happen in the initial instance afterwhich we'd still check for any addtional realtionship (which are more likely in low coverage datasets).**

```bash
mkdir M2_FInal_minGT5
populations -P M2_final/ -M popmap_clean.txt  --vcf --structure --plink --treemix --max-obs-het 0.65 -R 0.8  --write-single-snp --min-gt-depth 5 -O M2_FInal_minGT5
```
My results:
```
Removed 4068568 loci that did not pass sample/population constraints from 4090567 loci.
Kept 21999 loci, composed of 2605505 sites; 280618 of those sites were filtered, 11611 variant sites remained.
Filtered 912593 genotypes that fell below the minimum genotype depth threshold.
Mean genotyped sites per locus: 111.27bp (stderr 0.04).

Population summary statistics (more detail in populations.sumstats_summary.tsv):
  pop: 70.486 samples per locus; pi: 0.029416; all/variant/polymorphic sites: 2447914/11611/11611; private alleles: 0
Populations is done.
```
!!!!! Ineeed to remove this section, actually this will not impact Stairway Plot except in the negative sense because have so so many less loci .......
## Final Datasets
```bash
cd /home/mulha552/uoo04306/frogs_gbs
cp M2_FInal_minGT5/populations.snps.vcf HamFrogR08maxsnps1minGT5.vcf
module load VCFtools
vcftools --vcf HamFrogR08maxsnps1minGT5.vcf --het --out Heterozygosity
vcftools --vcf HamFrogR08maxsnps1minGT5.vcf --depth --out Depth
```
