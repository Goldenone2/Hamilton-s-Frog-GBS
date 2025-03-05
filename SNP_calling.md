# SNP calling

First, I create a folder for the source files and folders to run Stacks, raw (which will contain the trimmed data) and samples. For my future reference, my source files are called NS0229_S1_L004_R1_001.fastq.gz & NS0229_S1_L004_R2_001.fastq.gz

```
mkdir source_files
```

Once the raw data is in source_files, I go in there and take a subset of 1'00000 lines, 250'000 reads. I will use the subset to trial the code on the firest run through, and check there are no issues with my sequences before I submit a full job (which runs for ages).

```
cd source_files
zcat NS0229_S1_L004_R1_001.fastq.gz | head -n 1000000 > testR1.fastq 
zcat NS0229_S1_L004_R2_001.fastq.gz | head -n 1000000 > testR2.fastq
```
Let's check quality with FastQC
```
module load FastQC 
fastqc *.fastq # The star means everything that ends with fastq.
```
I then look at the generated testR1.html and testR2.html after downloading them (I just had a look in the Jupyter borwser). There's quite a few adapter's in there, especially in high read positions.

## Adapter trimming

Trimming off adapters and removing reads shorter than 125p with cutadapt. We choose 125bp based on the our FastQC quality check, stacks requeires everything to be the same length: if we cut too high the we will loose too many loose reads after the adaptors are removed, but if we make it too low then for those (in our case ~90%) of reads with no adaptors will be shortened too much loosing information. Also, it's good for you to remeber Ilumina reads are only 150bp long.

First on the sample files. The files with trimmed are the output files. Use `cutadapt --help` after loading the module to learn what all the parameters are:

-a specifies the adaptor sequence for forward reads (R1), and -A for reverse (R2) reads, -q specifies trimming low quaility bases below a QC score of 25 from the 3' end, -o (-p) specifies the output file for forward reads (and for reverse reads), minimum length specifies that reads must be 125bp long or should be discarded, length shortens all reads to the requeired 125 length.

```
cd source_files
module load cutadapt # 
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC  -q 25 -o trimmed_testR1.fastq  --minimum-length 125:125 --length 125  -p  trimmed_testR2.fastq testR1.fastq  testR2.fastq
cd ..
```
Let's check  the adapters are gone, we need to make this descicion based on the output. Specifically, looks for "bases preceeeding adaptors" these should be random, if a single nucleotide was overrepresented it indicates we haven't removed the whole adaptor. 
```
fastqc trimmed*fastq # the star make it that it runs on anything that start with trimmed and end with fastq
```
Checking the .html files, that worked. We do see there is an issue with 'per base sequence content' but this is caused by the presence of a limited the number of barcodes in our sequence, so we can ignore this for now.

Now let's now run that on all reads, with a 125bp reads limit, so that we have one common length for all reads. (the only code difference is 'j' which specifies the number of cores).

Now, we will submit this as a job because I have too much data to sit around and analyse on a Jupyter Session.

_Submitting a Slurm Job_

Create a new file, and open it (nano is a text editor)
```
nano Frogtrim.sl
```
We then specificy arguments with #SBATCH, and afterwards include our code. Let's copy in the below text, then save and exit the text editor with 'ctrl + x'. You can use 'sbatch --help' to see what we can specify.

```
#!/bin/bash -e
#SBATCH --job-name=Frogtrim # job name (shows up in the queue)
#SBATCH --time=48:00:00      # Walltime (HH:MM:SS), if our job finishes before this no worries but we can give ample time in case
#SBATCH --mem=8G          # Memory in G, minimum 1G per core.
#SBATCH --cpus-per-task=8 #number of cores for our job

cd /home/mulha552/uoo04306/frogs_gbs/source_files
module load cutadapt
cutadapt -j 8 -a AGATCGGAAGAGC -A AGATCGGAAGAGC  -q 25  -o trimmed_NS0229_Hamiltons_S1_R1_001.fastq  --minimum-length 125:125  --length 125  -p trimmed_NS0229_Hamiltons_S1_R2_001.fastq  NS0229_S1_L004_R1_001.fastq.gz  NS0229_S1_L004_R2_001.fastq.gz
cd ..

```
Now, beofre we submit this we want to check this is working. It may still breka later but we want no errors in the code. 
```
sh Frogtrim.sl
```
Once we see cutadapt start working we can kill it with ctrl+c, and submit the job:
```
sbatch Frogtrim.sl
```
Lets check we are in the queue. We can also cancel by ' scancel '<jobid>" '. Use squeue to get the ID. 
```
squeue -u mulha552 # -u specifies just my user's jobs
```
Once the Job is complete use 'Cat' or 'Less' to check the output. 

Finally, lets check that trimmed reads:
```
fastqc trimmed_NS0229_Hamiltons_S1_R1_001.fastq trimmed_NS0229_Hamiltons_S1_R2_001.fastq
```


## Demultiplexing

*Breif Information on barcode.txt and popmap.txt*

Barcode.txt was created by Ludo, which contains all three plates (remeber PstI-1a,2a,3a) with the unique barcode which identify the wells on each palte. I have changed the names of my frogs in each well to give an informative name on population etc; remeber this must simple to make future analysis easy. 

Popmap.txt is created/will be created by me and has a first column as the frogIDs taken from Barcode.txt, with the second column just all filled with 'Pop.' Basically this allows us to just analyse all of the frogs, ignore the skinks, snails etc. verything is labeled 'pop' as we do not want to run stacks based on population.  

Copy trimmed data to raw folder.
```
cd .. #assuming I'm in source files
mkdir raw
mkdir # I recived an error trying process_radtags without this directory.
```

```
cd raw
cp /home/mulha552/uoo04306/frogs_gbs/source_files/trimmed_NS0229_Hamiltons_S1_R1_001.fastq /home/mulha552/uoo04306/frogs_gbs/source_files/trimmed_NS0229_Hamiltons_S1_R2_001.fastq /home/mulha552/uoo04306/frogs_gbs/raw
cd ..
```
Run demultiplexing, here is some [info](https://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php) on files and run. I have run the code below as a SLURM job, as above in adaptor trimming. 

```
cd /home/mulha552/uoo04306/frogs_gbs/source_files/raw
module load Stacks #2.61
process_radtags -P   -p ../raw/ -o ../samples/ -b ../barcodes.txt -e pstI -r -c  --inline-inline # NO -q often used for process-radtags gives me an error because of it, but no worries, cutadapatalready took care of this
 ```

My Results (maybe to ask Ludo about at a later date..)
 ```
3666129448 total sequences
 489768784 barcode not found drops (13.4%)
    889709 low quality read drops (0.0%)
         0 poly-G run drops (0.0%)
  43575443 RAD cutsite not found drops (1.2%)
3131895512 retained reads (85.4%)
```
## Concatenate reads

The goal is to have one file per sample, that we care about, (remeber there are skinks, snails etc in this sequencing run) inside the folder samples_concat. The file popmap is the stacks population map [info](http://catchenlab.life.illinois.edu/stacks/manual/#popmap); Create popmap.txt:
```
nano popmap.txt #copy in the frogs sample prefixes and keep all in a single 'pop'
# mkdir samples_concat
```
As normal we will be running concatenating as a SLURM job, but because this is in a Python Enviornment we must do this a little differently....
```
nano frogconcat.py #create python script
```
We use a python shebang for this file (rather than a shell/bash one) and we cannot specify requirments using #SBATCH, this is now done later. Copy the code below into 'frogconcat.py'

Note: popmap.txt needs to be in the current working directoy (/frogs_gbs) from which it will look for the subdirectory samples, and create a new directory as output samples_concat

```
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
I'll check the code before submitting the job
```
module load Python #need to load a python environment//ipython console
python3 frogconcat.py
```
Finally, once this is working we can submit a job, and ofr python we specify arguments here. -A specifies the account to 'charge' the job, -t specifies time, -c specifies the number of cpus per task, --mem specifies memory, --job-name specifies a job name
```
sbatch -A uoo04306 -t 48:00:00 -c 8 --mem 8G --job-name=frogconcat frogconcat.py
squeue -u mulha552 #check its working
```
## Subsampling
*Hamilton's frog have a large genome so I have special code*
We will subsample from the concatenated reads to take the first 5 million reads (code reads 20,000,000 lines because each read is 4 lines in a .fastq) so we are working with an actually workable amount of sequence. We still run the previous code, because for those samples where we may have <5million reads we need to take all of it.

 We're also running as a pyhton slurm job, as above...
```
mkdir samples_subsampled
```
```
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

!!! Don't run all different numbers, but maybe just 2, with all the samples not a subsampled popmap; I can't remeber what excatly was meant by this, so I'll ask Ludo later... but, from now on I've run everything on my *subsampled data.....*

At this stage, I'll make a popmap and exclude all samples with less than 100'000 reads (combining forward and reverse, i.e. 50k retained reads). First check number of reads (remember .fq files have 4 lines for each read), then create popmap.

```
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

### Parameter optimisation

Run on a random 30 samples.
```
#mkdir para_opti
cd para_opti
shuf ../popmap.txt | head -n 30 > popmap_opti.txt # not done if no optimisation
```
This code is a "for loop." This loop automates the creation of directories, scripts, and slurm job for processing the specified random 30 samples using the denovo_map.pl tool. Each job is submitted with different parameters (based on the loop variable i), and each job will run with different configurations of parameters (see where variable i is included in the code). 

I believe this method of paramater optimisation, roughly, follow this [paper](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12775).

```
for i in 2 3 4 5 6 7 8
do mkdir -p M$i; echo '#!/bin/sh' > runM$i.sh
echo cd /home/mulha552/uoo04306/frogs_gbs/source_files/para_opti
echo "module load Stacks/2.61-gimkl-2022a" >> runM$i.sh
echo "denovo_map.pl --samples ../samples_subsampled/ --popmap popmap_opti.txt  -o M$i  -M $i -n $i -m 3 -T 8" >> runM$i.sh
sbatch -A uoo04306 -t 2-00:00:00 -J M$i -c 8 --mem=64G runM$i.sh
done
```
```
squeue -u mulha552
```
The *populations* command takes the the SNP data from *denovo_map.pl* and, can compute a variety of useful population genetics statistics, See [here](http://catchenlab.life.illinois.edu/stacks/comp/populations.php) for more information. 
Focus on the number of loci at -R 0.8 (loci covered in 80% of inds):
```
module load Stacks/2.61-gimkl-2022a
```
```
for i in 2 3 4 5 6 7 8 
do
populations -P M$i -R 0.8 -M popmap_opti.txt --vcf
done
```
See how many loci are left by looking in the log of populations:

cut -f 1, is looking only a column one which lists fragment ID, and uniq means we are only counting unique fragment IDs....
```
for i in 2 3 4 5 6 7 8
do
echo $i
cat M$i/populations.snps.vcf | grep -v '#' | cut -f 1 | uniq | wc -l
done
```
We have the most polymoprhic loci with M=n of 2, so that is now what we'll run for the enitre study!!


## Run on the whole dataset
I will run this as SLURM job, we have a lot of data...
```
mkdir M2_final
```
```
#!/bin/bash
cd /home/mulha552/uoo04306/frogs_gbs
module load Stacks/2.61-gimkl-2022a
denovo_map.pl --samples samples_subsampled/ --popmap popmap.txt  -o M2_final  -M 2 -n 2 -m 3 -T 8
```
```
sbatch -A uoo04306 -t 5-00:00:00 -J M2Final -c 8 --mem=300G -p hugemem M2Final.sl
```

I run populations again to obtain a VCF and check for low quality samples.

```
populations -P output_refmap/ -M popmap.txt  --vcf --structure --plink --treemix --max-obs-het 0.65 -R 0.5  -O M3_final
```


```
module load VCFtools 
vcftools --vcf M3_final/populations.snps.vcf --missing-indv
sort -k 4n out.imiss
``` 

The following samples have more than 60% missing data:

```
3308_AG_marg
3204_AT_ota
3442F_TE_ota
3424M_TE_ota
3427M_TE_ota
3425M_TE_ota
3060_TB_marg
3059_TB_marg
3256M_AL_ota
3454M_TW_ota
3424F_TE_ota
3051_TB_marg
3418F_AL_Notinmetadata_NA
3435M_TE_ota
3034_TR_marg
3438F_TE_ota
3259M_AL_ota
3456M_TW_marg
3254M_CD_ota
3258M_AL_ota
3264M_AL_ota
3303_AG_ota
3213_LP_marg
3237M_DQ_ota
3305_AG_ota
3299_AG_marg
3194_CD_ota
3309_AG_marg
3205_AT_marg
3041_TR_marg
3441M_TE_marg
3311_AG_marg
3298_AG_marg
3441F_TE_ota
3199_AT_ota
3304_AG_marg
3054_TB_marg
```

Remove them from the popmap.txt to create popmap_clean.txt . Added to the other samples we now miss 19 samples.

```
grep -Ev "^3308_AG_marg\\s|^3204_AT_ota\\s|^3442F_TE_ota\\s|^3424M_TE_ota\\s|^3427M_TE_ota\\s|^3425M_TE_ota\\s|^3060_TB_marg\\s|^3059_TB_marg\\s|^3256M_AL_ota\\s|^3454M_TW_ota\\s|^3424F_TE_ota\\s|^3051_TB_marg\\s|^3418F_AL_Notinmetadata_NA\\s|^3435M_TE_ota\\s|^3034_TR_marg\\s|^3438F_TE_ota\\s|^3259M_AL_ota\\s|^3456M_TW_marg\\s|^3254M_CD_ota\\s|^3258M_AL_ota\\s|^3264M_AL_ota\\s|^3303_AG_ota\\s|^3213_LP_marg\\s|^3237M_DQ_ota\\s|^3305_AG_ota\\s|^3299_AG_marg\\s|^3194_CD_ota\\s|^3309_AG_marg\\s|^3205_AT_marg\\s|^3041_TR_marg\\s|^3441M_TE_marg\\s|^3311_AG_marg\\s|^3298_AG_marg\\s|^3441F_TE_ota\\s|^3199_AT_ota\\s|^3304_AG_marg\\s|^3054_TB_marg\\s" popmap_100k.txt > popmap_clean.txt```
```


re run populations filtering positions found in less than 80% and keeping max one snp per locus.

```
populations -P M3_final/ -M popmap_clean.txt  --vcf --structure --plink --treemix --max-obs-het 0.65 -R 0.8  --write-single-snp -O M3_final
```


Only SNPs found in 80% of individuals are kept > 269k. 

Now let's have a look at how many SNP per locus.



23570 variant sites remained.

Save it with a meaningful name:

```
 cp M3_final/populations.snps.vcf phau190indsR08maxsnps1.recode.vcf
 ```
