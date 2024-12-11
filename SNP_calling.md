# SNP calling

First, I create a folder for the source files and folders to run Stacks, raw (which will contain the trimmed data) and samples. 

```
mkdir source_files raw samples 
```

Once the raw data is in source_files, I go in there and take a subset of 1'00000 lines, 250'000 reads.

```
cd source_files
zcat Grasshopper_GBS_S1_R1_001.fastq.gz | head -n 1000000 > testR1.fastq
zcat Grasshopper_GBS_S1_R2_001.fastq.gz | head -n 1000000 > testR2.fastq
```

```
module load FastQC # Anything after the # is ignored, # means comment, useful to document your code.
fastqc *.fastq # The star means everything that ends with fastq.
```
I then look at the generated testR1.html and testR2.html after downloading them. Quite a few adapter's in there.

## Adapter trimming

Trimming off adapters and removing reads shorter than 95bp with cutadapt.

First on the sample files. The files with trimmed are the output files. 
Use `cutadapt --help` after loading the module to learn what all the parameters are.

```
cd source_files
module load cutadapt #
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC  -q 25 -o trimmed_testR1.fastq  --minimum-length 95:95 --length 50  -p  trimmed_testR2.fastq testR1.fastq  testR2.fastq
cd ..
```
Let's check  the adapters are gone
```
fastqc trimmed*fastq # the star make it that it runs on anything that start with trimmed and end with fastq
```
Check the .html files. That worked!

Let's now run that on all reads, with a 95bp reads limit, so that we have one common length for all reads.
```
cd source_files
cutadapt -j 8 -a AGATCGGAAGAGC -A AGATCGGAAGAGC  -q 25 -o trimmed_Grasshopper_GBS_S1_R1_001.fastq   --minimum-length 95:95  --length 95  -p  trimmed_Grasshopper_GBS_S1_R2_001.fastq Grasshopper_GBS_S1_R1_001.fastq.gz  Grasshopper_GBS_S1_R2_001.fastq.gz
cd ..
```
Check that output:
```
fastqc trimmed_Grasshopper_GBS_S1_R2_001.fastq.gz Grasshopper_GBS_S1_R1_001.fastq.gz
```


## Demultiplexing




Copy trimmed data to raw folder.

```
cd raw
ln -s ../source_files/trimmed_Grasshopper* .
cd ..
```
Run demultiplexing ([info](https://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php) on files and run)

```
module load Stacks #2.61
process_radtags -P   -p raw/ -o ./samples/ -b barcodes.txt -e pstI -r -c  --inline-inline # NO -q often used for process-radtags gives me an error because of it, but no worries, cutadapatalready took care of this
 ```
 Good results:
 
 ```
1168510270 total sequences
 336643722 barcode not found drops (28.8%)
         0 low quality read drops (0.0%)
  17120888 RAD cutsite not found drops (1.5%)
 814745660 retained reads (69.7%)
```
### concatenate reads

The goal is to have one file per sample inside the folder samples_concat. The file popmap is the stacks population map. Google population map stacks and create popmap.txt.

Load any python module (module load Python) and then the code below after entering the ipython console.
```
import os
#os.mkdir("samples_concat")
allfiles_to_concat = os.listdir("samples/")
alltokeep=[]
with open("popmap.txt") as f:
	for line in f:
		print(line)
		samplename=line.split("\t")[0]
		#print (samplename)
		checkfiles=[filename for filename in allfiles_to_concat  if filename.startswith(samplename)]
		if len (checkfiles)!=4: # that was weirdly complicated because some sample name are contained in others different ways, but the vcheck above solve it uysing the rem file
			print(line)
			raise Exception
		else:
			os.system("zcat "+" ".join(["samples/"+checkfile for checkfile in checkfiles])+"|  gzip -c > samples_concat/"+samplename+".fq.gz" )

```

### Parameters optimisation

At this stage, I'll make a popmap and exclude all Celmisia samples as well as all samples with less than 100'000 reads (combining forward and reverse, i.e. 50k retained reads). 3450F_TW_marg and 3450M_TW_marg  are also found twice on the same plate. Something is off, I ignore them.



```
popmap_100k.txt
```


### Parameter optimisation

Run on a random 30 samples.
```
shuf popmap_100k.txt | head -n 30 > popmap_opti.txt # not done if no optimisation
```


```
for i in 2 3 4 5 6 7 8
do mkdir -p M$i; echo '#!/bin/sh' > runM$i.sh
echo "module load Stacks/2.61-gimkl-2022a" >> runM$i.sh
echo "denovo_map.pl --samples samples_concat/ --popmap popmap_opti.txt  -o M$i  -M $i -n $i -m 3 -T 8" >> runM$i.sh
sbatch -A uoo00116 -t 2-00:00:00 -J M$i -c 8 --mem=64G runM$i.sh
done
```

Focus on the number of loci at -R 0.8 (loci covered in 80% of inds):

```
for i in 2 3 4 5 6 7 8 
do
populations -P M$i -R 0.8 -M popmap.txt --vcf
done
```
see how many loci are left by looking in the log of populations:

```
for i in 2 3 4 5 6 7 8
do
echo $i
cat M$i/populations.snps.vcf | grep -v '#' | cut -f 1 | uniq | wc -l
done
```

We'll keep 3, it seems to have a lot of data.


## Run on the whole dataset


```
#!/bin/s
module load Stacks/2.61-gimkl-2022a
denovo_map.pl --samples samples_concat/ --popmap popmap_100k.txt  -o M3_final  -M 3 -n 3 -m 3 -T 8 

sbatch -A uoo00116 -t 5-00:00:00 -J M3 -c 8 --mem=300G -p hugemem runM3final.sh
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
