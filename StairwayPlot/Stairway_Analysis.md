# Stairway Plot 2
Stairway Plot is a method for inferring detailed population demographic history using the site frequency spectrum (SFS) from DNA sequence data. It does not need a pre-defined population model and can be applied to thousands of unphased sequences.

Our dataset is of a non-model organism, Hamilton's frog, without a reference genome i.e we are missing information on loci location within the genome etc. which other programme (like GONE) require. Stairway Plot 2 doesn't require any whole genome data or any information about the ancestral allele state. I used a folded SFS, avoiding all these issues!

Population sizes are inferred based on an understanding that large, expanding populations will contain an excess of rare alleles compared to a stable population and, declining, small populations will contain fewer. The SFS is related to change in effective population size by 1) mutation: novel changes in DNA sequence and 2) coalescence: the rate at which alleles ‘coalesce’ to a common ancestor. 

Note: with a reference genome superior inference can be made utilising *unfolded* data. Please see the [original GitHub Repository](https://github.com/xiaoming-liu/stairway-plot-v2) for all vignettes or tutorials relevant to running Stairway Plot 2. Mutation rate = rate/site/generation and generation time = time to first reproduction. 

## Site Frequency Spectra 
The site frequency spectra (SFS) is a summary of allele frequencies across my SNPs, it describes the proportion of samples carrying the minor alelle. I was interested in the history of the natural populations of Hamilton's frog on Stephen's Island and, Maud Island. I used [easySFS](https://github.com/isaacovercast/easySFS) to create an SFS of the polymorphic sites. It's important to note that I also needed the total number of sites from the "populations" output as Stairway Plot 2 requires information on the number of monomorphic loci. 

Stairway Plot will not incorporate information from founder effects or drift experienced by translocated populations will not impact these inferences. For natural populations it is important to understand the population structuring will confound results by showing artificial increases in Ne; Stairway Plot assumes your populations are discrete. 

```bash
mkdir stairway
cd stairway
git clone https://github.com/isaacovercast/easySFS.git
easySFS/easySFS.py -i ../HamFrogR08maxsnps1DP5.Under065.vcf -p pops_file_updated.txt --preview
```
Preview mode allows us to project segregating sites for multiple values of projecting down for each population. **We want to maximise the number of segregating sites** see easySFS help for more details, but depending on your data, downsampling will be a superior choice.
```    
Maud
(2, 447)        (3, 670)        (4, 873)        (5, 1067)       (6, 1255)       (7, 1438)       (8, 1619)       (9, 1797)       (10, 1972)      (11, 2147)      (12, 2319)      (13, 2491)      (14, 2661)      (15, 2830)      (16, 2998)      (17, 3164)      (18, 3330)      (19, 3490)      (20, 3654)     (21, 3799)       (22, 3961)      (23, 4063)      (24, 4221)      (25, 4261)      (26, 4414)      (27, 4298)      (28, 4440)      (29, 4042)      (30, 4166)      (31, 3547)      (32, 3649)      (33, 2808)      (34, 2884)      (35, 1912)      (36, 1962)      (37, 1028)      (38, 1053)      (39, 369)      (40, 378)
Takapourewa
(2, 596)        (3, 894)        (4, 1160)       (5, 1402)       (6, 1640)       (7, 1859)       (8, 2082)       (9, 2273)       (10, 2485)      (11, 2658)      (12, 2863)      (13, 3005)      (14, 3203)      (15, 3288)      (16, 3477)      (17, 3473)      (18, 3653)      (19, 3491)      (20, 3654)     (21, 3274)       (22, 3413)      (23, 2825)      (24, 2936)      (25, 2081)      (26, 2157)      (27, 1209)      (28, 1250)      (29, 466)       (30, 481)
```
Based on these results I will use 28, and 20. Note: this downsampling has important implications for calculating nseq for Stairway Plot v2 later.
```bash
easySFS/easySFS.py -i ../HamFrogR08maxsnps1DP5.Under065.vcf -p pops_file_updated.txt --proj 28,20
```
## Run Stairway Plot 2
I ran this for the two natural populations 'Maud' and 'Stephens.' Remember, that the Stairway plot is only using information from mutation and coalescence to understand historical demographics. Towards recent event it will become unreliable, what is more important to infer from this analysis is the overall trend we see in the data here. Other methods, based on recombinations or drift, such as GONE are more reliable for these more recent events, while PSMC utilising whole genome data can infer further into the past; Stairway sits in the 'middle'.

Note: I manually uploaded the .zip file from the [Stairway Plot Github](https://github.com/xiaoming-liu/stairway-plot-v2/blob/master/stairway_plot_v2.1.1.zip), I couldn't get wget or curl to work.

The folded SFS made by easySFS with *dadi* contains data for allele frequencies from 0 to N. For our purposes we removed the first column since data on invariant sites is included in the Stairway Plot blueprint file under "L" the total number of observed sites; the folded SFS is only describing the minor allele frequencies. 

```bash
mkdir stairway_plot_v2.1
unzip stairway_plot_v2.1.1.zip -d stairway_plot_v2.1
```
Navigate to the the folder with the directory *stairway_plot_es* and, create blueprint files for our two natural populations; pay attention to the parameters specified here always consult the [manual](https://github.com/xiaoming-liu/stairway-plot-v2/blob/master/READMEv2.1.pdf). Then, create the batch files (examples for Maud Island).

```bash
module load Java
java -cp stairway_plot_es Stairbuilder Maud_Island.blueprint
```
Finally, run Stairway Plot 2

```bash
bash Maud_Island.blueprint.sh
```
### SLURM Job
```
#!/bin/bash -e
#SBATCH --job-name=Frog_Staiway # job name (shows up in the queue)
#SBATCH --time=24:00:00      # Walltime (HH:MM:SS), if our job finishes before this no worries but we can give ample time in case
#SBATCH --mem=8G          # Memory in G, minimum 1G per core.
#SBATCH --cpus-per-task=8 #number of cores for our job

cd /home/mulha552/uoo04306/frogs_gbs/stairway/stairway_plot_v2.1/stairway_plot_v2.1.1
module load Java
bash Maud.blueprint.sh
bash Stephens.blueprint.sh
```
Check status
```
squeue -u mulha552
```
Before plotting it is important to understand that using to many breakpoints can 'overfit' your data. Depending on sample size select the fewer number of breakpoints see Lapierre et al. 2017.
