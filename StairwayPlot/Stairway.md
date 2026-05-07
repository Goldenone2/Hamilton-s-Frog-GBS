# Stariway Plot 2
Stairway Plot is a method for inferring detailed population demographic history using the site frequency spectrum (SFS) from DNA sequence data. It does not need a pre-defined population model and can be applied to thousands of unphased sequences; the key benefit here over other methods is that I can use this for my SNP dataset produced using genotyping-by-sequencing. My dataset is of a non-model organism, Hamilton's frog, without a reference genome i.e we are missing information on loci location within the genome etc. which other programme (like GONE) require. Stairway Plot 2 doesn't require any whole genome data or any information about the ancestral alelle state. I can use a folded SFS, avoiding all these isssues!

Population sizes are inferred based on an understanding that large, expanding populations will contain an excess of rare alleles compared to a stable population and, declining, small populations will contain fewer. The SFS is related to change in effective population size by 1) mutation: novel changes in DNA sequence and 2) coalescence: the rate at which alleles ‘coalesce’ to a common ancestor. 

Note: with a reference genome superiori inference can be made utilising *unfolded* data. Please see the [original GitHub Repository](https://github.com/xiaoming-liu/stairway-plot-v2) for all vignettes or tutorials relevant to running Stairway Plot 2. Mutation rate = rate/site/generation and generation time = time to first reproduction. 

### Site Frequencing Spectra 
The site frequnecy spectra (SFS) is a summary of allele frequencies across my SNPs, it describes the proportion of samples carrying the minor alelle. I am interested in the history of the natural populations of Hamilton's frog on Stephen's Island and, Maud Island. I will use [easySFS](https://github.com/isaacovercast/easySFS) to create an SFS of the polymorphic sites. It's important to note that I also need the total number of sites from the "populations" output as Stairway Plot 2 requires information on the number of monomorphic loci. 

Stairway Plot will not incorporate information from founder effects or drift experienced by translocated populations will not impact these inferences, however for natural populations it is important to understand the population structuring will confound results by showing artificial increases in Ne .... Stairway Plot assumes your populations are discrete. 

```bash
mkdir stairway
cd stairway
git clone https://github.com/isaacovercast/easySFS.git
easySFS/easySFS.py -i ../HamFrogR08maxsnps1DP5.recode.vcf -p pops_file_updated.txt --preview
```
Preview mode allows us to project segragating sites for multiple vlaues of projecting down for each populaiton. **We want to maximise the number of segregating sites** see easySFS help for more details, but depenindg on your data, downsampling will be a superoir choice; do not sample *more* than the real value of n.
```    
Maud
(2, 466)        (3, 700)        (4, 910)        (5, 1108)       (6, 1299)       (7, 1486)       (8, 1668)       (9, 1849)       (10, 2026)      (11, 2203)      (12, 2377)      (13, 2550)     (14, 2722)       (15, 2893)      (16, 3063)      (17, 3231)      (18, 3399)      (19, 3567)      (20, 3733)      (21, 3898)      (22, 4063)      (23, 4227)      (24, 4391)      (25, 4554)     (26, 4716)       (27, 4877)      (28, 5038)      (29, 5199)      (30, 5358)      (31, 5518)      (32, 5676)      (33, 5834)      (34, 5992)      (35, 6149)      (36, 6305)      (37, 6461)     (38, 6616)       (39, 6771)      (40, 6926)      (41, 7080)      (42, 7233)      (43, 7386)      (44, 7538)      (45, 7690)      (46, 7842)      (47, 7993)      (48, 8143)      (49, 8294)     (50, 8443)       (51, 8592)      (52, 8741)      (53, 8889)      (54, 9037)      (55, 9185)      (56, 9332)      (57, 9478)      (58, 9624)      (59, 9770)      (60, 9915)      (61, 10060)    (62, 10205)      (63, 10349)     (64, 10492)     (65, 10636)     (66, 10778)     (67, 10921)     (68, 11063)     (69, 11204)     (70, 11346)     (71, 11486)     (72, 11627)     (73, 11766)    (74, 11906)      (75, 12043)     (76, 12182)     (77, 12319)     (78, 12457)     (79, 12591)     (80, 12728)     (81, 12852)     (82, 12989)     (83, 13098)     (84, 13234)     (85, 13327)    (86, 13462)      (87, 13527)     (88, 13660)     (89, 13646)     (90, 13776)     (91, 13693)     (92, 13821)     (93, 13634)     (94, 13759)     (95, 13379)     (96, 13499)     (97, 12930)    (98, 13042)      (99, 12300)     (100, 12405)    (101, 11487)    (102, 11583)    (103, 10409)    (104, 10494)    (105, 9339)     (106, 9414)     (107, 8165)     (108, 8229)     (109, 7079)    (110, 7134)      (111, 5978)     (112, 6024)     (113, 5009)     (114, 5046)     (115, 4075)     (116, 4105)     (117, 3243)     (118, 3266)     (119, 2571)     (120, 2590)     (121, 1944)    (122, 1957)      (123, 1448)     (124, 1458)     (125, 1022)     (126, 1029)     (127, 702)      (128, 707)      (129, 420)      (130, 423)      (131, 232)      (132, 234)      (133, 72)      (134, 72)

Stephens
(2, 699)        (3, 1045)       (4, 1334)       (5, 1585)       (6, 1828)       (7, 2046)       (8, 2270)       (9, 2456)       (10, 2668)      (11, 2836)      (12, 3041)      (13, 3172)     (14, 3370)       (15, 3444)      (16, 3634)      (17, 3625)      (18, 3805)      (19, 3629)      (20, 3792)      (21, 3391)      (22, 3530)      (23, 2918)      (24, 3029)      (25, 2159)     (26, 2235)       (27, 1257)      (28, 1298)      (29, 487)       (30, 502)
```
Based on these results I will use 20, 20, 27 and, 15 i.e. the actual sample sizes.

```bash
easySFS/easySFS.py -i ../HamFrogR08maxsnps1minGT5.vcf -p pops_file_updated.txt --proj 67,15
```
### Run Stairway Plot 2
I'll run this for the two natural populations 'Maud' and 'Stephens.' Remember, that the stariway plot is only using information from mutation and coalescence to understand historical demogrphics. Towards recent event it will become unreliable, what is more important to infer from this analysis is the overall trend we see in the data here. Other methods, based on recombinations or drift, such as GONE are more reliable for these more recent events, while PSMC utilsing whole genome data can infer further into the past .... Stairway sits in the 'middle'.

Note: I manually uploaded the .zip file from the [Stairway Plot Github](https://github.com/xiaoming-liu/stairway-plot-v2/blob/master/stairway_plot_v2.1.1.zip), I couldn't get wget or curl to work.

The folded SFS made by easySFS contain data for alelle frequencies from 0 to (in our case) n .... for our purposes we will remove the first column. Data on invariant sites is included in the Stairway Plot blueprint file under "L" the total number of observed sites; the folded SFS is only describing the minor alelle frequencies. 

```bash
mkdir stairway_plot_v2.1
unzip stairway_plot_v2.1.1.zip -d stairway_plot_v2.1
```
Now, Navigate to the the folder with the directory *stairway_plot_es* and, create blueprint files for our two natural populations; pay attention to the paramters specified here always consult the [manual](https://github.com/xiaoming-liu/stairway-plot-v2/blob/master/READMEv2.1.pdf). Then, create the batch files (example for Maud Island).

```bash
module load Java
java -cp stairway_plot_es Stairbuilder Maud_Island.blueprint
```
Finally, run Stairway Plot 2

```bash
bash Maud_Island.blueprint.sh
```
### SLURM Job
Submitting all Maud Island & translocation samples seems to take quite a while to run; I will run it in a SLURM job. 

```bash
nano Maud67.sl
```
Copy & Paste
```
#!/bin/bash -e
#SBATCH --job-name=Maud67 # job name (shows up in the queue)
#SBATCH --time=24:00:00      # Walltime (HH:MM:SS), if our job finishes before this no worries but we can give ample time in case
#SBATCH --mem=8G          # Memory in G, minimum 1G per core.
#SBATCH --cpus-per-task=8 #number of cores for our job

cd /home/mulha552/uoo04306/frogs_gbs/stairway/stairway_plot_v2.1/stairway_plot_v2.1.1
module load Java
bash Maud67.blueprint.sh
```
Run
```
sbatch Maud67.sl
```
Check status
```
squeue -u mulha552
```

Before plotting it is important to understand that using to many breakpoints can 'overfit' your data. Depedning on sample size select the fewest number of breakpoints see Lapierre et al. 2017.
