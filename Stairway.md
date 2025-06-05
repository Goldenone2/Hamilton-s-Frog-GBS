# Stariway Plot 2
Stairway Plot is a method for inferring detailed population demographic history using the site frequency spectrum (SFS) from DNA sequence data. It does not need a pre-defined population model and can be applied to thousands of unphased sequences; the key benefit here over other methods is that I can use this for my SNP dataset produced using genotyping-by-sequencing. My dataset is of a non-model organism, Hamilton's frog, without a reference genome i.e we are missing information on loci location within the genome etc. which other programme (like GONE) require. Luckily Stairway Plot 2 doesn't require any whole genome data or any information about the ancestral alelle state. I can use a folded SFS, avoiding all these isssues!

### Site Frequencing Spectra 
The site frequnecy spectra (SFS) is a summary of allele frequencies across my SNPs, it describes the proportion of samples (or stacks ?) carrying the minor alelle; for Stairway Plot 2 the SFS *must* include monomorphic loci or else our data will be biased, overrepresenting rare alelles. But stacks & populations only outputs the polymorphic sites. 

I am interested in the history of the natural populations of Hamilton's frog on Stephen's Island and, Maud Island. I will use [easySFS](https://github.com/isaacovercast/easySFS) to create an SFS of the polymorphic sites, but then manually add a bin for the monomoprphic sites which *are* reported in my output from populations. 

```bash
mkdir stairway
cd stairway
git clone https://github.com/isaacovercast/easySFS.git
easySFS/easySFS.py -i ../HamFrogR08maxsnps1minGT5.vcf -p pops_file.txt --preview
```
```
  Processing 4 populations - ['Boat', 'Maud', 'Motuara', 'Stephens']
  Sampling one snp per locus (CHROM)

    Running preview mode. We will print out the results for # of segregating sites
    for multiple values of projecting down for each population. The dadi
    manual recommends maximizing the # of seg sites for projections, but also
    a balance must be struck between # of seg sites and sample size.

    For each population you should choose the value of the projection that looks
    best and then rerun easySFS with the `--proj` flag.
    
Boat
(2, 258)        (3, 387)        (4, 502)        (5, 611)        (6, 714)        (7, 815)        (8, 914)        (9, 1011)       (10, 1107)       (11, 1201)      (12, 1295)      (13, 1387)      (14, 1479)      (15, 1570)      (16, 1660)      (17, 1750)      (18, 1839)      (19, 1928)       (20, 2016)      (21, 2104)      (22, 2191)      (23, 2275)      (24, 2362)      (25, 2432)      (26, 2517)      (27, 2553)       (28, 2635)      (29, 2537)      (30, 2613)      (31, 2312)      (32, 2377)      (33, 1932)      (34, 1983)      (35, 1363)      (36, 1398)       (37, 817)       (38, 837)       (39, 316)       (40, 323)

Maud
(2, 252)        (3, 378)        (4, 490)        (5, 594)        (6, 694)        (7, 791)        (8, 886)        (9, 979)        (10, 1071)       (11, 1162)      (12, 1252)      (13, 1341)      (14, 1429)      (15, 1517)      (16, 1604)      (17, 1691)      (18, 1777)      (19, 1863)       (20, 1948)      (21, 2033)      (22, 2117)      (23, 2199)      (24, 2283)      (25, 2355)      (26, 2438)      (27, 2489)       (28, 2570)      (29, 2544)      (30, 2621)      (31, 2463)      (32, 2533)      (33, 2153)      (34, 2210)      (35, 1609)      (36, 1650)       (37, 953)       (38, 976)       (39, 374)       (40, 383)

Motuara
(2, 250)        (3, 376)        (4, 487)        (5, 591)        (6, 691)        (7, 788)        (8, 882)        (9, 975)        (10, 1067)       (11, 1157)      (12, 1247)      (13, 1335)      (14, 1423)      (15, 1511)      (16, 1597)      (17, 1683)      (18, 1769)      (19, 1854)       (20, 1938)      (21, 2022)      (22, 2106)      (23, 2189)      (24, 2272)      (25, 2355)      (26, 2437)      (27, 2519)       (28, 2600)      (29, 2675)      (30, 2756)      (31, 2825)      (32, 2904)      (33, 2960)      (34, 3038)      (35, 3056)      (36, 3132)       (37, 3123)      (38, 3197)      (39, 3124)      (40, 3194)      (41, 3031)      (42, 3095)      (43, 2775)      (44, 2831)       (45, 2355)      (46, 2401)      (47, 1811)      (48, 1845)      (49, 1163)      (50, 1184)      (51, 621)       (52, 632)       (53, 215)        (54, 219)

Stephens
(2, 432)        (3, 644)        (4, 816)        (5, 959)        (6, 1097)       (7, 1217)       (8, 1341)       (9, 1442)       (10, 1559)       (11, 1646)      (12, 1758)      (13, 1824)      (14, 1931)      (15, 1979)      (16, 2082)      (17, 2105)      (18, 2205)      (19, 2164)       (20, 2258)      (21, 2143)      (22, 2228)      (23, 1975)      (24, 2048)      (25, 1602)      (26, 1657)      (27, 1042)       (28, 1076)      (29, 438)       (30, 451)
```

We want to *maximise* the number of segregating sites based on our downsampling. So, based on these results I will use 20, 20, 27 and, 15 i.e. the actual sample sizes (obviously you cannot sample any higher than the real value of n).

```bash
easySFS/easySFS.py -i ../HamFrogR08maxsnps1minGT5.vcf -p pops_file.txt --proj 20,20,27,15
```

### Run Stairway Plot 2
I'll run this seperately for each natural population; it's key to remeber that the Stairway Plot method is only using information from mutation to understand historical demographics. So will become unreliable over the most recent ~100 generations, compared to methods which are based on recombination or drift like GONE. Nevertheless, it's super useful for our situation where we have such limited genomic resources for Hamilton's frog. I manually uploaded the .zip file from the [Stairway Plot Github](https://github.com/xiaoming-liu/stairway-plot-v2/blob/master/stairway_plot_v2.1.1.zip), I couldn't get wget or curl to work.
```bash
mkdir stairway_plot_v2.1
unzip stairway_plot_v2.1.1.zip -d stairway_plot_v2.1
```
Now, navigate to the the folder with the directory *stairway_plot_es* and create a new blueprint file, for Maud Island an example is below.
```
#example blueprint file
#input setting
popid: Maud # id of the population (no white space)
nseq: 40 # number of sequences
L: 2447914  # total number of observed nucleic sites, including polymorphic and monomorphic
whether_folded: true # whethr the SFS is folded (true or false)
SFS:9661.055806599743 1704.375413890788 83.63644452072168 29.10547129405674 23.85557294590284 23.52054970321189 22.94194565168617 20.48442680527145 17.35117731143766 15.32602919721953 7.347162079966721 0 0 0 0 0 0 0 0 0
smallest_size_of_SFS_bin_used_for_estimation: 1 # default is 1; to ignore singletons, uncomment this line and change this number to 2
largest_size_of_SFS_bin_used_for_estimation: 20 # default is nseq/2 for folded SFS
pct_training: 0.67 # percentage of sites for training
nrand: 9      19      28    38 # number of random break points for each try (separated by white space) integers only
project_dir: Muad_Island # project directory
stairway_plot_dir: stairway_plot_es # directory to the stairway plot files
ninput: 200 # number of input files to be created for each estimation
random_seed: 77
#output setting
mu: 2.7e-6 # assumed mutation rate per site per generation
year_per_generation: 7 # assumed generation time (in years)
#plot setting
plot_title: two-epoch_fold # title of the plot
xrange: 0.1,10000 # Time (1k year) range; format: xmin,xmax; "0,0" for default
yrange: 0,0 # Ne (1k individual) range; format: xmin,xmax; "0,0" for default
xspacing: 2 # X axis spacing
yspacing: 2 # Y axis spacing
fontsize: 12 # Font size
```
Create batch file.
```bash
module load Java
java -cp stairway_plot_es Stairbuilder Maud.blueprint
```
Run Stairway Plot 2
```bash
bash Maud.blueprint.sh
```
