# Stariway Plot 2
Stairway Plot is a method for inferring detailed population demographic history using the site frequency spectrum (SFS) from DNA sequence data. It does not need a pre-defined population model and can be applied to thousands of unphased sequences; the key benefit here over other methods is that I can use this for my SNP dataset produced using genotyping-by-sequencing. My dataset is of a non-model organism, Hamilton's frog, without a reference genome i.e we are missing information on loci location within the genome etc. which other programme (like GONE) require. Luckily Stairway Plot 2 doesn't require any whole genome data or any information about the ancestral alelle state. I can use a folded SFS, avoiding all these isssues!

### Site Frequencing Spectra 
The site frequnecy spectra (SFS) is a summary of allele frequencies across my SNPs, it describes the proportion of samples (or stacks ?) carrying the minor alelle; for Stairway Plot 2 the SFS *must* include monomorphic loci or else our data will be biased, overrepresenting rare alelles. But stacks & populations only outputs the polymorphic sites. 

I am interested in the history of the natural populations of Hamilton's frog on Stephen's Island and, Maud Island. I will use [easySFS](https://github.com/isaacovercast/easySFS) to create an SFS of the polymorphic sites, but then manually add a bin for the monomoprphic sites which *are* reported in my output from populations. 

```bash
mkdir stairway
cd stairway
git clone https://github.com/isaacovercast/easySFS.git
easySFS/easySFS.py -i ../HamFrogR08maxsnps1DP5.recode.vcf -p pops_file.txt --preview
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
(2, 484)        (3, 726)        (4, 945)        (5, 1152)       (6, 1352)       (7, 1548)       (8, 1740)       (9, 1929)   (10, 2115)       (11, 2300)      (12, 2483)      (13, 2664)      (14, 2844)      (15, 3023)      (16, 3201)      (17, 3373)  (18, 3549)       (19, 3714)      (20, 3887)      (21, 4017)      (22, 4186)      (23, 4252)      (24, 4416)      (25, 4326)  (26, 4479)       (27, 4208)      (28, 4346)      (29, 3858)      (30, 3976)      (31, 3197)      (32, 3288)      (33, 2428)  (34, 2493)       (35, 1573)      (36, 1613)      (37, 878)       (38, 899)       (39, 320)       (40, 328)

Maud
(2, 467)        (3, 700)        (4, 909)        (5, 1107)       (6, 1297)       (7, 1483)       (8, 1665)       (9, 1844)   (10, 2021)       (11, 2196)      (12, 2370)      (13, 2542)      (14, 2713)      (15, 2883)      (16, 3052)      (17, 3219)  (18, 3386)       (19, 3546)      (20, 3711)      (21, 3857)      (22, 4020)      (23, 4122)      (24, 4281)      (25, 4322)  (26, 4475)       (27, 4358)      (28, 4501)      (29, 4101)      (30, 4226)      (31, 3601)      (32, 3704)      (33, 2853)  (34, 2930)       (35, 1950)      (36, 2000)      (37, 1054)      (38, 1079)      (39, 387)       (40, 396)

Motuara
(2, 455)        (3, 683)        (4, 887)        (5, 1080)       (6, 1266)       (7, 1446)       (8, 1624)       (9, 1798)   (10, 1970)       (11, 2140)      (12, 2309)      (13, 2476)      (14, 2642)      (15, 2807)      (16, 2971)      (17, 3134)  (18, 3296)       (19, 3457)      (20, 3617)      (21, 3772)      (22, 3930)      (23, 4082)      (24, 4239)      (25, 4383)  (26, 4538)       (27, 4666)      (28, 4819)      (29, 4922)      (30, 5072)      (31, 5123)      (32, 5269)      (33, 5287)  (34, 5429)       (35, 5328)      (36, 5463)      (37, 5285)      (38, 5412)      (39, 5061)      (40, 5176)      (41, 4684)  (42, 4786)       (43, 4053)      (44, 4136)      (45, 3197)      (46, 3260)      (47, 2259)      (48, 2301)      (49, 1363)  (50, 1388)       (51, 676)       (52, 687)       (53, 223)       (54, 227)

Stephens
(2, 699)        (3, 1045)       (4, 1334)       (5, 1585)       (6, 1828)       (7, 2046)       (8, 2270)       (9, 2456)   (10, 2668)       (11, 2836)      (12, 3041)      (13, 3172)      (14, 3370)      (15, 3444)      (16, 3634)      (17, 3625)  (18, 3805)       (19, 3629)      (20, 3792)      (21, 3391)      (22, 3530)      (23, 2918)      (24, 3029)      (25, 2159)  (26, 2235)       (27, 1257)      (28, 1298)      (29, 487)       (30, 502)
```

We want to *maximise* the number of segregating sites based on our downsampling. So, based on these results I will use 20, 20, 27 and, 15 i.e. the actual sample sizes.

```bash
easySFS/easySFS.py -i ../HamFrogR08maxsnps1DP5.recode.vcf -p pops_file.txt --proj 20,20,27,15
```

### Run Stairway Plot 2
I'll run this seperately for each natural population; it's key to remeber that the Stairway Plot method is only using information from mutation to understand historical demographics. So will become unreliable over the most recent ~100 generations, compared to methods which are based on recombination or drift like GONE. Nevertheless, it's super useful for our situation where we have such limited genomic resources for Hamilton's frog. I manually uploaded the .zip file from the [Stairway Plot Github](https://github.com/xiaoming-liu/stairway-plot-v2/blob/master/stairway_plot_v2.1.1.zip), I couldn't get wget or curl to work.
```bash
mkdir stairway_plot_v2.1
unzip stairway_plot_v2.1.1.zip -d stairway_plot_v2.1
```
Now, navigate to the the folder with the directory *stairway_plot_es* and create a new blueprint file, for Maud Island an example is below. !!! To ask Ludo about number of sequences, and L, and mutation rate
```
#example blueprint file
#input setting
popid: Maud # id of the population (no white space)
nseq: 40 # number of sequences
L: 21250  # total number of observed nucleic sites, including polymorphic and monomorphic
whether_folded: true # whethr the SFS is folded (true or false)
SFS:17491.80837554265 3301.764618504911 151.7453490091901 46.28577101110631 37.24021196697853 35.629259603028 33.53631783139479 29.79094109313708 29.30350398686251 30.34035372805283 15.5552977227085 0 0 0 0 0 0 0 0 0 0
#smallest_size_of_SFS_bin_used_for_estimation: 1 # default is 1; to ignore singletons, uncomment this line and change this number to 2
#largest_size_of_SFS_bin_used_for_estimation: 10 # default is nseq/2 for folded SFS
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
java -cp stairway_plot_es Stairbuilder Maud.blueprint
```
Run Stairway Plot 2
```bash
bash Maud.blueprint.sh
```
