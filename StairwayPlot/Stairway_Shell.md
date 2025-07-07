# Stariway Plot 2
Stairway Plot is a method for inferring detailed population demographic history using the site frequency spectrum (SFS) from DNA sequence data. It does not need a pre-defined population model and can be applied to thousands of unphased sequences; the key benefit here over other methods is that I can use this for my SNP dataset produced using genotyping-by-sequencing. My dataset is of a non-model organism, Hamilton's frog, without a reference genome i.e we are missing information on loci location within the genome etc. which other programme (like GONE) require. Stairway Plot 2 doesn't require any whole genome data or any information about the ancestral alelle state. I can use a folded SFS, avoiding all these isssues!

### Site Frequencing Spectra 
The site frequnecy spectra (SFS) is a summary of allele frequencies across my SNPs, it describes the proportion of samples carrying the minor alelle. I am interested in the history of the natural populations of Hamilton's frog on Stephen's Island and, Maud Island. I will use [easySFS](https://github.com/isaacovercast/easySFS) to create an SFS of the polymorphic sites. It's important to note that I also need the total number of sites from the "populations" output as Stairway Plot 2 requires information on the number of monomorphic loci. 

```bash
mkdir stairway
cd stairway
git clone https://github.com/isaacovercast/easySFS.git
easySFS/easySFS.py -i ../HamFrogR08maxsnps1DP5.recode.vcf -p pops_file.txt --preview
```
Preview mode allows us to project segragating sites for multiple vlaues of projecting down for each populaiton; we want to maximise the number of segregating sites.
```    
Boat
(2, 484)        (3, 726)        (4, 945)        (5, 1152)       (6, 1352)       (7, 1548)       (8, 1740)       (9, 1929)       (10, 2115)      (11, 2300)      (12, 2483)      (13, 2664)    (14, 2844)      (15, 3023)      (16, 3201)      (17, 3373)      (18, 3549)      (19, 3714)      (20, 3887)      (21, 4017)      (22, 4186)      (23, 4252)      (24, 4416)   (25, 4326)       (26, 4479)      (27, 4208)      (28, 4346)      (29, 3858)      (30, 3976)      (31, 3197)      (32, 3288)      (33, 2428)      (34, 2493)      (35, 1573)      (36, 1613)    (37, 878)       (38, 899)       (39, 320)       (40, 328)

Maud
(2, 467)        (3, 700)        (4, 909)        (5, 1107)       (6, 1297)       (7, 1483)       (8, 1665)       (9, 1844)       (10, 2021)      (11, 2196)      (12, 2370)      (13, 2542)    (14, 2713)      (15, 2883)      (16, 3052)      (17, 3219)      (18, 3386)      (19, 3546)      (20, 3711)      (21, 3857)      (22, 4020)      (23, 4122)      (24, 4281)   (25, 4322)       (26, 4475)      (27, 4358)      (28, 4501)      (29, 4101)      (30, 4226)      (31, 3601)      (32, 3704)      (33, 2853)      (34, 2930)      (35, 1950)      (36, 2000)    (37, 1054)      (38, 1079)      (39, 387)       (40, 396)

Motuara
(2, 455)        (3, 683)        (4, 887)        (5, 1080)       (6, 1266)       (7, 1446)       (8, 1624)       (9, 1798)       (10, 1970)      (11, 2140)      (12, 2309)      (13, 2476)    (14, 2642)      (15, 2807)      (16, 2971)      (17, 3134)      (18, 3296)      (19, 3457)      (20, 3617)      (21, 3772)      (22, 3930)      (23, 4082)      (24, 4239)   (25, 4383)       (26, 4538)      (27, 4666)      (28, 4819)      (29, 4922)      (30, 5072)      (31, 5123)      (32, 5269)      (33, 5287)      (34, 5429)      (35, 5328)      (36, 5463)    (37, 5285)      (38, 5412)      (39, 5061)      (40, 5176)      (41, 4684)      (42, 4786)      (43, 4053)      (44, 4136)      (45, 3197)      (46, 3260)      (47, 2259)   (48, 2301)       (49, 1363)      (50, 1388)      (51, 676)       (52, 687)       (53, 223)       (54, 227)

Stephens
(2, 699)        (3, 1045)       (4, 1334)       (5, 1585)       (6, 1828)       (7, 2046)       (8, 2270)       (9, 2456)       (10, 2668)      (11, 2836)      (12, 3041)      (13, 3172)    (14, 3370)      (15, 3444)      (16, 3634)      (17, 3625)      (18, 3805)      (19, 3629)      (20, 3792)      (21, 3391)      (22, 3530)      (23, 2918)      (24, 3029)   (25, 2159)       (26, 2235)      (27, 1257)      (28, 1298)      (29, 487)       (30, 502)

```
Based on these results I will use 20, 20, 27 and, 15 i.e. the actual sample sizes. Obviously you cannot sample any higher than the real value of n.
```bash
easySFS/easySFS.py -i ../HamFrogR08maxsnps1minGT5.vcf -p pops_file.txt --proj 20,20,27,15
```
### Run Stairway Plot 2
I'll run this for the two natural populations 'Maud' and 'Stephens.' Remember, that the stariway plot is only using information from mutation to understand historical demogrphics. So, will become unreliable over the most recent ~100 generations, what is more important to infer from this analysis is the overall trend we see in the data here. Other methods, based on recombinations or drift, such as GONE are more reliable for these more recent events. Nevertheless, it's useful for our situaiton where we have limited geneitc resources for Hamilton's forg. 

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
