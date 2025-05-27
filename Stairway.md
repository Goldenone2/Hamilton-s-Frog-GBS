# Stariway Plot 2
Stairway Plot is a method for inferring detailed population demographic history using the site frequency spectrum (SFS) from DNA sequence data. It does not need a pre-defined population model and can be applied to thousands of unphased sequences; the key benefit here over other methods is that I can use this for my SNP dataset produced using genotyping-by-sequencing. My dataset is of a non-model organism, Hamilton's frog, without a reference genome i.e we are missing information on loci location within the genome etc. which other programme (like GONE) require. Luckily Stairway Plot 2 doesn't require any whole genome data or any information about the ancestral alelle state. I can use a folded SFS, avoiding all these isssues!

### Site Frequencing Spectra 
The site frequnecy spectra (SFS) is a summary of allele frequencies across my SNPs, it describes the proportion of samples (or stacks ?) carrying the minor alelle; for Stairway Plot 2 the SFS *must* include monomorphic loci or else our data will be biased, overrepresenting rare alelles. But stacks & populations only outputs the polymorphic sites. 

I am interested in the history of the natural populations of Hamilton's frog on Stephen's Island and, Maud Island. I will use [easySFS](https://github.com/isaacovercast/easySFS) to create an SFS of the polymorphic sites, but then manually add a bin for the monomoprphic sites which *are* reported in my output from populations. 


### Run Stairway Plot 2
I'll run this seperately for each natural population; it's key to remeber that the Stairway Plot method is only using information from mutation to understand historical demogrpahics. So will become unreliable over the most recent ~100 generations, compared to methods which are based on recombination or drift like GONE. Nevertheless, it's super useful for our situation where we have such limited genomic resources for Hamilton's frog.
