# Analysis Overview
I have measured inbreeding as F<sub>H</sub>, this estimate is inherently related to estimates of genetic diversity and, hence, you can see it effectively follows the opposite pattern to observed heterozygosity. Our inbreeding coefficient may be calculated to be ≤0 even when there is mating between relatives, because the estimate is based of a comparison to an "ideal" random mating population under Hardy Weinberg.In small populations random mating doesn't exclude relatives breeding, especially in our case where historical demographic factors (see Stairway Plots) are leading to excess heterozygosity in the smallest population, Takapourewa. 

### Setting up related
An annoying first time step...

First, download the .tar.gz file from the [github](https://github.com/timothyfrasier/related) into *this* working directory.

Second, install [Rtools](https://cran.r-project.org/bin/windows/Rtools/), which is required to compile a package direct from source. You must ensure that the environment variable "PATH" has been updated:
```{r}
Sys.which("gcc")
Sys.which("make")
```
If these both return a path then Rtools has installed correctly, if not you maust manually add Rtools to "PATH."

Now install related
```{r}
# install.packages("related_1.0.tar.gz", repos=NULL, type="source")
```

### Data Import and Setup
```{r}
# Clear workspace
rm(list = ls())

# Load Packages
library("vcfR")
library("ggplot2")
library("related")
library("adegenet")

```
Now, related is an R environment for an actual program called [COANCESTRY](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1755-0998.2010.02885.x) and we need to reformat our vcf file.....
```{r}
vcf <- read.vcfR("HamFrogR08maxsnps1DP5.recode.vcf")

# Pull the genotype calls form the vcf, into a matrix
gt <- extract.gt(vcf, element = "GT")

# Add the population data, this must be a two letter code affixed to sample names
indiv_ids <- colnames(gt)
        # vector to hold our pop codes
populations <- character(length(indiv_ids))
        # a for loop to assign pop codes used by related

for (i in seq_along(indiv_ids)) { # loops over all individuals
  id <- indiv_ids[i]  # takes the current individual ID
  
  if (grepl("^B_", id)) { # grepl is then used for pattern matching 
    populations[i] <- "BB"
    
  } else if (grepl("^G1_", id)) {
    populations[i] <- "MA"   # Special case, my formatting error
    
  } else if (grepl("^M_", id)) {
    populations[i] <- "MA" 
    
  } else if (grepl("^MT_", id)) {
    populations[i] <- "MT"
    
  } else if (grepl("^T_", id)) {
    populations[i] <- "TA"
    
  } 
}

        # then add population codes to individual IDs
new_id <- paste0(populations, "_", indiv_ids)
colnames(gt) <- new_id

# convert genotype string e.g. 0/0 to two separate columns required by related
convert_gt_to_alleles <- function(gt_string) {
  if (is.na(gt_string) || gt_string == "./.") return(c(NA, NA))  # handle missing data
  alleles <- unlist(strsplit(gt_string, "[/|]"))  # split genotype string on '/' or '|'
  as.numeric(alleles) # makes each string a numberic vector
}

# Finally, we make a data frame!
n_ind <- ncol(gt)
n_snp <- nrow(gt)

        # make an empty dataframe
allele_df <- data.frame(matrix(ncol = n_snp * 2, nrow = n_ind)) # 2 alleles * SNP

        # name columns 
colnames(allele_df) <-unlist(lapply(1:n_snp, function(i) c(paste0("SNP", i, "_A1"), paste0("SNP", i, "_A2"))))

        # for loop to convert genotypes
for (i in 1:n_ind) { # loop over individuals
            # take all SNP string for i, and run our function from above
            # returns a list the length of all SNPs
            # each element is a numeric vector of the two alleles
alleles_list <- lapply(gt[, i], convert_gt_to_alleles)
            # bind all these alleles vectors into a matrix by rows
            # i.e. rows = SNPs and column alleles (2)
alleles_matrix <- do.call(rbind, alleles_list)   
            # transpose the matrix, then flatten this into a single vector
            # flatten i.e. go through columns and collect allele 1 & 2
            # assign this vector to the i-th row of alleles df
allele_df[i, ] <- as.vector(t(alleles_matrix))
}

# add our individuals as the first column 
allele_df <- cbind(IndividualID = colnames(gt), allele_df)

# Finally, I discovered that we cannot encode our alleles as 0 and 1 because zero is treated as NA here
# we will change it to 1 and 2
allele_df[,-1] <- lapply(allele_df[,-1], function(col) {
  is_missing <- is.na(col)
  col_new <- col + 1
  col_new[is_missing] <- 0
  return(as.integer(col_new))
})

```

That was so hard to understand, and is confusing as hell, chatgpt is my friend here when it is *not* involved in real analysis. I'd recommend saving this as a .txt file for future use, we don't need to make this odd binary format again & again .... but it is important to understand what is going on with the the data formatting. 
```{r}
write.table(allele_df, "HamGeno.txt", quote=FALSE, row.names=FALSE, col.names=TRUE)
```

### Calculate Relatedness
I will calculate relatedness using the [triadic likelihood method](https://www.cambridge.org/core/journals/genetics-research/article/triadic-ibd-coefficients-and-applications-to-estimating-pairwise-relatedness/19C27DCC0F90870C52B5040132922281).

I must specify how many reference individuals are selected from my data. Reference individuals serve to estimate population allele frequencies and genotype probabilities, but this doesn't stop them from being used as focal individuals; the default is 100, so I will nominate my entire data set.

```{r}
results <- coancestry("HamGeno.txt", trioml = 1, trioml.num.reference = 82)
#RDS means we can save the R object for plotting later :) 
saveRDS(results, file = "TrioML_results.rds")
```
I ran this as a SLURM job, because it ran all weekend on my laptop and had not finished!

```
#!/bin/bash -e
#SBATCH --job-name= TrioML # job name (shows up in the queue)
#SBATCH --time=48:00:00      # Walltime (HH:MM:SS), if our job finishes before this no worries but we can give ample time in case
#SBATCH --mem= 8G          # Memory in G, minimum 1G per core.
#SBATCH --cpus-per-task= 8 #number of cores for our job

cd /home/mulha552/uoo04306/frogs_gbs/TrioML
module load R
Rscript --vanilla TrioML_SLURM.R
echo done
```
