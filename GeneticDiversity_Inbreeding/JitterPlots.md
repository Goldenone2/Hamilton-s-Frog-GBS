# Genetic Diversity and Inbreeding
I'd like to compare heterozygosity and, inbreeding across the four populations I've sequenced. 

## Data Import and Setup
```r
# Clear workspace
rm(list = ls())

# Load required Packages
library(ggplot2)
library(dplyr)
library(patchwork)

# Load data
Het <- read.table("Heterozygosity.het", header = TRUE)

#Add sample information to the dataframe
Het <- Het %>%
  mutate(Pop = case_when(
    grepl("^B", INDV) ~ "Boat Bay",   
    grepl("^T", INDV) ~ "Stephen's Island",    
    grepl("^M_", INDV) ~ "Maud Island",
    grepl("^MT", INDV) ~ "Motuara",   
    grepl("^G", INDV) ~ "Maud Island"
    ))

#Calculate the proportion of heterozygotes
Het$Prop.Het <- (Het$N_SITES - Het$O.HOM.) / Het$N_SITES
```
## Simple ANOVA
First for heterozygosity.
```r
m1 <- aov(Prop.Het ~ Pop, data = Het)
summary(m1)
TukeyHSD(m1)
```
```
> m1 <- aov(Prop.Het ~ Pop, data = Het)
> summary(m1)
            Df    Sum Sq   Mean Sq F value Pr(>F)    
Pop          3 0.0014831 0.0004944   89.35 <2e-16 ***
Residuals   78 0.0004316 0.0000055                   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> TukeyHSD(m1)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = Prop.Het ~ Pop, data = Het)

$Pop
                                      diff          lwr          upr     p adj
Maud Island-Boat Bay         -0.0004765059 -0.002429307 0.0014762949 0.9185007
Motuara-Boat Bay             -0.0016438475 -0.003465687 0.0001779915 0.0917234
Stephen's Island-Boat Bay     0.0100685285  0.007959263 0.0121777946 0.0000000
Motuara-Maud Island          -0.0011673416 -0.002989181 0.0006544975 0.3399711
Stephen's Island-Maud Island  0.0105450345  0.008435768 0.0126543005 0.0000000
Stephen's Island-Motuara      0.0117123761  0.009723741 0.0137010112 0.0000000
```
And for inbreeding
```r
m2 <- aov(F ~ Pop, data = Het)
summary(m2)
TukeyHSD(m2)
```
```
> m2 <- aov(F ~ Pop, data = Het)
> summary(m2)
            Df Sum Sq Mean Sq F value Pr(>F)    
Pop          3 1.7041  0.5680   90.68 <2e-16 ***
Residuals   78 0.4886  0.0063                   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> TukeyHSD(m2)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = F ~ Pop, data = Het)

$Pop
                                    diff          lwr         upr     p adj
Maud Island-Boat Bay          0.01873000 -0.046976492  0.08443649 0.8770989
Motuara-Boat Bay              0.05704978 -0.004250202  0.11834976 0.0773821
Stephen's Island-Boat Bay    -0.33987533 -0.410846456 -0.26890421 0.0000000
Motuara-Maud Island           0.03831978 -0.022980202  0.09961976 0.3619248
Stephen's Island-Maud Island -0.35860533 -0.429576456 -0.28763421 0.0000000
Stephen's Island-Motuara     -0.39692511 -0.463837327 -0.33001290 0.0000000
```
## Create Jitterplots
```r

#Total Sites
total_sites <- ggplot(Het, aes(x = Pop, y = Prop.Het)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, colour = "cornflowerblue", alpha = 0.3) +
  theme_light() +
  labs(y= "Observed Heterozygosity", x = "Population")


#Inbreeding F(is)
Inbreeding <- ggplot(Het, aes(x = Pop, y = F)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, colour = "darkorange", alpha = 0.3) +
  theme_light() +
  ylab(expression("F"[is])) +
  xlab("Population")
  
```
```r
ggsave(filename="total_sites.png", plot = total_sites, dpi = 300 )
ggsave(filename = "inbreeding.png", plot = Inbreeding, dpi =300)
```
![](GeneticDiversity_Inbreeding/total_sites.png)<!-- -->
![](Hamilton-s-Frog-GBS/GeneticDiversity_Inbreeding/inbreeding.png)<!-- -->
