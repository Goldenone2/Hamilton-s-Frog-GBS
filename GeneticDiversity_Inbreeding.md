# Genetic Diversity and Inbreeding
I'd like to compare heterozygosity and, inbreeding across the four populations I've sequenced. 

## Data Import and Setup
```{r}
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
```{r}
m1 <- aov(Prop.Het ~ Pop, data = Het)
summary(m1)
TukeyHSD(m1)

```
And for inbreeding
```{r}
m2 <- aov(F ~ Pop, data = Het)
summary(m2)
TukeyHSD(m2)
```
## Create Jitterplots
```{r}

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
```{r}
ggsave(filename="total_sites.png", plot = total_sites, dpi = 300 )
ggsave(filename = "inbreeding.png", plot = Inbreeding, dpi =300)
```
