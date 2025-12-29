# Stairway Plot

### Data Import and Setup

``` r
# Clear workspace
rm(list = ls())

# Load required packages
library(ggplot2)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(scales)
```

\###Read STWP2 results

``` r
Takapourewa <- read.table("Folded_SFS Takapourewa.rand7.summary", header = T)

Takapourewa_scaled <- Takapourewa %>% mutate(
  yeark = year / 1000,
    Ne_medianK = Ne_median / 1000,
    Ne_2.5.K = Ne_2.5. / 1000,
    Ne_97.5.K  = Ne_97.5. / 1000) 

TePakeka <- read.table("Folded SFS Maud Island.rand9.summary", header = TRUE)

TePakeka_scaled <- TePakeka %>%
  mutate(
    yeark       = year / 1000,
    Ne_medianK  = Ne_median / 1000,
    Ne_2.5.K    = Ne_2.5. / 1000,
    Ne_97.5.K   = Ne_97.5. / 1000)

Together <- bind_rows(
  TePakeka_scaled %>% mutate(pop = "Te Pākeka"),
  Takapourewa_scaled %>% mutate(pop = "Takapourewa"))
```

# Plotting

Geom_ribbon is creating the 95% confidence interval;Geom_Step is
creating the Ne median line. I’ve borrowed the log scales limits from
Stairway Plots native plots.

``` r
Stairway <- ggplot(Together, aes(x = yeark)) +
  geom_ribbon(aes(ymin = Ne_2.5.K, ymax = Ne_97.5.K, fill = pop), alpha = 0.1) +
  geom_step(aes(y = Ne_medianK, colour = pop), linewidth = 0.8) +
  # I want to avoid scientific notation in this plot
  scale_x_log10(limits = c(0.001, 20),
                breaks = c(0.001, 0.01, 0.1, 1, 10),
                labels = c("0.001", "0.01", "0.1", "1", "10")) +
  scale_y_log10(limits = c(0.08, 20)) +
  labs(x = "Time (thousands of years ago)", y = "Effective Population Size (Thousands of individuals)",
       colour = "Population", fill = "Population") +
    theme_light() +
    scale_colour_manual(values = c("cornflowerblue", "darkorange")) +
    scale_fill_manual(values = c("cornflowerblue", "darkorange"))
```

### Chronology of the Marlborough Sounds

I want to add vertical lines to represent the key events in the
chronology of the Marlborough Sounds:

stabilization of sea levels (7 Kya), arrival of Maori (0.75 Kya),
colonisation of Europeans (0.204 Kya and, the introduction of cats to
Takapourewa (0.056 Kya).

``` r
# create data frame
events <- data.frame(
  label = c("A","B","C","D"),
  x = c(7,0.75,0.204,0.056))


Stairway2 <- Stairway +
  geom_vline(data = events, aes(xintercept = x), colour = "red",
              linewidth = 0.5, linetype = "dashed") +
  geom_text(data = events, aes(x =x*1.4, y = 12, label = label),
            colour = "red", size = 3)
              
Stairway2
```

    ## Warning: Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_step()`).

![](Staiway_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
 # ggsave(filename="Stairway.png", plot = Stairway2, dpi = 300, width = 8, height = 6 )
```

