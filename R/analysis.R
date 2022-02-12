# compare decreased sample size length freq to original

# functions

library(tidyverse)
library(data.table)
library(vroom)
library(here)

# globals ----
region = "GOA"

# Set desired max sex sample sizes
max_sx <- 150 # seq(50, 150, 25)

# number of iterations
iters <- 10

# data ----

vroom::vroom(here::here("data", "strata.csv")) %>%
  filter(SURVEY==region) %>%
  rename_all(tolower) %>%
  dplyr::select(stratum, area) -> strata

vroom::vroom(here::here("data", "lfreq.csv")) %>%
  filter(SURVEY==region, STRATUM %in% unique(strata$stratum)) %>%
  rename_all(tolower) %>%
  dplyr::select(year, species_code, stratum, hauljoin, sex, length, frequency) -> lfreq

vroom::vroom(here::here("data", "sizepop.csv")) %>%
  filter(SURVEY==region) %>%
  rename_all(tolower) %>%
  dplyr::select(-1) -> race_pop

vroom::vroom(here::here("data", "CPUE.csv")) %>%
  filter(SURVEY==region) %>%
  rename_all(tolower) %>%
  left_join(strata) %>%
  dplyr::select(year, species_code, catchjoin, hauljoin, stratum, numcpue, area) -> cpue

# complete cases of unique lengths by species and sex
lfreq %>%
  group_by(species_code) %>%
  distinct(length, year, stratum) %>%
  expand(length, year, stratum) -> lngs


# estimate expanded length comps with full dataset
# also added the number of hauls and total number of lengths by haul
og <- pop_est(lfreq, cpue, lngs)
s20 <- replicate(10, pop_est(lfreq, cpue, lngs, samples = 20), simplify = FALSE)
s40 <- replicate(10, pop_est(lfreq, cpue, lngs, samples = 40), simplify = FALSE)
s60 <- replicate(10, pop_est(lfreq, cpue, lngs, samples = 60), simplify = FALSE)
s80 <- replicate(10, pop_est(lfreq, cpue, lngs, samples = 80), simplify = FALSE)
s100 <- replicate(10, pop_est(lfreq, cpue, lngs, samples = 100), simplify = FALSE)
s120 <- replicate(10, pop_est(lfreq, cpue, lngs, samples = 120), simplify = FALSE)



vroom::vroom_write(s20, here::here("test", "s20.csv"))


s20 = data.frame(a = 2)
