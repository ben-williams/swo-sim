# compare decreased sample size length freq to original

# functions

library(tidyverse)
library(data.table)
library(vroom)
library(here)

source("R/functions.R")

# globals ----
region = "GOA"
iters = 10
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


# annual abundance by year & species
og <- sims(iters = 1, lfreq, cpue, save = "og")
s20 <- sims(iters = iters, lfreq, cpue, samples = 20, strata = NULL, save = "s2")
s40 <- sims(iters = iters, lfreq, cpue, samples = 40, strata = NULL, save = "s40")
s60 <- sims(iters = iters, lfreq, cpue, samples = 60, strata = NULL, save = "s60")
s80 <- sims(iters = iters, lfreq, cpue, samples = 80, strata = NULL, save = "s80")
s100 <- sims(iters = iters, lfreq, cpue, samples = 100, strata = NULL, save = "s100")
s120 <- sims(iters = iters, lfreq, cpue, samples = 120, strata = NULL, save = "s120")
s140 <- sims(iters = iters, lfreq, cpue, samples = 140, strata = NULL, save = "s140")


# annual abundance by year & strata & species
og_st <- sims(iters = 1, lfreq, cpue, strata = TRUE, save = "og_st")
s20_st <- sims(iters = iters, lfreq, cpue, samples = 20, strata = TRUE, save = "s20_st")
s40_st <- sims(iters = iters, lfreq, cpue, samples = 40, strata = TRUE, save = "s40_st")
s60_st <- sims(iters = iters, lfreq, cpue, samples = 60, strata = TRUE, save = "s60_st")
s80_st <- sims(iters = iters, lfreq, cpue, samples = 80, strata = TRUE, save = "s80_st")
s100_st <- sims(iters = iters, lfreq, cpue, samples = 100, strata = TRUE, save = "s100_st")
s120_st <- sims(iters = iters, lfreq, cpue, samples = 120, strata = TRUE, save = "s120_st")
s140_st <- sims(iters = iters, lfreq, cpue, samples = 140, strata = TRUE, save = "s140_st")


# different no. of iterations
s8020 <- sims(iters = 20, lfreq, cpue, samples = 80, strata = NULL, save = "s8020")
s8040 <- sims(iters = 40, lfreq, cpue, samples = 80, strata = NULL, save = "s8040")
s8060 <- sims(iters = 60, lfreq, cpue, samples = 80, strata = NULL, save = "s8060")
s8080 <- sims(iters = 80, lfreq, cpue, samples = 80, strata = NULL, save = "s8080")
s80100 <- sims(iters = 100, lfreq, cpue, samples = 80, strata = NULL, save = "s80100")
s80120 <- sims(iters = 120, lfreq, cpue, samples = 80, strata = NULL, save = "s80120")


# more years ----
s80yr <- sims(iters = 10, lfreq, cpue, samples = 80, strata = NULL, save = "s80yr", yrs = 10)
