# compare decreased sample size length freq to original

# functions

library(tidyverse)
library(data.table)
library(vroom)
library(here)

source("R/functions.R")

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



# estimate expanded length comps with full dataset
# also added the number of hauls and total number of lengths by haul
og <- sims(iters = 1, lfreq, cpue)


# annual abundance by year & species
ptm = proc.time()
s20 <- sims(iters = 40, lfreq, cpue, samples = 20, strata = NULL)
proc.time() - ptm
# ~ 4 min
s40 <- sims(iters = 40, lfreq, cpue, samples = 40, strata = NULL)
s60 <- sims(iters = 40, lfreq, cpue, samples = 60, strata = NULL)
s80 <- sims(iters = 40, lfreq, cpue, samples = 80, strata = NULL)
s100 <- sims(iters = 40, lfreq, cpue, samples = 100, strata = NULL)
s120 <- sims(iters = 40, lfreq, cpue, samples = 120, strata = NULL)

getouts(og, save = "og.csv")
getouts(s20, save = "s20.csv")
getouts(s40, save = "s40.csv")
getouts(s60, save = "s60.csv")
getouts(s80, save = "s80.csv")
getouts(s100, save = "s100.csv")
getouts(s120, save = "s120.csv")

# annual abundance by year & strata & species
og_st <- sims(iters = 1, lfreq, cpue, samples = NULL, strata = TRUE)
ptm = proc.time()
s20_st <- sims(iters = 40, lfreq, cpue, samples = 20, strata = TRUE)
proc.time() - ptm
s40_st <- sims(iters = 40, lfreq, cpue, samples = 40, strata = TRUE)
s60_st <- sims(iters = 40, lfreq, cpue, samples = 60, strata = TRUE)
s80_st <- sims(iters = 40, lfreq, cpue, samples = 80, strata = TRUE)
s100_st <- sims(iters = 40, lfreq, cpue, samples = 100, strata = TRUE)
s120_st <- sims(iters = 40, lfreq, cpue, samples = 120, strata = TRUE)

# save strata-based comps
getouts(og_st, save = "og_st.csv")
getouts(s20_st, samples = TRUE, save = "s20_comp.csv")
getouts(s40_st, samples = TRUE, save = "s40_comp.csv")
getouts(s60_st, samples = TRUE, save = "s60_comp.csv")
getouts(s80_st, samples = TRUE, save = "s80_comp.csv")
getouts(s100_st, samples = TRUE, save = "s100_comp.csv")
getouts(s120_st, samples = TRUE, save = "s120_comp.csv")

# save strata-
getouts(s20_st, type = "removed", samples = TRUE, save = "s20_removed.csv")
getouts(s40_st, type = "removed", samples = TRUE, save = "s40_removed.csv")
getouts(s60_st, type = "removed", samples = TRUE, save = "s60_removed.csv")
getouts(s80_st, type = "removed", samples = TRUE, save = "s80_removed.csv")
getouts(s100_st, type = "removed", samples = TRUE, save = "s100_removed.csv")
getouts(s120_st, type = "removed", samples = TRUE, save = "s120_removed.csv")
