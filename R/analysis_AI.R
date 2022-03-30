# compare decreased sample size length freq to original

# functions

library(tidyverse)
library(tidytable)
library(vroom)
library(here)
library(purrr)
library(rsample)
library(data.table)

# globals ----
region = "AI"
iters = 100
boot = 10

# data ----

vroom::vroom(here::here("data", "strata.csv")) %>%
  filter(SURVEY == region) %>%
  rename_all(tolower) %>%
  dplyr::select(stratum, area) -> strata

vroom::vroom(here::here("data", "lfreq.csv")) %>%
  filter(SURVEY == region, STRATUM %in% unique(strata$stratum)) %>%
  rename_all(tolower) %>%
  dplyr::select(year, species_code, stratum, hauljoin, sex, length, frequency) -> lfreq # nolint

vroom::vroom(here::here("data", "sizepop.csv")) %>%
  filter(SURVEY == region) %>%
  rename_all(tolower) %>%
  dplyr::select(-1) -> race_pop

vroom::vroom(here::here("data", "CPUE.csv")) %>%
  filter(SURVEY == region) %>%
  rename_all(tolower) %>%
  left_join(strata) %>%
  dplyr::select(year, species_code, catchjoin, hauljoin, stratum, numcpue, area) -> cpue # nolint

vroom::vroom(here::here("data", "specimen.csv")) %>%
  filter(REGION == region, STRATUM %in% unique(strata$stratum)) %>%
  rename_all(tolower) %>%
  dplyr::select(year, species_code, stratum, hauljoin, sex, length, age) %>%
  dplyr::filter(!is.na(age)) -> specimen


###############################################################################################################
# Bootstrap & sub-sampling simulation function calls
source("R/functions.R")

rerun(1, size_pop_est(lfreq, cpue, samples = 10000, yrs = 2018)) %>%
  map_df.(., ~as.data.frame(.x), .id = "sim")  -> og_size
vroom::vroom_write(og_size, here::here("output", "bootstrap", region, "og_size.csv"), delim = ",")

base <- size_sims_boot(lfreq, cpue, samples = 10000, yrs = 2018, strata = NULL, og_data = og_size, boot = boot, write_comp = TRUE, save = "base", region = tolower(region))
s40 <- size_sims_boot(lfreq, cpue, samples = 40, yrs = 2018, strata = NULL, og_data = og_size, boot = boot, write_comp = TRUE, save = "s40", region = tolower(region))
s60 <- size_sims_boot(lfreq, cpue, samples = 60, yrs = 2018, strata = NULL, og_data = og_size, boot = boot, write_comp = TRUE, save = "s60", region = tolower(region))
s80 <- size_sims_boot(lfreq, cpue, samples = 80, yrs = 2018, strata = NULL, og_data = og_size, boot = boot, write_comp = TRUE, save = "s80", region = tolower(region))
s100 <- size_sims_boot(lfreq, cpue, samples = 100, yrs = 2018, strata = NULL, og_data = og_size, boot = boot, write_comp = TRUE, save = "s100", region = tolower(region))
s120 <- size_sims_boot(lfreq, cpue, samples = 120, yrs = 2018, strata = NULL, og_data = og_size, boot = boot, write_comp = TRUE, save = "s120", region = tolower(region))
s140 <- size_sims_boot(lfreq, cpue, samples = 140, yrs = 2018, strata = NULL, og_data = og_size, boot = boot, write_comp = TRUE, save = "s140", region = tolower(region))


###############################################################################################################
# Sub-sampling function calls

# annual abundance at size by year & species
#og <- sims(iters = 1, lfreq, cpue, yrs = 2006, save = "og_size")
#s20 <- sims(iters = iters, lfreq, cpue, samples = 20, yrs = 2006, strata = NULL, save = "s20", removed_summ=TRUE)
#s40 <- sims(iters = iters, lfreq, cpue, samples = 40, yrs = 2006, strata = NULL, save = "s40", removed_summ=TRUE)
#s60 <- sims(iters = iters, lfreq, cpue, samples = 60, yrs = 2006, strata = NULL, save = "s60", removed_summ=TRUE)
#s80 <- sims(iters = iters, lfreq, cpue, samples = 80, yrs = 2006, strata = NULL, save = "s80", removed_summ=TRUE)
#s100 <- sims(iters = iters, lfreq, cpue, samples = 100, yrs = 2006, strata = NULL, save = "s100", removed_summ=TRUE)
#s120 <- sims(iters = iters, lfreq, cpue, samples = 120, yrs = 2006, strata = NULL, save = "s120", removed_summ=TRUE)
#s140 <- sims(iters = iters, lfreq, cpue, samples = 140, yrs = 2006, strata = NULL, save = "s140", removed_summ=TRUE)

# annual abundance at age by year & species
#og_age <- age_pop_est(specimen, og, yrs = 2006, sim = NULL,  save = "og_age")
#s20_age <- age_pop_est(specimen, s20, yrs = 2006, sim = TRUE,  save = "s20_age")
#s40_age <- age_pop_est(specimen, s40, yrs = 2006, sim = TRUE,  save = "s40_age")
#s60_age <- age_pop_est(specimen, s60, yrs = 2006, sim = TRUE,  save = "s60_age")
#s80_age <- age_pop_est(specimen, s80, yrs = 2006, sim = TRUE,  save = "s80_age")
#s100_age <- age_pop_est(specimen, s100, yrs = 2006, sim = TRUE,  save = "s100_age")
#s120_age <- age_pop_est(specimen, s120, yrs = 2006, sim = TRUE,  save = "s120_age")
#s140_age <- age_pop_est(specimen, s140, yrs = 2006, sim = TRUE,  save = "s140_age")

# annual abundance at size by year & strata & species
#og_st <- sims(iters = 1, lfreq, cpue, strata = TRUE, yrs = 2006, save = "og_st")
#s20_st <- sims(iters = iters, lfreq, cpue, samples = 20, yrs = 2006, strata = TRUE, save = "s20_st", removed_summ=TRUE)
#s40_st <- sims(iters = iters, lfreq, cpue, samples = 40, yrs = 2006, strata = TRUE, save = "s40_st", removed_summ=TRUE)
#s60_st <- sims(iters = iters, lfreq, cpue, samples = 60, yrs = 2006, strata = TRUE, save = "s60_st", removed_summ=TRUE)
#s80_st <- sims(iters = iters, lfreq, cpue, samples = 80, yrs = 2006, strata = TRUE, save = "s80_st", removed_summ=TRUE)
#s100_st <- sims(iters = iters, lfreq, cpue, samples = 100, yrs = 2006, strata = TRUE, save = "s100_st", removed_summ=TRUE)
#s120_st <- sims(iters = iters, lfreq, cpue, samples = 120, yrs = 2006, strata = TRUE, save = "s120_st", removed_summ=TRUE)
#s140_st <- sims(iters = iters, lfreq, cpue, samples = 140, yrs = 2006, strata = TRUE, save = "s140_st", removed_summ=TRUE)

# Compute size comp ess at global and strata level
#ess_s20_size <- ess_size(s20, og, strata = NULL, "s20_ess_size")
#ess_s40_size <- ess_size(s40, og, strata = NULL, "s40_ess_size")
#ess_s60_size <- ess_size(s60, og, strata = NULL, "s60_ess_size")
#ess_s80_size <- ess_size(s80, og, strata = NULL, "s80_ess_size")
#ess_s100_size <- ess_size(s100, og, strata = NULL, "s100_ess_size")
#ess_s120_size <- ess_size(s120, og, strata = NULL, "s120_ess_size")
#ess_s140_size <- ess_size(s140, og, strata = NULL, "s140_ess_size")
#ess_s20_st_size <- ess_size(s20_st, og_st, strata = TRUE, "s20_st_ess_size")
#ess_s40_st_size <- ess_size(s40_st, og_st, strata = TRUE, "s40_st_ess_size")
#ess_s60_st_size <- ess_size(s60_st, og_st, strata = TRUE, "s60_st_ess_size")
#ess_s80_st_size <- ess_size(s80_st, og_st, strata = TRUE, "s80_st_ess_size")
#ess_s100_st_size <- ess_size(s100_st, og_st, strata = TRUE, "s100_st_ess_size")
#ess_s120_st_size <- ess_size(s120_st, og_st, strata = TRUE, "s120_st_ess_size")
#ess_s140_st_size <- ess_size(s140_st, og_st, strata = TRUE, "s140_st_ess_size")

# Compute age comp ess at global level
#ess_s20_age <- ess_age(s20_age, og_age, save = "s20_ess_age")
#ess_s40_age <- ess_age(s40_age, og_age, save = "s40_ess_age")
#ess_s60_age <- ess_age(s60_age, og_age, save = "s60_ess_age")
#ess_s80_age <- ess_age(s80_age, og_age, save = "s80_ess_age")
#ess_s100_age <- ess_age(s100_age, og_age, save = "s100_ess_age")
#ess_s120_age <- ess_age(s120_age, og_age, save = "s120_ess_age")
#ess_s140_age <- ess_age(s140_age, og_age, save = "s140_ess_age")

# Compute proportion of females to males from size pop'n data
#prop_fm_s20 <- prop_fm(s20, og, strata = NULL, "s20_prop_fm")
#prop_fm_s40 <- prop_fm(s40, og, strata = NULL, "s40_prop_fm")
#prop_fm_s60 <- prop_fm(s60, og, strata = NULL, "s60_prop_fm")
#prop_fm_s80 <- prop_fm(s80, og, strata = NULL, "s80_prop_fm")
#prop_fm_s100 <- prop_fm(s100, og, strata = NULL, "s100_prop_fm")
#prop_fm_s120 <- prop_fm(s120, og, strata = NULL, "s120_prop_fm")
#prop_fm_s140 <- prop_fm(s140, og, strata = NULL, "s140_prop_fm")
#prop_fm_s20_st <- prop_fm(s20_st, og_st, strata = TRUE, "s20_st_prop_fm")
#prop_fm_s40_st <- prop_fm(s40_st, og_st, strata = TRUE, "s40_st_prop_fm")
#prop_fm_s60_st <- prop_fm(s60_st, og_st, strata = TRUE, "s60_st_prop_fm")
#prop_fm_s80_st <- prop_fm(s80_st, og_st, strata = TRUE, "s80_st_prop_fm")
#prop_fm_s100_st <- prop_fm(s100_st, og_st, strata = TRUE, "s100_st_prop_fm")
#prop_fm_s120_st <- prop_fm(s120_st, og_st, strata = TRUE, "s120_st_prop_fm")
#prop_fm_s140_st <- prop_fm(s140_st, og_st, strata = TRUE, "s140_st_prop_fm")

# different no. of iterations
#s8020 <- sims(iters = 20, lfreq, cpue, samples = 80, strata = NULL, save = "s8020")
#s8040 <- sims(iters = 40, lfreq, cpue, samples = 80, strata = NULL, save = "s8040")
#s8060 <- sims(iters = 60, lfreq, cpue, samples = 80, strata = NULL, save = "s8060")
#s8080 <- sims(iters = 80, lfreq, cpue, samples = 80, strata = NULL, save = "s8080")
#s80100 <- sims(iters = 100, lfreq, cpue, samples = 80, strata = NULL, save = "s80100")
#s80120 <- sims(iters = 120, lfreq, cpue, samples = 80, strata = NULL, save = "s80120")

# more years ----
#s80yr <- sims(iters = 10, lfreq, cpue, samples = 80, strata = NULL, save = "s80yr", yrs = 10)
