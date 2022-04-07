# compare decreased sample size length freq to original

# functions

library(tidyverse)
library(tidytable)
library(vroom)
library(here)
library(purrr)
library(rsample)
library(data.table)

source("R/functions_dev.R")

# globals ----
region = "AI"
boot = 500

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
# Bootstrap & sub-sampling simulation function

boot_res <- rerun(1, age_pop_est_boot(lfreq, specimen, cpue, sx_subsample = 100000, yrs = 2014, sz_resample = NULL, age_resample = NULL))

do.call(mapply, c(list, boot_res, SIMPLIFY = FALSE))$age %>%
  map_df.(., ~as.data.frame(.x), .id = "sim") -> og_age

do.call(mapply, c(list, boot_res, SIMPLIFY = FALSE))$length %>%
  map_df.(., ~as.data.frame(.x), .id = "sim") -> og_size

vroom::vroom_write(og_size, here::here("output", region, "og_size.csv"), delim = ",")
vroom::vroom_write(og_age, here::here("output", region, "og_age.csv"), delim = ",")

st <- Sys.time()
base <- age_sims_boot(lfreq, specimen, cpue, og_age, og_size, sx_subsample = 10000, yrs = 2014, sz_resample = TRUE, age_resample = TRUE, boot = boot, write_comp = TRUE, save = "base", region = tolower(region))
s40 <- age_sims_boot(lfreq, specimen, cpue, og_age, og_size, sx_subsample = 10000, yrs = 2014, sz_resample = TRUE, age_resample = TRUE, boot = boot, write_comp = TRUE, save = "s40", region = tolower(region))
s60 <- age_sims_boot(lfreq, specimen, cpue, og_age, og_size, sx_subsample = 10000, yrs = 2014, sz_resample = TRUE, age_resample = TRUE, boot = boot, write_comp = TRUE, save = "s60", region = tolower(region))
s80 <- age_sims_boot(lfreq, specimen, cpue, og_age, og_size, sx_subsample = 10000, yrs = 2014, sz_resample = TRUE, age_resample = TRUE, boot = boot, write_comp = TRUE, save = "s80", region = tolower(region))
s100 <- age_sims_boot(lfreq, specimen, cpue, og_age, og_size, sx_subsample = 10000, yrs = 2014, sz_resample = TRUE, age_resample = TRUE, boot = boot, write_comp = TRUE, save = "s100", region = tolower(region))
s120 <- age_sims_boot(lfreq, specimen, cpue, og_age, og_size, sx_subsample = 10000, yrs = 2014, sz_resample = TRUE, age_resample = TRUE, boot = boot, write_comp = TRUE, save = "s120", region = tolower(region))
s140 <- age_sims_boot(lfreq, specimen, cpue, og_age, og_size, sx_subsample = 10000, yrs = 2014, sz_resample = TRUE, age_resample = TRUE, boot = boot, write_comp = TRUE, save = "s140", region = tolower(region))
end <- Sys.time()

run_time <- end-st

