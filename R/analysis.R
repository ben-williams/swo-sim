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
og_st <- pop_est(lfreq, cpue, strata, by_strata = TRUE)

s20 <- sims(iters = 10, lfreq, cpue, strata, samples = 20, by_strata = TRUE)
s50 <- sims(iters = 10, lfreq, cpue, strata, samples = 50, by_strata = TRUE)


og_st %>%
  # dplyr::filter(species_code == "10110") %>%
  group_by(year, species_code, length) %>%
  summarise(males = mean(males),
            id = "og") -> m

s50 %>%
  # dplyr::filter(species_code == "10110") %>%
  group_by(year, species_code, length) %>%
  summarise(males = mean(males),
            id = "50") %>%
  bind_rows(m) %>%
  pivot_wider(values_from = males, names_from = id) %>%
  mutate(diff = og - `50`,
         clr = ifelse(diff < 0, "negative", "positive")) %>%
  ggplot(aes(year, length, size = diff, color = clr)) +
  geom_point(alpha = 0.2) +
  scale_size_area() +
  facet_wrap(~species_code, scales = "free") +
  funcr::theme_report()
