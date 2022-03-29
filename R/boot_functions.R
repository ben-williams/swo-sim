library(tidyverse)
library(tidytable)
library(vroom)
library(here)
library(purrr)
library(rsample)
library(data.table)
ggplot2::theme_set(
  ggplot2::theme_light() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      # axis.ticks.length = grid::unit(base_ / 2.2, "pt"),
      strip.background = ggplot2::element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(colour = "black"),
      strip.text.y = element_text(colour = "black"),
      panel.border = element_rect(fill = NA),
      legend.key.size = grid::unit(0.9, "lines"),
      legend.key = ggplot2::element_rect(colour = NA, fill = NA),
      legend.background = ggplot2::element_rect(colour = NA, fill = NA)
    )
)


# globals ----
region <- "AI"

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

  pop_est <- function(lfreq, cpue, samples = 10000, yrs = NULL, strata = NULL){
    
    # year switch
    if (is.null(yrs)) yrs <- 0
    
    # complete cases of unique lengths by species and sex (for all years)
    
    lfreq %>%
      group_by(species_code) %>%
      distinct(length, year, stratum) %>%
      expand(length, year, stratum) %>%
      filter(year >= yrs) -> lngs
    
    # if no strata remove from
    if (is.null(strata)) {
      lngs %>%
        dplyr::select(-stratum) %>%
        group_by(species_code) %>%
        distinct(length, year) %>%
        expand(length, year) -> lngs
    }
    
    # first pass of filtering
    setDT(cpue) %>%
      filter.(year >= yrs) -> .cpue
    
    setDT(lfreq) %>%
      filter.(year >= yrs) %>% 
      uncount.(., frequency) -> .lfreq
    
    # pull out males and females, make single row for each, add an id
    
    .lfreq %>% 
      filter.(sex != 3) %>% 
      mutate.(id = .I) %>% 
      mutate.(n = .N, .by = c(year, species_code, stratum, hauljoin)) -> .inter
    
    # sample by sample size 
    .inter %>%
      group_by(year, species_code, stratum, hauljoin) %>%
      sample_n(if(n > samples) samples else n) -> .new_sexed 
    
    # find new unsexed, that were previously sexed
    .inter %>%
      anti_join.(.new_sexed, by = "id") %>%
      mutate.(sex = 3) %>% 
      count.(c(year, species_code, stratum, hauljoin,   sex, length), name = "frequency") -> .new_unsexed
    
    # rejoin original unsexed to the new_sexed samples
    .lfreq %>% 
      filter.(sex == 3) %>% 
      mutate.(id = .I,
              n = .N, .by = c(year, species_code, stratum, hauljoin))  %>% 
      bind_rows.(.new_sexed, .new_unsexed) %>% 
      summarise.(frequency = .N, 
                 .by = c(year, species_code, stratum, hauljoin, sex, length))  %>% 
      mutate.(nhauls = uniqueN(hauljoin), 
              .by = c(year, species_code, stratum)) %>% 
      mutate.(tot := sum(frequency), 
              .by = c(year, species_code, stratum, hauljoin)) %>% 
      summarise.(comp = sum(frequency) / mean(tot),
                 nhauls = mean(nhauls), 
                 .by = c(year, species_code, stratum, hauljoin, sex, length))  -> .lcomp
    
    
    # estimate for hauls w/o length samples
    .lcomp %>%
      summarise.(comp = sum(comp) / mean(nhauls), .by = c(year, species_code, stratum, sex, length)) -> .unk
    
    # id hauls without lengths
    .cpue %>%
      filter.(!is.na(catchjoin), !(hauljoin %in% .lcomp$hauljoin)) -> .no_length
    
    .cpue %>% 
      mutate.(st_num = mean(numcpue) * area,
              tot = sum(numcpue), .by = c(year, species_code, stratum)) %>% 
      summarise.(abund = mean(numcpue) / tot * st_num, 
                 .by = c(year, species_code, stratum, hauljoin)) -> .pop
    
    # if there are any samples w/o lengths rejoin them
    if(nrow(.no_length) == 0){
      .lcomp %>%
        left_join.(.pop) %>%
        mutate.(sz_pop = round(comp * abund, 0)) -> .temp
    } else {
      .no_length %>%
        left_join.(.unk) %>%
        select.(year, species_code, stratum, hauljoin, sex, length, comp) %>%
        bind_rows.(.lcomp) %>%
        left_join.(.pop) %>%
        mutate.(sz_pop = round(comp * abund, 0)) -> .temp
    }
    
    if(!is.null(strata)){
      .temp %>%
        summarise.(abund = sum(sz_pop, na.rm = T), .by = c(year, species_code, stratum, length, sex)) %>%
        pivot_wider.(names_from = sex, values_from = abund) %>%
        left_join.(lngs, .) %>%
        mutate.(across.(.cols = c(`1`, `2`, `3`), ~replace_na.(.x, 0))) %>%
        select.(year, species_code, stratum, length, males = `1`, females = `2`, unsexed = `3`) -> .out
      
    } else {
      .temp %>%
        select.(-stratum) %>% 
        summarise.(abund = sum(sz_pop, na.rm = T), .by = c(year, species_code, length, sex)) %>%
        pivot_wider.(names_from = sex, values_from = abund) %>%
        left_join.(lngs, .) %>%
        mutate.(across.(.cols = c(`1`, `2`, `3`), ~replace_na.(.x, 0))) %>%
        select.(year, species_code, length, males = `1`, females = `2`, unsexed = `3`) -> .out
    }
    
    if(!is.null(samples)){
      # list(new = .out, removed = .new_unsexed)
      .out
    } else {
      .out
    }
    
  }
pop_est_boot <- function(lfreq, cpue, samples = 10000, yrs = NULL, strata = NULL){

  # updated pop estimate with haul/length bootstraps (w/resampling)

  # year switch
  if (is.null(yrs)) yrs <- 0

  # complete cases of unique lengths by species and sex (for all years)
  lfreq %>%
    group_by(species_code) %>%
    distinct(length, year, stratum) %>%
    expand(length, year, stratum) %>%
    filter(year >= yrs) -> lngs

  # if no strata remove from
  if (is.null(strata)) {
    lngs %>%
      dplyr::select(-stratum) %>%
      group_by(species_code) %>%
      distinct(length, year) %>%
      expand(length, year) -> lngs
  }

  # first pass of filtering
  setDT(cpue) %>%
    filter.(year >= yrs) -> .cpue

  # bootstrap the hauls and the lengths
  setDT(lfreq) %>%
    filter.(year >= yrs) %>%
    uncount.(., frequency) %>%
    .[, hauljoin := hauljoin[sample.int(.N, .N, replace = TRUE)],
      by = c('year', 'species_code', 'stratum') ] %>%
    .[, length := length[sample.int(.N, .N, replace = TRUE)],
      by = c('year', 'species_code', 'stratum', 'hauljoin') ] -> .lfreq


  # pull out males and females, make single row for each, add an id

  .lfreq %>%
    filter.(sex != 3) %>%
    mutate.(id = .I) %>%
    mutate.(n = .N, .by = c(year, species_code, stratum, hauljoin)) -> .inter

  # sample by sample size
  .inter %>%
    group_by(year, species_code, stratum, hauljoin) %>%
    sample_n(if(n > samples) samples else n) -> .new_sexed

  # find new unsexed, that were previously sexed
  .inter %>%
    anti_join.(.new_sexed, by = "id") %>%
    mutate(sex = 3) -> .new_unsexed

  # rejoin original unsexed to the new_sexed samples
  .lfreq %>%
    filter.(sex == 3) %>%
    mutate.(id = .I,
            n = .N, .by = c(year, species_code, stratum, hauljoin))  %>%
    bind_rows.(.new_sexed, .new_unsexed) %>%
    summarise.(frequency = .N,
               .by = c(year, species_code, stratum, hauljoin, sex, length))  %>%
    mutate.(nhauls = uniqueN(hauljoin),
            .by = c(year, species_code, stratum)) %>%
    mutate.(tot := sum(frequency),
            .by = c(year, species_code, stratum, hauljoin)) %>%
    summarise.(comp = sum(frequency) / mean(tot),
               nhauls = mean(nhauls),
               .by = c(year, species_code, stratum, hauljoin, sex, length))  -> .lcomp


  # estimate for hauls w/o length samples
  .lcomp %>%
    summarise.(comp = sum(comp) / mean(nhauls), .by = c(year, species_code, stratum, sex, length)) -> .unk

  # id hauls without lengths
  .cpue %>%
    filter.(!is.na(catchjoin), !(hauljoin %in% .lcomp$hauljoin)) -> .no_length

  .cpue %>%
    mutate.(st_num = mean(numcpue) * area,
            tot = sum(numcpue), .by = c(year, species_code, stratum)) %>%
    summarise.(abund = mean(numcpue) / tot * st_num,
               .by = c(year, species_code, stratum, hauljoin)) -> .pop

  # if there are any samples w/o lengths rejoin them
  if(nrow(.no_length) == 0){
    .lcomp %>%
      left_join.(.pop) %>%
      mutate.(sz_pop = round(comp * abund, 0)) -> .temp
  } else {
    .no_length %>%
      left_join.(.unk) %>%
      select.(year, species_code, stratum, hauljoin, sex, length, comp) %>%
      bind_rows.(.lcomp) %>%
      left_join.(.pop) %>%
      mutate.(sz_pop = round(comp * abund, 0)) -> .temp
  }

  if(!is.null(strata)){
    .temp %>%
      summarise.(abund = sum(sz_pop, na.rm = T), .by = c(year, species_code, stratum, length, sex)) %>%
      pivot_wider.(names_from = sex, values_from = abund) %>%
      left_join.(lngs, .) %>%
      mutate.(across.(.cols = c(`1`, `2`, `3`), ~replace_na.(.x, 0))) %>%
      select.(year, species_code, stratum, length, males = `1`, females = `2`, unsexed = `3`) -> .out

  } else {
    .temp %>%
      select.(-stratum) %>%
      summarise.(abund = sum(sz_pop, na.rm = T), .by = c(year, species_code, length, sex)) %>%
      pivot_wider.(names_from = sex, values_from = abund) %>%
      left_join.(lngs, .) %>%
      mutate.(across.(.cols = c(`1`, `2`, `3`), ~replace_na.(.x, 0))) %>%
      select.(year, species_code, length, males = `1`, females = `2`, unsexed = `3`) -> .out
  }

  if(!is.null(samples)){
    # list(new = .out, removed = .new_unsexed)
    .out
  } else {
    .out
  }

}
ess_size <- function(sim_data, og_data, strata = NULL) {

  if ("stratum" %in% names(og_data) & is.null(strata) |
      "stratum" %in% names(sim_data) & is.null(strata)) {
    stop("check your strata")
  }

  if (!is.null(strata)) {
    og_data %>%
      mutate.(og_m = males / sum(males),
              og_f = females / sum(females),
              og_t = (males + females + unsexed)/ (sum(males) + sum(females) + sum(unsexed)),
              .by = c(year, species_code, stratum)) %>%
      filter.(og_f >= 0 &  og_m >= 0 & og_t >= 0) %>%
      select.(year, species_code, stratum, length, og_m, og_f, og_t) -> og_prop


    sim_data %>%
      mutate.(prop_m = males / sum(males),
              prop_f = females / sum(females),
              prop_t = (males + females + unsexed)/(sum(males) + sum(females) + sum(unsexed)),
              .by = c(year, species_code, stratum)) %>%
      filter.(prop_m >= 0 & prop_f >= 0 & prop_t >= 0) %>%
      left_join.(og_prop) %>%
      summarise.(ess_f = sum(prop_f * (1 - prop_f)) / sum((prop_f - og_f)^2),
                 ess_m = sum(prop_m * (1 - prop_m)) / sum((prop_m - og_m)^2),
                 ess_t = sum(prop_t * (1 - prop_t)) / sum((prop_t - og_t)^2)) %>%
      drop_na() %>%
      pivot_longer.(cols = c(ess_f, ess_m, ess_t), names_to = "ess") %>%
      mutate.(in_out = ifelse(is.infinite(value), "out", "in")) %>%
      group_by(year, species_code, stratum, ess, value, in_out) %>%
      distinct(value)

  } else {

    og_data %>%
      mutate.(og_m = males / sum(males),
              og_f = females / sum(females),
              og_t = (males + females + unsexed)/ (sum(males) + sum(females) + sum(unsexed)),
              .by = c(year, species_code)) %>%
      filter.(og_f >= 0 &  og_m >= 0 & og_t >= 0) %>%
      select.(year, species_code, length, og_m, og_f, og_t) -> og_prop

    sim_data %>%
      mutate.(prop_m = males / sum(males),
              prop_f = females / sum(females),
              prop_t = (males + females + unsexed)/(sum(males) + sum(females) + sum(unsexed)),
              .by = c(year, species_code)) %>%
      filter.(prop_m >= 0 & prop_f >= 0 & prop_t >= 0) %>%
      left_join.(og_prop) %>%
      mutate.(ess_f = sum(prop_f * (1 - prop_f)) / sum((prop_f - og_f)^2),
              ess_m = sum(prop_m * (1 - prop_m)) / sum((prop_m - og_m)^2),
              ess_t = sum(prop_t * (1 - prop_t)) / sum((prop_t - og_t)^2),
              .by = c(year, species_code)) %>%
      drop_na() %>%
      pivot_longer.(cols = c(ess_f, ess_m, ess_t), names_to = "ess") %>%
      mutate.(in_out = ifelse(is.infinite(value), "out", "in")) %>%
      group_by(year, species_code, ess, value, in_out) %>%
      distinct(value)

  }

}
sims_boot <- function(lfreq, cpue, samples = NULL, yrs = NULL, strata = NULL, og_data, boot = 500){

  rerun(boot, pop_est_boot(lfreq, cpue, samples = samples, yrs = yrs)) %>%
    map.(., ~ess_size(sim_data = .x, og_data = og_data)) %>%
    map_df.(., ~as.data.frame(.x), .id = "sim")

}

# compute data to compare, note samples = 10000 eliminates any of the simulation component
#
rerun(1, pop_est(lfreq, cpue, samples = 10000, yrs = 2006)) %>%
  map_df.(., ~as.data.frame(.x), .id = "sim")  -> og

sims_boot(lfreq, cpue, samples = 10000, yrs = 2006, strata = NULL, og_data = og, boot = 100) -> base_ai
# vroom::vroom_write(s40_ai, here::here("output", "base_ai.csv"), delim = ",")

sims_boot(lfreq, cpue, samples = 40, yrs = 2006, strata = NULL, og_data = og, boot = 100) -> s40_aib
# vroom::vroom_write(s40_ai, here::here("output", "s40_aib.csv"), delim = ",")
sims_boot(lfreq, cpue, samples = 60, yrs = 2006, strata = NULL, og_data = og, boot = 100) -> s60_aib
# vroom::vroom_write(s40_ai, here::here("output", "s46_aib.csv"), delim = ",")
sims_boot(lfreq, cpue, samples = 80, yrs = 2006, strata = NULL, og_data = og, boot = 100) -> s80_aib
# vroom::vroom_write(s40_ai, here::here("output", "s8_aib.csv"), delim = ",")
sims_boot(lfreq, cpue, samples = 100, yrs = 2006, strata = NULL, og_data = og, boot = 100) -> s100_aib
# vroom::vroom_write(s40_ai, here::here("output", "s100_aib.csv"), delim = ",")
sims_boot(lfreq, cpue, samples = 120, yrs = 2006, strata = NULL, og_data = og, boot = 100) -> s120_aib
# vroom::vroom_write(s40_ai, here::here("output", "s1200_aib.csv"), delim = ",")
sims_boot(lfreq, cpue, samples = 140, yrs = 2006, strata = NULL, og_data = og, boot = 100) -> s140_aib
# vroom::vroom_write(s40_ai, here::here("output", "s140_aib.csv"), delim = ",")



bind_rows(
  s40_aib %>%
    mutate(type = "s040",
           ss = 40)
) %>%
  bind_rows(
    s60_aib %>%
      mutate(type = "s060",
             ss = 60)
  ) %>%
  bind_rows(
    s80_aib %>%
      mutate(type = "s080",
             ss = 80)
  )%>%
  bind_rows(
    s100_aib %>%
      mutate(type = "s100",
             ss = 100)
  ) %>%
  bind_rows(
    s120_aib %>%
      mutate(type = "s120",
             ss = 120)
  ) %>%
  bind_rows(
    s140_aib %>%
      mutate(type = "s140",
             ss = 140)
  ) %>%
  bind_rows(
    base_ai %>%
      mutate(type = "base",
             sim = as.numeric(sim),
             ss = 0)
  ) -> ai

# boxplot
ai %>%
  filter(ess!= "ess_t", species_code %in% c(21921, 10110, 21720)) %>%
  ggplot(aes(type, 1/value)) +
  geom_boxplot(fill=4, alpha = 0.3) +
  facet_wrap(species_code~ess, scales = "free_y", ncol = 2)

# same fig but violin plot
ai %>%
  filter(ess!= "ess_t", species_code %in% c(21921, 10110, 21720)) %>%
  ggplot(aes(type, 1/value)) +
  geom_violin(draw_quantiles = 0.5, fill=4, alpha = 0.3) +
  facet_wrap(species_code~ess, scales = "free_y", ncol = 2)

# same fig but toggle by samplesize and year
ai %>%
  filter(ess!= "ess_t", species_code %in% c(21921, 10110, 21720)) %>%
  # group_by(year, species_code, ess, type, ss) %>%
  group_by(species_code, ess, type, ss) %>%
  summarise(ss = mean(ss),
            mean = 1/mean(value),
            se = 1/sd(value) / sqrt(n())) %>%
  mutate(lci = mean - se * 1.96,
         uci = mean + se * 1.96,
         lci = ifelse(lci <0, 0, lci)) %>%
  ggplot(aes(type, mean, color = type, group = ess)) +
  geom_pointrange(aes(ymin = lci, ymax = uci)) +
  facet_wrap(species_code~ess, scales = "free_y", ncol = 2)  +
  scico::scale_color_scico_d(palette = "roma")

# pointrange by year and sample size
ai %>%
  filter(ess!= "ess_t", species_code %in% c(21921, 10110, 21720)) %>%
  group_by(year, species_code, ess, type) %>%
  summarise(mean = 1/mean(value),
            se = 1/sd(value) / sqrt(n())) %>%
  mutate(lci = mean - se * 1.96,
         uci = mean + se * 1.96,
         lci = ifelse(lci <0, 0, lci)) %>%
  ggplot(aes(year, mean, color = type, group = ess)) +
  geom_pointrange(aes(ymin = lci, ymax = uci),
                  position=position_jitter(width=0.5)) +
  facet_wrap(species_code~ess, scales = "free_y", ncol = 2) +
  scale_x_continuous(breaks = seq(2006, 2018, 2)) +
  scico::scale_color_scico_d(palette = "roma")
