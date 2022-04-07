library(tidyverse)
library(tidytable)
library(vroom)
library(here)
library(purrr)
library(rsample)
library(data.table)

############################################################################################################
# Function to run bootstrap simulations for age & size
age_sims_boot <- function(lfreq, specimen, cpue, og_data_age, og_data_size, sx_subsample = NULL, yrs = NULL, sz_resample = NULL, age_resample = NULL, boot = 500, write_comp = NULL, save = NULL, region = NULL){

  boot_res <- rerun(boot, age_pop_est_boot(lfreq, specimen, cpue, sx_subsample = sx_subsample, yrs = yrs, sz_resample = sz_resample, age_resample = age_resample))

  age_res <- do.call(mapply, c(list, boot_res, SIMPLIFY = FALSE))$age
  length_res <- do.call(mapply, c(list, boot_res, SIMPLIFY = FALSE))$length

  # bootstrapped pop'n at age/length
  age_res %>%
    map_df.(., ~as.data.frame(.x), .id = "sim") -> boot_comp_age
  length_res %>%
    map_df.(., ~as.data.frame(.x), .id = "sim") -> boot_comp_length

  # ess of bootstrapped age/length
  age_res %>%
    map.(., ~ess_age(sim_data = .x, og_data = og_data_age)) %>%
    map_df.(., ~as.data.frame(.x), .id = "sim") -> boot_ess_age
  length_res %>%
    map.(., ~ess_size(sim_data = .x, og_data = og_data_size)) %>%
    map_df.(., ~as.data.frame(.x), .id = "sim") -> boot_ess_size

  # write out comp data
  if(!is.null(write_comp)){
    vroom::vroom_write(boot_comp_age, here::here("output", region, paste0("boot_", save, "_comp_ag.csv")), delim = ",")
    vroom::vroom_write(boot_comp_length, here::here("output", region, paste0("boot_", save, "_comp_sz.csv")), delim = ",")
  }

  # write out ess results
  if(!is.null(save)){
    vroom::vroom_write(boot_ess_age, here::here("output", region, paste0("boot_", save, "_ess_ag.csv")), delim = ",")
    vroom::vroom_write(boot_ess_size, here::here("output", region, paste0("boot_", save, "_ess_sz.csv")), delim = ",")
  }
  list(age = boot_ess_age, size = boot_ess_size)
}


############################################################################################################
# Function to bootstrap hauls -> lengths & ages for size & age pop'n estimates
age_pop_est_boot <- function(lfreq, specimen, cpue, sx_subsample = 10000, yrs = NULL, sz_resample = NULL, age_resample = NULL){

  # year switch
  if (is.null(yrs)) yrs <- 0

  # complete cases of unique lengths by species and sex (for all years)
  lfreq %>%
    group_by(species_code) %>%
    distinct(length, year, stratum) %>%
    expand(length, year, stratum) %>%
    filter(year >= yrs) %>%
    dplyr::select(-stratum) %>%
    group_by(species_code) %>%
    distinct(length, year) %>%
    expand(length, year) -> lngs

  # first pass of filtering
  setDT(cpue) %>%
    filter.(year >= yrs) -> .cpue


  if(!is.null(sz_resample)){
    # bootstrap the hauls and then lengths
    setDT(.cpue) %>%
      filter.(year >= yrs) %>%
      group_by(year, species_code) %>%
      distinct(hauljoin) %>%
      as_tidytable(.) %>%
      .[, hauljoin := hauljoin[sample.int(.N, .N, replace = TRUE)],
        by = c('year', 'species_code')] -> hls

    setDT(hls) %>%
      left_join.(.cpue) -> .cpue

    setDT(hls) %>%
      left_join.(lfreq) %>%
      filter.(!is.na(frequency)) %>%
      uncount.(., frequency) %>%
      mutate.(sex_ln = paste0(sex, "-", length)) %>%
      .[, sex_ln := sex_ln[sample.int(.N, .N, replace = TRUE)],
        by = c('year', 'species_code', 'hauljoin') ]  %>%
      separate.(., sex_ln, c('sex_res', 'len_res'), sep = '-') %>%
      mutate.(sex = as.numeric(sex_res), length = as.numeric(len_res))  %>%
      select.(-sex_res, -len_res) -> .lfreq
  } else{
    setDT(lfreq) %>%
      uncount.(., frequency) -> .lfreq
  }

  ## esimate size pop for males and females
  # pull out males and females, make single row for each, add an id
  .lfreq %>%
    filter.(sex != 3) %>%
    mutate.(id = .I) %>%
    mutate.(n = .N, .by = c(year, species_code, stratum, hauljoin)) -> .inter

  # sample by sample size
  .inter %>%
    group_by(year, species_code, stratum, hauljoin) %>%
    sample_n(if(n > sx_subsample) sx_subsample else n) -> .new_sexed

  .inter %>%
    anti_join.(.new_sexed, by = "id") %>%
    mutate(sex = 3) -> .new_unsexed

  # rejoin to original unsexed
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

  .temp %>%
    group_by(year, species_code) %>%
    count(sex) %>%
    filter.(sex == 3) %>%
    select.(year, species_code, n) -> .sex_cnt_sz

  if(length(.sex_cnt_sz$n)>0){
    .temp %>%
      select.(-stratum) %>%
      summarise.(abund = sum(sz_pop, na.rm = T), .by = c(year, species_code, length, sex)) %>%
      pivot_wider.(names_from = sex, values_from = abund) %>%
      left_join.(lngs, .) %>%
      mutate.(across.(.cols = c(`1`, `2`, `3`), ~replace_na.(.x, 0))) %>%
      select.(year, species_code, length, males = `1`, females = `2`, unsexed = `3`) -> .lenout
  } else{
    .temp %>%
      select.(-stratum) %>%
      summarise.(abund = sum(sz_pop, na.rm = T), .by = c(year, species_code, length, sex)) %>%
      pivot_wider.(names_from = sex, values_from = abund) %>%
      left_join.(lngs, .) %>%
      mutate.(across.(.cols = c(`1`, `2`), ~replace_na.(.x, 0))) %>%
      select.(year, species_code, length, males = `1`, females = `2`) %>%
      mutate.(unsexed = 0) -> .lenout
  }

  if(!is.null(age_resample)){
    # bootstrap ages within bootstrapped hauls
    setDT(hls) %>%
      left_join.(specimen) %>%
      mutate.(sex_ln_ag = paste0(sex, "-", length, "-", age)) %>%
      .[, sex_ln_ag := sex_ln_ag[sample.int(.N, .N, replace = TRUE)],
        by = c('year', 'species_code', 'hauljoin') ] %>%
      separate.(., sex_ln_ag, c('sex_res', 'len_res', "age_res"), sep = '-') %>%
      mutate.(sex = as.numeric(sex_res), length = as.numeric(len_res), age = as.numeric(age_res)) %>%
      select.(-sex_res, -len_res, -age_res) -> .agdat
  } else{
    # filter specimen data to year
    setDT(specimen) %>%
      filter.(year >= yrs) -> .agdat
  }

  # reformat size pop'n data
  .lenout %>%
    pivot_longer(cols = c(males, females, unsexed), names_to = "sex") %>%
    mutate.(sex = replace(sex, sex == 'males', 1),
           sex = replace(sex, sex == 'females', 2),
           sex = replace(sex, sex == 'unsexed', 3)) %>%
    mutate.(sex = as.numeric(sex)) %>%
    rename(sizepop = value) -> .sizepop_long

  # compute resampled age pop'n for females & males
  .agdat %>%
    group_by(year, species_code, sex, length, age) %>%
    summarise(age_num = length(age))  %>%
    mutate.(age_frac = age_num/sum(age_num), .by = c(year, species_code, sex, length)) %>%
    left_join.(.sizepop_long) %>%
    mutate.(agepop = age_frac * sizepop, .by = c(year, species_code, sex, length)) %>%
    mutate.(agepop = sum(agepop), .by = c(year, species_code, sex, age)) %>%
    select.(year, species_code, sex, age, agepop) %>%
    group_by(year, species_code, sex, age, agepop) %>%
    distinct(age) %>%
    filter.(sex <= 2) -> .agepop_mf

  # compute resampled age pop'n for unsexed (og rule is if you have a year with unsexed specimen data you use all the specimen data)
  .agdat %>%
    group_by(year, species_code) %>%
    count(sex) %>%
    filter.(sex == 3) %>%
    select.(year, species_code, n) -> .sex_cnt_ag

  if(length(.sex_cnt_ag$n)>0){
    .sizepop_long %>%
      filter.(sex == 3) -> .sizepop_long_un

    .agdat %>%
      left_join.(.sex_cnt_ag) %>%
      filter.(n > 0) %>%
      group_by(year, species_code, length, age) %>%
      summarise(age_num = length(age)) %>%
      mutate.(age_frac = age_num/sum(age_num), .by = c(year, species_code, length)) %>%
      left_join.(.sizepop_long_un) %>%
      mutate.(agepop = age_frac * sizepop, .by = c(year, species_code, length)) %>%
      mutate.(agepop = sum(agepop), .by = c(year, species_code, age)) %>%
      select.(year, species_code, sex, age, agepop) %>%
      group_by(year, species_code, sex, age, agepop) %>%
      distinct(age) %>%
      bind_rows.(.agepop_mf) %>%
      pivot_wider(names_from = sex, values_from = agepop, values_fill = 0) %>%
      rename.(unsexed = '3', males = '1', females = '2') -> .ageout
  } else{
    .agepop_mf %>%
      pivot_wider(names_from = sex, values_from = agepop, values_fill = 0) %>%
      rename.(males = '1', females = '2') %>%
      mutate.(unsexed = 0) -> .ageout
  }


  list(age = .ageout, length = .lenout)
}


############################################################################################################
# Function to estimate effective sample size for size comps
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
      select.(year, species_code, stratum, length, og_m, og_f, og_t) -> og_prop


    sim_data %>%
      mutate.(prop_m = males / sum(males),
              prop_f = females / sum(females),
              prop_t = (males + females + unsexed)/(sum(males) + sum(females) + sum(unsexed)),
              .by = c(year, species_code, stratum)) %>%
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
      select.(year, species_code, length, og_m, og_f, og_t) -> og_prop

    sim_data %>%
      mutate.(prop_m = males / sum(males),
              prop_f = females / sum(females),
              prop_t = (males + females + unsexed)/(sum(males) + sum(females) + sum(unsexed)),
              .by = c(year, species_code)) %>%
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

############################################################################################################
# Function to estimate effective sample size for age comps
ess_age <- function(sim_data, og_data){

  og_data %>%
    mutate.(og_m = males / sum(males),
           og_f = females / sum(females),
           og_t = (males + females + unsexed)/(sum(males) + sum(females) + sum(unsexed)),
           .by = c(year, species_code)) %>%
    select.(year, species_code, age, og_m, og_f, og_t) -> og_prop

  sim_data %>%
    mutate.(prop_m = males / sum(males),
           prop_f = females / sum(females),
           prop_t = (males + females + unsexed)/(sum(males) + sum(females) + sum(unsexed)),
           .by = c(year, species_code)) %>%
    left_join.(og_prop) %>%
    mutate.(ess_f = sum(prop_f * (1 - prop_f)) / sum((prop_f - og_f)^2),
           ess_m = sum(prop_m * (1 - prop_m)) / sum((prop_m - og_m)^2),
           ess_t = sum(prop_t * (1 - prop_t)) / sum((prop_t - og_t)^2),
           .by = c(year, species_code)) %>%
    pivot_longer.(cols = c(ess_f, ess_m, ess_t), names_to = "ess") %>%
    group_by(year, species_code, ess, value) %>%
    distinct(value)

}
