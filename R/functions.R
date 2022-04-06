library(tidyverse)
library(tidytable)
library(vroom)
library(here)
library(purrr)
library(rsample)
library(data.table)

############################################################################################################
# Function to estimate pop'n at size
size_pop_est <- function(lfreq, cpue, samples = 10000, yrs = NULL, strata = NULL){

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

############################################################################################################
# Function to estimate pop'n at age
age_pop_est <- function(specimen, sizepop, yrs = 2017, sim = NULL,  save = NULL){

  if(is.null(sim)){

    # filter specimen data to year
    setDT(specimen) %>%
      filter.(year >= yrs) -> .agdat

    # reformat size pop'n data
    setDT(sizepop) %>%
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
      select.(year, species_code, n) -> .sex_cnt

    if(length(.sex_cnt$n)>0){
      .sizepop_long %>%
        filter.(sex == 3) -> .sizepop_long_un

      .agdat %>%
        left_join.(.sex_cnt) %>%
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
        rename.(unsexed = '3', males = '1', females = '2') -> .agepop
    } else{
      .agepop_mf %>%
        pivot_wider(names_from = sex, values_from = agepop, values_fill = 0) %>%
        rename.(males = '1', females = '2') %>%
        mutate.(unsexed = 0) -> .agepop
    }

  } else{

    # filter specimen data to year
    specimen %>%
      group_by(species_code) %>%
      filter(year >= yrs) -> agdat

    # reformat size pop'n data
    sizepop %>%
      pivot_longer(cols = c(males, females, unsexed), names_to = "sex") %>%
      mutate(sex = replace(sex, sex == 'males', 1),
             sex = replace(sex, sex == 'females', 2),
             sex = replace(sex, sex == 'unsexed', 3)) %>%
      mutate(sex = as.numeric(sex)) %>%
      rename(sizepop = value) -> sizepop_long

    # compute resampled age pop'n for females & males
    agdat %>%
      group_by(year, species_code, sex, length, age) %>%
      summarise(age_num = length(age)) %>%
      group_by(year, species_code, sex, length) %>%
      mutate(age_frac = age_num/sum(age_num)) %>%
      left_join(sizepop_long) %>%
      mutate(agepop = age_frac * sizepop) %>%
      group_by(year, species_code, sex, sim, age) %>%
      mutate(agepop = sum(agepop)) %>%
      select(year, species_code, sex, age, sim, agepop) %>%
      group_by(year, species_code, sex, age, sim, agepop) %>%
      distinct(age) %>%
      filter(sex <= 2) -> agepop_mf

    # compute resampled age pop'n for unsexed (og rule is if you have a year with unsexed specimen data you use all the specimen data)
    agdat %>%
      group_by(year, species_code) %>%
      count(sex) %>%
      filter(sex == 3) %>%
      select(year, species_code, n) -> sex_cnt

    sizepop_long %>%
      filter(sex == 3) -> sizepop_long_un

    agdat %>%
      left_join(sex_cnt) %>%
      filter(n > 0) %>%
      group_by(year, species_code, length, age) %>%
      summarise(age_num = length(age)) %>%
      group_by(year, species_code, length) %>%
      mutate(age_frac = age_num/sum(age_num)) %>%
      left_join(sizepop_long_un) %>%
      mutate(agepop = age_frac * sizepop) %>%
      group_by(year, species_code, age, sim) %>%
      mutate(agepop = sum(agepop)) %>%
      select(year, species_code, sex, age, sim, agepop) %>%
      group_by(year, species_code, sex, age, sim, agepop) %>%
      distinct(age) %>%
      bind_rows(agepop_mf) %>%
      pivot_wider(names_from = sex, values_from = agepop, values_fill = 0) %>%
      rename(unsexed = '3', males = '1', females = '2') -> agepop

  }

  #vroom::vroom_write(agepop, here::here("output", region, paste0(save, ".csv")), delim = ",")
  if(is.null(sim)){
    .agepop} else{
      agepop
    }
}

############################################################################################################
# Function to run bootstrap simulations for size
size_sims_boot <- function(lfreq, cpue, samples = NULL, yrs = NULL, strata = NULL, og_data, boot = 500, write_comp = NULL, save = NULL, region = NULL){

  boot_res <- rerun(boot, size_pop_est_boot(lfreq, cpue, samples = samples, yrs = yrs))

  if(!is.null(write_comp) & boot>=100){
    boot_res %>%
      map_df.(., ~as.data.frame(.x), .id = "sim") %>%
      group_by(year, species_code, length) %>%
      sample_n(., 100) -> boot_comp
    vroom::vroom_write(boot_comp, here::here("output", region, paste0("boot_", save, "_comp_sz.csv")), delim = ",")
  }

  boot_res %>%
    map.(., ~ess_size(sim_data = .x, og_data = og_data)) %>%
    map_df.(., ~as.data.frame(.x), .id = "sim") -> boot_ess

  if(!is.null(save)){
    vroom::vroom_write(boot_ess, here::here("output", region, paste0("boot_", save, "_ess_sz.csv")), delim = ",")
  }
  boot_ess
}

############################################################################################################
# Function to run bootstrap simulations for age & size
age_sims_boot <- function(lfreq, specimen, cpue, samples = NULL, yrs = NULL, og_data_age, og_data_size, boot = 500, write_comp = NULL, save = NULL, region = NULL){

  boot_res <- rerun(boot, age_pop_est_boot(lfreq, specimen, cpue, samples = samples, yrs = yrs))

  age_res <- do.call(mapply, c(list, boot_res, SIMPLIFY = FALSE))$age
  length_res <- do.call(mapply, c(list, boot_res, SIMPLIFY = FALSE))$length

  # write out comp data
  age_res %>%
    map_df.(., ~as.data.frame(.x), .id = "sim") -> boot_comp_age
  vroom::vroom_write(boot_comp_age, here::here("output", region, paste0("boot_", save, "_comp_ag.csv")), delim = ",")
  length_res %>%
    map_df.(., ~as.data.frame(.x), .id = "sim") -> boot_comp_length
  vroom::vroom_write(boot_comp_length, here::here("output", region, paste0("boot_", save, "_comp_sz.csv")), delim = ",")

  # # sammple 100 simulations to write out comp data
  # if(!is.null(write_comp) & boot>=100){
  #   age_res %>%
  #     map_df.(., ~as.data.frame(.x), .id = "sim") %>%
  #     filter.(sim %in% sample.int(boot, 100)) -> boot_comp_age
  #   vroom::vroom_write(boot_comp_age, here::here("output", region, paste0("boot_", save, "_comp_ag.csv")), delim = ",")
  #   length_res %>%
  #     map_df.(., ~as.data.frame(.x), .id = "sim") %>%
  #     filter.(sim %in% sample.int(boot, 100)) -> boot_comp_length
  #   vroom::vroom_write(boot_comp_length, here::here("output", region, paste0("boot_", save, "_comp_sz.csv")), delim = ",")
  # }

  age_res %>%
    map.(., ~ess_age(sim_data = .x, og_data = og_data_age)) %>%
    map_df.(., ~as.data.frame(.x), .id = "sim") -> boot_ess_age

  length_res %>%
    map.(., ~ess_size(sim_data = .x, og_data = og_data_size)) %>%
    map_df.(., ~as.data.frame(.x), .id = "sim") -> boot_ess_size

  if(!is.null(save)){
    vroom::vroom_write(boot_ess_age, here::here("output", region, paste0("boot_", save, "_ess_ag.csv")), delim = ",")
    vroom::vroom_write(boot_ess_size, here::here("output", region, paste0("boot_", save, "_ess_sz.csv")), delim = ",")
  }
  list(age = boot_ess_age, size = boot_ess_size)
}

############################################################################################################
# Function to bootstrap hauls -> lengths for size pop'n estimates
size_pop_est_boot <- function(lfreq, cpue, samples = 10000, yrs = NULL, strata = NULL){

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

  # bootstrap the hauls and then lengths
  setDT(lfreq) %>%
    filter.(year >= yrs) %>%
    group_by(year, species_code) %>%
    distinct(hauljoin) %>%
    as_tidytable(.) %>%
    .[, hauljoin := hauljoin[sample.int(.N, .N, replace = TRUE)],
      by = c('year', 'species_code')] -> hls

  setDT(hls) %>%
    left_join.(lfreq) %>%
    uncount.(., frequency) %>%
    mutate.(sex_ln = paste0(sex, "-", length)) %>%
    .[, sex_ln := sex_ln[sample.int(.N, .N, replace = TRUE)],
      by = c('year', 'species_code', 'hauljoin') ] %>%
    separate.(., sex_ln, c('sex_res', 'len_res'), sep = '-') %>%
    mutate.(sex = as.numeric(sex_res), length = as.numeric(len_res)) %>%
    select.(-sex_res, -len_res) -> .lfreq

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

############################################################################################################
# Function to bootstrap hauls -> lengths & ages for size & age pop'n estimates
age_pop_est_boot <- function(lfreq, specimen, cpue, samples = 10000, yrs = NULL){

  # updated pop estimate with haul/length bootstraps (w/resampling)

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

  ## esimate size pop for males and females
  # pull out males and females, make single row for each, add an id
  .lfreq %>%
    filter.(sex != 3) %>%
    mutate.(id = .I) %>%
    mutate.(n = .N, .by = c(year, species_code, stratum, hauljoin)) -> .inter

  # sample by sample size
  .inter %>%
    group_by(year, species_code, stratum, hauljoin) %>%
    sample_n(if(n > samples) samples else n) -> .new_sexed

  # rejoin to original unsexed
  .lfreq %>%
    filter.(sex == 3) %>%
    mutate.(id = .I,
            n = .N, .by = c(year, species_code, stratum, hauljoin))  %>%
    bind_rows.(.new_sexed) %>%
    # bind_rows.(.new_sexed, .new_unsexed) %>%
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
      select.(year, species_code, length, males = `1`, females = `2`, unsexed = `3`) -> .lenout_mf
      #select.(-unsexed) -> .lenout_mf
  } else{
    .temp %>%
      select.(-stratum) %>%
      summarise.(abund = sum(sz_pop, na.rm = T), .by = c(year, species_code, length, sex)) %>%
      pivot_wider.(names_from = sex, values_from = abund) %>%
      left_join.(lngs, .) %>%
      mutate.(across.(.cols = c(`1`, `2`), ~replace_na.(.x, 0))) %>%
      select.(year, species_code, length, males = `1`, females = `2`) %>%
      mutate.(unsexed = 0) -> .lenout_mf
  }

  ## estimate size pop for unsexed
  # find new unsexed, that were previously sexed
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
      select.(year, species_code, length, males = `1`, females = `2`, unsexed = `3`) %>%
      select.(-males, -females) -> .lenout_u
  } else{
    .temp %>%
      select.(-stratum) %>%
      summarise.(abund = sum(sz_pop, na.rm = T), .by = c(year, species_code, length, sex)) %>%
      pivot_wider.(names_from = sex, values_from = abund) %>%
      left_join.(lngs, .) %>%
      mutate.(across.(.cols = c(`1`, `2`), ~replace_na.(.x, 0))) %>%
      select.(year, species_code, length, males = `1`, females = `2`) %>%
      mutate.(unsexed = 0) %>%
      select.(-males, -females) -> .lenout_u
  }

  .lenout_mf %>%
    select.(-unsexed) %>%
    left_join.(.lenout_u) -> .lenout

  # bootstrap ages within bootstrapped hauls
  setDT(hls) %>%
    left_join.(specimen) %>%
    mutate.(sex_ln_ag = paste0(sex, "-", length, "-", age)) %>%
    .[, sex_ln_ag := sex_ln_ag[sample.int(.N, .N, replace = TRUE)],
      by = c('year', 'species_code', 'hauljoin') ] %>%
    separate.(., sex_ln_ag, c('sex_res', 'len_res', "age_res"), sep = '-') %>%
    mutate.(sex = as.numeric(sex_res), length = as.numeric(len_res), age = as.numeric(age_res)) %>%
    select.(-sex_res, -len_res, -age_res) -> .agdat

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
# Function to simulate data at set sub-sampling level for lengths that are sexed
sims <- function(iters = 1, lfreq, cpue, strata = NULL, samples = NULL, yrs = 2017,  save = NULL, removed_summ=NULL){

  if(!is.null(samples) & is.null(save)){
    stop("you have to save the 'newly unsexed' samples - aka save = 's20'" )
  }

  .reps <- replicate(iters, size_pop_est(lfreq, cpue, samples, yrs, strata), simplify = FALSE)


  if(!is.null(samples) & !is.null(save) ){
    do.call(mapply, c(list, .reps, SIMPLIFY = FALSE))$new %>%
      purrr::map_df(., ~as.data.frame(.x), .id = "sim") -> .new
    do.call(mapply, c(list, .reps, SIMPLIFY = FALSE))$removed %>%
      purrr::map_df(., ~as.data.frame(.x), .id = "sim") -> .removed

    .removed %>%
      group_by(year, species_code, stratum, hauljoin) %>%
      summarise(ss_removed = length(id) / iters) -> .removed_summ

    vroom::vroom_write(.new, here::here("output", region, paste0(save, "_size.csv")), delim = ",")
    if(is.null(removed_summ)){
      vroom::vroom_write(.removed, here::here("output", region, paste0(save, "_removed.csv")), delim = ",")
    } else{
      vroom::vroom_write(.removed_summ, here::here("output", region, paste0(save, "_removed_summ.csv")), delim = ",")
    }
    .new
  } else if(is.null(samples) & is.null(save)){
    .reps  %>%
      purrr::map_df(., ~as.data.frame(.x), .id = "sim")

  } else {
    .reps  %>%
      purrr::map_df(., ~as.data.frame(.x), .id = "sim") -> .out
    vroom::vroom_write(.out, here::here("output", region, paste0(save, ".csv")), delim = ",")
    .out

  }


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

############################################################################################################
# Function to compute the proportion of females to males in sub-sampled data
prop_fm <- function(sim_data, og_data, strata = NULL, save){

  if("stratum" %in% names(og_data) & is.null(strata) | "stratum" %in% names(sim_data) & is.null(strata)){
    stop("check your strata")
  }

  if(!is.null(strata)){

    sim_data %>%
      group_by(sim, year, species_code, stratum) %>%
      mutate(prop_fm = sum(females)/sum(males)) %>%
      group_by(sim, year, species_code, stratum, prop_fm) %>%
      distinct(prop_fm) -> sim_st_prop

    og_data %>%
      group_by(year, species_code, stratum) %>%
      mutate(prop_fm_og = sum(females)/sum(males)) %>%
      group_by(year, species_code, stratum, prop_fm_og) %>%
      distinct(prop_fm_og) -> og_st_prop

    sim_st_prop %>%
      left_join(og_st_prop) %>%
      dplyr::filter(prop_fm >= 0 &  prop_fm_og >= 0) -> .out

  } else {

    sim_data %>%
      group_by(sim, year, species_code) %>%
      mutate(prop_fm = sum(females)/sum(males)) %>%
      group_by(sim, year, species_code, prop_fm) %>%
      distinct(prop_fm) -> sim_prop

    og_data %>%
      group_by(year, species_code) %>%
      mutate(prop_fm_og = sum(females)/sum(males)) %>%
      group_by(year, species_code, prop_fm_og) %>%
      distinct(prop_fm_og) -> og_prop

    sim_prop %>%
      left_join(og_prop) %>%
      dplyr::filter(prop_fm >= 0 &  prop_fm_og >= 0) -> .out

  }

  vroom::vroom_write(.out, here::here("output", region, paste0(save,".csv")), delim = ",")
  .out

}


