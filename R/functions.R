library(tidyverse)
library(vroom)
library(here)
library(purrr)

pop_est <- function(lfreq, cpue, samples = NULL, yrs = 2017, strata = NULL){

   # year switch
  if(is.null(yrs)) yrs = 0

  # complete cases of unique lengths by species and sex (for all years)

  lfreq %>%
    group_by(species_code) %>%
    distinct(length, year, stratum) %>%
    expand(length, year, stratum) %>%
    filter(year >= yrs) -> lngs

  # if no strata remove from
  if(is.null(strata)){
    lngs %>%
      dplyr::select(-stratum) %>%
      group_by(species_code) %>%
      distinct(length, year) %>%
      expand(length, year) -> lngs
  }

  # first pass of filtering
  cpue %>%
    filter(year >= yrs) -> .cpue

  lfreq %>%
    filter(year >= yrs) -> .lfreq

  if(is.null(samples)){

    .lfreq %>%
      group_by(year, species_code, stratum) %>%
      mutate(nhauls = length(unique(hauljoin))) %>%
      group_by(year, species_code, stratum, hauljoin) %>%
      mutate(tot = sum(frequency)) %>%
      group_by(year, species_code, stratum, hauljoin, sex, length) %>%
      summarise(comp = sum(frequency) / mean(tot),
                nhauls = mean(nhauls)) -> .lcomp
  } else {

    # pull out males and females, make single row for each, add an id
    .lfreq %>%
      filter(sex != 3) %>%
      uncount(frequency) %>%
      mutate(id = 1:n()) %>%
      group_by(year, species_code, stratum, hauljoin) %>%
      mutate(n = n()) -> .inter

    # sample by sample size or # of rows (it is still a grouped dataframe)
    .inter %>%
      sample_n(if(n() > samples) samples else n()) %>%
      mutate(new_samp = n()) -> .new_sexed

    # find new unsexed, that were previously sexed
    .inter %>%
      anti_join(.new_sexed) -> .new_unsexed

    # rejoin original unsexed to the new_sexed samples
    .lfreq %>%
      filter(sex == 3) %>%
      uncount(frequency) %>%
      group_by(year, species_code, stratum, hauljoin) %>%
      mutate(n = n())  %>%
      bind_rows(.new_sexed) %>%
      group_by(year, species_code, stratum, hauljoin, sex, length) %>%
      summarise(frequency = n()) %>%
      group_by(year, species_code, stratum) %>%
      mutate(nhauls = length(unique(hauljoin))) %>%
      group_by(year, species_code, stratum, hauljoin) %>%
      mutate(tot = sum(frequency)) %>%
      group_by(year, species_code, stratum, hauljoin, sex, length) %>%
      summarise(comp = sum(frequency) / mean(tot),
                nhauls = mean(nhauls)) -> .lcomp
  }

  # estimate for hauls w/o length samples
  .lcomp %>%
    group_by(year, species_code, stratum, sex, length) %>%
    summarise(comp = sum(comp) / mean(nhauls)) -> .unk

  # id hauls without lengths
  .cpue %>%
    filter(!is.na(catchjoin), !(hauljoin %in% .lcomp$hauljoin)) -> .no_length

  .cpue %>%
    group_by(year, species_code, stratum) %>%
    mutate(st_num = mean(numcpue) * area,
           tot = sum(numcpue)) %>%
    group_by(year, species_code, stratum, hauljoin) %>%
    summarise(abund = mean(numcpue) / tot * st_num) -> .pop

  # if there are any samples w/o lengths rejoin them
  if(nrow(.no_length == 0)){
    .lcomp %>%
      left_join(.pop) %>%
      mutate(sz_pop = round(comp * abund, 0)) -> .temp
  } else {
    .no_length %>%
      left_join(.unk) %>%
      dplyr::select(year, species_code, stratum, hauljoin, sex, length, comp) %>%
      bind_rows(.lcomp) %>%
      left_join(.pop) %>%
      mutate(sz_pop = round(comp * abund, 0)) -> .temp
  }


  if(!is.null(strata)){
    .temp %>%
      group_by(year, species_code, stratum, length, sex) %>%
      summarise(abund = sum(sz_pop), na.rm = T) %>%
      pivot_wider(names_from = sex, values_from = abund) %>%
      left_join(lngs, .) %>%
      mutate(across(tidyr::everything(), ~replace_na(.x, 0))) %>%
      dplyr::select(year, species_code, stratum, length, males = `1`, females = `2`, unsexed = `3`) -> .out

    if(!is.null(samples)){
      list(new = .out, removed = .new_unsexed)
      } else {
        .out
      }
  } else {
    .temp %>%
      dplyr::select(-stratum) %>%
      group_by(year, species_code, length, sex) %>%
      summarise(abund = sum(sz_pop, na.rm = T)) %>%
      pivot_wider(names_from = sex, values_from = abund) %>%
      left_join(lngs, .) %>%
      mutate(across(tidyr::everything(), ~replace_na(.x, 0))) %>%
      dplyr::select(year, species_code, length, males = `1`, females = `2`, unsexed = `3`)
  }

}


sims <- function(iters = 1, lfreq, cpue, strata = NULL, samples = NULL, yrs = 2017){
  replicate(iters, pop_est(lfreq, cpue, samples, yrs, strata), simplify = FALSE)
}

getouts <- function(data, type = "comp", samples = NULL, save = NULL){

  if(!is.null(samples)){
  if(type == "comp" ){
  do.call(mapply, c(list, data, SIMPLIFY = FALSE))$new %>%
      purrr::map_df(., ~as.data.frame(.x), .id = "sim") -> .out
  } else {
    do.call(mapply, c(list, data, SIMPLIFY = FALSE))$removed %>%
      purrr::map_df(., ~as.data.frame(.x), .id = "sim") -> .out
  }
  } else {
    data  %>%
      purrr::map_df(., ~as.data.frame(.x), .id = "sim") -> .out
  }

  if(!is.null(save)){
    vroom::vroom_write(.out, here::here("output", save), delim = ",")
  } else {
    .out
  }


}

