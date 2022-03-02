library(tidyverse)
library(vroom)
library(here)
library(purrr)
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

  } else {
    .temp %>%
      dplyr::select(-stratum) %>%
      group_by(year, species_code, length, sex) %>%
      summarise(abund = sum(sz_pop, na.rm = T)) %>%
      pivot_wider(names_from = sex, values_from = abund) %>%
      left_join(lngs, .) %>%
      mutate(across(tidyr::everything(), ~replace_na(.x, 0))) %>%
      dplyr::select(year, species_code, length, males = `1`, females = `2`, unsexed = `3`) -> .out
  }

  if(!is.null(samples)){
    list(new = .out, removed = .new_unsexed)
  } else {
    .out
  }

}


sims <- function(iters = 1, lfreq, cpue, strata = NULL, samples = NULL, yrs = 2017,  save = NULL){

  if(!is.null(samples) & is.null(save)){
    stop("you have to save the 'newly unsexed' samples - aka save = 's20'" )
  }

  .reps <- replicate(iters, pop_est(lfreq, cpue, samples, yrs, strata), simplify = FALSE)


  if(!is.null(samples) & !is.null(save) ){
    do.call(mapply, c(list, .reps, SIMPLIFY = FALSE))$new %>%
      purrr::map_df(., ~as.data.frame(.x), .id = "sim") -> .new
    do.call(mapply, c(list, .reps, SIMPLIFY = FALSE))$removed %>%
      purrr::map_df(., ~as.data.frame(.x), .id = "sim") -> .removed

    vroom::vroom_write(.new, here::here("output", paste0(save, "_comp.csv")), delim = ",")
    vroom::vroom_write(.removed, here::here("output", paste0(save, "_removed.csv")), delim = ",")
    .new
  } else if(is.null(samples) & is.null(save)){
    .reps  %>%
      purrr::map_df(., ~as.data.frame(.x), .id = "sim")

  } else {
    .reps  %>%
      purrr::map_df(., ~as.data.frame(.x), .id = "sim") -> .out
    vroom::vroom_write(.out, here::here("output", paste0(save, ".csv")), delim = ",")
    .out

  }


}


get_data <- function(data, id, strata = NULL, species = NULL, yrs = NULL){


  if(!is.null(strata)){
    data %>%
      pivot_longer(cols = c(males, females, unsexed)) %>%
      group_by(year, species_code, stratum, name, length) %>%
      summarise(mean = mean(value),
                se = sd(value) / sqrt(n()),
                .group = "drop") %>%
      mutate(lci = mean - se * 1.96,
             uci = mean + se * 1.96,
             id = id) -> .data
  } else {
  data %>%
    pivot_longer(cols = c(males, females, unsexed)) %>%
    group_by(year, species_code, name, length) %>%
    summarise(mean = mean(value),
              se = sd(value) / sqrt(n()),
              .group = "drop") %>%
    mutate(lci = mean - se * 1.96,
           uci = mean + se * 1.96,
           id = id) -> .data
  }

  if(!is.null(species) & !is.null(yrs)){
    .data %>%
      dplyr::filter(species_code %in% species, year %in% yrs)

  } else  if(!is.null(species) & is.null(yrs)){
    .data %>%
      dplyr::filter(species_code %in% species)

  } else  if(is.null(species) & !is.null(yrs)){

    .data %>%
      dplyr::filter(year %in% yrs)

  } else {
    .data
  }
}


plot_comp <- function(base_data, sim_data, species = NULL, yrs = NULL){

  if("stratum" %in% names(base_data) | "stratum" %in% names(sim_data)){
    stop("yeah, we aren't making all those plots \n use data that don't include strata...")
  }
  id1 = deparse(substitute(base_data))
  id2 = deparse(substitute(sim_data))

  get_data(base_data, id = id1, species, yrs) %>%
    bind_rows(get_data(sim_data, id = id2, species, yrs)) -> .data


  for (var in unique(.data$species_code)) {


   print(ggplot(.data[.data$species_code==var,], aes(length, mean, color = id)) +
             geom_line() +
             geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2, linetype = 0) +
             facet_wrap(name~year, scales = "free") +
             scale_color_manual(values = c(4, 1)) +
             scale_fill_manual(values = c(4, 1)) +
             guides(color = guide_legend(override.aes = list(fill=NA))) +
             ggtitle(paste0("species = ", var))
   )

   }
}


table_comp <- function(base_data, sim_data, strata = NULL, species = NULL, yrs = NULL){

  if("stratum" %in% names(base_data) & !("stratum" %in% names(sim_data)) |
     !("stratum" %in% names(base_data)) & "stratum" %in% names(sim_data)){
  stop("both files need to have stratum or not, can't be mixing them")
  }

  id1 = deparse(substitute(base_data))
  id2 = deparse(substitute(sim_data))

  get_data(base_data, id = id1, strata, species, yrs) %>%
    bind_rows(get_data(sim_data, id = id2, strata, species, yrs)) -> .data

  id1 = enquo(id1)
  id2 = enquo(id2)

  if("stratum" %in% names(base_data)){

    .data %>%
      group_by(year, species_code, stratum, name, id) %>%
      summarise(mean = mean(mean),
                .groups = "drop") %>%
      pivot_wider(names_from = id, values_from = mean) %>%
      mutate(diff = .data[[id1]] -.data[[id2]],
             pd = diff / ((.data[[id1]] + .data[[id2]]) / 2) * 100)

  } else {
    .data %>%
      group_by(year, species_code, name, id) %>%
      summarise(mean = mean(mean),
                .groups = "drop") %>%
      pivot_wider(names_from = id, values_from = mean) %>%
      mutate(diff = .data[[id1]] -.data[[id2]],
             pd = diff / ((.data[[id1]] + .data[[id2]]) / 2) * 100)
  }


}


plot_comp2 <- function(base_data, sim_data, species = NULL, yrs = NULL){


  table_comp(base_data, sim_data, species, yrs) %>%
    pivot_longer(cols = -c(year, species_code, name, diff, pd), names_to = "data") %>%
    ggplot(aes(value/1000, group = data, fill = data)) +
    geom_density(alpha = 0.3) +
    facet_wrap(species_code~name, scales = "free", dir = "v")

}


ess <- function(sim_data, og_data, strata = NULL, save){

  if(!is.null(strata)){
    og_data %>%
      group_by(year, species_code, stratum) %>%
      mutate(og_male = males / sum(males),
             og_female = females / sum(females)) -> og


    sim_data %>%
      group_by(sim, year, species_code, stratum) %>%
      mutate(prop_m = males / sum(males),
             prop_f = females / sum(females)) %>%
      distinct(prop_m, prop_f) %>%
      left_join(og) %>%
      summarise(ess_female = sum(prop_f * (1 - prop_f)) / sum((prop_f - og_female)^2),
                ess_male = sum(prop_m * (1 - prop_m)) / sum((prop_m - og_male)^2)) -> .out
  } else {

    og_data %>%
      group_by(year, species_code) %>%
      mutate(og_male = males / sum(males),
             og_female = females / sum(females)) %>%
      dplyr::select(year, species_code, length, og_male, og_female) -> og

    sim_data %>%
      group_by(sim, year, species_code) %>%
      mutate(prop_m = males / sum(males),
             prop_f = females / sum(females)) %>%
      left_join(og) %>%
      mutate(ess_female = sum(prop_f * (1 - prop_f)) / sum((prop_f - og_female)^2),
             ess_male = sum(prop_m * (1 - prop_m)) / sum((prop_m - og_male)^2)) -> .out

  }

  vroom::vroom_write(.out, here::here("output", paste0(save,".csv")), delim = ",")
  .out
}
