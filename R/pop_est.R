pop_est <- function(lfreq, cpue, lengths, samples = NULL, yrs = NULL){

  if(!is.null(samples)){
    lfreq %>%
      filter(sex != 3, year >= 2017) %>%
      uncount(frequency) %>%
      mutate(id = 1:n()) %>%
      group_by(year, species_code, stratum, hauljoin) %>%
      mutate(n = n()) -> inter

      inter %>%
        sample_n(if(n() > samples) samples else n()) %>%
        mutate(new_samp = n()) -> sexed

      inter %>%
        filter(!(id %in% sexed$id)) -> unsexed

    lfreq %>%
      filter(sex == 3, year >= 2017) %>%
      uncount(frequency) %>%
      group_by(year, species_code, stratum, hauljoin) %>%
      mutate(n = n())  %>%
      bind_rows(sexed) %>%
      group_by(year, species_code, stratum, hauljoin, sex, length) %>%
      summarise(frequency = n()) %>%
      group_by(year, species_code, stratum) %>%
      mutate(nhauls = length(unique(hauljoin))) %>%
      group_by(year, species_code, stratum, hauljoin) %>%
      mutate(tot = sum(frequency)) %>%
      group_by(year, species_code, stratum, hauljoin, sex, length) %>%
      summarise(comp = sum(frequency) / mean(tot),
                nhauls = mean(nhauls)) -> lcomp

    # estimate for hauls w/o length samples
    lcomp %>%
      group_by(year, species_code, stratum, sex, length) %>%
      summarise(comp = sum(comp) / mean(nhauls)) -> unk

    cpue %>%
      # filter(year == 2003, species_code==10110, stratum == 50,
      filter(year >= 2017, !is.na(catchjoin), !(hauljoin %in% lcomp$hauljoin)) -> no_length

    cpue %>%
      # filter(year == 2003, species_code==10110, stratum == 50) %>%
      group_by(year, species_code, stratum) %>%
      mutate(st_num = mean(numcpue) * area,
             tot = sum(numcpue)) %>%
      group_by(year, species_code, stratum, hauljoin) %>%
      summarise(abund = mean(numcpue) / tot * st_num) -> pop

    no_length %>%
      left_join(unk) %>%
      dplyr::select(year, species_code, stratum, hauljoin, sex, length, comp) %>%
      bind_rows(lcomp) %>%
      left_join(pop) %>%
      mutate(sz_pop = round(comp * abund, 0)) %>%
      group_by(year, species_code, stratum, length, sex) %>%
      summarise(abund = sum(sz_pop)) %>%
      pivot_wider(names_from = sex, values_from = abund) %>%
      left_join(lngs, .) %>%
      mutate(across(tidyr::everything(), ~replace_na(.x, 0))) %>%
      dplyr::select(year, species_code, stratum, length, males = `1`, females = `2`, unsexed = `3`)

  } else {

    lfreq %>%
      filter(year >= 2017) %>%
      group_by(year, species_code, stratum) %>%
      mutate(nhauls = length(unique(hauljoin))) %>%
      group_by(year, species_code, stratum, hauljoin) %>%
      mutate(tot = sum(frequency)) %>%
      group_by(year, species_code, stratum, hauljoin, sex, length) %>%
      summarise(comp = sum(frequency) / mean(tot),
                nhauls = mean(nhauls)) -> lcomp

    # estimate for hauls w/o length samples
    lcomp %>%
      group_by(year, species_code, stratum, sex, length) %>%
      summarise(comp = sum(comp) / mean(nhauls)) -> unk

    cpue %>%
      # filter(year == 2003, species_code==10110, stratum == 50,
      filter(year >= 2017, !is.na(catchjoin), !(hauljoin %in% lcomp$hauljoin)) -> no_length

    cpue %>%
      # filter(year == 2003, species_code==10110, stratum == 50) %>%
      group_by(year, species_code, stratum) %>%
      mutate(st_num = mean(numcpue) * area,
             tot = sum(numcpue)) %>%
      group_by(year, species_code, stratum, hauljoin) %>%
      summarise(abund = mean(numcpue) / tot * st_num) -> pop

    no_length %>%
      left_join(unk) %>%
      dplyr::select(year, species_code, stratum, hauljoin, sex, length, comp) %>%
      bind_rows(lcomp) %>%
      left_join(pop) %>%
      mutate(sz_pop = round(comp * abund, 0)) %>%
      group_by(year, species_code, stratum, length, sex) %>%
      summarise(abund = sum(sz_pop)) %>%
      pivot_wider(names_from = sex, values_from = abund) %>%
      left_join(lngs, .) %>%
      mutate(across(tidyr::everything(), ~replace_na(.x, 0))) %>%
      dplyr::select(year, species_code, stratum, length, males = `1`, females = `2`, unsexed = `3`)

  }
}
