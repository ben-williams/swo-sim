#' Analysis to examine survey history compared to Gerritsen and McGrath length
#' subsample suggested protocols of 10 samples per range of length bins

library(sumfish)
EBSdata <- getRacebase(2005:2021, 'EBS_SHELF')
speciesFilter <- 21720
juvFilter <- 21721
haulFilter <- select(EBSdata$haul, HAULJOIN) %>%
  unlist()

lengthTaken <- EBSdata$raw_length %>%
  filter(SPECIES_CODE %in% c(speciesFilter, juvFilter),
         HAULJOIN %in% haulFilter) %>%
  group_by(HAULJOIN,SPECIES_CODE) %>%
  summarize(num_lengths = sum(FREQUENCY)
  ) %>%
  mutate(juv_sampled = ifelse(SPECIES_CODE==juvFilter,1,0) ) %>%
  group_by(HAULJOIN) %>%
  summarize(sum_lengths = sum(num_lengths),
            juv_sampled = sum(juv_sampled)
  ) %>%
  inner_join(select(filter(EBSdata$catch, SPECIES_CODE==speciesFilter), HAULJOIN, NUMBER_FISH), by ='HAULJOIN')

lengthRange <- EBSdata$length %>%
  mutate(year = substr(CRUISE,1,4) ) %>%
  filter(SPECIES_CODE==speciesFilter) %>%
  group_by(HAULJOIN,year) %>%
  summarize(min_length = min(LENGTH),
            max_length = max(LENGTH)
  ) %>%
  mutate(ideal_subsample = (max_length - min_length)
  ) %>%
  inner_join(lengthTaken, by = 'HAULJOIN') %>%
  mutate(subsampled = ifelse(NUMBER_FISH-sum_lengths > 0, 1, 0) ) %>%
  filter(subsampled == 1)

ggplot() +
  geom_point(data = lengthRange,
             mapping= aes(x=ideal_subsample, y= sum_lengths,
                          color=year)) +
  geom_abline(color="purple") +
  geom_hline(yintercept = 200, color="red") +
  geom_hline(yintercept = 100, color="blue") +
  labs(title = tools::toTitleCase(EBSdata$species[EBSdata$species$SPECIES_CODE==speciesFilter,"COMMON_NAME"]),
       x = "Expected Subsample",
       y = "Actual Subsample",
       caption = "Actual EBS length subsample size plotted against expected subsample from Gerritsen & McGrath (2007). The purple
                  line is the one-to-one relationship, red is the historical minimum sample size, blue is the target sample
                  size recommended by SWO in 2022. Data are from valid EBS tows from 2005-2021 where a length subsample was taken
                  and extrapolated to estimated total numbers. This is a simplification, since some species were sampled separately
                  for juveniles and adults and this nuance in length range is not represented in the plot."
  ) +
  theme_bw() +
  theme(plot.caption = element_text(hjust=0)
  )
