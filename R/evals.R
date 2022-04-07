
# compare decreased sample size length freq to original

# functions ----
source("R/plot_functions.R")

# read in AI data ----
region = "ai"

og_size_ai <- vroom::vroom(here::here("output", region, "og_size.csv"))
ess_s40_size_ai <- vroom::vroom(here::here("output", region, "boot_s40_ess_sz.csv"))
ess_s60_size_ai <- vroom::vroom(here::here("output", region, "boot_s60_ess_sz.csv"))
ess_s80_size_ai <- vroom::vroom(here::here("output", region, "boot_s80_ess_sz.csv"))
ess_s100_size_ai <- vroom::vroom(here::here("output", region, "boot_s100_ess_sz.csv"))
ess_s120_size_ai <- vroom::vroom(here::here("output", region, "boot_s120_ess_sz.csv"))
ess_s140_size_ai <- vroom::vroom(here::here("output", region, "boot_s140_ess_sz.csv"))
ess_base_size_ai <- vroom::vroom(here::here("output", region, "boot_base_ess_sz.csv"))

ai_size <- collate4plot(ess_s40_size_ai, ess_s60_size_ai, ess_s80_size_ai, ess_s100_size_ai, ess_s120_size_ai, ess_s140_size_ai, ess_base_size_ai)

og_age_ai <- vroom::vroom(here::here("output", region, "og_age.csv"))
ess_s40_age_ai <- vroom::vroom(here::here("output", region, "boot_s40_ess_ag.csv"))
ess_s60_age_ai <- vroom::vroom(here::here("output", region, "boot_s60_ess_ag.csv"))
ess_s80_age_ai <- vroom::vroom(here::here("output", region, "boot_s80_ess_ag.csv"))
ess_s100_age_ai <- vroom::vroom(here::here("output", region, "boot_s100_ess_ag.csv"))
ess_s120_age_ai <- vroom::vroom(here::here("output", region, "boot_s120_ess_ag.csv"))
ess_s140_age_ai <- vroom::vroom(here::here("output", region, "boot_s140_ess_ag.csv"))
ess_base_age_ai <- vroom::vroom(here::here("output", region, "boot_base_ess_ag.csv"))

ai_age <- collate4plot(ess_s40_age_ai, ess_s60_age_ai, ess_s80_age_ai, ess_s100_age_ai, ess_s120_age_ai, ess_s140_age_ai, ess_base_age_ai)


comp_s40_age_ai <- vroom::vroom(here::here("output", region, "boot_s40_comp_ag.csv"))
comp_s60_age_ai <- vroom::vroom(here::here("output", region, "boot_s60_comp_ag.csv"))
comp_s80_age_ai <- vroom::vroom(here::here("output", region, "boot_s80_comp_ag.csv"))
comp_s100_age_ai <- vroom::vroom(here::here("output", region, "boot_s100_comp_ag.csv"))
comp_s120_age_ai <- vroom::vroom(here::here("output", region, "boot_s120_comp_ag.csv"))
comp_s140_age_ai <- vroom::vroom(here::here("output", region, "boot_s140_comp_ag.csv"))
comp_base_age_ai <- vroom::vroom(here::here("output", region, "boot_base_comp_ag.csv"))

ai_age_comp <- collate4plot(comp_s40_age_ai, comp_s60_age_ai, comp_s80_age_ai, comp_s100_age_ai, comp_s120_age_ai, comp_s140_age_ai, comp_base_age_ai)

comp_s40_size_ai <- vroom::vroom(here::here("output", region, "boot_s40_comp_sz.csv"))
comp_s60_size_ai <- vroom::vroom(here::here("output", region, "boot_s60_comp_sz.csv"))
comp_s80_size_ai <- vroom::vroom(here::here("output", region, "boot_s80_comp_sz.csv"))
comp_s100_size_ai <- vroom::vroom(here::here("output", region, "boot_s100_comp_sz.csv"))
comp_s120_size_ai <- vroom::vroom(here::here("output", region, "boot_s120_comp_sz.csv"))
comp_s140_size_ai <- vroom::vroom(here::here("output", region, "boot_s140_comp_sz.csv"))
comp_base_size_ai <- vroom::vroom(here::here("output", region, "boot_base_comp_sz.csv"))

ai_size_comp <- collate4plot(comp_s40_size_ai, comp_s60_size_ai, comp_s80_size_ai, comp_s100_size_ai, comp_s120_size_ai, comp_s140_size_ai, comp_base_size_ai)

# read in EBS data ----
# region = "ebs"
#
# og_size_ebs <- vroom::vroom(here::here("output", region, "og_size.csv"))
#
# ess_s40_size_ebs <- vroom::vroom(here::here("output", region, "boot_s40_ess_sz.csv"))
# ess_s60_size_ebs <- vroom::vroom(here::here("output", region, "boot_s60_ess_sz.csv"))
# ess_s80_size_ebs <- vroom::vroom(here::here("output", region, "boot_s80_ess_sz.csv"))
# ess_s100_size_ebs <- vroom::vroom(here::here("output", region, "boot_s100_ess_sz.csv"))
# ess_s120_size_ebs <- vroom::vroom(here::here("output", region, "boot_s120_ess_sz.csv"))
# ess_s140_size_ebs <- vroom::vroom(here::here("output", region, "boot_s140_ess_sz.csv"))
# ess_base_size_ebs <- vroom::vroom(here::here("output", region, "boot_base_ess_sz.csv"))
#
# # collate plot data
# ebs <- collate4plot(ess_s40_size_ebs, ess_s60_size_ebs, ess_s80_size_ebs, ess_s100_size_ebs, ess_s120_size_ebs, ess_s140_size_ebs, ess_base_size_ebs)
#

# Scatter plot of og proportions to mean simulated proportions

og_age_ai %>%
  mutate.(fem_comp_og = females/sum(females), .by = c(year, species_code)) %>%
  select.(-sim, -males, -females, -unsexed) -> fem_comp_og


ai_age_comp %>%
  mutate.(fem_comp_sim = females/sum(females), .by = c(year, species_code, type, sim)) %>%
  summarise.(fem_comp_sim = mean(fem_comp_sim), .by = c(year, species_code, age, type)) %>%
  left_join.(fem_comp_og) %>%
  filter.(species_code != 21740) %>%
  ggplot(aes(fem_comp_og, fem_comp_sim, color = species_code)) +
  geom_point() +
  ylim(0, 0.7) +
  xlim(0, 0.7) +
  geom_segment(aes(x=0, y=0, xend=0.7, yend=0.7))

# Scatter plot of og proportions to mean simulated proportions

og_size_ai %>%
  mutate.(fem_comp_og = females/sum(females), .by = c(year, species_code)) %>%
  select.(-sim, -males, -females, -unsexed) -> fem_comp_og

ai_size_comp %>%
  mutate.(fem_comp_sim = females/sum(females), .by = c(year, species_code, type, sim)) %>%
  summarise.(fem_comp_sim = mean(fem_comp_sim), .by = c(year, species_code, length, type)) %>%
  left_join.(fem_comp_og) %>%
  filter.(species_code != 21740) %>%
  ggplot(aes(fem_comp_og, fem_comp_sim, color = species_code)) +
  geom_point() +
  ylim(0, 0.25) +
  xlim(0, 0.25) +
  geom_segment(aes(x=0, y=0, xend=0.25, yend=0.25))


# line plot

ai_size %>%
  #filter(species_code %in% c(30060, 30420)) %>% # rockfish
  #filter(species_code %in% c(10110, 10130)) %>% # flatfish
  #filter(species_code %in% c(21740, 21720, 21921)) %>% # roundfish
  filter(ess != "ess_t", species_code %in% c(21720, 10110, 30060)) %>% # example species
  mutate(comm_name = as.character(species_code)) %>%
  mutate(comm_name = replace(comm_name, comm_name == "10110", "arrowtooth"),
         comm_name = replace(comm_name, comm_name == "10130", "flathead"),
         comm_name = replace(comm_name, comm_name == "21720", "Pcod"),
         comm_name = replace(comm_name, comm_name == "21740", "pollock"),
         comm_name = replace(comm_name, comm_name == "21921", "atka"),
         comm_name = replace(comm_name, comm_name == "30020", "shortspine"),
         comm_name = replace(comm_name, comm_name == "30060", "POP"),
         comm_name = replace(comm_name, comm_name == "30420", "n rockfish")) %>%
  mutate.(ess = replace(ess, ess== "ess_f", "females"), ess = replace(ess, ess == "ess_m", "males"), ess = replace(ess, ess == "ess_t", "total")) %>%
  group_by(species_code, comm_name, ess, type, ss) %>%
  mutate(hmean = length(unique(sim))*length(unique(year))/sum(value^(-1)),
         lci = quantile(value, 0.025),
         uci = quantile(value, 0.975)) %>%
  summarise(ss = mean(ss),
            mean = mean(value),
            hmean = mean(hmean),
            lci = mean(lci),
            uci = mean(uci)) %>%
  ggplot(aes(type, mean, color = type, group = ess)) +
  geom_pointrange(aes(ymin = lci, ymax = uci))  +
  facet_grid(comm_name~ess, scales = "free_y")  +
  labs(title = "AI example",y = "Effective sample size of resampled length comps", x = "Subsampling case", color = "SS case") +
  scico::scale_color_scico_d(palette = "roma")

ai_age %>%
  #filter(species_code %in% c(30060, 30420)) %>% # rockfish
  #filter(species_code %in% c(10110, 10130)) %>% # flatfish
  #filter(species_code %in% c(21740, 21720, 21921)) %>% # roundfish
  filter(species_code %in% c(21720, 10110, 30060)) %>% # roundfish
  mutate(comm_name = as.character(species_code)) %>%
  mutate(comm_name = replace(comm_name, comm_name == "10110", "arrowtooth"),
         comm_name = replace(comm_name, comm_name == "10130", "flathead"),
         comm_name = replace(comm_name, comm_name == "21720", "Pcod"),
         comm_name = replace(comm_name, comm_name == "21740", "pollock"),
         comm_name = replace(comm_name, comm_name == "21921", "atka"),
         comm_name = replace(comm_name, comm_name == "30020", "shortspine"),
         comm_name = replace(comm_name, comm_name == "30060", "POP"),
         comm_name = replace(comm_name, comm_name == "30420", "n rockfish")) %>%
  group_by(species_code, comm_name, ess, type, ss) %>%
  mutate(hmean = length(unique(sim))*length(unique(year))/sum(value^(-1)),
         lci = quantile(value, 0.025),
         uci = quantile(value, 0.975)) %>%
  summarise(ss = mean(ss),
            mean = mean(value),
            hmean = mean(hmean),
            lci = mean(lci),
            uci = mean(uci)) %>%
  mutate.(ess = replace(ess, ess== "ess_f", "females"), ess = replace(ess, ess == "ess_m", "males"), ess = replace(ess, ess == "ess_t", "total")) %>%
  ggplot(aes(type, mean, color = type, group = ess)) +
  geom_pointrange(aes(ymin = lci, ymax = uci))  +
  facet_grid(comm_name~ess, scales = "free_y")  +
  labs(title = "AI example", y = "Effective sample size of resampled age comps", x = "Subsampling case", color = "SS case") +
  scico::scale_color_scico_d(palette = "roma")


# bar plot of harmonic mean
ai_size %>%
  #filter(ess!= "ess_t", species_code %in% c(30060)) %>% # rockfish example
  #filter(ess!= "ess_t",species_code %in% c(10110)) %>% # flatfish example
  #filter(ess!= "ess_t",species_code %in% c(21720)) %>% # roundfish example
  filter(ess!= "ess_t",year == 2018, species_code %in% c(10110, 21720, 30060)) %>% # species example
  mutate(comm_name = as.character(species_code)) %>%
  mutate(comm_name = replace(comm_name, comm_name == "10110", "arrowtooth"),
         comm_name = replace(comm_name, comm_name == "10130", "flathead"),
         comm_name = replace(comm_name, comm_name == "21720", "Pcod"),
         comm_name = replace(comm_name, comm_name == "21740", "pollock"),
         comm_name = replace(comm_name, comm_name == "21921", "atka"),
         comm_name = replace(comm_name, comm_name == "30020", "shortspine"),
         comm_name = replace(comm_name, comm_name == "30060", "POP"),
         comm_name = replace(comm_name, comm_name == "30420", "n rockfish")) %>%
  mutate.(ess = replace(ess, ess== "ess_f", "females"), ess = replace(ess, ess == "ess_m", "males")) %>%
  group_by(year, species_code, comm_name, ess, type, ss) %>%
  mutate(hmean = length(unique(sim))/sum(value^(-1)),
         lci = quantile(value, 0.025),
         uci = quantile(value, 0.975)) %>%
  summarise(ss = mean(ss),
            mean = mean(value),
            hmean = mean(hmean),
            lci = mean(lci),
            uci = mean(uci)) %>%
  ggplot(aes(type, hmean, fill = type, group = ess)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_grid(comm_name~ess, scales = "free_y")  +
  #labs(title = "Pacific ocean perch", y = "Survey length comp input sample size", x = "Subsampling case", fill = "SS case") + # rockfish example
  #labs(title = "Arrowtooth flounder", y = "Survey length comp input sample size", x = "Subsampling case", fill = "SS case") + # flatfish example
  #labs(title = "Pacific cod", y = "Survey length comp input sample size", x = "Subsampling case", fill = "SS case") + # roundfish example
  labs(title = "AI example", y = "Survey length comp input sample size (2018)", x = "Subsampling case", fill = "SS case") + # flatfish example
  scico::scale_fill_scico_d(palette = "roma")

ai_age %>%
  #filter(ess!= "ess_t", species_code %in% c(30060)) %>% # rockfish example
  #filter(ess!= "ess_t",species_code %in% c(10110)) %>% # flatfish example
  #filter(ess!= "ess_t",species_code %in% c(21720)) %>% # roundfish example
  filter(species_code %in% c(21720, 10110, 30060), year == 2018) %>% # example speies
  mutate(comm_name = as.character(species_code)) %>%
  mutate(comm_name = replace(comm_name, comm_name == "10110", "arrowtooth"),
         comm_name = replace(comm_name, comm_name == "10130", "flathead"),
         comm_name = replace(comm_name, comm_name == "21720", "Pcod"),
         comm_name = replace(comm_name, comm_name == "21740", "pollock"),
         comm_name = replace(comm_name, comm_name == "21921", "atka"),
         comm_name = replace(comm_name, comm_name == "30020", "shortspine"),
         comm_name = replace(comm_name, comm_name == "30060", "POP"),
         comm_name = replace(comm_name, comm_name == "30420", "n rockfish")) %>%
  group_by(year, species_code, comm_name, ess, type, ss) %>%
  mutate(hmean = length(unique(sim))/sum(value^(-1)),
         lci = quantile(value, 0.025),
         uci = quantile(value, 0.975)) %>%
  summarise(ss = mean(ss),
            mean = mean(value),
            hmean = mean(hmean),
            lci = mean(lci),
            uci = mean(uci)) %>%
  mutate.(ess = replace(ess, ess== "ess_f", "females"), ess = replace(ess, ess == "ess_m", "males"), ess = replace(ess, ess == "ess_t", "total")) %>%
  ggplot(aes(type, hmean, fill = type, group = ess)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_grid(comm_name~ess, scales = "free_y")  +
  #labs(title = "Pacific ocean perch", y = "Survey length comp input sample size", x = "Subsampling case", fill = "SS case") + # rockfish example
  #labs(title = "Arrowtooth flounder", y = "Survey age comp input sample size", x = "Subsampling case", fill = "SS case") + # flatfish example
  #labs(title = "Pacific cod", y = "Survey age comp input sample size", x = "Subsampling case", fill = "SS case") + # roundfish example
  labs(title = "AI example", y = "Survey age comp input sample size (2018)", x = "Subsampling case", fill = "SS case") + # flatfish example
  scico::scale_fill_scico_d(palette = "roma")



# Percent reduction plots

ai_size %>%
  filter(ess!= "ess_t") %>%
  mutate.(comm_name = as.character(species_code)) %>%
  mutate.(comm_name = replace(comm_name, comm_name == "10110", "arrowtooth"),
         comm_name = replace(comm_name, comm_name == "10130", "flathead"),
         comm_name = replace(comm_name, comm_name == "21720", "Pcod"),
         comm_name = replace(comm_name, comm_name == "21740", "pollock"),
         comm_name = replace(comm_name, comm_name == "21921", "atka"),
         comm_name = replace(comm_name, comm_name == "30020", "shortspine"),
         comm_name = replace(comm_name, comm_name == "30060", "POP"),
         comm_name = replace(comm_name, comm_name == "30420", "n rockfish")) %>%
  mutate.(hmean = length(unique(sim))/sum(value^(-1)), .by = c(year, species_code, comm_name, ess, type, ss)) %>%
  group_by(year, species_code, comm_name, ess, type, ss) %>%
  summarise(ss = mean(ss), hmean = mean(hmean)) -> .ai_size

.ai_size %>%
  filter.(type == "base") %>%
  mutate.(hmean_base = hmean, .by = c(year, species_code, comm_name)) %>%
  select.(-type, -ss, -hmean) -> .base_size

.ai_size %>%
  filter.(type != "base") %>%
  left_join.(.base_size) %>%
  mutate.(p_redux = 1-(hmean_base - hmean)/hmean_base, .by = c(year, species_code, comm_name, ess)) %>%
  #mutate.(p_redux = ifelse(p_redux > 1, 1, p_redux), .by = c(year, species_code, comm_name, ess)) %>%
  #filter(species_code %in% c(30060, 30420)) %>% # rockfish
  #filter(species_code %in% c(10110, 10130)) %>% # flatfish
  #filter(species_code %in% c(21740, 21720, 21921)) %>% # roundfish
  filter(species_code %in% c(21720, 10110, 30060)) %>% # example species
  mutate.(ess = replace(ess, ess== "ess_f", "females"), ess = replace(ess, ess == "ess_m", "males"), ess = replace(ess, ess == "ess_t", "total")) %>%
  ggplot(aes(type, p_redux, group = ess)) +
  geom_point(aes(shape = as.factor(year), color = as.factor(year))) +
  ylim(0, 1) +
  facet_grid(comm_name~ess, scales = "free_y") +
  labs(y = "Percent of base length comp input sample size", x = "Subsampling case", shape = "Year", color = "Year") +
  scico::scale_color_scico_d(palette = "roma")



ai_age %>%
  #filter(ess!= "ess_t") %>%
  mutate.(comm_name = as.character(species_code)) %>%
  mutate.(comm_name = replace(comm_name, comm_name == "10110", "arrowtooth"),
          comm_name = replace(comm_name, comm_name == "10130", "flathead"),
          comm_name = replace(comm_name, comm_name == "21720", "Pcod"),
          comm_name = replace(comm_name, comm_name == "21740", "pollock"),
          comm_name = replace(comm_name, comm_name == "21921", "atka"),
          comm_name = replace(comm_name, comm_name == "30020", "shortspine"),
          comm_name = replace(comm_name, comm_name == "30060", "POP"),
          comm_name = replace(comm_name, comm_name == "30420", "n rockfish")) %>%
  mutate.(hmean = length(unique(sim))/sum(value^(-1)), .by = c(year, species_code, comm_name, ess, type, ss)) %>%
  group_by(year, species_code, comm_name, ess, type, ss) %>%
  summarise(ss = mean(ss), hmean = mean(hmean)) -> .ai_age

.ai_age %>%
  filter.(type == "base") %>%
  mutate.(hmean_base = hmean, .by = c(year, species_code, comm_name)) %>%
  select.(-type, -ss, -hmean) -> .base_age

.ai_age %>%
  filter.(type != "base") %>%
  left_join.(.base_age) %>%
  mutate.(p_redux = 1-(hmean_base - hmean)/hmean_base, .by = c(year, species_code, comm_name, ess)) %>%
  mutate.(p_redux = ifelse(p_redux > 1, 1, p_redux), .by = c(year, species_code, comm_name, ess)) %>%
  mutate.(ess = replace(ess, ess== "ess_f", "females"), ess = replace(ess, ess == "ess_m", "males"), ess = replace(ess, ess == "ess_t", "total")) %>%
  #filter(species_code %in% c(30060, 30420)) %>% # rockfish
  #filter(species_code %in% c(10110, 10130)) %>% # flatfish
  #filter(species_code %in% c(21740, 21720, 21921)) %>% # roundfish
  filter(species_code %in% c(10110, 21720, 30060)) %>% # example species
  ggplot(aes(type, p_redux, group = ess)) +
  geom_point(aes(shape = as.factor(year), color = as.factor(year))) +
  ylim(0, 1) +
  facet_grid(comm_name~ess, scales = "free_y") +
  labs(y = "Percent of base age comp input sample size", x = "Subsampling case", shape = "Year", color = "Year") +
  scico::scale_color_scico_d(palette = "roma")







# boxplot
# ai %>%
#   filter(ess!= "ess_t", species_code %in% c(21921, 10110, 21720)) %>%
#   ggplot(aes(type, 1/value)) +
#   geom_boxplot(fill=4, alpha = 0.3) +
#   facet_wrap(species_code~ess, scales = "free_y", ncol = 2)

# same fig but violin plot
# ai %>%
#   filter(ess!= "ess_t", species_code %in% c(21921, 10110, 21720)) %>%
#   ggplot(aes(type, 1/value)) +
#   geom_violin(draw_quantiles = 0.5, fill=4, alpha = 0.3) +
#   facet_wrap(species_code~ess, scales = "free_y", ncol = 2)


# # same fig but toggle by sample size and year
# ai %>%
#   filter(ess!= "ess_t", species_code %in% c(21921, 10110, 21740)) %>%
#   # group_by(year, species_code, ess, type, ss) %>%
#   group_by(species_code, ess, type, ss) %>%
#   summarise(ss = mean(ss),
#             mean = 1/mean(value),
#             se = 1/sd(value) / sqrt(n())) %>%
#   mutate(lci = mean - se * 1.96,
#          uci = mean + se * 1.96,
#          lci = ifelse(lci <0, 0, lci)) %>%
#   ggplot(aes(type, mean, color = type, group = ess)) +
#   geom_pointrange(aes(ymin = lci, ymax = uci)) +
#   facet_wrap(species_code~ess, scales = "free_y", ncol = 2)  +
#   scico::scale_color_scico_d(palette = "roma")
#
#
# # ess toggled by sample size and year, uci/lci based on percentiles
# ai %>%
#   filter(ess!= "ess_t", species_code %in% c(21921, 10110, 21740)) %>%
#   group_by(year, species_code, ess, type, ss) %>%
#   #group_by(species_code, ess, type, ss) %>%
#   mutate(hmean = length(unique(sim))/sum(value^(-1)),
#          lci = quantile(value, 0.125),
#          uci = quantile(value, 0.875)) %>%
#   summarise(ss = mean(ss),
#             mean = mean(value),
#             hmean = mean(hmean),
#             lci = mean(lci),
#             uci = mean(uci)) %>%
#   ggplot(aes(type, hmean, color = type, group = ess)) +
#   geom_pointrange(aes(ymin = lci, ymax = uci),
#                   position=position_jitter(width=0.15))  +
#   facet_wrap(species_code~ess, scales = "free_y", ncol = 2)  +
#   scico::scale_color_scico_d(palette = "roma")
#









