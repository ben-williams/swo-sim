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

# read in EBS data ----
region = "ebs"

og_size_ebs <- vroom::vroom(here::here("output", region, "og_size.csv"))

ess_s40_size_ebs <- vroom::vroom(here::here("output", region, "boot_s40_ess_sz.csv"))
ess_s60_size_ebs <- vroom::vroom(here::here("output", region, "boot_s60_ess_sz.csv"))
ess_s80_size_ebs <- vroom::vroom(here::here("output", region, "boot_s80_ess_sz.csv"))
ess_s100_size_ebs <- vroom::vroom(here::here("output", region, "boot_s100_ess_sz.csv"))
ess_s120_size_ebs <- vroom::vroom(here::here("output", region, "boot_s120_ess_sz.csv"))
ess_s140_size_ebs <- vroom::vroom(here::here("output", region, "boot_s140_ess_sz.csv"))
ess_base_size_ebs <- vroom::vroom(here::here("output", region, "boot_base_ess_sz.csv"))

# collate plot data
ebs <- collate4plot(ess_s40_size_ebs, ess_s60_size_ebs, ess_s80_size_ebs, ess_s100_size_ebs, ess_s120_size_ebs, ess_s140_size_ebs, ess_base_size_ebs)


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

# same fig but toggle by sample size and year
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

