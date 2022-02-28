# compare decreased sample size length freq to original

# functions ----
source("R/functions.R")

# read in data ----
og <- vroom::vroom(here::here("output", "og.csv"))
og_st <- vroom::vroom(here::here("output", "og_st.csv"))

s20 <- vroom::vroom(here::here("output", "s20_comp.csv"))
s40 <- vroom::vroom(here::here("output", "s40_comp.csv"))
s60 <- vroom::vroom(here::here("output", "s60_comp.csv"))
s80 <- vroom::vroom(here::here("output", "s80_comp.csv"))
s100 <- vroom::vroom(here::here("output", "s100_comp.csv"))
s120 <- vroom::vroom(here::here("output", "s120_comp.csv"))
s140 <- vroom::vroom(here::here("output", "s140_comp.csv"))

# stratified data
s20_st <- vroom::vroom(here::here("output", "s20_comp.csv"))
s40_st <- vroom::vroom(here::here("output", "s40_comp.csv"))
s60_st <- vroom::vroom(here::here("output", "s60_comp.csv"))
s80_st <- vroom::vroom(here::here("output", "s80_comp.csv"))
s100_st <- vroom::vroom(here::here("output", "s100_comp.csv"))
s120_st <- vroom::vroom(here::here("output", "s120_comp.csv"))
s140_st <- vroom::vroom(here::here("output", "s140_comp.csv"))

# all years
s80yr <- vroom::vroom(here::here("output", "s80yr.csv"))


# figures - do not plot strata...
plot_comp(og, s20, species = "10110")
plot_comp(og, s40, species = "10110")
# plot_comp(og, s60, species = "10210")
plot_comp(og, s80, species = "10110")
# plot_comp(og, s100, species = "10210")
plot_comp(og, s120, species = "10110")
plot_comp(og, s140, species = "10110")

plot_comp(og, s80yr, species = "21740", yrs = 2015:2019)

plot_comp2(og, s20, species = c("10110", "10120"))
plot_comp2(og, s40, species = c("10110", "10120"))
plot_comp2(og, s60, species = c("10110", "10120"))
plot_comp2(og, s80, species = c("10110", "10120"))
plot_comp2(og, s100, species = c("10110", "10120"))
plot_comp2(og, s120, species = c("10110", "10120"))
plot_comp2(og, s140, species = c("10110", "10120"))

plot_comp2(og, s80yr, species = c("10110", "10120"))
