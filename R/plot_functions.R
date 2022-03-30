# put your plotting code in here

# Get data together for plotting ess
collate4plot <- function (s40, s60, s80, s100, s120, s140, base){

  bind_rows(
    s40 %>%
      mutate(type = "s040",
             ss = 40)
  ) %>%
    bind_rows(
      s60 %>%
        mutate(type = "s060",
               ss = 60)
    ) %>%
    bind_rows(
      s80 %>%
        mutate(type = "s080",
               ss = 80)
    )%>%
    bind_rows(
      s100 %>%
        mutate(type = "s100",
               ss = 100)
    ) %>%
    bind_rows(
      s120 %>%
        mutate(type = "s120",
               ss = 120)
    ) %>%
    bind_rows(
      s140 %>%
        mutate(type = "s140",
               ss = 140)
    ) %>%
    bind_rows(
      base %>%
        mutate(type = "base",
               sim = as.numeric(sim),
               ss = 0)
    ) -> dat
  dat
}

# plots here
