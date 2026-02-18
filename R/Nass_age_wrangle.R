# script for wrangling raw age data from LGL 
library(tidyverse)
library(here)
library(ggsidekick)

# read in and wrangle data ----
nass_ages_raw <- read.csv(here("data/Upper lower coastal Nass Chinook ages 1989-2025.csv"))

nass_ages <- nass_ages_raw |>
  select(Year, CU, Location, NFL, MEF, GR_age, Sex) |>
  mutate(salt_age = case_when(GR_age == "21" ~ 1, 
                              GR_age == "31" ~ 2, 
                              GR_age == "41" ~ 3,
                              GR_age == "51" ~ 4, 
                              GR_age == "32" ~ 1, 
                              GR_age == "42" ~ 2,
                              GR_age == "52" ~ 3, 
                              GR_age == "62" ~ 4, 
                              GR_age == "72" ~ 5,
                              GR_age == "43" ~ 1, 
                              GR_age == "53" ~ 2, 
                              GR_age == "63" ~ 3,
                              GR_age == "73" ~ 4, 
                              GR_age == "83" ~ 5),
         fresh_age = case_when(GR_age == "21" ~ 1, 
                              GR_age == "31" ~ 1, 
                              GR_age == "41" ~ 1,
                              GR_age == "51" ~ 1, 
                              GR_age == "32" ~ 2, 
                              GR_age == "42" ~ 2,
                              GR_age == "52" ~ 2, 
                              GR_age == "62" ~ 2, 
                              GR_age == "72" ~ 2,
                              GR_age == "43" ~ 3, 
                              GR_age == "53" ~ 3, 
                              GR_age == "63" ~ 3,
                              GR_age == "73" ~ 3, 
                              GR_age == "83" ~ 4),
         Sex_full = case_when(Sex == "M" ~ "Male",
                              Sex == "F" ~ "Female"),
         age = fresh_age + salt_age)

# exploratory plots ----
# age comps by CU
ggplot(nass_ages, aes(x = Year, fill = as.factor(age))) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~CU, nrow  = 2) +
  scale_fill_brewer(palette = "Dark2") +
  theme_sleek() +
  ylab("") + 
  labs(fill = "Age")
  
# age comps by sex for upper CU
ggplot(nass_ages |> filter(CU == "Upper Nass", Sex %in% c("M","F")), aes(x = Year, fill = as.factor(age))) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~Sex_full, nrow  = 2) +
  scale_fill_brewer(palette = "Dark2") +
  theme_sleek() +
  ylab("") + 
  labs(fill = "Age")

# age comps by location
ggplot(nass_ages, aes(x = Year, fill = as.factor(age))) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~Location, ncol  = 3) +
  scale_fill_brewer(palette = "Dark2") +
  theme_sleek() +
  ylab("") + 
  labs(fill = "Age")

# female length at age
ggplot(nass_age_length_raw, aes(x = Year, y = MEF)) +
  stat_summary(fun = mean, geom = "point", size = 3,col="dark grey") + # Add points for the mean
  stat_summary(
    fun.data = mean_cl_boot, # Calculate mean and standard error
    geom = "errorbar",
    width = 0, # Control the width of the error bar whiskers
    col="dark grey"
  ) +
  labs(x = "Year", y = "Female length (MEF)")  +
  facet_wrap(~age) +
  geom_smooth(method = "lm", level = 0.95, color = "dark grey", fill = "lightblue") +
  theme_sleek()

# plot age comps excluding jacks and extremely rare older age classes
ggplot(nass_ages |> filter(Location == "GW fishwheels", between(age,4,6)), aes(x = Year, fill = as.factor(age))) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~CU, nrow  = 2) +
  scale_fill_brewer(palette = "Dark2") +
  theme_sleek() +
  ylab("") + 
  labs(fill = "Age")

# age comps for spawner recruitment modelling (use GW fishwheels as only source of data across years)
SR_age_comps <- nass_ages |>
  filter(Location == "GW fishwheels",
         between(age,4,6)) |>
  count(Year, age) |>
  group_by(Year) |>
  mutate(year_total = sum(n),
         prop_count = (n/year_total)*100) |>
  select(Year, age, prop_count) |>
  pivot_wider(names_from = age, values_from = prop_count) |>
  rename(a4 = '4', 
         a5 = '5', 
         a6 = '6')
write.csv(SR_age_comps, here("data/NassAobs.16Feb2026.csv"),row.names = FALSE)
