library(tidyverse)
library(here)
library(ggsidekick)

nass_ages_raw <- read.csv(here("data/Upper lower coastal Nass Chinook ages 1989-2025.csv"))

nass_ages <- nass_ages_raw |>
  select(Year, CU, Location, MEF, GR_age, Sex) |>
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


ggplot(nass_ages, aes(x = Year, fill = as.factor(age))) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~CU, nrow  = 2) +
  scale_fill_brewer(palette = "Dark2") +
  theme_sleek() +
  ylab("") + 
  labs(fill = "Age")
  
ggplot(nass_ages |> filter(CU == "Upper Nass", Sex %in% c("M","F")), aes(x = Year, fill = as.factor(age))) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~Sex_full, nrow  = 2) +
  scale_fill_brewer(palette = "Dark2") +
  theme_sleek() +
  ylab("") + 
  labs(fill = "Age")
