# script for wrangling raw age data from LGL 
library(tidyverse)
library(here)
library(ggsidekick)

# read in and wrangle data ----
kit_mort_dist <- read.csv(here("data/Kitsumkalum_Mortality_distribution_data_2026-10-27.csv")) |>
  filter(FisheryTypeGrouped != "ESCAP") |>
  group_by(Year)


ggplot(kit_mort_dist, aes(x = Year, y = PercentMortality, fill = FisheryTypeGrouped)) +
  geom_col(position = "fill") + 
  scale_y_continuous(labels = scales::percent) + # Formats axis from 0.0-1.0 to 0%-100%
  labs(y = "Contribution to total fishing mortality") +
  scale_fill_brewer(palette = "Dark2") +
  theme_sleek()  
ggsave(here("plots/Skeena/kitsum-fishing-mortality-relative.PNG"), height=4, width=9)


ggplot(kit_mort_dist, aes(x = Year, y = PercentMortality, fill = FisheryTypeGrouped)) +
  geom_col() + 
  labs(y = "Fishing mortality") +
  scale_fill_brewer(palette = "Dark2") +
  theme_sleek()  
ggsave(here("plots/Skeena/kitsum-fishing-mortality-contribution.PNG"), height=4, width=9)
