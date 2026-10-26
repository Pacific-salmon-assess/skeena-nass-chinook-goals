# script for plotting tributary escapements 
library(tidyverse)
library(here)
library(ggsidekick)

# read in and wrangle data ----
skeena_tribs <- read.csv(here("data/Area_4_Chinook_QAQC_cm.csv"))

tribs <- skeena_tribs |>
  mutate(spawn = as.numeric(Recommended.Estimate)) |>
  filter(Year>1974)


ggplot(tribs, aes(x = Year, y =spawn )) +
  geom_line(size=1, color="grey30")+
  facet_wrap(~combined, scales = "free_y") +
  xlab("Year") +
  ylab("Spawners") +
  theme_sleek() +
  coord_cartesian(xlim = c(1975,2025)) 
  
ggsave("plots/skeena/tribs.jpeg",width=10,height=5)
