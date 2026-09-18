# visualize realized HCR for Skeena aggregate over time

library(tidyverse)
library(here)
library(ggrepel)
source(here("R/data_functions.R"))
benchmarks <- read.csv(here("data/generated/bench_par_table_forTWG.4Aug2026.csv"))

# Skeena -----
agg_RR <- sp_har |>
  filter(CU == "Skeena Aggregate") |>
  mutate(Run = spwn + harv,
           ER = harv/Run)


ggplot(data = agg_RR, aes(x = Run, y = ER, col = year)) +
  geom_point(size=3) +
  labs(x = "Run-size", y = "Exploitation rate") +
  coord_cartesian(ylim=c(0,1), xlim=c(0,250000)) +
  theme_sleek() +
  geom_text_repel(aes(label = year), size=3)


ggsave("plots/Skeena/realized-HCR.jpeg", width = 6, height=4,units="in", dpi=600)

agg_Umsy <- benchmarks |>
  filter(CU == "Skeena Aggregate",
         par == "Umsy")
agg_Smsy <- benchmarks |>
  filter(CU == "Skeena Aggregate",
         par == "Smsy")

ggplot(data = agg_RR, aes(x = spwn, y = ER, col = year)) +
  geom_point() +
  labs(x = "Spawners", y = "Exploitation rate") +
  coord_cartesian(ylim=c(0,1), xlim=c(0,150000)) +
  theme_sleek() +
  geom_hline(yintercept = agg_Umsy$median, lty=2, col = "grey") +
  geom_vline(xintercept = agg_Smsy$median, lty=2, col = "grey") 

ggplot(df, aes(x = x, y = y)) +
  geom_point(size = 3) +
  

# Nass -----
agg_RR <- sp_har |>
  filter(CU == "Nass Aggregate") |>
  mutate(Run = spwn + harv,
         ER = harv/Run)


ggplot(data = agg_RR, aes(x = Run, y = ER, col = year)) +
  geom_point(size=3) +
  labs(x = "Run-size", y = "Exploitation rate") +
  coord_cartesian(ylim=c(0,1), xlim=c(0,60000)) +
  theme_sleek()+
  geom_text_repel(aes(label = year), size=3)

ggsave("plots/Nass/realized-HCR.jpeg", width = 6, height=4,units="in", dpi=600)


agg_Umsy <- benchmarks |>
  filter(CU == "Nass Aggregate",
         par == "Umsy")
agg_Smsy <- benchmarks |>
  filter(CU == "Nass Aggregate",
         par == "Smsy")

ggplot(data = agg_RR, aes(x = spwn, y = ER, col = year)) +
  geom_point() +
  labs(x = "Spawners", y = "Exploitation rate") +
  coord_cartesian(ylim=c(0,1), xlim=c(0,50000)) +
  theme_sleek() +
  geom_hline(yintercept = agg_Umsy$median, lty=2, col = "grey") +
  geom_vline(xintercept = agg_Smsy$median, lty=2, col = "grey") 

