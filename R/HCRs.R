# visualize realized HCR for Skeena and Nass aggregates over time

library(tidyverse)
library(here)
library(ggrepel)
source(here("R/data_functions.R"))
benchmarks <- read.csv(here("data/generated/bench_par_table_forTWG.4Aug2026.csv"))
theme_sleek <- function(base_size = 11, base_family = "") {
  half_line <- base_size/2
  t <- theme_light(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(colour = "grey30"),
      strip.text.y = element_text(colour = "grey30"),
      axis.text = element_text(colour = "grey30"),
      axis.title = element_text(colour = "grey30"),
      legend.title = element_text(colour = "grey30", size = rel(0.9)),
      panel.border = element_rect(fill = NA, colour = "grey70", linewidth = 1),
      legend.key.size = unit(0.9, "lines"),
      legend.text = element_text(size = rel(0.7), colour = "grey30"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA),
      plot.title = element_text(colour = "grey30", size = rel(1)),
      plot.subtitle = element_text(colour = "grey30", size = rel(.85))
    )
}

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

# visualize hypothetical HRCs ----

visualize_HCR <- function(EG=NULL, ER.cap=NULL, ER.floor=NULL, Max.Run=200000, int=1000) {
  out <- NULL
  
  for(i in seq(1,Max.Run, int)) {
    run.size = i
    catch <- ifelse(run.size<=EG, 0, run.size-EG)
    HR.all <- ifelse(run.size==0, 0, catch/run.size)
    
    itn <- data.frame(run_size=run.size,
                      HR = HR.all)
    out <- rbind(out, itn)
    
  }
  if(is.numeric(ER.floor) ){
    out <- out |>
      mutate(HR = if_else(HR < ER.floor, ER.floor, HR))
  }
  
  if(is.numeric(ER.cap) ){
    out <- out |>
      mutate(HR = if_else(HR > ER.cap, ER.cap, HR))
  }
  
  
  return(out)
}


simple.eg <- visualize_HCR(EG=56000)
eg.floor <- visualize_HCR(EG=56000, ER.floor = 0.1)
eg.floor.cap <- visualize_HCR(EG=56000, ER.cap=0.60,ER.floor = 0.1)
eg.floor.cap.lower <- visualize_HCR(EG=26000, ER.cap=0.60,ER.floor = 0.1)
eg.floor.cap.lower$HR2 <- eg.floor.cap.lower$HR
eg.floor.cap.lower$HR2[c(30:141)] <- seq(0.1,0.6, length.out=112)

ggplot(simple.eg, aes(x = run_size, y = HR*100)) +
  geom_line(lwd=1.5, col="dark grey")+
  labs(x="Run size", y="Exploitation rate (%)") +
  scale_y_continuous(limits = c(0, 100)) +
  theme_sleek()

ggplot(eg.floor, aes(x = run_size, y = HR*100)) +
  geom_line(lwd=1.5, col="dark grey")+
  labs(x="Run size", y="Exploitation rate (%)") +
  scale_y_continuous(limits = c(0, 100)) +
  theme_sleek()

ggplot(eg.floor.cap, aes(x = run_size, y = HR*100)) +
  geom_line(lwd=1.5, col="dark grey")+
  labs(x="Run size", y="Exploitation rate (%)") +
  scale_y_continuous(limits = c(0, 100)) +
  theme_sleek()

ggplot(eg.floor.cap.lower, aes(x = run_size, y = HR2*100)) +
  geom_line(lwd=1.5, col="dark grey")+
  labs(x="Run size", y="Exploitation rate (%)") +
  scale_y_continuous(limits = c(0, 100)) +
  theme_sleek()

