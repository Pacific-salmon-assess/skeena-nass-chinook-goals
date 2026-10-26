# visualize realized HCR for Skeena and Nass aggregates over time

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

# Skeena & Nass -----
agg_RR <- sp_har |>
  filter(CU %in% c("Skeena Aggregate","Nass Aggregate")) |>
  mutate(Run = spwn + harv,
         ER = harv/Run)


ggplot(data = agg_RR, aes(x = Run, y = ER, col = year)) +
  geom_point(size=3) +
  labs(x = "Run-size", y = "Exploitation rate") +
  coord_cartesian(ylim=c(0,1), xlim=c(0,NA)) +
  facet_wrap(~CU, scales = "free_x") +
  theme_sleek() +
  geom_text_repel(aes(label = year), size=3)


ggsave("plots/realized-HCR.jpeg", width = 10, height=4,units="in", dpi=600)

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
  theme_minimal()
ggsave("plots/eg-HCR.jpeg",width=5,height=3.5)

ggplot(eg.floor, aes(x = run_size, y = HR*100)) +
  geom_line(lwd=1.5, col="dark grey")+
  labs(x="Run size", y="Exploitation rate (%)") +
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal()
ggsave("plots/eg-floor-HCR.jpeg",width=5,height=3.5)

ggplot(eg.floor.cap, aes(x = run_size, y = HR*100)) +
  geom_line(lwd=1.5, col="dark grey")+
  labs(x="Run size", y="Exploitation rate (%)") +
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal()
ggsave("plots/eg-floor-ceiling-HCR.jpeg",width=5,height=3.5)

ggplot(eg.floor.cap.lower, aes(x = run_size, y = HR2*100)) +
  geom_line(lwd=1.5, col="dark grey")+
  labs(x="Run size", y="Exploitation rate (%)") +
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal()
ggsave("plots/ramp-HCR.jpeg",width=5,height=3.5)
