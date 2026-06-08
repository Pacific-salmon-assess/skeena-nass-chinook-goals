# script for wrangling Tseax genetics and reconstructing Tseax spawners over time 
library(tidyverse)
library(here)
library(ggsidekick)

# read in GSI data ----
Nass_gsi <- read.csv(here("data/nass_gsi.csv"))
Middle_nass_escape <- read.csv(here("data/Nass_SpHar_16Feb2026.csv"))

tseax <- matrix(NA,32,1000)

for(i in 1:1000){
  Tseax_gsi_comps <- Nass_gsi |>
    filter(Area == "MIDDLE") |>
    mutate(beta_dist_phi = gsi*(1-gsi)/(sd*sd)-1,
           beta_dist_alpha = gsi*beta_dist_phi,
           beta_dist_beta = (1-gsi)*beta_dist_phi)|>
    group_by(Year,stock) |>
    mutate(random_prop = rbeta(1,beta_dist_alpha,beta_dist_beta))|>
    group_by(Year) %>%
    mutate(proportion = random_prop / sum(random_prop, na.rm = TRUE)) %>%
    ungroup()|>
    filter(stock == "Tseax") |>
    rename(year=Year) |>
    select(year,proportion )
  
  avg_Tseax_prop = mean(Tseax_gsi_comps$proportion)
  sd_Tseax_prop = sd(Tseax_gsi_comps$proportion)
  beta_Tseax_phi = avg_Tseax_prop*(1-avg_Tseax_prop)/(sd_Tseax_prop*sd_Tseax_prop)-1
  beta_Tseax_alpha = avg_Tseax_prop*beta_Tseax_phi
  beta_Tseax_beta = (1-avg_Tseax_prop)*beta_Tseax_phi
  rbeta(1,beta_Tseax_alpha,beta_Tseax_beta)
  
  Tseax_esc <-  Middle_nass_escape |> 
    filter(CU == "Middle Nass") |>
    left_join(Tseax_gsi_comps, by = "year") |>
    mutate(proportion_infill = replace_na(proportion,rbeta(1,beta_Tseax_alpha,beta_Tseax_beta))) |>
    group_by(year) |>
    mutate(esc = rnorm(1,spwn,spwn*spwn_cv)*proportion_infill) |>
    select(year,esc)
  
  tseax[,i] <- as.vector(Tseax_esc$esc)
}

tseax_error <- t(apply(tseax,1,quantile, probs=(c(0.5,0.1,0.9))))
tseax_error <- as.data.frame(tseax_error)
tseax_error$year <- seq(1994,2025,by=1)

ggplot(tseax_error, aes(x = year, y = `50%`)) +
  annotate(geom = "rect", xmin = 2009.5, xmax = 2025, ymin = -Inf, ymax = Inf,
           fill = "lightblue", alpha = 0.5) +
  annotate(geom = "rect", xmin = 2006.5, xmax = 2007.5, ymin = -Inf, ymax = Inf,
           fill = "lightblue", alpha = 0.5) +  geom_ribbon(aes(ymin = `10%`, ymax = `90%`), fill = "light grey", alpha=0.85) +
  geom_line(size=1.5, color="grey20")+
  xlab("Year") +
  ylab("Tseax reconstructed spawner abundance") +
  theme_sleek()   

ggsave("plots/nass/tseax-reconstruction.jpeg",width=7,height=5)

mean(tseax_error$`50%`[27:32])/mean(tseax_error$`50%`)
quantile(tseax_error$`50%`, probs=seq(0.79,2.33))

# remove Tseax from Upper CU ---
upper_recon <- Middle_nass_escape |>
  filter(CU != "Lower Nass") |>
  group_by(year) |>
  summarize(spwn = sum(spwn),
            harv = sum(harv),
            spwn_cv = max(spwn_cv),
            CU="Upper Nass",
            SMU="Nass") |>
  left_join(tseax_error, by="year") |>
  mutate(spwn_corr = spwn-`50%`,
         prop_tseax=`50%`/spwn,
         harv_corr=harv-(harv*prop_tseax)) |>
  select(CU,year,spwn_corr,spwn_cv,harv_corr,SMU) |>
  rename(spwn=spwn_corr,
         harv=harv_corr)

# write escapement time series to file ----

write.csv(tseax_error,"data/Tseax-esc.csv")
write.csv(upper_recon,"data/Upper-esc-harv.csv",row.names = F)
