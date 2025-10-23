library(tidyverse)
library(here)

#-----------------------------------------------------------------------------------------
# read in data for SRR fit. Skeena data comes from Luke Warkentin (DFO) and Nass data from
  # Ian Beveridge & Richard Alexander (LGL)

skeena <- readRDS(here("data/twg_data_for_examples_Nov2025.rds"))

#NEED TO check on and account for `complete_brood_year` in table below
sp_har <- as.data.frame(skeena$run_recon) |>
  select(population, return_year, wild_spawners, total_harvest_estimate) |>
  rename(CU = population, 
         year = return_year) |>
  group_by(CU, year) |>
  summarise(spwn = round(sum(wild_spawners),0), #sum across ages
            harv = round(sum(total_harvest_estimate),0)) |>
  mutate(SMU = "Skeena") |>
  arrange(CU, year) |>
  filter(CU != "Skeena") #remove aggregate

A_obs <- NULL
for(i in 1:dim(skeena$age_observations)[1]){
  a_obs <- as.data.frame(skeena$age_comps_arr[i,,]) |>
    mutate(CU = dimnames(skeena$age_comps_arr)$i[i])
  a_obs <- mutate(a_obs, year = rownames(a_obs))
  A_obs <- bind_rows(A_obs, a_obs)
}

A_obs <- A_obs |>
  rename(a4 = `4`, 
         a5 = `5`, 
         a6 = `6`) |>
  select(CU, year, a4, a5, a6) |>
  arrange(CU, year) |>
  mutate(SMU = "Skeena")  |>
  filter(CU != "Skeena") #remove aggregate
rownames(A_obs) <- NULL

a_min <- 4
a_max <- 6
nyrs <- max(sp_har$year)-min(sp_har$year)+1 #number of years of observations
A <- a_max - a_min + 1 #total age classes
nRyrs <- nyrs + A - 1 #number of recruitment years: unobserved age classes at ? to predict last year of spawners
