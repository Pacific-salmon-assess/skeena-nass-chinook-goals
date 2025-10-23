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
  summarise(spwn = sum(wild_spawners), #sum across ages
            harv = sum(total_harvest_estimate)) |>
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
  mutate(SMU = "Skeena", 
         a4 = a4*100, #* all by 100 to have a 100 fish balanced sample size for now 
         a5 = a5*100, 
         a6 = a6*100, 
         year = as.numeric(year))  |>
  filter(CU != "Skeena") #remove aggregate
rownames(A_obs) <- NULL

# read in Nass data - this data was processed by copy & pasting LGL's data and calculating 
  # new vars (i.e. combining age groups and calc'ing harvest) all in excel as a placeholder

nass_Aobs <- read.csv(here("data/Nass_Aobs.csv"))
nass_sp_har <- read.csv(here("data/Nass_SpHar.csv"))

#bind so we have 2 dfs for model fitting
A_obs <- bind_rows(A_obs, nass_Aobs) |>
  mutate_if(is.numeric, round, 0)
  
sp_har <- bind_rows(sp_har, nass_sp_har) |>
  filter(!is.na(spwn))|>
  mutate_if(is.numeric, round, 0)

#check timeseries of complete dataset
A_obs |>
  group_by(CU) |>
  summarise(min(year))

sp_har |>
  group_by(CU) |>
  summarise(min(year))

#filter dfs so they start in the oldest complete year of data
  # could make this robust depending on our final data... 
A_obs <- filter(A_obs, year>=1989)
sp_har <- filter(sp_har, year>=1989)

#calc other indicies
a_min <- 4
a_max <- 6
nyrs <- max(sp_har$year)-min(sp_har$year)+1 #number of years of observations
A <- a_max - a_min + 1 #total age classes
nRyrs <- nyrs + A - 1 #number of recruitment years: unobserved age classes at ? to predict last year of spawners


rm(a_obs, nass_Aobs, nass_sp_har, skeena, i)
