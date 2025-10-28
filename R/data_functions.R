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
  filter(CU != "Skeena") |> #remove aggregate
  mutate(harv = abs(harv)) ##FIX NEGATIVE Harvest! Need to ask luke about this... 

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

#bind skeena and nass
A_obs <- bind_rows(A_obs, nass_Aobs) |>
  filter(!is.na(a4))|>        #hope a4 alone is representative of NA rows
  mutate_if(is.numeric, round, 0)
  
sp_har <- bind_rows(sp_har, nass_sp_har) |>
  filter(!is.na(spwn))|>
  mutate_if(is.numeric, round, 0)


#check timeseries of complete dataset
A_obs |>
  group_by(CU) |>
  summarise(min(year), max(year))

sp_har |>
  group_by(CU) |>
  summarise(min(year), max(year))

#merge and filter dfs so they have complete years of data among CUs
  # could make this robust depending on our final data... 
sp_har <- left_join(sp_har, A_obs) |> 
  filter(!is.na(a4)) |> #remove NAs once joined
#  filter(year>=1992 & year<= 2019) |> #turn on if only doing even years
  mutate(spwn_cv = 0.5, #assumed CVs for now
         harv_cv = 0.3) |>
  as.data.frame()


#calc other indicies
a_min <- 4
a_max <- 6
A <- a_max - a_min + 1 #total age classes

#nyrs <- max(sp_har$year)-min(sp_har$year)+1 #number of years of observations
#nRyrs <- nyrs + A - 1 #number of recruitment years: unobserved age classes at ? to predict last year of spawners

yrs <- sp_har |>
  group_by(CU) |>
  summarise(nyrs = max(year)-min(year)+1, 
            nRyrs = nyrs + A - 1)

rm(A_obs, a_obs, nass_Aobs, nass_sp_har, skeena, i)



# functions ------------------------------------------------------------------------------
load <- base::load # make sure renv::load() does not mask base::load()

my.ggsave <- function(filename = default_name(plot), plot = last_plot(),
                      width= 9, height = 5.562, dpi= 180){
  ggsave(filename=filename, plot = last_plot(), width=width, height=height, dpi=dpi, bg="white")
}

# benchmark functions ---
get_Smsy <- function(a, b){
  Smsy <- (1-lambert_W0(exp(1-a)))/b
  if(Smsy <0){Smsy <- 0.001} #dumb hack for low draws so Smsy doesnt go negative
  return(Smsy)
}

get_Sgen <- function(a, b, int_lower, int_upper, Smsy){
  fun_Sgen <- function(Sgen, a, b, Smsy){Sgen*a*exp(-b*Sgen)-Smsy}
  Sgen <- uniroot(fun_Sgen, interval=c(int_lower, int_upper), a=a, b=b, Smsy=Smsy)$root
  return(Sgen)
}

# below funs adapted from BC's kusko code (https://github.com/brendanmichaelconnors/Kusko-harvest-diversity-tradeoffs/blob/master/functions.R#L237)

#------------------------------------------------------------------------------#
# Status function (to estimate whether stock is overfished or predicted to go
#  extinct at a given harvest rate, over the long-term)
#------------------------------------------------------------------------------#
# U <- harvest rate
# a <- productivity (Ricker a parameter)
# b <- density dependence (Ricker beta parameter)
SC.eq <- function(U,a,b){ ## Think we can delete this?
  a <- log(a)
  S.eq <- max(0,(a-(-log(1-U)))/b)
  C.eq <- max(0,((a-(-log(1-U)))/b)*exp(a-b*((a-(-log(1-U)))/b))-((a-(-log(1-U)))/b))
  OF <- ifelse(U>0.5*a-0.07*a^2,1,0)
  EX <- ifelse(S.eq==0,1,0)
  return(c(S.eq,C.eq,OF,EX))
}

#take a slice of the posterior ---
process.iteration = function(samp) {
  # 1.) extract names
  nms = names(samp)
  A = 4
  ns = length(unique(sp_har$CU))
  
  # 2.) extract elements according to the names and put them into the appropriate data structure
  # parameters
  alpha = unname(samp[substr(nms, 1, 5) == "alpha"])
  beta = unname(samp[substr(nms, 1, 4) == "beta"])
  last_resid = unname(samp[substr(nms, 1, 10) == "last_resid"])
  #phi = unname(samp["phi"])
  #Sigma_R = matrix(samp[substr(nms, 1, 7) == "Sigma_R"], ns, ns)
  pis = c(as.numeric(samp["pi_1"]), as.numeric(samp["pi_2"]), as.numeric(samp["pi_3"]), as.numeric(samp["pi_4"]))
  
  # states
  S = matrix(samp[substr(nms, 1, 2) == "S_"], A, ns)
  R = matrix(samp[substr(nms, 1, 2) == "R_"], A - 1, ns)
  
  # 3.) create output list
  output = list(
    alpha = alpha,
    beta = beta,
    last_resid = last_resid,
    S = S,
    R = R,
    pis = pis
  )
  # 4.) return output
  return(output)
}
