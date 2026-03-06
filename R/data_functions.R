library(tidyverse)
library(here)
#library(devtools) #comment out below if you need to install skrunchy
#install_github("Pacific-salmon-assess/skrunchy")
#install_github("lukewarkentin/skrunchy2025")
library(skrunchy2025)

#-----------------------------------------------------------------------------------------
# read in data for SRR fit. Skeena data comes from Luke Warkentin (DFO)'s package
  #and Nass data from Ian Beveridge & Richard Alexander (LGL)

sp_har <- skrunchy2025::run_reconstruction_table |>
  select(i_population, y_return_year, W_star_wild_spawners, total_harvest_estimate) |>
  rename(CU = i_population, 
         year = y_return_year) |>
  group_by(CU, year) |>
  summarise(spwn = round(sum(W_star_wild_spawners, na.rm = TRUE)), #sum across ages
            harv = round(sum(total_harvest_estimate, na.rm = TRUE))) |>
  mutate(SMU = "Skeena", 
         CU = ifelse(CU == "Skeena", "Skeena Aggregate", CU)) |>
  arrange(CU, year) |>
  as.data.frame()

A_obs <- NULL #empty object to bind to
for(i in 1:dim(skrunchy2025::n_age_observations)[1]){
  a_obs <- as.data.frame(skrunchy2025::omega$omega[i,,]) |>
    mutate(CU = dimnames(skrunchy2025::omega$omega)$i[i])
  a_obs <- mutate(a_obs, year = rownames(a_obs))
  A_obs <- bind_rows(A_obs, a_obs)
}
 
A_obs <- A_obs |>
  mutate(`6` = `6` +`7`) |>
  rename(a4 = `4`, 
         a5 = `5`, 
         a6 = `6`) |>
  select(CU, year, a4, a5, a6) |>
  arrange(CU, year) |>
  mutate(SMU = "Skeena", 
         a4 = a4*100, #* all by 100 to have a 100 fish balanced sample size for now 
         a5 = a5*100, 
         a6 = a6*100, 
         year = as.numeric(year), 
         CU = ifelse(CU == "Skeena", "Skeena Aggregate", CU))

  rownames(A_obs) <- NULL

# read in Nass data - this data was processed by copy & pasting LGL's data and calculating 
  # new vars (i.e. combining age groups and calc'ing harvest) all in excel as a placeholder

# age data is aggregate for nass - no CU level data
nass_A_obs <- read.csv(here("data/Nass_Aobs_16Feb2026.csv"))

#make dataset long and add CU names since they are shared among all CUs 
  #hacky fix so it lines up for how skeena code is written
nass_A_obs <- bind_rows(nass_A_obs, nass_A_obs, nass_A_obs, nass_A_obs) |> 
  mutate(CU = c(rep("Upper Nass", nrow(nass_A_obs)), 
                rep("Middle Nass", nrow(nass_A_obs)), 
                rep("Lower Nass", nrow(nass_A_obs)), 
                rep("Nass Aggregate", nrow(nass_A_obs))), 
         SMU = "Nass") |>
  rename(year = Year)

nass_sp_har <- read.csv(here("data/Nass_SpHar_16Feb2026.csv"))

#check years so we know to put only fully observed in aggregate
nass_sp_har |>
  group_by(CU) |>
  summarise(min(year), 
            max(year))

nass_sp_har_agg <- nass_sp_har |>
  filter(year >= 1994) |> #based on above
  ##ASSUMPTION TO ADDRESS BELOW - infill NAs, then take mean?
  mutate(spwn_cv = ifelse(is.na(spwn_cv), 0.5, spwn_cv)) |>  #0.5: upper end of empirical CVs
  group_by(SMU, year) |>
  summarise(spwn = sum(spwn), 
            harv = sum(harv), 
            spwn_cv = mean(spwn_cv)) |>
  mutate(CU = "Nass Aggregate")

nass_sp_har <- bind_rows(nass_sp_har, nass_sp_har_agg)

#bind skeena and nass
prop_age <- bind_rows(A_obs, nass_A_obs) |> #for plotting later
  mutate(a4 = a4/100, 
         a5 = a5/100, 
         a6 = a6/100) |>
  relocate(SMU, .after = 1)

if(FALSE){ #NEED TO FIX THIS LATER - IT'S DROPPING AGES
  # correct Skeena for uneven ESS - STILL NEED TO DO FOR NASS! 
  ESS_Skeena <- t(skrunchy2025::n_age_observations[,,1] + skrunchy2025::n_age_observations[,,2] + skrunchy2025::n_age_observations[,,2]) |>
    as.data.frame() |>
    mutate(year = as.numeric(rownames(t(skrunchy2025::n_age_observations[,,1])))) |>
    pivot_longer(-year, names_to = "CU") |>
    rename(ESS = value)
  
  A_obs <- A_obs |>
    left_join(ESS_Skeena, by = c("CU", "year")) |>
    mutate(a4 = (a4/100)*ESS, #revet fromm scaled to 100 to calculated ESS. 
           a5 = (a5/100)*ESS, 
           a6 = (a6/100)*ESS) |>
    select(-ESS)
}


A_obs <- bind_rows(A_obs, nass_A_obs) |>
  mutate_if(is.numeric, round, 0)

sp_har <- bind_rows(sp_har, nass_sp_har) |>
  filter(!is.na(spwn)) |>
  relocate(SMU, .after=1)

#check timeseries of complete dataset
A_obs |>
  group_by(CU) |>
  summarise(min(year), max(year))

sp_har |>
  group_by(CU) |>
  summarise(min(year), max(year))

#merge and filter dfs so they have complete years of data among CUs
  # could make this robust depending on our final data... 
sp_har <- left_join(sp_har, A_obs, by = c("CU", "year", "SMU")) |> 
  filter(!is.na(a4)) |> #remove years of NA age obs once joined
  mutate(harv_cv = 0.3, #assumed harvest CV
         spwn_cv = replace_na(spwn_cv, 0.5)) |> #replace unknown CVs with assumed
  arrange(CU, year) |> 
  relocate(harv_cv, .after = spwn_cv) |>
  as.data.frame() 

#calc other indicies
a_min <- 4
a_max <- 6
A <- a_max - a_min + 1 #total age classes

yrs <- sp_har |>
  group_by(CU) |>
  summarise(nyrs = max(year)-min(year)+1, 
            nRyrs = nyrs + A - 1) ## want to be certain this is A not a_min

rm(A_obs, a_obs, nass_A_obs, nass_sp_har, nass_sp_har_agg, i)


# functions ------------------------------------------------------------------------------
load <- base::load # make sure renv::load() does not mask base::load()

my.ggsave <- function(filename = default_name(plot), plot = last_plot(),
                      width= 9, height = 5.562, dpi= 180){
  ggsave(filename=filename, plot = last_plot(), width=width, height=height, dpi=dpi, bg="white")
}


# theme from https://github.com/seananderson/ggsidekick/blob/master/R/theme_sleek.R
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
  # Only set tagger theme element if tagger package is available
  if (requireNamespace("tagger", quietly = TRUE)) {
    t <- t + theme(tagger.panel.tag.text = element_text(colour = "grey30"))
  }
  t
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
