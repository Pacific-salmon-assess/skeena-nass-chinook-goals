library(tidyverse)
library(here)
library(rstan)

source(here("R/data_functions_new.R")) #try with informed A-obs

#refit <- FALSE #toggle T/F if you want to refit models
refit <- TRUE

# fit AR1 and time varying productivity (TVA) models--------------------------------------
if(refit == TRUE){
  for(i in unique(sp_har$CU)){
    
    sp_har1 <- filter(sp_har, CU == i)
    
    yrs1 <- filter(yrs, CU == i)
    
    A_obs1 <- select(sp_har1, a4:a6) |>
      as.matrix()
    
    stan.data <- list("nyrs" = yrs1$nyrs,
                      "a_min" = a_min,
                      "a_max" = a_max,
                      "A" = A,
                      "nRyrs" = yrs1$nRyrs,
                      "A_obs" = A_obs1,
                      "S_obs" = sp_har1$spwn,
                      "H_obs" = sp_har1$harv,
                      "S_cv" = sp_har1$spwn_cv,
                      "H_cv" = sp_har1$harv_cv,
                      "prior_S_max" = max(sp_har1$spwn)) #data for Smax priors
    
    cat("Fitting AR1 to", i, "===================================\n")

    AR1.fit <- stan(file = here("Stan/SS-SR-AR1.stan"),
                    data = stan.data,
                    seed = 1,
                    iter = 4000,
                    refresh = 0, #hide long list of timing updates
                    control = list(adapt_delta = 0.995,
                                   max_treedepth = 20))
    
    check_divergences(AR1.fit)
    
    saveRDS(AR1.fit, here("data/generated/model_fits/AR1",
                          paste0(i, "_AR1.rds")))
    
    cat("Fitting TVA to", i, "===================================\n")

    TV.fit <- stan(file = here("Stan/SS-SR-TVA.stan"),
                   data = stan.data,
                   seed = 1,
                   iter = 4000,
                   refresh = 0,
                   control = list(adapt_delta = 0.995,
                                  max_treedepth = 20))
    check_divergences(TV.fit)
    
    saveRDS(TV.fit, here("data/generated/model_fits/TVA/",
                         paste0(i, "_TVA.rds")))
  }
}

#fit IWAM prior models--------------------------------------------------------------------

#only going to do this with an AR1 version for now
IWAM_sp_har <- filter(sp_har, CU %in% IWAM_priors$CU) |>
  arrange(CU)

#check these CUs line up
cbind(unique(IWAM_sp_har$CU),
      IWAM_priors$CU)

if(TRUE){
  for(i in unique(IWAM_sp_har$CU)){
    cat("Fitting IWAM AR1 to", i, "===================================\n")

    sp_har1 <- filter(IWAM_sp_har, CU == i)
    
    yrs1 <- filter(yrs, CU == i)
    
    A_obs1 <- select(sp_har1, a4:a6) |>
      as.matrix()
    
    IWAM_prior <- filter(IWAM_priors, CU == i)
    
    IWAM.stan.data <- list("nyrs" = yrs1$nyrs,
                      "a_min" = a_min,
                      "a_max" = a_max,
                      "A" = A,
                      "nRyrs" = yrs1$nRyrs,
                      "A_obs" = A_obs1,
                      "S_obs" = sp_har1$spwn,
                      "H_obs" = sp_har1$harv,
                      "S_cv" = sp_har1$spwn_cv, 
                      "H_cv" = sp_har1$harv_cv,
                      "IWAM_prior" = IWAM_prior$SMAX_mean, #data for Smax priors
                      "IWAM_prior_sd" = IWAM_prior$sd) 
    
    IWAM.AR1.fit <- stan(file = here("Stan/SS-SR-AR1-IWAM-prior.stan"),
                    data = IWAM.stan.data,
                    seed = 1,
                    iter = 4000,
                    refresh = 0,
                    control = list(adapt_delta = 0.995,
                                   max_treedepth = 20))
    
    check_divergences(IWAM.AR1.fit)
    
    saveRDS(IWAM.AR1.fit, here("data/generated/model_fits/IWAM",
                          paste0(i, "_IWAM_AR1.rds")))
  }
}
