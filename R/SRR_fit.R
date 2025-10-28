library(tidyverse)
library(here)
library(rstan)

source(here("R/data_functions.R"))

set.seed(2)

#refit <- FALSE #toggle T/F if you want to refit models
refit <- TRUE

# fit AR1 and time varying productivity (TVA) models--------------------------------------
if(refit == TRUE){
  for(i in unique(sp_har$CU)){
    
    sp_har1 <- filter(sp_har, CU == i)
    
    A_obs1 <- select(sp_har1, a4:a6) |>
      as.matrix()
    
    stan.data <- list("nyrs" = nyrs,
                      "a_min" = a_min,
                      "a_max" = a_max,
                      "A" = A,
                      "nRyrs" = nRyrs,
                      "A_obs" = A_obs1,
                      "S_obs" = sp_har1$spwn,
                      "H_obs" = sp_har1$harv,
                      "S_cv" = sp_har1$spwn_cv,
                      "H_cv" = sp_har1$harv_cv,
                      "Smax_p" = 0.75*max(sp_har1$spwn), #data for priors in semi_inform models, can tinker based on what assumed Smax is
                      "Smax_p_sig" = 0.75*max(sp_har1$spwn))
    
    #AR1.fit <- stan(file = here("Stan/SS-SR_AR1.stan"),
    #                data = stan.data,
    #                cores = 4,
    #                seed = 2,
    #                iter = 4000)
    
    #saveRDS(AR1.fit, here("data/generated/model_fits/AR1",
    #                      paste0(i, "_AR1.rds")))
    
    TV.fit <- stan(file = here("Stan/SS-SR_TVA.stan"),
                   data = stan.data,
                   cores = 4,
                   seed = 2,
                   iter = 4000,
                   control = list(adapt_delta = 0.999,
                                  max_treedepth = 20))
    
    saveRDS(TV.fit, here("data/generated/model_fits/TVA/",
                         paste0(i, "_TVA.rds")))
  }
}
