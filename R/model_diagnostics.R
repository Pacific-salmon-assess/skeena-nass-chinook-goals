library(tidyverse)
library(here)
library(rstan)

source(here("R/data_functions.R"))

# read in model fits
AR1.fits <- lapply(list.files(here("data/generated/model_fits/AR1"),
                              full.names = T), readRDS)
names(AR1.fits) <- unique(sp_har$CU)

TVA.fits <- lapply(list.files(here("data/generated/model_fits/TVA"),
                              full.names = T), readRDS)
names(TVA.fits) <- unique(sp_har$CU)

# make big summary df of pars
AR1.summaries <- NULL
TVA.summaries <- NULL

for(i in 1:length(unique(sp_har$CU))){
  AR1.summaries <- bind_rows(AR1.summaries, as.data.frame(rstan::summary(AR1.fits[[i]])$summary) |>
                               mutate(CU = unique(sp_har$CU)[i]))
  TVA.summaries <- bind_rows(TVA.summaries, as.data.frame(rstan::summary(TVA.fits[[i]])$summary) |>
                               mutate(CU = unique(sp_har$CU)[i]))
}
