#processing model fits
library(here)
library(tidyverse)
library(gsl)
source(here("R/data_functions.R"))

# read in data ---------------------------------------------------------------------------
# model fits ---
AR1.fits <- lapply(list.files(here("data/generated/model_fits/AR1"),
                              full.names = T), readRDS)
names(AR1.fits) <- unique(sp_har$CU)

TVA.fits <- lapply(list.files(here("data/generated/model_fits/TVA"),
                              full.names = T), readRDS)
names(TVA.fits) <- unique(sp_har$CU)

# process data and fits to make plots later ----------------------------------------------
bench.par.table <- NULL #empty objects to rbind CU's outputs to
bench.posts <- NULL
par.posts <- NULL
par.posts.tva <- NULL ##need to reshape to account for changing years
tva.par.summary.out <- NULL
SR.preds <- NULL
AR1.spwn <- NULL
AR1.harv <- NULL
AR1.ages <- NULL
AR1.resids <- NULL
brood.all <- NULL
brood.all.long <- NULL
a.yrs.all <- NULL
TV.resids <- NULL
TV.SR.preds <- NULL
TV.spwn <- NULL
TV.harv <- NULL

#quantiles for reporting
probs_80 <- c(0.1,0.5,0.9)
probs_50 <- c(0.25,0.5,0.75)
probs_80_50 <- c(0.1,0.25,0.5,0.75,0.9)

for(i in unique(sp_har$CU)){ # Loop over CUs to process model outputs
  
  # AR1 (base S-R) models ----------------------------------------------------------------
  sub_dat <- filter(sp_har, CU==i)
  nyrs <- filter(yrs, CU==i)$nyrs
  nRyrs <- filter(yrs, CU==i)$nRyrs
  
  sub_pars <- rstan::extract(AR1.fits[[i]])
  # latent states
  AR1.spwn <- rbind(AR1.spwn, bind_cols(t(apply(sub_pars$S, 2, quantile, probs = probs_50)),
                                        unique(sub_dat$year),
                                        i))
  AR1.harv <- rbind(AR1.harv, bind_cols(t(apply(sub_pars$H, 2, quantile, probs = probs_50)),
                                        unique(sub_dat$year),
                                        i))
  
  #latent states of spawners and recruits---
  spwn.quant <- apply(sub_pars$S, 2, quantile, probs = probs_80)[,1:(nyrs-a_min)]
  rec.quant <- apply(sub_pars$R, 2, quantile, probs = probs_80)[,(a_max+1):nRyrs]
  
  brood_t <- as.data.frame(cbind(sub_dat$year[1:(nyrs-a_min)], t(spwn.quant), t(rec.quant))) |>
    round(2) 
  
  colnames(brood_t) <- c("BroodYear","S_lwr","S_med","S_upr","R_lwr","R_med","R_upr")
  
  brood_t <- mutate(brood_t, CU = i)
  
  brood.all <- rbind(brood.all, brood_t)
  
  spwn.quant.long <- apply(sub_pars$S, 2, quantile, probs = probs_80)[,1:(nyrs-a_min)]
  rec.quant.long <- apply(sub_pars$R, 2, quantile, probs = probs_80)[,(A+a_min):nRyrs]
  
  brood_t.long <- as.data.frame(cbind(sub_dat$year[1:(nyrs-a_min)],t(spwn.quant.long), t(rec.quant.long))) |>
    round(2)
  colnames(brood_t.long) <- c("BroodYear","S_lwr","S_med","S_upr","R_lwr","R_med","R_upr")
  
  brood_t.long <- mutate(brood_t.long, CU = i)
  
  brood.all.long <- rbind(brood.all.long, brood_t.long)
  
  #SR relationship based on full posterior---
  spw <- seq(0,max(brood_t$S_upr),length.out=100)
  SR.pred <- matrix(NA,length(spw), length(sub_pars$lnalpha))
  bench <- matrix(NA,length(sub_pars$lnalpha),9,
                  dimnames = list(seq(1:length(sub_pars$lnalpha)), c("Sgen", "Smsy", "Umsy", 
                                                                     "Seq", "Smsr", "S.recent",
                                                                     "Smsr.20","Smsr.40", "Smsy.80")))
  
  par <- matrix(NA,length(sub_pars$lnalpha),4,
                dimnames = list(seq(1:length(sub_pars$lnalpha)), c("sample","ln_a","beta", "S_max")))
  
  #latent states of proportions at age --- 
  if(TRUE){
    sub_age <- NULL
    for(j in 1:A){
      sub_age <- bind_rows(sub_age, 
                           data.frame(t(apply(sub_pars$p[,,j], 2, quantile, probs = probs_50)), 
                                      age = j+a_min-1, 
                                      year = (min(sub_dat$year)-2):max(sub_dat$year), ##NEED TO CHECK INDEX, WHY 2?
                                      CU = i))
    }
    
    AR1.ages <- bind_rows(AR1.ages, sub_age)
  }
  # get benchmarks & pars (AR1 model)------------------------------------------------------------------
  
  for(j in 1:length(sub_pars$lnalpha)){
    ln_a <- sub_pars$lnalpha[j]
    b <- sub_pars$beta[j]
    S_max <- sub_pars$S_max[j]
    
    SR.pred[,j] <- (exp(ln_a)*spw*exp(-b*spw))
    
    bench[j,2] <- get_Smsy(ln_a, b)
    bench[j,1] <- get_Sgen(exp(ln_a),b,-1,1/b*2, bench[j,2])
    bench[j,3] <- (1 - lambert_W0(exp(1 - ln_a))) #U_MSY
    bench[j,4] <- ln_a/b #S_eq
    bench[j,5] <- 1/b #S_msr
    bench[j,6] <- exp(mean(log(sub_pars$S[j, (nyrs-5):nyrs]))) #S recent - mean spawners in last generation
    bench[j,7] <- (1/b)*0.2
    bench[j,8] <- (1/b)*0.4
    bench[j,9] <- get_Smsy(ln_a, b)*0.8 # 80% Smsy
    
    par[j,1] <- j
    par[j,2] <- ln_a
    par[j,3] <- b
    par[j,4] <- S_max
    
  }
  
  SR.pred <- as.data.frame(cbind(spw,t(apply(SR.pred, 1, quantile, probs = probs_80, na.rm=T))))|>
    round(2) |>
    mutate(CU = i, 
           SMU = unique(sub_dat$SMU))
  
  SR.preds <- rbind(SR.preds, SR.pred)
  
  bench.posts <- rbind(bench.posts, as.data.frame(bench) |> mutate(CU = i))
  
  par.posts <- rbind(par.posts, as.data.frame(par) |> mutate(CU = i))
  
  bench.quant <- apply(bench[,1:9], 2, quantile, probs = probs_80, na.rm=T) |>
    t()
  
  mean <- apply(bench[,1:9],2,mean, na.rm=T) #get means of each
  
  sub_benchmarks <- cbind(bench.quant, mean) |>
    as.data.frame() |>
    mutate(CU = i) |>
    relocate('50%', 1)
  
  #other pars to report
  alpha <- quantile(exp(sub_pars$lnalpha), probs = probs_80)
  beta <- quantile(sub_pars$beta, probs = probs_80)
  sigma <- quantile(sub_pars$sigma_R, probs = probs_80)
  phi <- quantile(sub_pars$phi, probs = probs_80)
  
  par.quants <- rbind(alpha, beta, sigma, phi)
  
  #make big table of bench and pars
  par.summary <- as.data.frame(rstan::summary(AR1.fits[[i]])$summary) |>
    select(mean, n_eff, Rhat)
  
  #summarise not other pars...
  par.summary <- filter(par.summary, row.names(par.summary) %in% c('lnalpha', 'beta',
                                                                   'sigma_R', 'phi')) |>
    mutate(CU = i)
  par.summary[1,1] <- exp(par.summary[1,1]) #exp ln_alpha
  
  pars <- cbind(par.quants, par.summary)
  
  sub.bench.par.table <- bind_rows(sub_benchmarks, pars) |>
    mutate(n_eff = round(n_eff, 0),
           Rhat = round(Rhat, 4))
  
  sub.bench.par.table <- mutate(sub.bench.par.table, bench.par = rownames(sub.bench.par.table))
  
  bench.par.table <- bind_rows(bench.par.table, sub.bench.par.table)
  
  # then residuals---
  resid.quant <- apply(sub_pars$lnresid, 2, quantile, probs = probs_80_50)[,(A):nRyrs]
  
  resids <- as.data.frame(cbind(sub_dat$year, t(resid.quant))) |>
    mutate(CU = i)
  colnames(resids) <- c("year","lwr","midlwr","mid","midupr","upr", "CU")
  
  AR1.resids <- rbind(AR1.resids, resids)
  
  # TVA (time varying alpha) models ------------------------------------------------------
  sub_pars_TVA <- rstan::extract(TVA.fits[[i]])
  par_TVA <- matrix(NA,length(sub_pars_TVA$beta),nyrs+2, 
                    dimnames = list(seq(1:length(sub_pars_TVA$beta)),c("sample", seq(1:nyrs),"beta")))
  par_TVA[,1] <- seq(1:length(sub_pars_TVA$beta))
  par_TVA[,c(2:(nyrs+1))] <- sub_pars_TVA$ln_alpha
  par_TVA[,c(nyrs+2)] <- sub_pars_TVA$beta
  
  par.posts.tva <- bind_rows(par.posts.tva, as.data.frame(par_TVA) |> mutate(CU = i))
  
  tva.par.summary <- rstan::summary(TVA.fits[[i]], pars=c("ln_alpha", "beta", "pi", "alpha_dev"),
                                    probs = probs_80)$summary
  tva.par.summary.out <- rbind(tva.par.summary.out, as.data.frame(tva.par.summary) |> mutate(CU = i))
  a.yrs <- apply(exp(sub_pars_TVA$ln_alpha), 2, quantile, probs = probs_80)
  a.yrs <- as.data.frame(cbind(sub_dat$year, t(a.yrs)))
  
  colnames(a.yrs) <- c("brood_year", "lwr", "mid", "upr")
  
  a.yrs.all <- rbind(a.yrs.all, data.frame(a.yrs, CU = i))
  
  #time varying alpha residuals
  resid.quant <- apply(sub_pars_TVA$lnresid, 2, quantile,
                       probs = probs_80_50)[,(A):nRyrs]
  
  resids <- as.data.frame(cbind(sub_dat$year, t(resid.quant))) |>
    mutate(CU = i)
  TV.resids <- bind_rows(TV.resids, resids)
  
  # median SR pred fit by year
  for(j in 1:ncol(sub_pars_TVA$ln_alpha)){
    ln_alpha_yr <- sub_pars_TVA$ln_alpha[,j]
    pred <- data.frame(pred.R = (exp(median(ln_alpha_yr))*spw*exp(-median(sub_pars_TVA$beta)*spw))) |>
      mutate(CU = i,
             year = unique(sp_har$year)[j]) |>
      cbind(spw)
    TV.SR.preds <- rbind(TV.SR.preds, pred)
  }
  
  #storing spawner quantiles for fwd sim plot
  TV.spwn.quant <- data.frame(t(apply(sub_pars_TVA$S, 2, quantile, probs = probs_50))) |>
    mutate(CU =i,
           year = unique(sub_dat$year))
  
  TV.spwn <- bind_rows(TV.spwn, TV.spwn.quant)
  
  TV.harv.quant <- data.frame(t(apply(sub_pars_TVA$C, 2, quantile, probs = probs_50))) |>
    mutate(CU =i,
           year = unique(sub_dat$year))
  
  TV.harv <- bind_rows(TV.harv, TV.harv.quant)
}  # End data wrangling loop by CU


colnames(SR.preds) <- c("Spawn", "Rec_lwr","Rec_med","Rec_upr", "CU", "SMU")
colnames(AR1.resids) <- c("year","lwr","midlwr","mid","midupr","upr", "CU")
colnames(AR1.spwn) <- c("S.25", "S.50", "S.75", "year", "CU")
colnames(AR1.harv) <- c("H.25", "H.50", "H.75", "year", "CU")
colnames(AR1.ages) <- c("A.25", "A.50", "A.75", "age", "year", "CU")
colnames(TV.resids) <- c("year","lwr","midlwr","mid","midupr","upr", "CU")
colnames(TV.spwn) <- c("S.25", "S.50", "S.75", "CU", "year")
colnames(TV.harv) <- c("H.25", "H.50", "H.75", "CU", "year")

#add SMU to output
CU_SMU <- sp_har |>
  select(CU, SMU) |>
  distinct()

a.yrs.all <- left_join(a.yrs.all, CU_SMU, by = "CU") #coulda map()'d this i guess... 
AR1.ages <- left_join(AR1.ages, CU_SMU, by = "CU")
AR1.harv <- left_join(AR1.harv, CU_SMU, by = "CU")
AR1.resids <- left_join(AR1.resids, CU_SMU, by = "CU")
AR1.spwn <- left_join(AR1.spwn, CU_SMU, by = "CU")
bench.par.table <- left_join(bench.par.table, CU_SMU, by = "CU")
bench.posts <- left_join(bench.posts, CU_SMU, by = "CU")
brood.all <- left_join(brood.all, CU_SMU, by = "CU")
brood.all.long <- left_join(brood.all.long, CU_SMU, by = "CU")
par.posts <- left_join(par.posts, CU_SMU, by = "CU")
par.posts.tva <- left_join(par.posts.tva, CU_SMU, by = "CU")
TV.harv.quant <- left_join(TV.harv.quant, CU_SMU, by = "CU")
TV.resids <- left_join(TV.resids, CU_SMU, by = "CU")
TV.spwn <- left_join(TV.spwn, CU_SMU, by = "CU")
TV.SR.preds <- left_join(TV.SR.preds, CU_SMU, by = "CU")
tva.par.summary.out <- left_join(tva.par.summary.out, CU_SMU, by = "CU")

  # write important tables to repo ---------------------------------------------------------
bench.par.table <- bench.par.table |> # AR1 SR model pars and benchmarks
  relocate(CU, 1) |>
  relocate(SMU, .after = 1) |>
  relocate(bench.par, .after = 2) |>
  relocate(mean, .after = 3) |>
  mutate_at(4:7, ~round(.,5)) |>
  arrange(bench.par, CU)

colnames(bench.par.table) <- c("CU", "SMU", "par", "mean", "median", "lwr", "upr", "n_eff", "R_hat")

write.csv(bench.par.table, here("data/generated/bench_par_table.csv"),
          row.names = FALSE)

#remove helpers
rm(a.yrs, AR1.fits, bench, bench.quant, brood_t, brood_t.long, par, par_TVA, par.quants, 
   par.summary, pars, pred, rec.quant, rec.quant.long, resid.quant, resids, spwn.quant, 
   spwn.quant.long, SR.pred, sub_age, sub_dat, sub_pars, sub_pars_TVA, sub.bench.par.table,
   TV.harv, TV.spwn.quant, TVA.fits, tva.par.summary, alpha, b, beta, i, j, ln_a, 
   ln_alpha_yr, mean, phi, sigma, spw, sub_benchmarks)