#inference_figs code from yukon... will cut out stuff to adapt
library(here)
library(tidyverse)
library(gsl)
library(ggcorrplot)
library(ggpubr)
library(viridis)
source(here("R/data_functions.R"))

# read in data ---------------------------------------------------------------------------
# model fits ---
AR1.fits <- lapply(list.files(here("data/generated/model_fits/AR1"),
                              full.names = T), readRDS)
names(AR1.fits) <- unique(sp_har$CU)

TVA.fits <- lapply(list.files(here("data/generated/model_fits/TVA"),
                              full.names = T), readRDS)
names(TVA.fits) <- unique(sp_har$CU)

# escapement estimates ---
  # these must be the latent states and their error? 
#esc <- read.csv(here("analysis/data/generated/esc-data.csv")) |>
#  mutate_at(2:6, as.numeric)

# process data and fits to make plots later ----------------------------------------------
bench.par.table <- NULL #empty objects to rbind CU's outputs to
bench.posts <- NULL
par.posts <- NULL
par.posts.tva <- NULL ##need to reshape to account for changing years
tva.par.summary.out <- NULL
SR.preds <- NULL
AR1.spwn <- NULL
AR1.harv <- NULL
AR1.resids <- NULL
brood.all <- NULL
brood.all.long <- NULL
a.yrs.all <- NULL
TV.resids <- NULL
TV.SR.preds <- NULL
TV.spwn <- NULL
TV.harv <- NULL


for(i in unique(sp_har$CU)){ # Loop over CUs to process model outputs
  
  # AR1 (base S-R) models ----------------------------------------------------------------
  sub_dat <- filter(sp_har, CU==i)
  nyrs <- filter(yrs, CU==i)$nyrs
  nRyrs <- filter(yrs, CU==i)$nRyrs
  
  sub_pars <- rstan::extract(AR1.fits[[i]])
  
  AR1.spwn <- rbind(AR1.spwn, bind_cols(t(apply(sub_pars$S, 2, quantile, c(0.25, .5, .75))),
                                        unique(sub_dat$year),
                                        i))
  
  AR1.harv <- rbind(AR1.harv, bind_cols(t(apply(sub_pars$H, 2, quantile, c(0.25, .5, .75))),
                                        unique(sub_dat$year),
                                        i))
  
  #latent states of spawners and recruits---
  spwn.quant <- apply(sub_pars$S, 2, quantile, probs=c(0.1,0.5,0.9))[,1:(nyrs-a_min)]
  rec.quant <- apply(sub_pars$R, 2, quantile, probs=c(0.1,0.5,0.9))[,(A+a_min):nRyrs]
  
  brood_t <- as.data.frame(cbind(sub_dat$year[1:(nyrs-A-1)], t(spwn.quant), t(rec.quant))) |>
    round(2) ##Why doesnt this line up for skeena? added -1 as hacky fix
  colnames(brood_t) <- c("BroodYear","S_lwr","S_med","S_upr","R_lwr","R_med","R_upr")
  
  brood_t <- mutate(brood_t, CU = i)
  
  brood.all <- rbind(brood.all, brood_t)
  
  spwn.quant.long <- apply(sub_pars$S, 2, quantile, probs=c(0.1,0.5,0.9))[,1:(nyrs)]
  rec.quant.long <- apply(sub_pars$R, 2, quantile, probs=c(0.1,0.5,0.9))[,(A+a_min):nRyrs]
  
  brood_t.long <- as.data.frame(cbind(sub_dat$year[1:nyrs],t(spwn.quant.long), rbind(t(rec.quant.long),matrix(NA,4,3)))) |>
    round(2)
  colnames(brood_t.long) <- c("BroodYear","S_lwr","S_med","S_upr","R_lwr","R_med","R_upr")
  
  brood_t.long <- mutate(brood_t.long, CU = i)
  
  brood.all.long <- rbind(brood.all.long, brood_t.long)
  
  #SR relationship based on full posterior---
  spw <- seq(0,max(brood_t$S_upr),length.out=100)
  SR.pred <- matrix(NA,length(spw), length(sub_pars$lnalpha))
  bench <- matrix(NA,length(sub_pars$lnalpha),8,
                  dimnames = list(seq(1:length(sub_pars$lnalpha)), c("Sgen", "Smsy", "Umsy", "Seq", "Smsr", "S.recent","Smsr.20","Smsr.40")))
  
  par <- matrix(NA,length(sub_pars$lnalpha),3,
                dimnames = list(seq(1:length(sub_pars$lnalpha)), c("sample","ln_a","beta")))
  
  # get benchmarks & pars (AR1 model)------------------------------------------------------------------
  
  for(j in 1:length(sub_pars$lnalpha)){
    ln_a <- sub_pars$lnalpha[j]
    b <- sub_pars$beta[j]
    SR.pred[,j] <- (exp(ln_a)*spw*exp(-b*spw))
    
    bench[j,2] <- get_Smsy(ln_a, b) #S_MSY
    bench[j,1] <- get_Sgen(exp(ln_a),b,-1,1/b*2, bench[j,2]) #S_gen
    bench[j,3] <- (1 - lambert_W0(exp(1 - ln_a))) #U_MSY
    bench[j,4] <- ln_a/b #S_eq
    bench[j,5] <- 1/b #S_msr
    bench[j,6] <- exp(mean(log(sub_pars$S[j, (nyrs-5):nyrs]))) #S recent - mean spawners in last generation
    bench[j,7] <- (1/b)*0.2
    bench[j,8] <- (1/b)*0.4
    
    par[j,1] <- j
    par[j,2] <- ln_a
    par[j,3] <- b
  }
  
  SR.pred <- as.data.frame(cbind(spw,t(apply(SR.pred, 1, quantile,probs=c(0.1,0.5,0.9), na.rm=T))))|>
    round(2) |>
    mutate(CU = i)
  
  SR.preds <- rbind(SR.preds, SR.pred)
  
  bench.posts <- rbind(bench.posts, as.data.frame(bench) |> mutate(CU = i))
  
  par.posts <- rbind(par.posts, as.data.frame(par) |> mutate(CU = i))
  
  bench.quant <- apply(bench[,1:8], 2, quantile, probs=c(0.1,0.5,0.9), na.rm=T) |>
    t()
  
  mean <- apply(bench[,1:8],2,mean, na.rm=T) #get means of each
  
  sub_benchmarks <- cbind(bench.quant, mean) |>
    as.data.frame() |>
    mutate(CU = i) |>
    relocate('50%', 1)
  
  #other pars to report
  alpha <- quantile(exp(sub_pars$lnalpha), probs = c(.1, .5, .9))
  beta <- quantile(sub_pars$beta, probs = c(.1, .5, .9))
  sigma <- quantile(sub_pars$sigma_R, probs = c(.1, .5, .9))
  phi <- quantile(sub_pars$phi, probs = c(.1, .5, .9))
  
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
  resid.quant <- apply(sub_pars$lnresid, 2, quantile, probs=c(0.1,0.25,0.5,0.75,0.9))[,(A):nRyrs]
  
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
                                    probs=c(0.025, 0.5, 0.975))$summary
  tva.par.summary.out <- rbind(tva.par.summary.out, as.data.frame(tva.par.summary) |> mutate(CU = i))
  a.yrs <- apply(exp(sub_pars_TVA$ln_alpha), 2, quantile, probs=c(0.1,0.5,0.9))
  a.yrs <- as.data.frame(cbind(sub_dat$year, t(a.yrs)))
  
  colnames(a.yrs) <- c("brood_year", "lwr", "mid", "upr")
  
  a.yrs.all <- rbind(a.yrs.all, data.frame(a.yrs, CU = i))
  
  #time varying alpha residuals
  resid.quant <- apply(sub_pars_TVA$lnresid, 2, quantile,
                       probs=c(0.1,0.25,0.5,0.75,0.9))[,(A):nRyrs]
  
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
  TV.spwn.quant <- data.frame(t(apply(sub_pars_TVA$S, 2, quantile, probs=c(0.25,0.5,0.75)))) |>
    mutate(CU =i,
           year = unique(sub_dat$year))
  
  TV.spwn <- bind_rows(TV.spwn, TV.spwn.quant)
  
  TV.harv.quant <- data.frame(t(apply(sub_pars_TVA$C, 2, quantile, probs=c(0.25,0.5,0.75)))) |>
    mutate(CU =i,
           year = unique(sub_dat$year))
  
  TV.harv <- bind_rows(TV.harv, TV.harv.quant)
}  # End data wrangling loop by CU


colnames(SR.preds) <- c("Spawn", "Rec_lwr","Rec_med","Rec_upr", "CU")
colnames(AR1.resids) <- c("year","lwr","midlwr","mid","midupr","upr", "CU")
colnames(AR1.spwn) <- c("S.25", "S.50", "S.75", "year", "CU")
colnames(AR1.harv) <- c("H.25", "H.50", "H.75", "year", "CU")
colnames(TV.resids) <- c("year","lwr","midlwr","mid","midupr","upr", "CU")
colnames(TV.spwn) <- c("S.25", "S.50", "S.75", "CU", "year")
colnames(TV.harv) <- c("H.25", "H.50", "H.75", "CU", "year")

if(FALSE){
SR.preds$CU_f <- factor(SR.preds$CU, levels = CU_order)
AR1.spwn$CU_f <- factor(AR1.spwn$CU, levels = CU_order)
AR1.resids$CU_f <- factor(AR1.resids$CU, levels = CU_order)
AR1.harv$CU_f <- factor(AR1.harv$CU, levels = CU_order)
TV.resids$CU_f <- factor(TV.resids$CU, levels = CU_order)
TV.spwn$CU_f <- factor(TV.spwn$CU, levels = CU_order)
TV.harv$CU_f <- factor(TV.harv$CU, levels = CU_order)
brood.all$CU_f <- factor(brood.all$CU, levels = CU_order)
brood.all.long$CU_f <- factor(brood.all.long$CU, levels = CU_order)
TV.SR.preds$CU_f <- factor(TV.SR.preds$CU, levels = CU_order)
esc$CU_f <- factor(esc$stock, levels = CU_order)
}#turn off

# write important tables to repo ---------------------------------------------------------
bench.par.table.out <- bench.par.table |> # AR1 SR model pars and benchmarks
  relocate(CU, 1) |>
  relocate(bench.par, .after = 1) |>
  relocate(mean, .after = 2) |>
  mutate_at(3:7, ~round(.,5)) |>
  arrange(bench.par, CU)

bench.par.table.out[73:80,3:6] <- bench.par.table.out[73:80,3:6]*10000 ## transform beta - do it better?

write.csv(bench.par.table.out, here("data/generated/bench_par_table.csv"),
          row.names = FALSE)

write_rds(bench.posts, here("data/generated/benchmark_posteriors.rds"))

write.csv(par.posts, here("data/generated/AR1_posteriors.csv"))

write.csv(par.posts.tva, here("data/generated/TVA_posteriors.csv"))

write.csv(tva.par.summary.out, here("data/generated/TVA_par_summary.csv"))

write.csv(brood.all.long, here("data/generated/brood_table_long.csv"),
          row.names = FALSE)

a.yrs.all$model<-"spw"
write.csv(a.yrs.all, here("data/generated/spw_TVA.csv"),
          row.names = FALSE)


# make key plots for pub -----------------------------------------------------------------

# SR fits ----
#Skeena
ggplot() +
  geom_abline(intercept = 0, slope = 1,col="dark grey") +
  geom_ribbon(data = filter(SR.preds, !(CU %in% c("Lower Nass", "Upper Nass"))), 
              aes(x = Spawn/1000, ymin = Rec_lwr/1000, ymax = Rec_upr/1000),
              fill = "grey80", alpha=0.5, linetype=2, colour="gray46") +
  geom_errorbar(data = filter(brood.all, !(CU %in% c("Lower Nass", "Upper Nass"))), 
                aes(x= S_med/1000, y = R_med/1000, ymin = R_lwr/1000, ymax = R_upr/1000),
                colour="grey", width=0, linewidth=0.3) +
  geom_errorbarh(data = filter(brood.all, !(CU %in% c("Lower Nass", "Upper Nass"))), 
                 aes(y = R_med/1000, xmin = S_lwr/1000, xmax = S_upr/1000),
                 height=0, colour = "grey", linewidth = 0.3) +
  geom_point(data = filter(brood.all, !(CU %in% c("Lower Nass", "Upper Nass"))),
             aes(x = S_med/1000, y = R_med/1000, color=BroodYear), size = 1.5) +
  geom_line(data = filter(SR.preds, !(CU %in% c("Lower Nass", "Upper Nass"))), 
            aes(x = Spawn/1000, y = Rec_med/1000)) +
  facet_wrap(~CU, scales = "free") +
  scale_colour_viridis_c(name = "Brood Year")+
  labs(x = "Spawners (000s)", y = "Recruits (000s)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

my.ggsave(here("plots/Skeena/SR_fits_AR1.PNG"))

#Nass
ggplot() +
  geom_abline(intercept = 0, slope = 1,col="dark grey") +
  geom_ribbon(data = filter(SR.preds, CU %in% c("Lower Nass", "Upper Nass")), 
              aes(x = Spawn/1000, ymin = Rec_lwr/1000, ymax = Rec_upr/1000),
              fill = "grey80", alpha=0.5, linetype=2, colour="gray46") +
  geom_errorbar(data = filter(brood.all, CU %in% c("Lower Nass", "Upper Nass")), 
                aes(x= S_med/1000, y = R_med/1000, ymin = R_lwr/1000, ymax = R_upr/1000),
                colour="grey", width=0, linewidth=0.3) +
  geom_errorbarh(data = filter(brood.all, CU %in% c("Lower Nass", "Upper Nass")), 
                 aes(y = R_med/1000, xmin = S_lwr/1000, xmax = S_upr/1000),
                 height=0, colour = "grey", linewidth = 0.3) +
  geom_point(data = filter(brood.all, CU %in% c("Lower Nass", "Upper Nass")),
             aes(x = S_med/1000, y = R_med/1000, color=BroodYear), size = 1.5) +
  geom_line(data = filter(SR.preds, CU %in% c("Lower Nass", "Upper Nass")), 
            aes(x = Spawn/1000, y = Rec_med/1000)) +
  facet_wrap(~CU, scales = "free") +
  scale_colour_viridis_c(name = "Brood Year")+
  labs(x = "Spawners (000s)", y = "Recruits (000s)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

my.ggsave(here("plots/Nass/SR_fits_AR1.PNG"))

# AR1 resids ----
ggplot(AR1.resids, aes(x=year, y = mid)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),  fill = "darkgrey", alpha = 0.5) +
  geom_ribbon(aes(ymin = midlwr, ymax = midupr),  fill = "black", alpha=0.2) +
  geom_line(lwd = 1.1) +
  coord_cartesian(ylim=c(-2,2)) +
  labs(x = "Return year",
       y = "Recruitment residuals",
       title = "AR1 recruitment residuals") +
  facet_wrap(~CU) +
  theme(legend.position = "none",
        panel.grid = element_blank()) +
  geom_abline(intercept = 0, slope = 0, col = "dark grey", lty = 2)

my.ggsave(here("plots/AR1_resids.PNG"))

# TV resids ----
# Skeena
ggplot(TV.resids, aes(x=year, y = mid)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),  fill = "darkgrey", alpha = 0.5) +
  geom_ribbon(aes(ymin = midlwr, ymax = midupr),  fill = "black", alpha=0.2) +
  geom_line(lwd = 1.1) +
  labs(x = "Return year",
       y = "Recruitment residuals",
       title = "Time-varying recruitment residuals") +
  facet_wrap(~CU) +
  theme(legend.position = "none",
        panel.grid = element_blank()) +
  geom_abline(intercept = 0, slope = 0, col = "dark grey", lty = 2)

my.ggsave(here("plots/TV_resids.PNG"))

# TV alpha ----
#Skeena
a.yrs.all |>
  filter(brood_year < 2018) |> ## why?
  filter(!(CU %in% c("Lower Nass", "Upper Nass"))) |>
  ggplot(aes(color=CU)) +
  geom_line(aes(x = brood_year , y = mid), lwd = 1.5) +
  scale_color_viridis_d() +
  geom_hline(yintercept = 1, lty=2, col = "grey") +
  labs(y ="Productivity (\U03B1)", x = "Brood year")+
  guides(color=guide_legend(title="Conservation Unit"))

my.ggsave(here("plots/Skeena/changing_productivity.PNG"))

#Nass
a.yrs.all |>
  filter(brood_year < 2018) |> ## why?
  filter(CU %in% c("Lower Nass", "Upper Nass")) |>
  ggplot(aes(color=CU)) +
  geom_line(aes(x = brood_year , y = mid), lwd = 1.5) +
  scale_color_viridis_d() +
  geom_hline(yintercept = 1, lty=2, col = "grey") +
  labs(y ="Productivity (\U03B1)", x = "Brood year")+
  guides(color=guide_legend(title="Conservation Unit"))

my.ggsave(here("plots/Nass/changing_productivity.PNG"))

# TV SR fits ----
ggplot() +
  geom_point(data = brood.all,
             aes(x = S_med/1000,
                 y = R_med/1000),
             size = 1.5) +
  geom_line(data = TV.SR.preds, aes(x = spw/1000, y = pred.R/1000, color = year, group = year)) +
  facet_wrap(~CU, scales = "free") +
  scale_colour_viridis_c(name = "Brood Year")+
  geom_abline(intercept = 0, slope = 1,col="dark grey") +
  labs(x = "Spawners (000s)",
       y = "Recruits (000s)",
       title = "Time varying productivity spawner-recruit fits") +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size=9),
        legend.text = element_text(size=8))

tmy.ggsave(here("plots/TV_SR_fits.PNG"))

# deleted yukon code below this fig. To see more look here 
# (https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/R/inference_figs.R#L383)


#post hoc plots for skeena

#latent states ---
AR1.obs <- AR1.harv |>
  mutate(obs = "Harvest") |>
  rename(S.25 = H.25, 
         S.50 = H.50, 
         S.75 = H.75) |>
  bind_rows(AR1.spwn) |> 
  mutate(obs = ifelse(is.na(obs), "Spawners", obs))
colnames(AR1.obs) <- c("p.25", "p.50", "p.75", "year", "CU", "obs")

#Skeena
ggplot(filter(AR1.obs, !(CU %in% c("Lower Nass", "Upper Nass")))) +
  geom_ribbon(aes(x = year, ymin = p.25/1000, ymax = p.75/1000, fill = CU), alpha = 0.2) +
  geom_line(aes(x = year, y = p.50/1000, color = CU)) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  facet_wrap(~obs, nrow = 2, scales = "free") +
  labs(x = "Brood year", y = "Fish (thousands)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size=9),
        legend.text = element_text(size=8))
  
my.ggsave(here("plots/Skeena/sp_har_obs.PNG"))

#Nass
ggplot(filter(AR1.obs, CU %in% c("Lower Nass", "Upper Nass"))) +
  geom_line(aes(x = year, y = p.50/1000, color = CU), lwd=1) +
  geom_ribbon(aes(x = year, ymin = p.25/1000, ymax = p.75/1000, fill = CU), alpha = 0.2) +
  scale_color_viridis_d(end = 0.6) +
  scale_fill_viridis_d(end = 0.6) +
  facet_wrap(~obs, nrow = 2, scales = "free") +
  labs(x = "Brood year", y = "Fish (thousands)") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

my.ggsave(here("plots/Nass/sp_har_obs.PNG"))

#plot A_obs 
A_plot <- sp_har |>
  select(CU, year, a4, a5, a6) |>
  pivot_longer(cols = a4:a6)

ggplot(filter(A_plot, !(CU %in% c("Lower Nass", "Upper Nass"))), 
       aes(x=year, y=value/100)) +
  geom_line() +
  facet_grid(CU~name) +
  scale_color_viridis_d() +
  labs(x = "Return year", y = "Proportion at age") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#better option
#Skeena
ggplot(filter(A_plot, !(CU %in% c("Lower Nass", "Upper Nass"))), 
       aes(x=year, y=value/100, fill = name)) +
  geom_area(alpha = 0.8) +  
  theme_minimal() +
  facet_wrap(~CU) +
  scale_fill_viridis_d(name = "Age class") +
  labs(x = "Return year", y = "Proportion at age")

my.ggsave(here("plots/Skeena/A_obs.PNG"))

#Nass
ggplot(filter(A_plot, CU %in% c("Lower Nass", "Upper Nass")), 
       aes(x=year, y=value/100, fill = name)) +
  geom_area(alpha = 0.8) +  
  theme_minimal() +
  facet_wrap(~CU) +
  scale_fill_viridis_d(name = "Age class") +
  labs(x = "Return year", y = "Proportion at age")

my.ggsave(here("plots/Nass/A_obs.PNG"))

