# script to generate yield and recruitment profiles for Kitsumkalum, and Skeena and Nass aggregates ----
library(here)
library(tidyverse)
library(gsl)
library(viridis)
library(ggpubr)

source(here("R/data_functions.R"))

# read in data ---------------------------------------------------------------------------
# model fits ---
AR1.fits <- lapply(list.files(here("data/generated/model_fits/AR1"),
                              full.names = T), readRDS)
names(AR1.fits) <- unique(sp_har$CU)

# Kitsumkalum profiles ----
Kitsumkalum <- sp_har |> filter(CU == "Kitsumkalum")

pars <- rstan::extract(AR1.fits[["Kitsumkalum"]])
max_samples <- dim(pars$lnalpha)
spw <- seq(0,25000,length.out=100)

Profiles <- array(NA, dim = c(length(spw),max_samples,9))

for(i in 1:max_samples){
  r <- sample(seq(1,max_samples),1,replace=T)
  a <- pars$lnalpha[r]
  b <- pars$beta[r]
  S_msy_c <- get_Smsy(a, b)
  R_msy_c <- S_msy_c*exp(a-b*S_msy_c)
  msy_c <- R_msy_c - S_msy_c
  U_max_c <- 1-(1/exp(a))
  Umsy_Sch_corr <- (1 - lambert_W0(exp(1 - a)))
  
  for(j in 1:length(spw)){
    spw_star <- spw[j]
    R_star <- spw_star*exp(a-b*spw_star)
    R_star[is.nan(R_star)] <- 0
    sus_yield <- R_star - spw_star
    U <- ifelse(sus_yield/R_star < 0, 1, sus_yield/R_star)
    Profiles[j,i,1] <- ifelse(sus_yield > 0.9*msy_c, 1, 0) # probability yield > 90% of MSY
    Profiles[j,i,2] <- ifelse(sus_yield > 0.8*msy_c, 1, 0) # probability yield > 80% of MSY
    Profiles[j,i,3] <- ifelse(sus_yield > 0.7*msy_c, 1, 0) # probability yield > 70% of MSY
    Profiles[j,i,4] <- ifelse(R_star > 0.9*R_msy_c, 1, 0) # probability recruitment > 90% of MSR
    Profiles[j,i,5] <- ifelse(R_star > 0.8*R_msy_c, 1, 0) # probability recruitment > 80% of MSR
    Profiles[j,i,6] <- ifelse(R_star > 0.7*R_msy_c, 1, 0) # probability recruitment > 70% of MSR
    if(sus_yield < (0.9*msy_c) & spw_star < S_msy_c){
      Profiles[j,i,7] <- 1
    } else {
      Profiles[j,i,7] <- 0
    }
    if(sus_yield < (0.8*msy_c) & spw_star < S_msy_c){
      Profiles[j,i,8] <- 1
    } else {
      Profiles[j,i,8] <- 0
    }
    if(sus_yield < (0.7*msy_c) & spw_star < S_msy_c){
      Profiles[j,i,9] <- 1
    } else {
      Profiles[j,i,9] <- 0
    }    
  }
}

profiles_median <- cbind(spw,apply(Profiles,c(1,3),mean, na.rm=T))

colnames(profiles_median) <- c("Spawners","90% MSY","80% MSY","70% MSY", "90% MSR","80% MSR","70% MSR", "90% MSY ", "80% MSY ", "70% MSY ")
profiles_median <- as.data.frame(profiles_median)

long_yield_profiles <- profiles_median%>%
  select("Spawners","90% MSY","80% MSY","70% MSY")%>%
  pivot_longer(!Spawners, names_to = "Metric", values_to = "perc")

a <- ggplot(long_yield_profiles, aes(x = Spawners, y = perc, group = Metric, color = Metric)) +
  geom_line(size = 1.5) +
  scale_color_viridis(discrete = TRUE) +
  ylab("Probability") +
  theme(legend.position = "none") +
  theme_sleek() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank()) + 
  geom_rug(data = Kitsumkalum,
           aes(x = spwn),
           inherit.aes = FALSE,
           sides="b")

long_rec_profiles <- profiles_median%>%
  select("Spawners","90% MSR","80% MSR","70% MSR")%>%
  pivot_longer(!Spawners, names_to = "Metric", values_to = "perc")

b <- ggplot(long_rec_profiles, aes(x = Spawners, y = perc, group = Metric, color = Metric)) +
  geom_line(size = 1.5) +
  scale_color_viridis(discrete = TRUE) +
  ylab("Probability") +
  theme(legend.position = "none") +
  theme_sleek() +
  theme(legend.title = element_blank()) + 
  geom_rug(data = Kitsumkalum,
           aes(x = spwn),
           inherit.aes = FALSE,
           sides="b")

g <- ggarrange(a,b,nrow =2, labels = c("a", "b"), heights = c(0.75,0.9))
print(g)
my.ggsave(here("plots/Skeena/Kitsum-profiles.PNG"))


# Nass aggregate profiles ----
Nass_agg <- sp_har |> filter(CU == "Nass Aggregate")

pars <- rstan::extract(AR1.fits[["Nass Aggregate"]])
max_samples <- dim(pars$lnalpha)
spw <- seq(0,55000,length.out=100)

Profiles <- array(NA, dim = c(length(spw),max_samples,9))

for(i in 1:max_samples){
  r <- sample(seq(1,max_samples),1,replace=T)
  a <- pars$lnalpha[r]
  b <- pars$beta[r]
  S_msy_c <- get_Smsy(a, b)
  R_msy_c <- S_msy_c*exp(a-b*S_msy_c)
  msy_c <- R_msy_c - S_msy_c
  U_max_c <- 1-(1/exp(a))
  Umsy_Sch_corr <- (1 - lambert_W0(exp(1 - a)))
  
  for(j in 1:length(spw)){
    spw_star <- spw[j]
    R_star <- spw_star*exp(a-b*spw_star)
    R_star[is.nan(R_star)] <- 0
    sus_yield <- R_star - spw_star
    U <- ifelse(sus_yield/R_star < 0, 1, sus_yield/R_star)
    Profiles[j,i,1] <- ifelse(sus_yield > 0.9*msy_c, 1, 0) # probability yield > 90% of MSY
    Profiles[j,i,2] <- ifelse(sus_yield > 0.8*msy_c, 1, 0) # probability yield > 80% of MSY
    Profiles[j,i,3] <- ifelse(sus_yield > 0.7*msy_c, 1, 0) # probability yield > 70% of MSY
    Profiles[j,i,4] <- ifelse(R_star > 0.9*R_msy_c, 1, 0) # probability recruitment > 90% of MSR
    Profiles[j,i,5] <- ifelse(R_star > 0.8*R_msy_c, 1, 0) # probability recruitment > 80% of MSR
    Profiles[j,i,6] <- ifelse(R_star > 0.7*R_msy_c, 1, 0) # probability recruitment > 70% of MSR
    if(sus_yield < (0.9*msy_c) & spw_star < S_msy_c){
      Profiles[j,i,7] <- 1
    } else {
      Profiles[j,i,7] <- 0
    }
    if(sus_yield < (0.8*msy_c) & spw_star < S_msy_c){
      Profiles[j,i,8] <- 1
    } else {
      Profiles[j,i,8] <- 0
    }
    if(sus_yield < (0.7*msy_c) & spw_star < S_msy_c){
      Profiles[j,i,9] <- 1
    } else {
      Profiles[j,i,9] <- 0
    }    
  }
}

profiles_median <- cbind(spw,apply(Profiles,c(1,3),mean, na.rm=T))

colnames(profiles_median) <- c("Spawners","90% MSY","80% MSY","70% MSY", "90% MSR","80% MSR","70% MSR", "90% MSY ", "80% MSY ", "70% MSY ")
profiles_median <- as.data.frame(profiles_median)

long_yield_profiles <- profiles_median%>%
  select("Spawners","90% MSY","80% MSY","70% MSY")%>%
  pivot_longer(!Spawners, names_to = "Metric", values_to = "perc")

a <- ggplot(long_yield_profiles, aes(x = Spawners, y = perc, group = Metric, color = Metric)) +
  geom_line(size = 1.5) +
  scale_color_viridis(discrete = TRUE) +
  ylab("Probability") +
  theme(legend.position = "none") +
  theme_sleek() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank()) + 
  geom_rug(data = Nass_agg,
           aes(x = spwn),
           inherit.aes = FALSE,
           sides="b")

long_rec_profiles <- profiles_median%>%
  select("Spawners","90% MSR","80% MSR","70% MSR")%>%
  pivot_longer(!Spawners, names_to = "Metric", values_to = "perc")

b <- ggplot(long_rec_profiles, aes(x = Spawners, y = perc, group = Metric, color = Metric)) +
  geom_line(size = 1.5) +
  scale_color_viridis(discrete = TRUE) +
  ylab("Probability") +
  theme(legend.position = "none") +
  theme_sleek() +
  theme(legend.title = element_blank()) + 
  geom_rug(data = Nass_agg,
           aes(x = spwn),
           inherit.aes = FALSE,
           sides="b")

g <- ggarrange(a,b,nrow =2, labels = c("a", "b"), heights = c(0.75,0.9))
print(g)
my.ggsave(here("plots/Nass/Nass-profiles.PNG"))

# Skeena aggregate profiles ----
Skeena_agg <- sp_har |> filter(CU == "Skeena Aggregate")

pars <- rstan::extract(AR1.fits[["Skeena Aggregate"]])
max_samples <- dim(pars$lnalpha)
spw <- seq(0,200000,length.out=100)

Profiles <- array(NA, dim = c(length(spw),max_samples,9))

for(i in 1:max_samples){
  r <- sample(seq(1,max_samples),1,replace=T)
  a <- pars$lnalpha[r]
  b <- pars$beta[r]
  S_msy_c <- get_Smsy(a, b)
  R_msy_c <- S_msy_c*exp(a-b*S_msy_c)
  msy_c <- R_msy_c - S_msy_c
  U_max_c <- 1-(1/exp(a))
  Umsy_Sch_corr <- (1 - lambert_W0(exp(1 - a)))
  
  for(j in 1:length(spw)){
    spw_star <- spw[j]
    R_star <- spw_star*exp(a-b*spw_star)
    R_star[is.nan(R_star)] <- 0
    sus_yield <- R_star - spw_star
    U <- ifelse(sus_yield/R_star < 0, 1, sus_yield/R_star)
    Profiles[j,i,1] <- ifelse(sus_yield > 0.9*msy_c, 1, 0) # probability yield > 90% of MSY
    Profiles[j,i,2] <- ifelse(sus_yield > 0.8*msy_c, 1, 0) # probability yield > 80% of MSY
    Profiles[j,i,3] <- ifelse(sus_yield > 0.7*msy_c, 1, 0) # probability yield > 70% of MSY
    Profiles[j,i,4] <- ifelse(R_star > 0.9*R_msy_c, 1, 0) # probability recruitment > 90% of MSR
    Profiles[j,i,5] <- ifelse(R_star > 0.8*R_msy_c, 1, 0) # probability recruitment > 80% of MSR
    Profiles[j,i,6] <- ifelse(R_star > 0.7*R_msy_c, 1, 0) # probability recruitment > 70% of MSR
    if(sus_yield < (0.9*msy_c) & spw_star < S_msy_c){
      Profiles[j,i,7] <- 1
    } else {
      Profiles[j,i,7] <- 0
    }
    if(sus_yield < (0.8*msy_c) & spw_star < S_msy_c){
      Profiles[j,i,8] <- 1
    } else {
      Profiles[j,i,8] <- 0
    }
    if(sus_yield < (0.7*msy_c) & spw_star < S_msy_c){
      Profiles[j,i,9] <- 1
    } else {
      Profiles[j,i,9] <- 0
    }    
  }
}

profiles_median <- cbind(spw,apply(Profiles,c(1,3),mean, na.rm=T))

colnames(profiles_median) <- c("Spawners","90% MSY","80% MSY","70% MSY", "90% MSR","80% MSR","70% MSR", "90% MSY ", "80% MSY ", "70% MSY ")
profiles_median <- as.data.frame(profiles_median)

long_yield_profiles <- profiles_median%>%
  select("Spawners","90% MSY","80% MSY","70% MSY")%>%
  pivot_longer(!Spawners, names_to = "Metric", values_to = "perc")

a <- ggplot(long_yield_profiles, aes(x = Spawners, y = perc, group = Metric, color = Metric)) +
  geom_line(size = 1.5) +
  scale_color_viridis(discrete = TRUE) +
  ylab("Probability") +
  theme(legend.position = "none") +
  theme_sleek() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank()) + 
  geom_rug(data = Skeena_agg,
           aes(x = spwn),
           inherit.aes = FALSE,
           sides="b")

long_rec_profiles <- profiles_median%>%
  select("Spawners","90% MSR","80% MSR","70% MSR")%>%
  pivot_longer(!Spawners, names_to = "Metric", values_to = "perc")

b <- ggplot(long_rec_profiles, aes(x = Spawners, y = perc, group = Metric, color = Metric)) +
  geom_line(size = 1.5) +
  scale_color_viridis(discrete = TRUE) +
  ylab("Probability") +
  theme(legend.position = "none") +
  theme_sleek() +
  theme(legend.title = element_blank()) + 
  geom_rug(data = Skeena_agg,
           aes(x = spwn),
           inherit.aes = FALSE,
           sides="b")

g <- ggarrange(a,b,nrow =2, labels = c("a", "b"), heights = c(0.75,0.9))
print(g)
my.ggsave(here("plots/Skeena/Skeena-profiles.PNG"))
