# script to calculate reference points for demographic models
library(here)
library(tidyverse)
library(rstan)
library(ggsidekick)

source(here("R/data_functions.R"))

# Skeena aggregate ----
male_S_comps <- read.csv(here("data/Agg_male_age_comps.csv"))[,3:5] |> # male spawner age composition
  mutate(X6 = replace_na(X6, 0))

fem_S_comps <- read.csv(here("data/Agg_fem_age_comps.csv"))[,3:5] # female spawner age composition
fem_egg_mass_comps <- read.csv(here("data/Agg_fem_age_egg_mass.csv"))[,2:4] # female spawner egg mass by age

AR1.EM.fit <- readRDS("data/generated/model_fits/EM/Skeena Aggregate_EM.rds")

bench_AR1_eggs <- NULL
ER.pred <- NULL
brood.all <- NULL
a.yrs.all <- NULL
demo.par.summary.out <- NULL

sub_dat <- filter(sp_har, CU=="Skeena Aggregate")
a_min <- 4
a_max <- 6
nyrs <- max(sub_dat$year)-min(sub_dat$year)+1 #number of calendar years of observations
A <- a_max - a_min + 1 #total age classes
nRyrs <- nyrs + A - 1 #total number of fully and partially observed recruitments from observed spawners

#load AR1 model with egg mass parameters ---
EM_pars <- rstan::extract(AR1.EM.fit)

#eggmass-recruit relationship
egg.quant <- apply(exp(EM_pars$lnEM), 2, quantile, probs=c(0.1,0.5,0.9))[,1:(nyrs-a_min)]
rec.quant <- apply(EM_pars$R, 2, quantile, probs=c(0.1,0.5,0.9))[,(A+a_min):nRyrs]

brood_t <- as.data.frame(cbind(sub_dat$year[1:(nyrs-A)],t(egg.quant), t(rec.quant))) |>
  round(2)

colnames(brood_t) <- c("BroodYear","EM_lwr","EM_med","EM_upr","R_lwr","R_med","R_upr")


#SR relationship based on full posterior---
em <- seq(0,max(brood_t$EM_upr),length.out=100)
ER.pred <- matrix(NA,length(em), length(EM_pars$lnalpha))

for(j in 1:length(EM_pars$lnalpha)){
  ln_a <- EM_pars$lnalpha[j]
  b <- EM_pars$beta[j]
  ER.pred[,j] <- (exp(ln_a)*em*exp(-b*em))
}

ER.pred <- as.data.frame(cbind(em,t(apply(ER.pred, 1, quantile, probs=c(0.1,0.5,0.9), na.rm=T))))|>
  round(2) 

colnames(ER.pred) <- c("em","Rec_lwr","Rec_med","Rec_upr")

# plot EM recruitment relationship 

ggplot() +
  geom_ribbon(data = ER.pred, 
              aes(x = em/100000, ymin = Rec_lwr/1000, ymax = Rec_upr/1000),
              fill = "grey80", alpha=0.5, linetype=2, colour="gray46") +
  geom_errorbar(data = brood_t, 
                aes(x= EM_med/100000, y = R_med/1000, ymin = R_lwr/1000, ymax = R_upr/1000),
                colour="grey", width=0, linewidth=0.3) +
  geom_errorbar(data = brood_t, 
                aes(y = R_med/1000, xmin = EM_lwr/100000, xmax = EM_upr/100000),
                height=0, colour = "grey", linewidth = 0.3) +
  geom_point(data = brood_t,
             aes(x = EM_med/100000, y = R_med/1000, color=BroodYear), size = 1.5) +
  geom_line(data = ER.pred, 
            aes(x = em/100000, y = Rec_med/1000)) +
  scale_colour_viridis_c(name = "Brood Year")+
  labs(x = "Egg mass (00 kgs)", y = "Recruits (000s)") +
  theme_sleek()

ggsave("plots/skeena/Agg-EMR.jpeg",width=7,height=5)

#benchmarks based on early demographics
bench_early <- matrix(NA,length(EM_pars$lnalpha),2,
                      dimnames = list(seq(1:length(EM_pars$lnalpha)), c("Smsy","Smsr")))

for(j in 1:length(EM_pars$lnalpha)){
  ln_a <- EM_pars$lnalpha[j]
  beta <- EM_pars$beta[j]
  year_index = c(1:6) # first 6 brood years
  eta_as_F <- colMeans(fem_S_comps[year_index,]) # probability of returning by age for females
  eta_as_M <- colMeans(male_S_comps[year_index,]) # probability of returning by age for males
  eta_as <- c(eta_as_F,eta_as_M) # probability of returning by age and sex
  z_as <- c(colMeans(fem_egg_mass_comps[year_index,]),c(0,0,0)) #  age/sex class-specific reproductive output
  opt <- optimize(calcCeq, c(0,5), maximum=TRUE)
  bench_early[j,1] <- sum(Seq_as)
  opt <- optimize(calcNeq, c(0,5), maximum=TRUE)
  bench_early[j,2] <- sum(Seq_as)
}
bench_early <- as.data.frame(bench_early)
bench_early$Smsy[bench_early$Smsy <0] <- NA
bench_early$Smsr[bench_early$Smsr <0] <- NA
bench_early$period <- "early"
bench_early$CU <- "Skeena Aggregate"

#benchamrks based on long-term average demographics
bench_avg <- matrix(NA,length(EM_pars$lnalpha),2,
                    dimnames = list(seq(1:length(EM_pars$lnalpha)), c("Smsy","Smsr")))

for(j in 1:length(EM_pars$lnalpha)){
  ln_a <- EM_pars$lnalpha[j]
  beta <- EM_pars$beta[j]
  year_index = c(1:nyrs)
  eta_as_F <- colMeans(fem_S_comps[year_index,]) # probability of returning by age for females
  eta_as_M <- colMeans(male_S_comps[year_index,]) # probability of returning by age for males
  eta_as <- c(eta_as_F,eta_as_M) # probability of returning by age and sex
  z_as <- c(colMeans(fem_egg_mass_comps[year_index,]),c(0,0,0)) #  age/sex class-specific reproductive output
  opt <- optimize(calcCeq, c(0,5), maximum=TRUE)
  bench_avg[j,1] <- sum(Seq_as)
  opt <- optimize(calcNeq, c(0,5), maximum=TRUE)
  bench_avg[j,2] <- sum(Seq_as)
}
bench_avg <- as.data.frame(bench_avg)
bench_avg$Smsy[bench_avg$Smsy <0] <- NA
bench_avg$Smsr[bench_avg$Smsr <0] <- NA
bench_avg$period <- "avg"
bench_avg$CU <- "Skeena Aggregate"

#benchmarks based on recent demographics
bench_recent <- matrix(NA,length(EM_pars$lnalpha),2,
                       dimnames = list(seq(1:length(EM_pars$lnalpha)), c("Smsy","Smsr")))

for(j in 1:length(EM_pars$lnalpha)){
  ln_a <- EM_pars$lnalpha[j]
  beta <- EM_pars$beta[j]
  year_index = c((nyrs-7):nyrs) # last 6 brood years
  eta_as_F <- colMeans(fem_S_comps[year_index,]) # probability of returning by age for females
  eta_as_M <- colMeans(male_S_comps[year_index,]) # probability of returning by age for males
  eta_as <- c(eta_as_F,eta_as_M) # probability of returning by age and sex
  z_as <- c(colMeans(fem_egg_mass_comps[year_index,]),c(0,0,0)) #  age/sex class-specific reproductive output
  opt <- optimize(calcCeq, c(0,5), maximum=TRUE)
  bench_recent[j,1] <- sum(Seq_as)
  opt <- optimize(calcNeq, c(0,5), maximum=TRUE)
  bench_recent[j,2] <- sum(Seq_as)
}
bench_recent <- as.data.frame(bench_recent)
bench_recent$Smsy[bench_recent$Smsy <0] <- NA
bench_recent$Smsr[bench_recent$Smsr <0] <- NA
bench_recent$period <- "recent"
bench_recent$CU <- "Skeena Aggregate"

bench_AR1_eggmass <- rbind(bench_early, bench_avg, bench_recent)

colnames(ER.pred) <- c("EM", "Rec_lwr","Rec_med","Rec_upr")

bench_egg_mass_long <- pivot_longer(bench_AR1_eggmass, cols = c(Smsy, Smsr), names_to = "par")

bench_egg_mass_long_summary <- bench_egg_mass_long |>
  group_by(period, par) |>
  summarize(median = quantile(value,probs=c(0.5),na.rm=T),
            lower = quantile(value,probs=c(0.025),na.rm=T),
            upper = quantile(value,probs=c(0.975),na.rm=T))

ggplot(bench_egg_mass_long |> filter(par=="Smsy", value <150000) ,aes(value, fill = period, color = period)) +
  geom_density(alpha = 0.8, adjust = 5) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme_sleek() +
  coord_cartesian(xlim = c(NA,100000)) +
  xlab("Smsy") 

ggsave("plots/skeena/Aggregate-demographic-Smsy-post.jpeg",width=7,height=5)


# Kitsumkalum ----
male_S_comps <- read.csv(here("data/Kitsum_male_age_comps.csv"))[,3:5] |> # male spawner age composition
  mutate(X6 = replace_na(X6, 0))

fem_S_comps <- read.csv(here("data/Kitsum_fem_age_comps.csv"))[,3:5] # female spawner age composition
fem_egg_mass_comps <- read.csv(here("data/Kitsum_fem_age_egg_mass.csv"))[,2:4] # female spawner egg mass by age

AR1.EM.fit <- readRDS("data/generated/model_fits/EM/Kitsumkalum_EM.rds")

bench_AR1_eggs <- NULL
ER.pred <- NULL
brood.all <- NULL
a.yrs.all <- NULL
demo.par.summary.out <- NULL

sub_dat <- filter(sp_har, CU=="Kitsumkalum")
a_min <- 4
a_max <- 6
nyrs <- max(sub_dat$year)-min(sub_dat$year)+1 #number of calendar years of observations
A <- a_max - a_min + 1 #total age classes
nRyrs <- nyrs + A - 1 #total number of fully and partially observed recruitments from observed spawners

#load AR1 model with egg mass parameters ---
EM_pars <- rstan::extract(AR1.EM.fit)

#eggmass-recruit relationship
egg.quant <- apply(exp(EM_pars$lnEM), 2, quantile, probs=c(0.1,0.5,0.9))[,1:(nyrs-a_min)]
rec.quant <- apply(EM_pars$R, 2, quantile, probs=c(0.1,0.5,0.9))[,(A+a_min):nRyrs]

brood_t <- as.data.frame(cbind(sub_dat$year[1:(nyrs-A)],t(egg.quant), t(rec.quant))) |>
  round(2)

colnames(brood_t) <- c("BroodYear","EM_lwr","EM_med","EM_upr","R_lwr","R_med","R_upr")


#SR relationship based on full posterior---
em <- seq(0,max(brood_t$EM_upr),length.out=100)
ER.pred <- matrix(NA,length(em), length(EM_pars$lnalpha))

for(j in 1:length(EM_pars$lnalpha)){
  ln_a <- EM_pars$lnalpha[j]
  b <- EM_pars$beta[j]
  ER.pred[,j] <- (exp(ln_a)*em*exp(-b*em))
}

ER.pred <- as.data.frame(cbind(em,t(apply(ER.pred, 1, quantile, probs=c(0.1,0.5,0.9), na.rm=T))))|>
  round(2) 

colnames(ER.pred) <- c("em","Rec_lwr","Rec_med","Rec_upr")

# plot EM recruitment relationship 

ggplot() +
  geom_ribbon(data = ER.pred, 
              aes(x = em/100000, ymin = Rec_lwr/1000, ymax = Rec_upr/1000),
              fill = "grey80", alpha=0.5, linetype=2, colour="gray46") +
  geom_errorbar(data = brood_t, 
                aes(x= EM_med/100000, y = R_med/1000, ymin = R_lwr/1000, ymax = R_upr/1000),
                colour="grey", width=0, linewidth=0.3) +
  geom_errorbar(data = brood_t, 
                aes(y = R_med/1000, xmin = EM_lwr/100000, xmax = EM_upr/100000),
                height=0, colour = "grey", linewidth = 0.3) +
  geom_point(data = brood_t,
             aes(x = EM_med/100000, y = R_med/1000, color=BroodYear), size = 1.5) +
  geom_line(data = ER.pred, 
            aes(x = em/100000, y = Rec_med/1000)) +
  scale_colour_viridis_c(name = "Brood Year")+
  labs(x = "Egg mass (00 kgs)", y = "Recruits (000s)") +
  theme_sleek()

ggsave("plots/skeena/kitsumkalum-EMR.jpeg",width=7,height=5)

#benchmarks based on early demographics
bench_early <- matrix(NA,length(EM_pars$lnalpha),2,
                dimnames = list(seq(1:length(EM_pars$lnalpha)), c("Smsy","Smsr")))

for(j in 1:length(EM_pars$lnalpha)){
  ln_a <- EM_pars$lnalpha[j]
  beta <- EM_pars$beta[j]
  year_index = c(1:6) # first 6 brood years
  eta_as_F <- colMeans(fem_S_comps[year_index,]) # probability of returning by age for females
  eta_as_M <- colMeans(male_S_comps[year_index,]) # probability of returning by age for males
  eta_as <- c(eta_as_F,eta_as_M) # probability of returning by age and sex
  z_as <- c(colMeans(fem_egg_mass_comps[year_index,]),c(0,0,0)) #  age/sex class-specific reproductive output
  opt <- optimize(calcCeq, c(0,5), maximum=TRUE)
  bench_early[j,1] <- sum(Seq_as)
  opt <- optimize(calcNeq, c(0,5), maximum=TRUE)
  bench_early[j,2] <- sum(Seq_as)
}
bench_early <- as.data.frame(bench_early)
bench_early$Smsy[bench_early$Smsy <0] <- NA
bench_early$Smsr[bench_early$Smsr <0] <- NA
bench_early$period <- "early"
bench_early$CU <- "Kitsumkalum"

#benchamrks based on long-term average demographics
bench_avg <- matrix(NA,length(EM_pars$lnalpha),2,
                      dimnames = list(seq(1:length(EM_pars$lnalpha)), c("Smsy","Smsr")))

for(j in 1:length(EM_pars$lnalpha)){
  ln_a <- EM_pars$lnalpha[j]
  beta <- EM_pars$beta[j]
  year_index = c(1:nyrs)
  eta_as_F <- colMeans(fem_S_comps[year_index,]) # probability of returning by age for females
  eta_as_M <- colMeans(male_S_comps[year_index,]) # probability of returning by age for males
  eta_as <- c(eta_as_F,eta_as_M) # probability of returning by age and sex
  z_as <- c(colMeans(fem_egg_mass_comps[year_index,]),c(0,0,0)) #  age/sex class-specific reproductive output
  opt <- optimize(calcCeq, c(0,5), maximum=TRUE)
  bench_avg[j,1] <- sum(Seq_as)
  opt <- optimize(calcNeq, c(0,5), maximum=TRUE)
  bench_avg[j,2] <- sum(Seq_as)
}
bench_avg <- as.data.frame(bench_avg)
bench_avg$Smsy[bench_avg$Smsy <0] <- NA
bench_avg$Smsr[bench_avg$Smsr <0] <- NA
bench_avg$period <- "avg"
bench_avg$CU <- "Kitsumkalum"

#benchmarks based on recent demographics
bench_recent <- matrix(NA,length(EM_pars$lnalpha),2,
                    dimnames = list(seq(1:length(EM_pars$lnalpha)), c("Smsy","Smsr")))

for(j in 1:length(EM_pars$lnalpha)){
  ln_a <- EM_pars$lnalpha[j]
  beta <- EM_pars$beta[j]
  year_index = c((nyrs-7):nyrs) # last 6 brood years
  eta_as_F <- colMeans(fem_S_comps[year_index,]) # probability of returning by age for females
  eta_as_M <- colMeans(male_S_comps[year_index,]) # probability of returning by age for males
  eta_as <- c(eta_as_F,eta_as_M) # probability of returning by age and sex
  z_as <- c(colMeans(fem_egg_mass_comps[year_index,]),c(0,0,0)) #  age/sex class-specific reproductive output
  opt <- optimize(calcCeq, c(0,5), maximum=TRUE)
  bench_recent[j,1] <- sum(Seq_as)
  opt <- optimize(calcNeq, c(0,5), maximum=TRUE)
  bench_recent[j,2] <- sum(Seq_as)
}
bench_recent <- as.data.frame(bench_recent)
bench_recent$Smsy[bench_recent$Smsy <0] <- NA
bench_recent$Smsr[bench_recent$Smsr <0] <- NA
bench_recent$period <- "recent"
bench_recent$CU <- "Kitsumkalum"

bench_AR1_eggmass <- rbind(bench_early, bench_avg, bench_recent)

colnames(ER.pred) <- c("EM", "Rec_lwr","Rec_med","Rec_upr")

bench_egg_mass_long <- pivot_longer(bench_AR1_eggmass, cols = c(Smsy, Smsr), names_to = "par")

bench_egg_mass_long_summary <- bench_egg_mass_long |>
  group_by(period, par) |>
  summarize(median = quantile(value,probs=c(0.5),na.rm=T),
            lower = quantile(value,probs=c(0.025),na.rm=T),
            upper = quantile(value,probs=c(0.975),na.rm=T))

ggplot(bench_egg_mass_long |> filter(par=="Smsy") ,aes(value, fill = period, color = period)) +
  geom_density(alpha = 0.8, adjust = 4) +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  theme_sleek() +
  coord_cartesian(xlim = c(NA,12000)) +
  xlab("Smsy") 

ggsave("plots/skeena/kitsumkalum-demographic-Smsy-post.jpeg",width=7,height=5)


# Nass Aggregate----
male_S_comps <- read.csv(here("data/Nass_male_age_comps.csv"))[,3:5] |> # male spawner age composition
  mutate(X6 = replace_na(X6, 0))

fem_S_comps <- read.csv(here("data/Nass_fem_age_comps.csv"))[,3:5] # female spawner age composition
fem_egg_mass_comps <- read.csv(here("data/Nass_fem_age_eggg_mass.csv"))[,2:4] # female spawner egg mass by age

AR1.EM.fit <- readRDS("data/generated/model_fits/EM/Nass Aggregate_EM.rds")

bench_AR1_eggs <- NULL
ER.pred <- NULL
brood.all <- NULL
a.yrs.all <- NULL
demo.par.summary.out <- NULL

sub_dat <- filter(sp_har, CU=="Nass Aggregate")
a_min <- 4
a_max <- 6
nyrs <- max(sub_dat$year)-min(sub_dat$year)+1 #number of calendar years of observations
A <- a_max - a_min + 1 #total age classes
nRyrs <- nyrs + A - 1 #total number of fully and partially observed recruitments from observed spawners

#load AR1 model with egg mass parameters ---
EM_pars <- rstan::extract(AR1.EM.fit)

#eggmass-recruit relationship
egg.quant <- apply(exp(EM_pars$lnEM), 2, quantile, probs=c(0.1,0.5,0.9))[,1:(nyrs-a_min)]
rec.quant <- apply(EM_pars$R, 2, quantile, probs=c(0.1,0.5,0.9))[,(A+a_min):nRyrs]

brood_t <- as.data.frame(cbind(sub_dat$year[1:(nyrs-A)],t(egg.quant), t(rec.quant))) |>
  round(2)

colnames(brood_t) <- c("BroodYear","EM_lwr","EM_med","EM_upr","R_lwr","R_med","R_upr")


#SR relationship based on full posterior---
em <- seq(0,max(brood_t$EM_upr),length.out=100)
ER.pred <- matrix(NA,length(em), length(EM_pars$lnalpha))

for(j in 1:length(EM_pars$lnalpha)){
  ln_a <- EM_pars$lnalpha[j]
  b <- EM_pars$beta[j]
  ER.pred[,j] <- (exp(ln_a)*em*exp(-b*em))
}

ER.pred <- as.data.frame(cbind(em,t(apply(ER.pred, 1, quantile, probs=c(0.1,0.5,0.9), na.rm=T))))|>
  round(2) 

colnames(ER.pred) <- c("em","Rec_lwr","Rec_med","Rec_upr")

# plot EM recruitment relationship 

ggplot() +
  geom_ribbon(data = ER.pred, 
              aes(x = em/100000, ymin = Rec_lwr/1000, ymax = Rec_upr/1000),
              fill = "grey80", alpha=0.5, linetype=2, colour="gray46") +
  geom_errorbar(data = brood_t, 
                aes(x= EM_med/100000, y = R_med/1000, ymin = R_lwr/1000, ymax = R_upr/1000),
                colour="grey", width=0, linewidth=0.3) +
  geom_errorbar(data = brood_t, 
                aes(y = R_med/1000, xmin = EM_lwr/100000, xmax = EM_upr/100000),
                height=0, colour = "grey", linewidth = 0.3) +
  geom_point(data = brood_t,
             aes(x = EM_med/100000, y = R_med/1000, color=BroodYear), size = 1.5) +
  geom_line(data = ER.pred, 
            aes(x = em/100000, y = Rec_med/1000)) +
  scale_colour_viridis_c(name = "Brood Year")+
  labs(x = "Egg mass (00 kgs)", y = "Recruits (000s)") +
  theme_sleek()

ggsave("plots/nass/nass-aggregate-EMR.jpeg",width=7,height=5)

#benchmarks based on early demographics
bench_early <- matrix(NA,length(EM_pars$lnalpha),2,
                      dimnames = list(seq(1:length(EM_pars$lnalpha)), c("Smsy","Smsr")))

for(j in 1:length(EM_pars$lnalpha)){
  ln_a <- EM_pars$lnalpha[j]
  beta <- EM_pars$beta[j]
  year_index = c(1:6) # first 6 brood years
  eta_as_F <- colMeans(fem_S_comps[year_index,]) # probability of returning by age for females
  eta_as_M <- colMeans(male_S_comps[year_index,]) # probability of returning by age for males
  eta_as <- c(eta_as_F,eta_as_M) # probability of returning by age and sex
  z_as <- c(colMeans(fem_egg_mass_comps[year_index,]),c(0,0,0)) #  age/sex class-specific reproductive output
  opt <- optimize(calcCeq, c(0,5), maximum=TRUE)
  bench_early[j,1] <- sum(Seq_as)
  opt <- optimize(calcNeq, c(0,5), maximum=TRUE)
  bench_early[j,2] <- sum(Seq_as)
}
bench_early <- as.data.frame(bench_early)
bench_early$Smsy[bench_early$Smsy <0] <- NA
bench_early$Smsr[bench_early$Smsr <0] <- NA
bench_early$period <- "early"
bench_early$CU <- "Upper Nass"

#benchamrks based on long-term average demographics
bench_avg <- matrix(NA,length(EM_pars$lnalpha),2,
                    dimnames = list(seq(1:length(EM_pars$lnalpha)), c("Smsy","Smsr")))

for(j in 1:length(EM_pars$lnalpha)){
  ln_a <- EM_pars$lnalpha[j]
  beta <- EM_pars$beta[j]
  year_index = c(1:nyrs)
  eta_as_F <- colMeans(fem_S_comps[year_index,]) # probability of returning by age for females
  eta_as_M <- colMeans(male_S_comps[year_index,]) # probability of returning by age for males
  eta_as <- c(eta_as_F,eta_as_M) # probability of returning by age and sex
  z_as <- c(colMeans(fem_egg_mass_comps[year_index,]),c(0,0,0)) #  age/sex class-specific reproductive output
  opt <- optimize(calcCeq, c(0,5), maximum=TRUE)
  bench_avg[j,1] <- sum(Seq_as)
  opt <- optimize(calcNeq, c(0,5), maximum=TRUE)
  bench_avg[j,2] <- sum(Seq_as)
}
bench_avg <- as.data.frame(bench_avg)
bench_avg$Smsy[bench_avg$Smsy <0] <- NA
bench_avg$Smsr[bench_avg$Smsr <0] <- NA
bench_avg$period <- "avg"
bench_avg$CU <- "Upper Nass"

#benchmarks based on recent demographics
bench_recent <- matrix(NA,length(EM_pars$lnalpha),2,
                       dimnames = list(seq(1:length(EM_pars$lnalpha)), c("Smsy","Smsr")))

for(j in 1:length(EM_pars$lnalpha)){
  ln_a <- EM_pars$lnalpha[j]
  beta <- EM_pars$beta[j]
  year_index = c((nyrs-7):nyrs) # last 6 brood years
  eta_as_F <- colMeans(fem_S_comps[year_index,]) # probability of returning by age for females
  eta_as_M <- colMeans(male_S_comps[year_index,]) # probability of returning by age for males
  eta_as <- c(eta_as_F,eta_as_M) # probability of returning by age and sex
  z_as <- c(colMeans(fem_egg_mass_comps[year_index,]),c(0,0,0)) #  age/sex class-specific reproductive output
  opt <- optimize(calcCeq, c(0,5), maximum=TRUE)
  bench_recent[j,1] <- sum(Seq_as)
  opt <- optimize(calcNeq, c(0,5), maximum=TRUE)
  bench_recent[j,2] <- sum(Seq_as)
}
bench_recent <- as.data.frame(bench_recent)
bench_recent$Smsy[bench_recent$Smsy <0] <- NA
bench_recent$Smsr[bench_recent$Smsr <0] <- NA
bench_recent$period <- "recent"
bench_recent$CU <- "Upper Nass"

bench_AR1_eggmass <- rbind(bench_early, bench_avg, bench_recent)

colnames(ER.pred) <- c("EM", "Rec_lwr","Rec_med","Rec_upr")

bench_egg_mass_long <- pivot_longer(bench_AR1_eggmass, cols = c(Smsy, Smsr), names_to = "par")

bench_egg_mass_long_summary <- bench_egg_mass_long |>
  group_by(period, par) |>
  summarize(median = quantile(value,probs=c(0.5),na.rm=T),
            lower = quantile(value,probs=c(0.025),na.rm=T),
            upper = quantile(value,probs=c(0.975),na.rm=T))

ggplot(bench_egg_mass_long |> filter(par=="Smsy") ,aes(value, fill = period, color = period)) +
  geom_density(alpha = 0.8, adjust = 4) +
  coord_cartesian(xlim = c(4000,17000)) +
  xlab("Smsy") +
  scale_fill_brewer(palette = "Dark2",
                    labels = c(
                      "avg" = "1994-2025", 
                      "early" = "1994-1999", 
                      "recent" = "2020-2025"
                    )) +
  scale_color_brewer(palette = "Dark2",
                     labels = c(
    "avg" = "1994-2025", 
    "early" = "1994-1999", 
    "recent" = "2020-2025"
  )) +
  theme_sleek()

ggsave("plots/nass/nass-aggregate-demographic-Smsy-post.jpeg",width=7,height=5)



# plot regular spawner-recruit relationship

source(here("R/SR_posterior_processing.R")) #can take some time

ggplot() +
  geom_abline(intercept = 0, slope = 1,col="dark grey") +
  geom_ribbon(data = filter(SR.preds, CU == "Nass Aggregate"), 
              aes(x = Spawn/1000, ymin = Rec_lwr/1000, ymax = Rec_upr/1000),
              fill = "grey80", alpha=0.5, linetype=2, colour="gray46") +
  geom_errorbar(data = filter(brood.all, CU == "Upper Nass"), 
                aes(x= S_med/1000, y = R_med/1000, ymin = R_lwr/1000, ymax = R_upr/1000),
                colour="grey", width=0, linewidth=0.3) +
  geom_errorbar(data = filter(brood.all, CU == "Upper Nass"), 
                aes(y = R_med/1000, xmin = S_lwr/1000, xmax = S_upr/1000),
                height=0, colour = "grey", linewidth = 0.3) +
  geom_point(data = filter(brood.all, CU == "Upper Nass"),
             aes(x = S_med/1000, y = R_med/1000, color=BroodYear), size = 1.5) +
  geom_line(data = filter(SR.preds, CU == "Upper Nass"), 
            aes(x = Spawn/1000, y = Rec_med/1000)) +
  scale_colour_viridis_c(name = "Brood Year")+
  labs(x = "Spawners (000s)", y = "Recruits (000s)") +
  coord_cartesian(xlim = c(0,33)) +
  theme_sleek()

ggsave("plots/nass/nass-aggregate-SR.jpeg",width=7,height=5)

ggplot() +
  geom_abline(intercept = 0, slope = 1,col="dark grey") +
  geom_ribbon(data = filter(SR.preds, CU == "Kitsumkalum"), 
              aes(x = Spawn/1000, ymin = Rec_lwr/1000, ymax = Rec_upr/1000),
              fill = "grey80", alpha=0.5, linetype=2, colour="gray46") +
  geom_errorbar(data = filter(brood.all, CU == "Kitsumkalum"), 
                aes(x= S_med/1000, y = R_med/1000, ymin = R_lwr/1000, ymax = R_upr/1000),
                colour="grey", width=0, linewidth=0.3) +
  geom_errorbar(data = filter(brood.all, CU == "Kitsumkalum"), 
                aes(y = R_med/1000, xmin = S_lwr/1000, xmax = S_upr/1000),
                height=0, colour = "grey", linewidth = 0.3) +
  geom_point(data = filter(brood.all, CU == "Kitsumkalum"),
             aes(x = S_med/1000, y = R_med/1000, color=BroodYear), size = 1.5) +
  geom_line(data = filter(SR.preds, CU == "Kitsumkalum"), 
            aes(x = Spawn/1000, y = Rec_med/1000)) +
  scale_colour_viridis_c(name = "Brood Year")+
  labs(x = "Spawners (000s)", y = "Recruits (000s)") +
  theme_sleek()

ggsave("plots/skeena/kitsumkalum-SR.jpeg",width=7,height=5)
