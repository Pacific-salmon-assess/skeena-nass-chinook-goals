# script for wrangling raw age data from LGL 
library(tidyverse)
library(here)
library(ggsidekick)

# read in and wrangle data ----
nass_recon <- read.csv(here("data/Nass_SpHar_16Feb2026.csv"))
coast_nass_esc <- read.csv(here("data/nass_coastal_esc.csv"))

lower_nass_esc <-coast_nass_esc |>
  mutate(
    na_systems = rowSums(is.na(pick(2:5)))
  ) |>
  group_by(year) |>
  mutate(esc = sum(pick(2:5))
  )
  

df <- read.csv(here("data/nass_coastal_esc.csv"))

lower_nass <- matrix(NA,34,1000)

for(i in 1:1000){
  Iknouk_cont = rbeta(1,0.9324756, 7.7788844)
  Ishkeenickh_cont = rbeta(1,1.4512994 , 2.9874005)
  Kincolith_cont = rbeta(1,1.2995996 , 4.4590524)
  Kwinamass_cont = rbeta(1,5.828208 ,12.747912 )
  Kitsault_cont = rbeta(1,0.7336723, 4.1870655)
  
  sum_cont = Iknouk_cont + Ishkeenickh_cont + Kincolith_cont + Kwinamass_cont + Kitsault_cont
  
  Iknouk_cont_std = Iknouk_cont/sum_cont
  Ishkeenickh_cont_std = Ishkeenickh_cont/sum_cont
  Kincolith_cont_std = Kincolith_cont/sum_cont
  Kwinamass_cont_std = Kwinamass_cont/sum_cont
  Kitsault_cont_std = Kitsault_cont/sum_cont
  
  Khutze_cont = rbeta(1,2.1750142, 7.3460290)
  aerial_exp = rnorm(1,2.33,0.62)
  lower_cont = rbeta(1,1.3951783 ,13.5421138)
  
df2 <- df |>
    mutate(
      na_count = rowSums(is.na(across(2:6))),
      col_sum  = case_when(
        na_count == 0 ~ rowSums(across(2:6))*aerial_exp,
        is.na(Iknouk) & na_count == 1 ~ (rowSums(across(2:6), na.rm = TRUE) / (Ishkeenickh_cont_std+Kincolith_cont_std+Kwinamass_cont_std+Kitsault_cont_std))*aerial_exp,
        is.na(Ishkeenickh) & na_count == 1 ~ (rowSums(across(2:6), na.rm = TRUE) / (Iknouk_cont_std+Kincolith_cont_std+Kwinamass_cont_std+Kitsault_cont_std))*aerial_exp,
        is.na(Kincolith)  & na_count == 1 ~ (rowSums(across(2:6), na.rm = TRUE) / (Iknouk_cont_std+Ishkeenickh_cont_std+Kwinamass_cont_std+Kitsault_cont_std))*aerial_exp,
        is.na(Kwinamass) & na_count == 1 ~ (rowSums(across(2:6), na.rm = TRUE) / (Iknouk_cont_std+Ishkeenickh_cont_std+Kincolith_cont_std+Kitsault_cont_std))*aerial_exp,
        is.na(Kitsault)  & na_count == 1 ~ (rowSums(across(2:6), na.rm = TRUE) / (Iknouk_cont_std+Ishkeenickh_cont_std+Kincolith_cont_std+Kwinamass_cont_std))*aerial_exp,
        is.na(Iknouk) & is.na(Kincolith) & na_count == 2 ~ (rowSums(across(2:6), na.rm = TRUE) / (Ishkeenickh_cont_std+Kwinamass_cont_std+Kitsault_cont_std))*aerial_exp,
        is.na(Kwinamass) & is.na(Kitsault) & na_count == 2 ~ (rowSums(across(2:6), na.rm = TRUE) / (Iknouk_cont_std+Ishkeenickh_cont_std+Kincolith_cont_std))*aerial_exp,
        is.na(Ishkeenickh) & is.na(Kincolith) & na_count == 2 ~ (rowSums(across(2:6), na.rm = TRUE) / (Iknouk_cont_std+Kwinamass_cont_std+Kitsault_cont_std))*aerial_exp,
        is.na(Ishkeenickh) & is.na(Iknouk) & na_count == 2 ~ (rowSums(across(2:6), na.rm = TRUE) / (Kincolith_cont_std+Kwinamass_cont_std+Kitsault_cont_std))*aerial_exp,
        is.na(Iknouk) & is.na(Kitsault) & na_count == 2~ (rowSums(across(2:6), na.rm = TRUE) / (Ishkeenickh_cont_std+Kincolith_cont_std+Kwinamass_cont_std))*aerial_exp,
        is.na(Ishkeenickh) & is.na(Iknouk) & is.na(Kincolith) & na_count == 3 ~ (rowSums(across(2:6), na.rm = TRUE) / (Kwinamass_cont_std+Kitsault_cont_std))*aerial_exp,
        is.na(Iknouk) & is.na(Kincolith) & is.na(Kitsault) & na_count == 3 ~ (rowSums(across(2:6), na.rm = TRUE) / (Ishkeenickh_cont_std+Kwinamass_cont_std))*aerial_exp,
        is.na(Ishkeenickh) & is.na(Iknouk) & is.na(Kitsault) & na_count == 3 ~ (rowSums(across(2:6), na.rm = TRUE) / (Kincolith_cont_std+Kwinamass_cont_std))*aerial_exp,
        na_count > 3 ~ (Upper_Nass / (1-lower_cont)) - Upper_Nass,
        TRUE ~ NA_real_
      )
    )|> 
    mutate(
      lower_esc = case_when(
        !is.na(Khutzeeymateen) ~ col_sum + (Khutzeeymateen*aerial_exp),
        is.na(Khutzeeymateen) ~ col_sum/(1-Khutze_cont)
      )
    ) 

lower_nass[,i] <- as.vector(df2$lower_esc)
}

nass_sp_har <- read.csv(here("data/Nass_SpHar_16Feb2026.csv"))
lower_escapement_error <- t(apply(lower_nass,1,quantile, probs=(c(0.5,0.05,0.95))))
lower_escapement_origonal <- nass_sp_har |> filter(CU == "Lower Nass") |> select(spwn, year)
lower_escapement <- cbind(lower_escapement_error, lower_escapement_origonal)

ggplot(lower_escapement, aes(x = year, y = `50%`)) +
  geom_ribbon(aes(ymin = `5%`, ymax = `95%`), fill = "light grey") +
  geom_line(size=1.5, color="grey30")+
  xlab("Year") +
  ylab("Lower Nass reconstructed spawner abundance") +
  geom_line(aes(y=spwn),size=1, color="black", lty=2) +
  theme_sleek()
ggsave("plots/nass/lowerCU-reconstruction.jpeg",width=7,height=5)

Iknouk_cont = rbeta(10000,0.9324756, 7.7788844)
a <- ggplot(data.frame(x = Iknouk_cont), aes(x = Iknouk_cont)) +
  geom_density(fill = "dark grey", size = 0.5) +
  labs(title = "Xnukw (Iknouk)", x = "Contribution", y = "") +
  theme_sleek() +
  theme(axis.text.y = element_blank())

Ishkeenickh_cont = rbeta(10000,1.4512994 , 2.9874005)
b <- ggplot(data.frame(x = Ishkeenickh_cont), aes(x = Ishkeenickh_cont)) +
  geom_density(fill = "dark grey", size = 0.5) +
  labs(title = "Ksi Hlginx (Ishkeenickh)", x = "Contribution", y = "") +
  theme_sleek()+
  theme(axis.text.y = element_blank())

Kincolith_cont = rbeta(10000,1.2995996 , 4.4590524)
c <- ggplot(data.frame(x = Kincolith_cont), aes(x = Kincolith_cont)) +
  geom_density(fill = "dark grey", size = 0.5) +
  labs(title = "Gingolx (Kincolith)", x = "Contribution", y = "") +
  theme_sleek()+
  theme(axis.text.y = element_blank())

Kwinamass_cont = rbeta(10000,5.828208 ,12.747912 )
d <- ggplot(data.frame(x = Kwinamass_cont), aes(x = Kwinamass_cont)) +
  geom_density(fill = "dark grey", size = 0.5) +
  labs(title = "Ksi X'anmas (Kwinamass)", x = "Contribution", y = "") +
  theme_sleek()+
  theme(axis.text.y = element_blank())

Kitsault_cont = rbeta(10000,0.7336723, 4.1870655)
e <- ggplot(data.frame(x = Kitsault_cont), aes(x = Kitsault_cont)) +
  geom_density(fill = "dark grey", size = 0.5) +
  labs(title = "Kitsault", x = "Contribution", y = "") +
  theme_sleek()+
  theme(axis.text.y = element_blank())


Khutze_cont = rbeta(10000,2.1750142, 7.3460290)
f <- ggplot(data.frame(x = Khutze_cont), aes(x = Khutze_cont)) +
  geom_density(fill = "dark grey", size = 0.5) +
  labs(title = "Khutzamateen", x = "Contribution", y = "") +
  theme_sleek()+
  theme(axis.text.y = element_blank())

aerial_exp = rnorm(10000,2.33,0.62)
h <- ggplot(data.frame(x = aerial_exp), aes(x = aerial_exp)) +
  geom_density(fill = "dark grey", size = 0.5) +
  labs(title = "Aerial survey expansion", x = "Expansion factor", y = "") +
  theme_sleek()+
  theme(axis.text.y = element_blank())

lower_cont = rbeta(10000,1.3951783 ,13.5421138)
g <- ggplot(data.frame(x = lower_cont), aes(x = lower_cont)) +
  geom_density(fill = "dark grey", size = 0.5) +
  labs(title = "Lower CU vs agg. Nass ", x = "Contribution", y = "") +
  theme_sleek()+
  theme(axis.text.y = element_blank())

library(cowplot)
plot_grid(a,b,c,d,e,f,g,h, labels = "", ncol = 3)
ggsave("plots/nass/lowerCU-uncertianty-dist.jpeg",width=7,height=5)

# update spawners for lower CU and combine with updated upper CU---

update_upper_recon <- read.csv("data/Upper-esc-harv.csv")
update_lower_recon <- lower_escapement |>
  mutate(spwn_cv= ((`95%`-`5%`)/4)/spwn) |>
  left_join(nass_recon|>filter(CU=="Lower Nass"),by="year") |>
  rename(spwn=`50%`,
         spwn_cv=spwn_cv.x) |>
  select(CU,year,spwn,spwn_cv,harv,SMU)

update_agg_recon <- nass_recon |>
  filter(CU != "Lower Nass") |>
  group_by(year) |>
  summarize(spwn_agg = sum(spwn),
            harv_agg = sum(harv),
            spwn_cv_agg = max(spwn_cv)) |>
  select(year,spwn_agg,spwn_cv_agg,harv_agg) |>
  left_join(update_lower_recon, by = "year") |>
  mutate(spwn_all = spwn_agg + spwn,
         harv_all = harv_agg + harv,
         spwn_cv = spwn_cv_agg,
         CU = "Nass Aggregate",
         SMU = "Nass") |>
  select(CU,year,spwn_all,spwn_cv,harv_all,SMU) |>
  rename(spwn = spwn_all,
         harv = harv_all) |>
  select(CU,year,spwn,spwn_cv,harv,SMU)

update_nass_recon <-rbind(update_upper_recon,update_lower_recon,update_agg_recon)

# write escapement time series to file ----
write.csv(update_nass_recon,"data/Nass_SpHar_3June2026.csv")
