# script for wrangling raw age data from LGL 
library(tidyverse)
library(here)
source(here("R/data_functions.R"))

# read in and wrangle data ----
nass_ages_raw <- read.csv(here("data/Upper lower coastal Nass Chinook ages 1989-2025.csv"))

nass_ages <- nass_ages_raw |>
  select(Year, CU, Location, NFL, MEF, GR_age, Sex) |>
  mutate(salt_age = case_when(GR_age == "21" ~ 1, 
                              GR_age == "31" ~ 2, 
                              GR_age == "41" ~ 3,
                              GR_age == "51" ~ 4, 
                              GR_age == "32" ~ 1, 
                              GR_age == "42" ~ 2,
                              GR_age == "52" ~ 3, 
                              GR_age == "62" ~ 4, 
                              GR_age == "72" ~ 5,
                              GR_age == "43" ~ 1, 
                              GR_age == "53" ~ 2, 
                              GR_age == "63" ~ 3,
                              GR_age == "73" ~ 4, 
                              GR_age == "83" ~ 5),
         fresh_age = case_when(GR_age == "21" ~ 1, 
                              GR_age == "31" ~ 1, 
                              GR_age == "41" ~ 1,
                              GR_age == "51" ~ 1, 
                              GR_age == "32" ~ 2, 
                              GR_age == "42" ~ 2,
                              GR_age == "52" ~ 2, 
                              GR_age == "62" ~ 2, 
                              GR_age == "72" ~ 2,
                              GR_age == "43" ~ 3, 
                              GR_age == "53" ~ 3, 
                              GR_age == "63" ~ 3,
                              GR_age == "73" ~ 3, 
                              GR_age == "83" ~ 4),
         Sex_full = case_when(Sex == "M" ~ "Male",
                              Sex == "F" ~ "Female"),
         age = fresh_age + salt_age)

# exploratory plots ----
# age comps by CU ----
ggplot(nass_ages, aes(x = Year, fill = as.factor(age))) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~CU, nrow  = 2) +
  scale_fill_brewer(palette = "Dark2") +
  theme_sleek() +
  ylab("") + 
  labs(fill = "Age")
  
# age comps by sex for upper CU ----
ggplot(nass_ages |> filter(CU == "Upper Nass", Sex %in% c("M","F")), aes(x = Year, fill = as.factor(age))) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~Sex_full, nrow  = 2) +
  scale_fill_brewer(palette = "Dark2") +
  theme_sleek() +
  ylab("") + 
  labs(fill = "Age")

# age comps by location ----
ggplot(nass_ages, aes(x = Year, fill = as.factor(age))) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~Location, ncol  = 3) +
  scale_fill_brewer(palette = "Dark2") +
  theme_sleek() +
  ylab("") + 
  labs(fill = "Age")

# female length at age ----
ggplot(nass_age_length_raw, aes(x = Year, y = MEF)) +
  stat_summary(fun = mean, geom = "point", size = 3,col="dark grey") + # Add points for the mean
  stat_summary(
    fun.data = mean_cl_boot, # Calculate mean and standard error
    geom = "errorbar",
    width = 0, # Control the width of the error bar whiskers
    col="dark grey"
  ) +
  labs(x = "Year", y = "Female length (MEF)")  +
  facet_wrap(~age) +
  geom_smooth(method = "lm", level = 0.95, color = "dark grey", fill = "lightblue") +
  theme_sleek()

# plot age comps excluding jacks and extremely rare older age classes ----
ggplot(nass_ages |> filter(Location == "GW fishwheels", between(age,4,6)), aes(x = Year, fill = as.factor(age))) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~CU, nrow  = 2) +
  scale_fill_brewer(palette = "Dark2") +
  theme_sleek() +
  ylab("") + 
  labs(fill = "Age")

# age comps for spawner recruitment modelling (use GW fishwheels as only source of data across years) ----
SR_age_comps <- nass_ages |>
  filter(Location == "GW fishwheels",
         between(age,4,6)) |>
  count(Year, age) |>
  group_by(Year) |>
  mutate(year_total = sum(n),
         prop_count = (n/year_total)*100) |>
  select(Year, age, prop_count) |>
  pivot_wider(names_from = age, values_from = prop_count) |>
  rename(a4 = '4', 
         a5 = '5', 
         a6 = '6')
write.csv(SR_age_comps, here("data/NassAobs.16Feb2026.csv"),row.names = FALSE)

# reproductive potential multi-panel plot ----
nass_age_length_raw <- nass_ages |>
  filter(Location == "GW fishwheels",
         between(age,4,6))

# relationships between NFL and MEF to predict MEF in years with only NFL
nass_age_length_raw_fem <- nass_age_length_raw |> filter(Sex == "F")
plot(nass_age_length_raw_fem$MEF ~ nass_age_length_raw_fem$NFL)
lm(nass_age_length_raw_fem$MEF ~ nass_age_length_raw_fem$NFL)

nass_age_length_infill <- nass_age_length_raw |> 
  filter(Sex == "F") |>
  mutate(MEF = case_when(
    Year < 2011 ~ 2.225 + NFL*0.8907,
    Year > 2010 ~ MEF))
 
nas_fem_avg_len <- nass_age_length_infill|>
  group_by(Year, age) |>
  summarise(mean_length = mean(MEF, na.rm = T)) |>
  mutate(fecundity = 9.35e-4*(mean_length*10)^2.36,
         egg_mass = 8.71e-12*(mean_length*10)^4.83)

nass_sex_age_comps <- nass_age_length_raw |>
  filter(Location == "GW fishwheels",
         between(age,4,6)) |>
  count(Year, Sex, age) |>
  group_by(Year) |>
  mutate(year_total = sum(n),
         prop_count = (n/year_total)*100) |>
  select(Year, age, Sex, prop_count)

nass_sex_ratio <- nass_sex_age_comps |>
  group_by(Year, Sex) |>
  summarize(ratio = sum(prop_count))

nass_repro_potential <- nass_sex_age_comps |>
  filter(Sex == "F") |>
  left_join(nas_fem_avg_len, by=c("Year", "age")) |>
  mutate(sex_cor_fecundity = (prop_count/100)*fecundity, 
         sex_cor_egg_mass = (prop_count/100)*egg_mass) |>
  group_by(Year) |>
  summarize(eggsperspwn = sum(sex_cor_fecundity),
            eggmassperspwn = sum(sex_cor_egg_mass))


a <- ggplot(nass_sex_ratio, aes(x = Year, y= ratio/100, fill= as.factor(Sex))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Dark2") +
  theme_sleek() +
  ylab("Spawners of given sex") + 
  labs(fill = "Sex") +
  theme(legend.position = "top")

b <- ggplot(nass_age_length_infill |> filter(Sex == "F"), aes(x = Year, y = MEF)) +
   stat_summary(fun = mean, geom = "point", size = 2, col="dark grey") + # Add points for the mean
   stat_summary(
     fun.data = mean_cl_boot, # Calculate mean and standard error
     geom = "errorbar",
     width = 0, # Control the width of the error bar whiskers
     col="dark grey"
   ) +
   labs(x = "Year", y = "Female length (MEF)")  +
   facet_wrap(~age) +
   geom_smooth(method = "lm", level = 0.99, color = "dark grey", fill = "lightblue") +
   theme_sleek()+
   theme(plot.margin = margin(20,10,0.5,10))
 
c <- ggplot(nass_ages |> filter(Location == "GW fishwheels", between(age,4,6), Sex == "F"), aes(x = Year, fill = as.factor(age))) +
   geom_bar(position = "fill") +
   scale_y_continuous(labels = scales::percent) +
   scale_fill_brewer(palette = "Dark2") +
   theme_sleek() +
   ylab("Female age composition") + 
   labs(fill = "Age") +
   theme(legend.position = "top")
 
d <- ggplot(nass_repro_potential, aes(x = Year, y = eggsperspwn)) +
  geom_smooth(method="lm", color="grey", fill = "lightblue") +
  geom_point(size=2, color="dark grey")+
  xlab("Year") +
  ylab("Average reproductive output \n (total eggs per spawner)") +
  theme_sleek() +
  theme(plot.margin = margin(40,10,0.5,0.5))
 

g <- ggarrange(a,b,c,d,
               labels = c("a", "b","c", "d"),
               heights = c(1,1))


g
my.ggsave(here("plots/Nass/repro-potential.PNG"))

