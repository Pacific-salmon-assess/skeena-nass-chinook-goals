library(here)
library(tidyverse)
library(viridis)
library(ggpubr)

source(here("R/data_functions.R"))

kitsumkalum_ASL <- read.csv(here("data/kitsumkalum-female-brood-age-length-egg-data.csv")) 

kitsumkalum_ASL |>
  group_by(total_age) |>
  summarise(n()) #why so many NAs?, not many 4 either

kitsumkalum_ASL <- kitsumkalum_ASL |>
  filter(total_age %in% c(4,5,6), 
         total_n_green_eggs_taken_at_collection_sum != 0) |>
  select(year, poh_length_cm, total_age, total_n_green_eggs_taken_at_collection_sum, 
         egg_volume_ml)

colnames(kitsumkalum_ASL) <- c("year", "poh_length_cm", "age", "eggs", "egg_volume")

kitsumkalum_sex_ratio <- skrunchy2025::kitsumkalum |>
  mutate(f = female_escapement/(female_escapement+male_escapement), 
         m = male_escapement/(male_escapement+female_escapement)) |>
  select(year, f, m) |>
  pivot_longer(-year)

colnames(kitsumkalum_sex_ratio) <- c("Year", "Sex", "ratio")

# sex ratio
a <- ggplot(kitsumkalum_sex_ratio, aes(x = Year, y= ratio/100, fill= as.factor(Sex))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Dark2") +
  theme_sleek() +
  ylab("Spawners of given sex") + 
  labs(fill = "Sex") +
  theme(legend.position = "top")
  
#length at age 
b <- ggplot(kitsumkalum_ASL, aes(x = year, y = poh_length_cm)) +
  stat_summary(fun = mean, geom = "point", size = 2, col="dark grey") + # Add points for the mean
  stat_summary(fun.data = mean_cl_boot, # Calculate mean and standard error
               geom = "errorbar",
               width = 0, # Control the width of the error bar whiskers
               col="dark grey") +
  labs(x = "Year", y = "Female length (POH; cm)") +
  scale_x_continuous(n.breaks = 3) +
  facet_wrap(~age) +
  geom_smooth(method = "lm", level = 0.99, color = "dark grey", fill = "lightblue") +
  theme_sleek()+
  theme(plot.margin = margin(20,10,0.5,10))

c <- ggplot(kitsumkalum_ASL, aes(x = year, fill = as.factor(age))) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Dark2") +
  theme_sleek() +
  ylab("Female age composition") + 
  labs(fill = "Age") +
  theme(legend.position = "top")

d <- ggplot(kitsumkalum_ASL, aes(x = year, y = eggs)) +
  geom_smooth(method="lm", color="grey", fill = "lightblue") +
  geom_point(size=2, color="dark grey")+
  xlab("Year") +
  ylab("Reproductive output \n (total eggs per spawner)") +
  theme_sleek() +
  theme(plot.margin = margin(40,10,0.5,0.5))

g <- ggarrange(a,b,c,d,
               labels = c("a", "b","c", "d"),
               heights = c(1,1))
g

my.ggsave(here("plots/Skeena/kalum-repro-potential.PNG"))

# length-fecundity relationship ----------------------------------------------------------
mod <- lm(log(kitsumkalum_ASL$egg_volume) ~ log(kitsumkalum_ASL$poh_length_cm)) #fit log-log model

ggplot(kitsumkalum_ASL, aes(poh_length_cm, egg_volume)) +
  geom_point() +
  stat_function( fun = function(x) exp(coef(mod)[1]) * x^coef(mod)[2],
                 color = "black", lwd = 2) +
  labs(x = "POH length (cm)", y = "egg volume", 
       title = "Kitsumkalum length-egg relationship") +
  theme_sleek()
  
