library(here)
library(tidyverse)
library(viridis)

source(here("R/data_functions.R"))

skeena_ASL <- skrunchy2025::tyee_biodata_age_sex_length_CU |>
  arrange(i, y) |>
  as.data.frame() |>
  rename(CU = i, 
         Year = y, 
         Sex = sex, 
         age = a) |>
  filter(!is.na(CU), !is.na(age))


length(which(is.na(skeena_ASL$nose_fork_length_mm)))
length(which(is.na(skeena_ASL$hypural_length_mm)))

# exploratory plots ----
# age comps by CU ----
ggplot(skeena_ASL, aes(x = Year, fill = as.factor(age))) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~CU, nrow  = 2) +
  scale_fill_brewer(palette = "Dark2") +
  theme_sleek() +
  ylab("") + 
  labs(fill = "Age")

# female length at age ----
ggplot(skeena_ASL, aes(x = Year, y = hypural_length_mm)) +
  stat_summary(fun = mean, geom = "point", size = 3,col="dark grey") + # Add points for the mean
  stat_summary(
    fun.data = mean_cl_boot, # Calculate mean and standard error
    geom = "errorbar",
    width = 0, # Control the width of the error bar whiskers
    col="dark grey"
  ) +
  labs(x = "Year", y = "Female length \n (nose fork length, mm)")  +
  facet_wrap(~age) +
  geom_smooth(method = "lm", level = 0.95, color = "dark grey", fill = "lightblue") +
  theme_sleek()

# process data for multipanel plot -------------------------------------------------------
kitsumkalum_ASL <- skeena_ASL |>
  filter(CU == "Kitsumkalum", 
         age %in% c(4,5,6), 
         Sex != "u") 
  
kitsumkalum_sex_age_comps <- kitsumkalum_ASL |>
  count(Year, Sex, age) |>
  group_by(Year) |>
  mutate(year_total = sum(n),
         prop_count = (n/year_total)*100) |>
  select(Year, age, Sex, prop_count)

kitsumkalum_sex_ratio <- kitsumkalum_sex_age_comps |>
  group_by(Year, Sex) |>
  summarize(ratio = sum(prop_count))

# wrangle avg fecundity ---
plot(filter(kitsumkalum_ASL, Sex== "f")$nose_fork_length_mm ~ filter(kitsumkalum_ASL, Sex== "f")$hypural_length_mm)
lm(filter(kitsumkalum_ASL, Sex== "f")$nose_fork_length_mm ~ filter(kitsumkalum_ASL, Sex== "f")$hypural_length_mm)

kitsumkalum_ASL_infill <- kitsumkalum_ASL |>
  filter(Sex == "f") |>
  # Convert hypural length to nose fork length
  mutate(nose_fork_length_mm = ifelse(is.na(nose_fork_length_mm), (23.04 + hypural_length_mm*1.19), 
                                      nose_fork_length_mm), 
         MEF = 2.225 + (nose_fork_length_mm/10)* 0.8907) #then convert NFL (in cms) to MEF

kitsumkalum_fem_avg_len <- kitsumkalum_ASL_infill|>
  group_by(Year, age) |>
  summarise(mean_length = mean(MEF, na.rm = T)) |>
  mutate(fecundity = 9.35e-4*(mean_length*10)^2.36, #coefficients from elsewhere (from BC)
         egg_mass = 8.71e-12*(mean_length*10)^4.83) 

kitsumkalum_repro_potential <- kitsumkalum_sex_age_comps |>
  filter(Sex == "f") |>
  left_join(kitsumkalum_fem_avg_len, by=c("Year", "age")) |>
  mutate(sex_cor_fecundity = (prop_count/100)*fecundity, 
         sex_cor_egg_mass = (prop_count/100)*egg_mass) |>
  group_by(Year) |>
  summarize(eggsperspwn = sum(sex_cor_fecundity),
            eggmassperspwn = sum(sex_cor_egg_mass))

a <- ggplot(kitsumkalum_sex_ratio, aes(x = Year, y= ratio, fill= as.factor(Sex))) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Dark2") +
  theme_sleek() +
  ylab("Spawners of given sex") + 
  labs(fill = "Sex") +
  theme(legend.position = "top")

b <- ggplot(kitsumkalum_ASL |> filter(Sex == "f"), aes(x = Year, y = hypural_length_mm)) +
  stat_summary(fun = mean, geom = "point", size = 2, col="dark grey") + # Add points for the mean
  stat_summary(fun.data = mean_cl_boot, # Calculate mean and standard error
               geom = "errorbar",
               width = 0, # Control the width of the error bar whiskers
               col="dark grey") +
  labs(x = "Year", y = "Female length (hypural; mm)") +
  scale_x_continuous(n.breaks = 3) +
  facet_wrap(~age) +
  geom_smooth(method = "lm", level = 0.99, color = "dark grey", fill = "lightblue") +
  theme_sleek()+
  theme(plot.margin = margin(20,10,0.5,10))

c <- ggplot(kitsumkalum_ASL, aes(x = Year, fill = as.factor(age))) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Dark2") +
  theme_sleek() +
  ylab("Female age composition") + 
  labs(fill = "Age") +
  theme(legend.position = "top")

d <- ggplot(kitsumkalum_repro_potential, aes(x = Year, y = eggsperspwn)) +
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
my.ggsave(here("plots/Skeena/repro-potential.PNG"))
