library(here)
library(tidyverse)
library(viridis)

source(here("R/data_functions"))

skeena_ages <- skrunchy2025::tyee_biodata_age_sex_length_CU |>
  arrange(i, y) |>
  as.data.frame() |>
  rename(CU = i, 
         Year = y, 
         Sex = sex, 
         age = a) |>
  filter(!is.na(CU), !is.na(age))

# age comps by CU ----
ggplot(skeena_ages, aes(x = Year, fill = as.factor(age))) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~CU, nrow  = 2) +
  scale_fill_brewer(palette = "Dark2") +
  theme_sleek() +
  ylab("") + 
  labs(fill = "Age")

# female length at age ----
ggplot(skeena_ages, aes(x = Year, y = nose_fork_length_mm)) +
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

#  sex ratio
skeena_sex_ratio <- skeena_ages |>
  group_by(Year, Sex) |>
  summarize(ratio = sum(prop_count)) #hmm just calc proportions based on this df? aggregate it?

# length at age by year

# female age comp

# avg reproductive output 