

# make key plots for pub -----------------------------------------------------------------

# SR fits ----
#Skeena
ggplot() +
  geom_abline(intercept = 0, slope = 1,col="dark grey") +
  geom_ribbon(data = filter(SR.preds, SMU == "Skeena"), 
              aes(x = Spawn/1000, ymin = Rec_lwr/1000, ymax = Rec_upr/1000),
              fill = "grey80", alpha=0.5, linetype=2, colour="gray46") +
  geom_errorbar(data = filter(brood.all, SMU == "Skeena"), 
                aes(x= S_med/1000, y = R_med/1000, ymin = R_lwr/1000, ymax = R_upr/1000),
                colour="grey", width=0, linewidth=0.3) +
  geom_errorbarh(data = filter(brood.all, SMU == "Skeena"), 
                 aes(y = R_med/1000, xmin = S_lwr/1000, xmax = S_upr/1000),
                 height=0, colour = "grey", linewidth = 0.3) +
  geom_point(data = filter(brood.all, SMU == "Skeena"),
             aes(x = S_med/1000, y = R_med/1000, color=BroodYear), size = 1.5) +
  geom_line(data = filter(SR.preds, SMU == "Skeena"), 
            aes(x = Spawn/1000, y = Rec_med/1000)) +
  facet_wrap(~CU, scales = "free") +
  scale_colour_viridis_c(name = "Brood Year")+
  labs(x = "Spawners (000s)", y = "Recruits (000s)") +
  theme_sleek()

my.ggsave(here("plots/Skeena/SR_fits_AR1.PNG"))

#Nass
ggplot() +
  geom_abline(intercept = 0, slope = 1,col="dark grey") +
  geom_ribbon(data = filter(SR.preds, SMU == "Nass"), 
              aes(x = Spawn/1000, ymin = Rec_lwr/1000, ymax = Rec_upr/1000),
              fill = "grey80", alpha=0.5, linetype=2, colour="gray46") +
  geom_errorbar(data = filter(brood.all, SMU == "Nass"), 
                aes(x= S_med/1000, y = R_med/1000, ymin = R_lwr/1000, ymax = R_upr/1000),
                colour="grey", width=0, linewidth=0.3) +
  geom_errorbarh(data = filter(brood.all, SMU == "Nass"), 
                 aes(y = R_med/1000, xmin = S_lwr/1000, xmax = S_upr/1000),
                 height=0, colour = "grey", linewidth = 0.3) +
  geom_point(data = filter(brood.all, SMU == "Nass"),
             aes(x = S_med/1000, y = R_med/1000, color=BroodYear), size = 1.5) +
  geom_line(data = filter(SR.preds, SMU == "Nass"), 
            aes(x = Spawn/1000, y = Rec_med/1000)) +
  facet_wrap(~CU, nrow=2,scales = "free") +
  scale_colour_viridis_c(name = "Brood Year")+
  labs(x = "Spawners (000s)", y = "Recruits (000s)") +
  theme_sleek()

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
  geom_abline(intercept = 0, slope = 0, col = "dark grey", lty = 2) +
  theme_sleek()

my.ggsave(here("plots/AR1_resids.PNG"))

# TV resids ----
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
  geom_abline(intercept = 0, slope = 0, col = "dark grey", lty = 2) +
  theme_sleek()

my.ggsave(here("plots/TV_resids.PNG"))

# TV alpha ----
#Skeena

a.yrs.all |>
  filter(SMU == "Skeena", 
         brood_year < max(a.yrs.all$brood_year)-a_max) |> # index fully observed rec
  ggplot() +
  geom_ribbon(aes(x = brood_year, ymin = lwr, ymax = upr), fill = "lightgrey") +
  geom_line(aes(x = brood_year , y = mid)) +
  geom_hline(yintercept = 1, lty=2, col = "grey") +
  labs(y ="Productivity (\U03B1) and 80th percentile", x = "Brood year")+
  facet_wrap(~factor(CU, levels = skeena_order),scales = "free_y") +
  theme_sleek() +
  theme(legend.position = "bottom") 

my.ggsave(here("plots/Skeena/changing_productivity.PNG"))

#Nass
a.yrs.all |>
  filter(SMU == "Nass", 
         brood_year < max(a.yrs.all$brood_year)-a_max) |> # index fully observed rec
  ggplot() +
  geom_ribbon(aes(x = brood_year, ymin = lwr, ymax = upr), fill = "lightgrey") +
  geom_line(aes(x = brood_year , y = mid)) +
  geom_hline(yintercept = 1, lty=2, col = "grey") +
  labs(y ="Productivity (\U03B1) and 80th percentile", x = "Brood year")+
  facet_wrap(~factor(CU, levels = nass_order),scales = "free_y") +
  theme_sleek() +
  theme(legend.position = "bottom") 


my.ggsave(here("plots/Nass/changing_productivity.PNG"), height=3.5)

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
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size=9),
        legend.text = element_text(size=8)) +
  theme_sleek()

my.ggsave(here("plots/TV_SR_fits.PNG"))

# deleted yukon code below this fig. To see more look here 
# (https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/R/inference_figs.R#L383)


#post hoc plots for skeena-nass ----------------------------------------------------------

# the data ---
plot_data <- sp_har |>
  select(CU, year, spwn, harv) |>
  pivot_longer(cols = c(spwn, harv), names_to = "obs") |>
  mutate(obs = ifelse(obs == "spwn", "Spawners", "Harvest"))

#Skeena
ggplot(filter(plot_data, !(CU %in% c("Lower Nass", "Middle Nass", "Upper Nass"))), 
       aes(year, value/1000, color = CU)) +
  geom_point() +
  geom_line() +
  scale_color_viridis_d() +
  facet_wrap(~obs, nrow = 2) +
  labs(x = "Brood year", y = "Fish (thousands)") +
  theme_sleek()
  
my.ggsave(here("plots/Skeena/sp_har_data.PNG"))

#Nass
ggplot(filter(plot_data, CU %in% c("Lower Nass", "Middle Nass", "Upper Nass")), 
       aes(year, value/1000, color = CU)) +
  geom_point() +
  geom_line() +
  scale_color_viridis_d(end = 0.6) +
  facet_wrap(~obs, nrow = 2) +
  labs(x = "Brood year", y = "Fish (thousands)") +
  theme_sleek()

my.ggsave(here("plots/Nass/sp_har_data.PNG"))

#latent states and data by CU ---
AR1.obs <- AR1.harv |>
  mutate(obs = "Harvest") |>
  rename(S.25 = H.25, 
         S.50 = H.50, 
         S.75 = H.75) |>
  bind_rows(AR1.spwn) |> 
  mutate(obs = ifelse(is.na(obs), "Spawners", obs))

colnames(AR1.obs) <- c("p.25", "p.50", "p.75", "year", "CU", "obs")

# Harvest obs and data ---

#Skeena
ggplot(filter(AR1.obs, !(CU %in% c("Lower Nass", "Middle Nass", "Upper Nass")), obs == "Harvest")) +
  geom_ribbon(aes(x = year, ymin = p.25/1000, ymax = p.75/1000), alpha = 0.2) +
  geom_line(aes(x = year, y = p.50/1000)) +
  geom_point(data=filter(plot_data, !(CU %in% c("Lower Nass", "Middle Nass", "Upper Nass")), obs == "Harvest"), 
                         aes(x=year, y=value/1000)) +
  facet_wrap(~CU,scales = "free_y", ncol = 2) +
  labs(x = "Brood year", y = "Harvest (thousands)") +
  theme_sleek()
  
my.ggsave(here("plots/Skeena/har_obs.PNG"), height = 10, width = 9)

#Nass
ggplot(filter(AR1.obs, CU %in% c("Lower Nass", "Middle Nass", "Upper Nass"), obs == "Harvest")) +
  geom_ribbon(aes(x = year, ymin = p.25/1000, ymax = p.75/1000), alpha = 0.2) +
  geom_line(aes(x = year, y = p.50/1000)) +
  geom_point(data=filter(plot_data, CU %in% c("Lower Nass", "Middle Nass", "Upper Nass"), obs == "Harvest"), 
             aes(x=year, y=value/1000)) +
  facet_wrap(~CU,scales = "free_y", ncol = 2) +
  labs(x = "Brood year", y = "Harvest (thousands)") +
  theme_sleek()

my.ggsave(here("plots/Nass/har_obs.PNG"), height = 10, width = 9)

# Spawner obs and data ---

#Skeena
ggplot(filter(AR1.obs, !(CU %in% c("Lower Nass", "Middle Nass", "Upper Nass")), obs == "Spawners")) +
  geom_ribbon(aes(x = year, ymin = p.25/1000, ymax = p.75/1000), alpha = 0.2) +
  geom_line(aes(x = year, y = p.50/1000)) +
  geom_point(data=filter(plot_data, !(CU %in% c("Lower Nass", "Middle Nass", "Upper Nass")), obs == "Spawners"), 
             aes(x=year, y=value/1000)) +
  facet_wrap(~CU,scales = "free_y", ncol = 2) +
  labs(x = "Brood year", y = "Spawners (thousands)") +
  theme_sleek()

my.ggsave(here("plots/Skeena/spwn_obs.PNG"), height = 10, width = 9)

#Nass
ggplot(filter(AR1.obs, CU %in% c("Lower Nass", "Middle Nass", "Upper Nass"), obs == "Spawners")) +
  geom_ribbon(aes(x = year, ymin = p.25/1000, ymax = p.75/1000), alpha = 0.2) +
  geom_line(aes(x = year, y = p.50/1000)) +
  geom_point(data=filter(plot_data, CU %in% c("Lower Nass", "Middle Nass", "Upper Nass"), obs == "Spawners"), 
             aes(x=year, y=value/1000)) +
  facet_wrap(~CU,scales = "free_y", ncol = 1) +
  labs(x = "Brood year", y = "Spawners (thousands)")+
  theme_sleek()

my.ggsave(here("plots/Nass/spwn_obs.PNG"), height = 10, width = 9)

#plot A_obs 
A_plot <- sp_har |>
  select(CU, year, a4, a5, a6) |>
  pivot_longer(cols = a4:a6)

#Skeena
ggplot(filter(A_plot, !(CU %in% c("Lower Nass", "Middle Nass", "Upper Nass"))), 
       aes(x=year, y=value/100, fill = name)) +
  geom_area(alpha = 0.8) +  
  facet_wrap(~CU) +
  scale_fill_viridis_d(name = "Age class") +
  labs(x = "Return year", y = "Proportion at age") +
  theme_sleek()

my.ggsave(here("plots/Skeena/A_obs.PNG"))

#Nass
ggplot(filter(A_plot, CU %in% c("Lower Nass", "Middle Nass", "Upper Nass")), 
       aes(x=year, y=value/100, fill = name)) +
  geom_area(alpha = 0.8) +  
  facet_wrap(~CU, nrow=2) +
  scale_fill_viridis_d(name = "Age class") +
  labs(x = "Return year", y = "Proportion at age")+
  theme_sleek()

my.ggsave(here("plots/Nass/A_obs.PNG"))

# plot latent states of prop at age ---
A_plot <- A_plot |>
  mutate(age = case_when(name == "a4" ~ 4, 
                         name == "a5" ~ 5, 
                         name == "a6" ~ 6))

#individual fit EXAMPLE
ggplot(filter(AR1.ages, CU == "Kitsumkalum"), aes(year, A.50)) +
  geom_ribbon(aes(ymin=A.25, ymax=A.75), alpha = 0.4) +
  geom_point(data = filter(A_plot, CU == "Kitsumkalum"), aes(x=year, y= value/100)) +
  geom_line() +
  facet_wrap(~age, nrow=3) +
  labs(y= "Proportion at age", x = "Return year")+
  theme_sleek()

my.ggsave(here("plots/Skeena/prop_age_example.PNG"), height = 10, width = 9)

ggplot(filter(AR1.ages, CU == "Lower Nass"), aes(year, A.50)) +
  geom_ribbon(aes(ymin=A.25, ymax=A.75), alpha = 0.4) +
  geom_point(data = filter(A_plot, CU == "Lower Nass"), aes(x=year, y= value/100)) +
  geom_line() +
  facet_wrap(~age, nrow=3) +
  labs(y= "Proportion at age", x = "Return year")+
  theme_sleek()

my.ggsave(here("plots/Nass/prop_age_example.PNG"), height = 10, width = 9)

# posterior distribution of reference points ----
bench.long <- pivot_longer(bench.posts, cols = c(Smsr.20, Smsr.40, S.recent), names_to = "par") |>
  arrange(CU, par, value) |>
  filter(value <= 10000)

#Skeena 
b <- ggplot(filter(bench.long, !(CU %in% c("Lower Nass", "Middle Nass", "Upper Nass"))),  
            aes(Smsr/1000, fill = CU, color = CU)) +
  geom_density(alpha = 0.3,bw=0.6) +
  theme(legend.position = "bottom") +
  labs(x = expression(italic(S[MSR])), y = "Posterior density") +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme_sleek()   +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title=element_blank(),
        legend.position=c(0.68, 0.65)) +
  guides(fill=guide_legend(ncol=2, theme = theme(legend.byrow = TRUE)),
         colour=guide_legend(ncol=2, theme = theme(legend.byrow = TRUE))) +
  scale_x_continuous(limits = c(0, 30))

c <- ggplot(filter(bench.long, !(CU %in% c("Lower Nass", "Middle Nass", "Upper Nass"))), 
            aes(Umsy, fill = CU, color = CU)) +
  geom_density(alpha = 0.3,bw=0.03) +
  theme(legend.position = "bottom") +
  labs(x = expression(italic(U[MSY])), y = "Posterior density") +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme_sleek()   +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title=element_blank(),
        legend.position="none") +
  scale_x_continuous(limits = c(0, 1))

par.long <- par.posts |>
  mutate(alpha = exp(ln_a))

a <- ggplot(filter(par.long, !(CU %in% c("Lower Nass", "Middle Nass", "Upper Nass"))), 
            aes(alpha, fill = CU, color = CU)) +
  geom_density(alpha = 0.3,bw=0.4) +
  labs(x = "Intrinsic productivity", y = "Posterior density") +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme_sleek()   +
  theme(legend.position = "none",
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title=element_blank()) +
  scale_x_continuous(limits = c(0, 12))

cowplot::plot_grid(a, b, c, labels="auto", ncol=1)
my.ggsave(here("plots/Skeena/ref_pts.PNG"))

#Nass
b <- ggplot(filter(bench.long, CU %in% c("Lower Nass", "Middle Nass", "Upper Nass")),  
            aes(Smsr/1000, fill = CU, color = CU)) +
  geom_density(alpha = 0.3,bw=0.6) +
  theme(legend.position = "bottom") +
  labs(x = expression(italic(S[MSR])), y = "Posterior density") +
  scale_color_viridis_d(end = 0.6) +
  scale_fill_viridis_d(end = 0.6) +
  theme_sleek()   +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title=element_blank(),
        legend.position=c(0.68, 0.65)) +
  guides(fill=guide_legend(ncol=2, theme = theme(legend.byrow = TRUE)),
         colour=guide_legend(ncol=2, theme = theme(legend.byrow = TRUE))) +
  scale_x_continuous(limits = c(0, 30))

c <- ggplot(filter(bench.long, CU %in% c("Lower Nass", "Middle Nass", "Upper Nass")), 
            aes(Umsy, fill = CU, color = CU)) +
  geom_density(alpha = 0.3,bw=0.03) +
  theme(legend.position = "bottom") +
  labs(x = expression(italic(U[MSY])), y = "Posterior density") +
  scale_color_viridis_d(end = 0.6) +
  scale_fill_viridis_d(end = 0.6) +
  theme_sleek()   +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title=element_blank(),
        legend.position="none") +
  scale_x_continuous(limits = c(0, 1))

a <- ggplot(filter(par.long, CU %in% c("Lower Nass", "Middle Nass", "Upper Nass")), 
            aes(alpha, fill = CU, color = CU)) +
  geom_density(alpha = 0.3,bw=0.4) +
  labs(x = "Intrinsic productivity", y = "Posterior density") +
  scale_color_viridis_d(end = 0.6) +
  scale_fill_viridis_d(end = 0.6) +
  theme_sleek()   +
  theme(legend.position = "none",
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title=element_blank()) +
  scale_x_continuous(limits = c(0, 12))

cowplot::plot_grid(a, b, c, labels="auto", ncol=1)
my.ggsave(here("plots/Nass/ref_pts.PNG"))


#spawners through time relative to ref points
# Skeena
bench_plot <- filter(bench.par.table, SMU == "Skeena", 
                     par %in% c("Sgen", "Smsy.80")) |>
  select(CU, par, median)

ggplot(filter(AR1.spwn, SMU == "Skeena"), aes(year, S.50/1000)) +
  geom_line() +
  geom_ribbon(aes(ymin = S.25/1000, ymax = S.75/1000), fill = "darkgrey", alpha = 0.5) +
  geom_hline(data = bench_plot, 
             aes(yintercept = median/1000, col = par), lty = "dashed") +
  geom_vline(xintercept = max(AR1.spwn$year)- a_max, lwd = 0.5, lty = 2) +
  scale_colour_manual(values=c(Smsy.80 = "forestgreen", 
                               Sgen = "darkred"),
                      labels=c(Smsy.80 = expression(paste("80%", S[MSY])),
                               Sgen = expression(S[Gen]))) +
  facet_wrap(~factor(CU, levels = skeena_order), nrow = 3, scales = "free_y") +
  theme_sleek() +
  coord_cartesian(ylim = c(0,NA)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(x = "Retun Year", y = "Spawners (000s)", col = "")

my.ggsave(here("plots/Skeena/spwn-benchmarks-time.PNG"))
# Nass
bench_plot <- filter(bench.par.table, SMU == "Nass", 
                     par %in% c("Sgen", "Smsy.80")) |>
  select(CU, par, median)

ggplot(filter(AR1.spwn, SMU == "Nass"), aes(year, S.50/1000)) +
  geom_line() +
  geom_ribbon(aes(ymin = S.25/1000, ymax = S.75/1000), fill = "darkgrey", alpha = 0.5) +
  geom_hline(data = bench_plot, 
             aes(yintercept = median/1000, col = par), lty = "dashed") +
  geom_vline(xintercept = max(AR1.spwn$year)- a_max, lwd = 0.5, lty = 2) +
  scale_colour_manual(values=c(Smsy.80 = "forestgreen", 
                               Sgen = "darkred"),
                      labels=c(Smsy.80 = expression(paste("80%", S[MSY])),
                               Sgen = expression(S[Gen]))) +
  facet_wrap(~factor(CU, levels = nass_order), nrow = 2, scales = "free_y") +
  theme_sleek() +
  coord_cartesian(ylim = c(0,NA)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(x = "Retun Year", y = "Spawners (000s)", col = "")
my.ggsave(here("plots/Nass/spwn-benchmarks-time.PNG"))
