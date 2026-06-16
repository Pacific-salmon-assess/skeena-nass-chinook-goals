# read data
library(here)
source(here("R/SR_posterior_processing.R")) #can take some time
nass_order <- c("Lower Nass", "Upper Nass", "Nass Aggregate")
nass_order_spw <- c("Lower Nass", "Tseax","Upper Nass", "Nass Aggregate")
all_order <- c("Lower Nass", "Upper Nass","Lower Skeena", "Kitsumkalum", "Zymoetz", "Middle Skeena", "Large Lakes", 
               "Upper Skeena")
#skeena_order
tseax<-read.csv("data/Tseax-esc.csv") |>
  rename(spwn = `X50.`,
         spwn_LCI = `X5.`,
         spwn_UCI = `X95.`) |>
  mutate(CU = "Tseax",
         SMU = "Nass",
         harv = NA,    
         X = NA, 
         spwn_cv  = NA,
         harv_cv = NA,    
         a4 = NA,    
         a5  = NA,    
         a6 = NA,
         harv_LCI = NA,
         harv_UCI = NA) |>
  select(CU,SMU,year,spwn,harv,X, spwn_cv, harv_cv,a4 ,a5,a6, spwn_LCI, spwn_UCI, harv_LCI,harv_UCI)


sp_har_fig <- sp_har |>
  group_by(CU, year) |>
  #calculate random observation error, 4000 = iter
  mutate(spwn_LCI = quantile(rnorm(4000, spwn, spwn_cv*spwn), 0.25),  
         spwn_UCI = quantile(rnorm(4000, spwn, spwn_cv*spwn), 0.75), 
         harv_LCI = quantile(rnorm(4000, harv, harv_cv*harv), 0.25),  
         harv_UCI = quantile(rnorm(4000, harv, harv_cv*harv), 0.75),
         er = harv/(harv+spwn))

sp_fig <- rbind(sp_har_fig,as.data.frame(tseax))

ggplot(filter(sp_fig, SMU == "Nass"), aes(year, spwn/1000)) +
  geom_ribbon(aes(xintercept = year, ymin = spwn_LCI/1000, ymax = spwn_UCI/1000), 
              width = 0, fill = "grey") +
  geom_line(size=1) +
  facet_wrap(~factor(CU, levels = nass_order_spw), scales = "free_y") +
  theme_sleek() +
  labs(x = "Return year", y = "Spawners (000s) and 50th percentile")

ggsave("plots/nass/spawner-reconstructions.jpeg",width=6.5,height=4)


ggplot(filter(sp_fig, SMU == "Nass" & CU !="Tseax"), aes(year, er)) +
  geom_line(size=0.75) +
  facet_wrap(~factor(CU, levels = nass_order_spw), ncol=2,scales = "free_y") +
  theme_sleek() +
  labs(x = "Return year", y = "Harvest rate (%)")

ggsave("plots/nass/er-reconstructions.jpeg",width=6.5,height=4)

ggplot(filter(sp_har_fig, SMU == "Nass"), aes(year, harv/1000)) +
  geom_ribbon(aes(xintercept = year, ymin = harv_LCI/1000, ymax = harv_UCI/1000), 
              width = 0, fill = "grey") +
  geom_line() +
  facet_wrap(~factor(CU, levels = nass_order), ncol=2, scales = "free_y") +
  theme_sleek() +
  labs(x = "Return year", y = "Harvest (000s) and 50th percentile")

ggsave("plots/nass/harvest-reconstructions.jpeg",width=8,height=5)

prop_age |>
  pivot_longer(4:6, names_to = "age") |>
  filter(SMU == "Nass",
         CU == "Upper Nass") |>
  mutate(age = gsub("a", "", age)) |>
  ggplot(aes(year, value, fill = age)) +
  geom_col(position = position_stack(reverse = TRUE), width = 1) +
  scale_fill_brewer(palette = "Dark2")+
  theme_sleek() +
  labs(x = "Return year", y = "Proportion-at-age", fill = "Age")
ggsave("plots/nass/agecomp-reconstructions.jpeg",width=6,height=4)

# SR fits ----
ggplot() +
  geom_abline(intercept = 0, slope = 1,col="dark grey") +
  geom_ribbon(data = filter(SR.preds, SMU == "Nass"), 
              aes(x = Spawn/1000, ymin = Rec_lwr/1000, ymax = Rec_upr/1000),
              fill = "grey80", alpha=0.5, linetype=2, colour="gray46") +
  geom_errorbar(data = filter(brood.all, SMU == "Nass"), 
                aes(x= S_med/1000, y = R_med/1000, ymin = R_lwr/1000, ymax = R_upr/1000),
                colour="grey", width=0, linewidth=0.3) +
  geom_errorbar(data = filter(brood.all, SMU == "Nass"), 
                aes(y = R_med/1000, xmin = S_lwr/1000, xmax = S_upr/1000),
                height=0, colour = "grey", linewidth = 0.3) +
  geom_point(data = filter(brood.all, SMU == "Nass"),
             aes(x = S_med/1000, y = R_med/1000, color=BroodYear), size = 1.5) +
  geom_line(data = filter(SR.preds, SMU == "Nass"), 
            aes(x = Spawn/1000, y = Rec_med/1000)) +
  facet_wrap(~factor(CU, levels = nass_order), scales = "free",ncol=2) +
  scale_colour_viridis_c(name = "Brood Year")+
  labs(x = "Spawners (000s)", y = "Recruits (000s)") +
  theme_sleek()
ggsave("plots/nass/SR-relationships.jpeg",width=8,height=5)

# spawn wqith benchamrks

bench_plot <- filter(bench.par.table, SMU == "Nass", CU != "Nass Aggregate", 
                     par %in% c("Sgen", "Smsy.80")) |>
  select(CU, par, median)

bench_plot_Tseax<-rbind(c("Tseax",'Sgen',1012),
                        c("Tseax",'Smsy.80',2985))
bench_plot_Tseax <- as.data.frame(bench_plot_Tseax)
colnames(bench_plot_Tseax) <- colnames(bench_plot)
bench_plot_all <- rbind(bench_plot,bench_plot_Tseax)
bench_plot_all$median <- as.numeric(bench_plot_all$median)

ggplot(filter(sp_fig, SMU == "Nass", CU!= "Nass Aggregate"), aes(year, spwn/1000)) +
  geom_line() +
  geom_ribbon(aes(ymin = spwn_LCI/1000, ymax = spwn_UCI/1000), fill = "darkgrey", alpha = 0.5) +
  geom_hline(data = bench_plot_all, 
             aes(yintercept = median/1000, col = par), lty = "dashed") +
  scale_colour_manual(values=c(Smsy.80 = "forestgreen", 
                               Sgen = "darkred"),
                      labels=c(Smsy.80 = expression(paste("80%", S[MSY])),
                               Sgen = expression(S[Gen]))) +
  facet_wrap(~factor(CU, levels = nass_order_spw), nrow = 2, scales = "free_y") +
  theme_sleek() +
  coord_cartesian(ylim = c(0,NA)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(x = "Retun Year", y = "Spawners (000s)", col = "")


ggsave("plots/nass/spawner-status.jpeg",width=8,height=5)

# time varying alpha estimates ---
data <- a.yrs.all |>
  filter(SMU == "Nass", 
         brood_year < max(a.yrs.all$brood_year)-a_max) 

# index fully observed rec
  ggplot(data) +
  geom_ribbon(aes(x = brood_year, ymin = lwr, ymax = upr), fill = "lightgrey") +
  geom_line(aes(x = brood_year , y = mid)) +
  geom_hline(yintercept = 1, lty=2, col = "grey") +
  labs(y ="Productivity (\U03B1) and 80th percentile", x = "Brood year")+
  facet_wrap(~factor(CU, levels = nass_order), ncol=2) +
  scale_y_continuous(limits = c(0,8)) +
  theme_sleek() +
  theme(legend.position = "bottom") 
  ggsave("plots/nass/TV-productivity.jpeg",width=8,height=5)
  
  
# time varying productivity for all CUs ----  
  a.yrs.all |>
    filter(brood_year < max(a.yrs.all$brood_year)-a_max, # index fully observed rec
           !(CU %in% c("Skeena Aggregate", "Nass Aggregate"))) |> 
    ggplot(aes(color=CU)) +
    geom_line(aes(x = brood_year , y = mid), lwd = 1.5) +
    scale_color_viridis_d() +
    geom_hline(yintercept = 1, lty=2, col = "grey") +
    labs(y ="Productivity (recruits-per-spawner)", x = "Brood year")+
    guides(color=guide_legend(title="Conservation Unit")) +
    scale_y_continuous(limits = c(0,8)) +
    annotate("text", x = 1990, y = 1.25, label = "replacement", color = "grey") +
    theme_sleek() +
    theme(legend.position = "right") 
  
ggsave("plots/time-varying-prod-all-CUs.jpeg",width=6.5,height=4)


# spawners and benchmarks plots all CUs ----

bench_plot <- filter(bench.par.table, !CU %in% c("Skeena Aggregate","Nass Aggregate"),
                     par %in% c("Sgen", "Smsy.80")) |>
  select(CU, par, median)

ggplot(filter(AR1.spwn, !CU %in% c("Skeena Aggregate","Nass Aggregate")), aes(year, S.50/1000)) +
  geom_line() +
  geom_ribbon(aes(ymin = S.25/1000, ymax = S.75/1000), fill = "darkgrey", alpha = 0.5) +
  geom_hline(data = bench_plot, 
             aes(yintercept = median/1000, col = par), lty = "dashed") +
  geom_vline(xintercept = max(AR1.spwn$year) - a_max, lwd  = 0.3, lty = 2) +
  scale_colour_manual(values=c(Smsy.80 = "forestgreen", 
                               Sgen = "darkred"),
                      labels=c(Smsy.80 = expression(paste("80%", S[MSY])),
                               Sgen = expression(S[Gen]))) +
  facet_wrap(~factor(CU, levels = all_order), scales = "free_y") +
  theme_sleek() +
  labs(x = "Retun Year", y = "Spawners (000s)", col = "")
ggsave("plots/spawn-benchmarks-all-cus.jpeg",width=7.5,height=4.5)

  