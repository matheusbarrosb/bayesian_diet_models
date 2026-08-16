#rm(list = ls())
set.seed(444)
# Packages ---------------------------------------------------------------------
pack_list = c("dplyr", "rstan", "tidyr", "ggplot2", "purrr", "readr",
              "stringr", "here", "reshape2", "PNWColors", "ggpubr", "ggmcmc", "binom")
if (!all(pack_list %in% rownames(installed.packages()))) {
  install.packages(pack_list[!pack_list %in% rownames(installed.packages())])
}
lapply(pack_list, library, character.only = TRUE)

# Source functions -------------------------------------------------------------
fun_dir = here::here("R")
fun_files = list.files(fun_dir, pattern = "\\.R$", full.names = TRUE)
lapply(fun_files, source)

# Load data --------------------------------------------------------------------
data_dir = here::here("data")
raw_data = read.csv(file.path(data_dir, "rawData.csv"))

# Data wrangling ---------------------------------------------------------------
spps  = sort(c("LAGRHO"))
sites = sort(c("AM", "DR", "HWP", "LB", "NEPaP", "SA"))
prey  = sort(c("Amphipod", "Crustacean", "Fish", "Isopod", "Polychaete", "SAV", "Tanaidacea"))
RS_levels = sort(c("CT", "LS"))

wg_data =
  raw_data %>%
  filter(Species.code %in% spps) %>%
  filter(Site %in% sites) %>%
  filter(Treatment %in% RS_levels) %>%
  filter(Group %in% prey) %>%
  mutate(present = 1) %>%
  pivot_wider(
    id_cols = c(Fish.ID, Fish_ID_year, Species.code, Site, Treatment, Length, Wet.Weight, Gut.weight),
    names_from = Group,
    values_from = present,
    values_fn = ~1,
    values_fill = 0
  )

# Gather data for stan ---------------------------------------------------------
prey_mat = 
wg_data %>%
  select(Species.code, Site, Treatment, Amphipod, Crustacean, Fish, Isopod, Polychaete, SAV, Tanaidacea) %>%
  mutate(Species.code = as.numeric(factor(Species.code, levels = spps))) %>%
  mutate(Site = as.numeric(factor(Site, levels = sites))) %>%
  mutate(Treatment = as.numeric(factor(Treatment, levels = RS_levels))) %>%
  select(-Species.code, -Site, -Treatment) %>%
  as.matrix()

site   = as.numeric(factor(wg_data$Site, levels = sites))
status = as.numeric(factor(wg_data$Treatment, levels = RS_levels))
spp    = as.numeric(factor(wg_data$Species.code, levels = spps))
TL     = wg_data$Length; TL[is.na(TL)] = mean(TL, na.rm = TRUE)
meanTL = mean(TL, na.rm = TRUE)

stan_data = list(
  N = nrow(prey_mat),
  S = length(unique(sites)),
  K = length(unique(RS_levels)),
  G = length(unique(prey)),
  site     = site,
  status   = status,
  TL       = TL,
  prey_mat = prey_mat,
  meanTL   = meanTL
)


# fit the model ----------------------------------------------------------------
model_dir = here::here("stan/")
stanc(paste0(model_dir, "logit_int.stan"))
rstan_options(auto_write = TRUE)

fit_pinfish = stan(
  file   = paste0(model_dir, "logit_int.stan"),  
  data   = stan_data,
  chains = 3,        
  iter   = 10000,      
  warmup = 1000,       
  cores  = 3,        
  seed   = 444       
)

pinfish_diff = get_FO_diff(fit_pinfish, prey, RS_levels)
pinfish_diff$Species = "Pinfish"

# Plotting ---------------------------------------------------------------------

### 1. Probabilities of occurrence ###
post = rstan::extract(fit_pinfish, pars = "p_status")$p_status

df =
ggmcmc::ggs(fit_pinfish, family = "p_status")

prey_names = prey
status_names = c("C", "R")

df_parsed = df %>%
  mutate(
    indices = str_extract(Parameter, "\\[(.*?)\\]"),
    indices = str_remove_all(indices, "\\[|\\]"),
    prey_idx   = as.integer(str_split_fixed(indices, ",", 2)[,1]),
    status_idx = as.integer(str_split_fixed(indices, ",", 2)[,2]),
    Prey   = prey_names[prey_idx],
    Status = status_names[status_idx]
  ) %>%
  select(-indices, -prey_idx, -status_idx)

head(df_parsed)

pinfish_probs = 
df_parsed %>%
  group_by(Prey, Status) %>%
  summarise(
    mean = mean(value),
    lower = quantile(value, 0.1),
    upper = quantile(value, 0.9)
  ) %>%
  
  ggplot(aes(x = Prey, y = mean, color = Status)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0,
                position = position_dodge(width = 0.5)) +
  custom_theme() +
  coord_flip() +
  ylab("") +
  xlab("") +
  ggtitle("") +
  theme(legend.title = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  scale_color_manual(values = pnw_palette("Bay", 2)) 

## 2. beta_tl plot ##
post = rstan::extract(fit_pinfish, pars = "beta_TL")$beta_TL

df =
ggmcmc::ggs(fit_pinfish, family = "beta_TL") %>%
  mutate(
    # Extract indices from Parameter: beta_TL[prey]
    indices = str_extract(Parameter, "\\[(.*?)\\]"),
    indices = str_remove_all(indices, "\\[|\\]"),
    prey_idx   = as.integer(str_split_fixed(indices, ",", 1)[,1]),
    Prey   = prey_names[prey_idx]
  ) %>%
  select(-indices, -prey_idx) %>%
  
  group_by(Prey) %>%
  summarise(
    mean = mean(value),
    lower = quantile(value, 0.1),
    upper = quantile(value, 0.9)
  )

pinfish_beta_TL =
df %>%
  ggplot(aes(x = Prey, y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  custom_theme() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  coord_flip() +
  xlab("") +
  ggtitle("Pinfish") +
  ylab("")

# P(LS > CT) ------------------------------------------
pinfish_diff = pinfish_diff %>%
  mutate(sig = sign(Diff_lower) == sign(Diff_upper))

pinfish_FO_diff =
  pinfish_diff %>%
  ggplot(aes(x = Prey, y = Diff_mean)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point() +
  geom_errorbar(aes(ymin = Diff_lower, ymax = Diff_upper), width = 0) +
  geom_text(data = filter(pinfish_diff, sig),
            aes(y = Diff_upper + 0.12*sign(Diff_upper), label = paste0("P=", round(P_LS_greater, 2))),
            size = 2.8) +
  custom_theme() +
  coord_flip() +
  ggtitle("") +
  xlab("") +
  ylab("") +
  ylim(-0.5,0.8) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())

pinfish_logit_plot = 
  ggarrange(
    pinfish_beta_TL,
    pinfish_probs,
    pinfish_FO_diff,
    ncol = 3,
    widths = c(1, 0.75, 0.75)
  );pinfish_logit_plot

## 3. Posterior predictive check ##
obs_df = as.data.frame(prey_mat)
obs_df$site = site

obs_freq = obs_df %>%
  group_by(site) %>%
  summarise(
    across(everything(), list(
      Observed = ~mean(.), 
      Successes = ~sum(.), 
      Trials = ~n()
    )),
    .groups = "drop"
  ) %>%
  # Reshape to long format
  pivot_longer(
    -site, 
    names_to = c("Prey", ".value"), 
    names_pattern = "(.*)_(Observed|Successes|Trials)"
  ) %>%
  mutate(Site = sites[site]) %>%
  mutate(
    ci = binom::binom.confint(x = Successes, n = Trials, conf.level = 0.95, methods = "wilson")
  ) %>%
  mutate(
    Obs_Lower = ci$lower,
    Obs_Upper = ci$upper
  ) %>%
  select(Site, Prey, Observed, Obs_Lower, Obs_Upper)

post_pred = rstan::extract(fit_pinfish, pars = "p_site")$p_site  

post_pred_df = reshape2::melt(post_pred)
colnames(post_pred_df) = c("Iteration", "Prey_idx", "Site_idx", "Modelled")
post_pred_df$Prey = prey[post_pred_df$Prey_idx]
post_pred_df$Site = sites[post_pred_df$Site_idx]

ppc_summary = post_pred_df %>%
  group_by(Prey, Site) %>%
  summarise(
    Model_Mean  = mean(Modelled),
    Model_Lower = quantile(Modelled, 0.025),
    Model_Upper = quantile(Modelled, 0.975),
    .groups = "drop"
  )

pinfish_ppcheck =
  obs_freq %>%
  left_join(ppc_summary, by = c("Prey", "Site")) %>%
  ggplot(aes(x = Observed, y = Model_Mean, color = Prey)) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = Model_Lower, ymax = Model_Upper), width = 0.0, alpha = 0.6) +
  geom_errorbarh(aes(xmin = Obs_Lower, xmax = Obs_Upper), height = 0.0, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  labs(
    x = "Observed FO",
    y = "Predicted FO",
    title = "Pinfish"
  ) +
  theme(legend.title = element_blank()) +
  custom_theme() +
  xlim(0, 1) + ylim(0, 1); print(pinfish_ppcheck)


