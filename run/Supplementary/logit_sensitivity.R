# Prior sensitivity analysis for logit model: vary SD for alpha and beta terms (all species, joint plot)
library(here)
library(rstan)
library(dplyr)
library(ggplot2)
library(purrr)
library(stringr)
library(tidyr)
library(ggmcmc)
library(reshape2)
options(error = NULL)

# Source custom functions
fun_dir = here::here("R")
fun_files = list.files(fun_dir, pattern = "\\.R$", full.names = TRUE)
lapply(fun_files, source)

# Species definitions
species_list = list(
  LAGRHO = list(
    sites = sort(c("AM", "DR", "HWP", "LB", "NEPaP", "SA")),
    prey = sort(c("Amphipod", "Crustacean", "Fish", "Isopod", "Polychaete", "SAV", "Tanaidacea"))
  ),
  MICUND = list(
    sites = sort(c("CI", "NEPaP", "SA")),
    prey = sort(c("Amphipod", "Crustacean", "Fish", "Isopod", "Polychaete", "Tanaidacea"))
  ),
  BAICHR = list(
    sites = sort(c("CI", "NEPaP", "SA")),
    prey = sort(c("Amphipod", "Crustacean", "Fish", "Isopod", "Polychaete", "Tanaidacea"))
  )
)
RS_levels = sort(c("CT", "LS"))

# SD grid for priors (skip 1)
prior_sds = c(seq(0.5, 0.75, by = 0.25), seq(1.25, 2, by = 0.25))

# Helper to change Stan model code for prior SDs
update_prior_sd = function(stan_file, out_file, sd_val) {
  stan_code = readLines(stan_file)
  stan_code = gsub("normal\\(0, 1\\)", sprintf("normal(0, %.2f)", sd_val), stan_code)
  writeLines(stan_code, out_file)
  invisible(out_file)
}

# Directory setup
model_dir = here::here("stan")
model_template = file.path(model_dir, "logit_int.stan")
output_dir = here::here("res", "logit_prior_sensitivity")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

results = list()
species_names = names(species_list)

for (sp in species_names) {
  message("Processing species: ", sp)
  
  # Load data
  data_dir = here::here("data")
  raw_data = read.csv(file.path(data_dir, "rawData.csv"))
  sites = species_list[[sp]]$sites
  prey = species_list[[sp]]$prey
  spps = sp
  
  wg_data =
    raw_data %>%
    filter(Species.code %in% spps) %>%
    filter(Site %in% sites) %>%
    filter(Treatment %in% RS_levels) %>%
    filter(Group %in% prey) %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = Group, values_from = present, values_fill = list(present = 0))
  
  prey_mat = 
    wg_data %>%
    select(Species.code, Site, Treatment, all_of(prey)) %>%
    mutate(Species.code = as.numeric(factor(Species.code, levels = spps))) %>%
    mutate(Site = as.numeric(factor(Site, levels = sites))) %>%
    mutate(Treatment = as.numeric(factor(Treatment, levels = RS_levels))) %>%
    select(-Species.code, -Site, -Treatment) %>%
    as.matrix()
  
  site   = as.numeric(as.factor(sort(wg_data$Site)))
  status = as.numeric(as.factor(sort(wg_data$Treatment)))
  TL     = wg_data$Length; TL[is.na(TL)] = mean(TL, na.rm = TRUE)
  
  stan_data = list(
    N = nrow(prey_mat),
    S = length(unique(sites)),
    K = length(unique(RS_levels)),
    G = length(unique(prey)),
    site     = site,
    status   = status,
    TL       = TL,
    prey_mat = prey_mat,
    meanTLs_site = tapply(TL, site, mean),
    meanTLs_status = tapply(TL, status, mean)
  )
  
  # Loop over prior SDs
  for (sd_val in prior_sds) {
    message("  Fitting model with prior SD: ", sd_val)
    model_file = file.path(output_dir, sprintf("logit_int_sd_%.2f.stan", sd_val))
    update_prior_sd(model_template, model_file, sd_val)
    
    fit = stan(
      file = model_file,
      data = stan_data,
      chains = 3,
      iter = 4000,
      warmup = 1000,
      cores = 3,
      seed = 444,
      control = list(adapt_delta = 0.95, max_treedepth = 15)
    )
    saveRDS(fit, file = file.path(output_dir, sprintf("fit_%s_sd_%.2f.rds", sp, sd_val)))
    results[[paste(sp, sd_val, sep = "_")]] <- list(fit = fit, species = sp, sd_val = sd_val, prey = prey)
  }
}

param_names = c("alpha", "beta_status", "beta_site_status", "beta_TL")
plot_df_list = list()

for (res_name in names(results)) {
  fit = results[[res_name]]$fit
  sp = results[[res_name]]$species
  sd_val = results[[res_name]]$sd_val
  prey = results[[res_name]]$prey
  
  for (param in param_names) {
    post = as.matrix(fit, pars = param)
    post_df = as.data.frame(post)
    colnames(post_df) = paste0(param, "[", seq_len(ncol(post_df)), "]")
    post_df_long = post_df %>%
      mutate(draw = 1:n()) %>%
      pivot_longer(-draw, names_to = "parameter", values_to = "value") %>%
      mutate(
        prior_sd = as.numeric(sd_val),
        family = param,
        species = sp,
        prey_idx = as.integer(str_match(parameter, "\\[(\\d+)\\]")[,2]),
        Prey = ifelse(is.na(prey_idx), NA, prey[prey_idx])
      )
    plot_df_list[[paste(res_name, param, sep = "_")]] <- post_df_long
  }
}

plot_df = bind_rows(plot_df_list) %>%
  mutate(family = as.character(family))

# Joint density plot: alpha and beta_TL
plot1 = 
plot_df %>%
  filter(species== "LAGRHO") %>%
  group_by(family, Prey, prior_sd) %>%
  summarize(
    mean_value = mean(value),
    sd_value = sd(value),
    .groups = "drop"
  ) %>%
  na.exclude() %>%
  ggplot(aes(x = Prey, y = mean_value, fill = factor(prior_sd))) +
  geom_point(position = position_dodge(width = 0.7), size = 2, shape = 21, color = "black") +
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), position = position_dodge(width = 0.5), width = 0.2) +
  facet_wrap(~ family, scales = "free") +
  xlab("") +
  ylab("") +
  coord_flip() +
  custom_theme() +
  ggtitle("Pinfish") +
  theme(legend.tile = element_blank())


plot2 =   
plot_df %>%
    filter(species== "MICUND") %>%
    group_by(family, Prey, prior_sd) %>%
    summarize(
      mean_value = mean(value),
      sd_value = sd(value),
      .groups = "drop"
    ) %>%
    na.exclude() %>%
  ggplot(aes(x = Prey, y = mean_value, fill = factor(prior_sd))) +
    geom_point(position = position_dodge(width = 0.6), size = 2, shape = 21, color = "black") +
    geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), position = position_dodge(width = 0.5), width = 0.2) +
    facet_wrap(~ family, scales = "free") +
    xlab("") +
    ylab("") +
    coord_flip() +
    custom_theme()+
  ggtitle("Croaker") +
  theme(legend.tile = element_blank())


plot3 = 
plot_df %>%
  filter(species== "BAICHR") %>%
  group_by(family, Prey, prior_sd) %>%
  summarize(
    mean_value = mean(value),
    sd_value = sd(value),
    .groups = "drop"
  ) %>%
  na.exclude() %>%
  ggplot(aes(x = Prey, y = mean_value, fill = factor(prior_sd))) +
  xlab("") +
  geom_point(position = position_dodge(width = 0.6), size = 2, shape = 21, color = "black") +
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), position = position_dodge(width = 0.5), width = 0.2) +
  facet_wrap(~ family, scales = "free") +
  coord_flip() +
  ggtitle("Silver perch") +
  custom_theme()
  
mainplot = ggarrange(
  plot1, plot2, plot3,
  ncol = 1,
  common.legend = TRUE,
  legend = "right"
  ); print(mainplot)

fig_dir = here::here("res", "figures", "supplementary")

ggsave(
  filename = file.path(fig_dir, "logit_sens.pdf"),
  width = 6, height = 15, units = "in"
)

# Optional: Save raw posterior data for further custom plots
saveRDS(plot_df, file = file.path(output_dir, "logit_sens.rds"))