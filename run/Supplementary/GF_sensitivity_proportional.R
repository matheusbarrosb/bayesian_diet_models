library(here)
library(rstan)
library(dplyr)
library(ggplot2)
library(purrr)
library(readr)
library(stringr)
library(reshape2)
library(tidyr)
library(PNWColors)
library(ggpubr)
options(error = NULL)

# Source custom functions
fun_dir = here::here("R")
fun_files = list.files(fun_dir, pattern = "\\.R$", full.names = TRUE)
lapply(fun_files, source)

# Load data
data_dir = here::here("data")
raw_data = read.csv(file.path(data_dir, "rawData.csv"))

spps      = sort(c("LAGRHO", "MICUND", "BAICHR"))
sites     = sort(c("AM", "DR", "HWP", "LB", "NEPaP", "SA", "CI"))
RS_levels = sort(c("CT", "LS"))

wg_data = raw_data %>%
  filter(Species.code %in% spps) %>%
  filter(Site %in% sites) %>%
  filter(Treatment %in% RS_levels) %>%
  group_by(Fish_ID_year)

spp   = as.numeric(factor(wg_data$Species.code, levels = spps))
site  = as.numeric(factor(wg_data$Site, levels = sites))
RS    = as.numeric(factor(wg_data$Treatment, levels = RS_levels))
GW    = as.numeric(wg_data$Gut.weight)
FW    = as.numeric(wg_data$Wet.Weight)
TL    = as.numeric(wg_data$Length)
group = as.character(wg_data$Group)

df = data.frame(
  spp   = spp,
  site  = site,
  RS    = RS,
  GW    = GW,
  FW    = FW,
  TL    = TL,
  Group = group
) %>%
  filter(!is.na(spp) & !is.na(site) & !is.na(RS) & !is.na(FW) & !is.na(TL)) %>%
  filter(FW > 0) %>%
  mutate(
    Z  = ifelse(Group == "Empty", 1, 0),
    GF = ifelse(Z == 1, 0.5, GW / FW)
  ) %>%
  filter(!is.na(GF))

stan_data = list(
  N       = nrow(df),
  S       = length(spps),         
  I       = length(sites),         
  TT      = length(RS_levels),    
  species = df$spp,
  site    = df$site,              
  status  = df$RS,                
  GF      = df$GF,
  Z       = df$Z,
  TL      = df$TL
)

# Sensitivity grid for prior scale
prior_scales = seq(0.25, 2, length.out = 8) 

# Helper to change Stan model code (replace sigma prior scale)
update_prior_scale = function(stan_file, out_file, prior_scale) {
  stan_code = readLines(stan_file)
  stan_code = gsub("cauchy\\(0, 1\\)", sprintf("cauchy(0, %.2f)", prior_scale), stan_code)
  writeLines(stan_code, out_file)
  invisible(out_file)
}

# Directory setup
model_dir = here::here("stan")
model_template = file.path(model_dir, "GF_beta_reg.stan")
output_dir = here::here("res", "gf_prior_sensitivity")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Loop over prior scales
sigma_posteriors = list(sigma_site = list(), sigma_RS = list(), sigma_int = list())

for (scale_val in prior_scales) {
  message("Fitting model with sigma prior scale: ", scale_val)
  model_file = file.path(output_dir, sprintf("GF_beta_reg_scale_%.1f.stan", scale_val))
  update_prior_scale(model_template, model_file, scale_val)
  
  fit = stan(
    file    = model_file,
    data    = stan_data,
    chains  = 1,
    iter    = 500,
    warmup  = 100,
    cores   = 1,
    seed    = 427,
    control = list(adapt_delta = 0.95, max_treedepth = 20)
  )
  
  saveRDS(fit, file = file.path(output_dir, sprintf("fit_scale_%.1f.rds", scale_val)))
  
  sigma_posteriors$sigma_site[[as.character(scale_val)]] = as.matrix(fit, pars = "sigma_site")
  sigma_posteriors$sigma_RS[[as.character(scale_val)]]   = as.matrix(fit, pars = "sigma_RS")
  sigma_posteriors$sigma_int[[as.character(scale_val)]]  = as.matrix(fit, pars = "sigma_int")
  
  rm(fit)
  gc()
}

# Prepare posteriors for plotting
plot_df_list = list()

for (param_name in names(sigma_posteriors)) {
  for (scale_val in names(sigma_posteriors[[param_name]])) {
    posterior = sigma_posteriors[[param_name]][[scale_val]] # Matrix of size [n_iter, S]
    
    df = as.data.frame(posterior)
    colnames(df) = spps # Assign species names to columns
    
    df_long = df %>%
      mutate(draw = 1:nrow(df), scale = as.numeric(scale_val), parameter = param_name) %>%
      pivot_longer(cols = all_of(spps), names_to = "species", values_to = "value")
    
    plot_df_list[[length(plot_df_list) + 1]] = df_long
  }
}

plot_df = bind_rows(plot_df_list)

# Clean up parameter labels for plotting
plot_df$parameter = factor(plot_df$parameter,
                           levels = c("sigma_site", "sigma_RS", "sigma_int"),
                           labels = c(expression(sigma[site]), expression(sigma[status]), expression(sigma[site %*% status]))
)

# Add descriptive species labels
spp_labs = c(LAGRHO = "Pinfish", MICUND = "Croaker", BAICHR = "Silver perch")
plot_df$species = recode(plot_df$species, !!!spp_labs)
plot_df$species = factor(plot_df$species, levels = c("Pinfish", "Croaker", "Silver perch"))

# Plot density for all sigma terms, faceted by parameter and species
ggplot(plot_df, aes(x = value, fill = factor(scale))) +
  geom_density(alpha = 0.3) +
  facet_grid(species ~ parameter, scales = "free") +
  labs(
    x = "Posterior standard deviation",
    fill = "Prior scale"
  ) +
  custom_theme() +
  theme(legend.position = "top", strip.text = element_text(size = 10, face = "bold")) +
  xlim(0, 3.5)

supp_fig_dir = file.path(here::here(), "res", "figures", "supplementary/")
dir.create(supp_fig_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = file.path(supp_fig_dir, "GF_sigma_sensitivity.pdf"),
  width = 7,
  height = 6,
  units = "in"
)