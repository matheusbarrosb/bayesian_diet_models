# Prior Sensitivity Analysis for tau hyperpriors in Stan model
library(here)
library(rstan)
library(dplyr)
library(ggplot2)
library(purrr)
library(readr)
library(stringr)
library(reshape2)
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

spps  = sort(c("LAGRHO", "MICUND", "BAICHR"))
sites = sort(c("AM", "DR", "HWP", "LB", "NEPaP", "SA", "CI"))
RS_levels = sort(c("CT", "LS"))

wg_data = raw_data %>%
  filter(Species.code %in% spps) %>%
  filter(Site %in% sites) %>%
  filter(Treatment %in% RS_levels) %>%
  group_by(Fish_ID_year)

spp  = as.numeric(factor(wg_data$Species.code, levels = spps))
site = as.numeric(factor(wg_data$Site, levels = sites))
RS   = as.numeric(factor(wg_data$Treatment, levels = RS_levels))
GW   = as.numeric(wg_data$Gut.weight)
FW   = as.numeric(wg_data$Wet.Weight)
TL   = as.numeric(wg_data$Length)

df = data.frame(
  spp  = spp,
  site = site,
  RS   = RS,
  GW   = GW,
  FW   = FW,
  TL   = TL
) %>%
  filter(!is.na(spp) & !is.na(site) & !is.na(RS) & !is.na(GW) & !is.na(FW) & !is.na(TL))

stan_data = list(
  N       = nrow(df),
  S       = length(spps),         
  I       = length(sites),        
  TT      = length(RS_levels),    
  species = df$spp,               
  site    = df$site,             
  status  = df$RS,                
  GW      = df$GW,
  FW      = df$FW,
  TL      = df$TL
)

# Sensitivity grid for prior scale
tau_scales = seq(1, 3.5, length.out = 6) # 1, 1.5, 2, 2.5, 3, 3.5

# Helper to change Stan model code (replace tau prior scale)
update_tau_scale <- function(stan_file, out_file, tau_scale) {
  stan_code <- readLines(stan_file)
  # Replace all cauchy(0, 2.5) --> cauchy(0, tau_scale)
  stan_code <- gsub("cauchy\\(0, 2.5\\)", sprintf("cauchy(0, %.2f)", tau_scale), stan_code)
  writeLines(stan_code, out_file)
  invisible(out_file)
}

# Directory setup
model_dir = here::here("stan")
model_template = file.path(model_dir, "GF_centered.stan")
output_dir = here::here("res", "gf_prior_sensitivity")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Loop over tau scales
results <- list()
for (tau_scale in tau_scales) {
  message("Fitting model with tau prior scale: ", tau_scale)
  model_file <- file.path(output_dir, sprintf("GF_centered_tau_%.1f.stan", tau_scale))
  update_tau_scale(model_template, model_file, tau_scale)
  
  fit <- stan(
    file = model_file,
    data = stan_data,
    chains = 3,
    iter = 5000,
    warmup = 1000,
    cores = 3,
    seed = 427,
    control = list(adapt_delta = 0.95, max_treedepth = 20)
  )
  
  saveRDS(fit, file = file.path(output_dir, sprintf("fit_tau_%.1f.rds", tau_scale)))
  results[[as.character(tau_scale)]] <- fit
}

tau_posteriors <- list(
  tau_alpha    = map(results, ~ as.matrix(.x, pars = "tau_alpha")),
  tau_TL       = map(results, ~ as.matrix(.x, pars = "tau_TL")),
  tau_FW       = map(results, ~ as.matrix(.x, pars = "tau_FW")),
  tau_site     = map(results, ~ as.matrix(.x, pars = "tau_site")),
  tau_RS       = map(results, ~ as.matrix(.x, pars = "tau_RS")),
  tau_site_RS  = map(results, ~ as.matrix(.x, pars = "tau_site_RS"))
)

reshape_tau_posterior <- function(posterior, scale, tau_name, S = NULL, I = NULL) {
  ndraws <- nrow(posterior)
  if (is.null(S)) {
    # scalar tau
    data.frame(value = posterior[,1], scale = scale, tau = tau_name, element = "")
  } else if (tau_name %in% c("tau_site", "tau_RS")) {
    # vector tau
    df <- as.data.frame(posterior)
    colnames(df) <- paste0("sp", seq_len(S))
    df_long <- df %>%
      mutate(draw = 1:ndraws, scale = scale, tau = tau_name) %>%
      pivot_longer(cols = starts_with("sp"), names_to = "element", values_to = "value")
    df_long
  } else if (tau_name == "tau_site_RS") {
    # matrix tau: n_draws x (S*(I-1))
    df <- as.data.frame(posterior)
    colnames(df) <- paste0("sp", rep(seq_len(S), each = (I - 1)), "_i", rep(seq_len(I - 1), times = S))
    df$draw <- 1:ndraws
    df$scale <- scale
    df$tau <- tau_name
    df_long <- df %>%
      pivot_longer(cols = starts_with("sp"), names_to = "element", values_to = "value")
    df_long
  }
}

# Prepare posteriors for plotting
S <- length(spps)
I <- length(sites)

plot_df_list <- list()

for (tau_name in names(tau_posteriors)) {
  for (scale in names(tau_posteriors[[tau_name]])) {
    posterior <- tau_posteriors[[tau_name]][[scale]]
    ndraws <- nrow(posterior)
    if (tau_name %in% c("tau_alpha", "tau_TL", "tau_FW")) {
      # scalar
      plot_df_list[[paste0(tau_name,"_",scale)]] <- data.frame(
        value = posterior[,1],
        scale = as.numeric(scale),
        tau = tau_name,
        element = ""
      )
    } else if (tau_name %in% c("tau_site", "tau_RS")) {
      plot_df_list[[paste0(tau_name,"_",scale)]] <- reshape_tau_posterior(posterior, as.numeric(scale), tau_name, S = S)
    } else if (tau_name == "tau_site_RS") {
      plot_df_list[[paste0(tau_name,"_",scale)]] <- reshape_tau_posterior(posterior, as.numeric(scale), tau_name, S = S, I = I)
    }
  }
}

plot_df <- bind_rows(plot_df_list)

# Clean up tau labels for plotting
plot_df$tau <- factor(plot_df$tau,
                      levels = c("tau_alpha", "tau_TL", "tau_FW", "tau_site", "tau_RS", "tau_site_RS"),
                      labels = c(expression(tau[alpha]), expression(tau[TL]), expression(tau[FW]),
                                 expression(tau[site]), expression(tau[RS]), expression(tau[site*RS]))
)

# Plot density for all tau terms
ggplot(plot_df, aes(x = value, fill = factor(scale))) +
  geom_density(alpha = 0.3) +
  facet_wrap(~ tau, scales = "free", ncol = 2,
             labeller = label_parsed) +
  labs(
    x = "Posterior value",
    fill = "Prior scale"
  ) +
  custom_theme() +
  theme(legend.position = "right") +
  xlim(0,10)

supp_fig_dir = file.path(here::here(), "res", "figures", "supplementary/")
ggsave(
  filename = file.path(supp_fig_dir, "GF_tau_sensitivity.pdf"),
  width = 5,
  height = 6,
  units = "in"
)



