library(ggmcmc)
library(ggplot2)
library(dplyr)
library(stringr)

data_dir = here::here("data")
raw_data = read.csv(file.path(data_dir, "rawData.csv"))

# plot size histograms ---------------------------------------------------------

spps  = sort(c("LAGRHO", "MICUND", "BAICHR"))
sites = sort(c("AM", "DR", "HWP", "LB", "NEPaP", "SA"))
RS_levels = sort(c("CT", "LS"))

wg_data =
  raw_data %>%
  filter(Species.code %in% spps) %>%
  filter(Site %in% sites) %>%
  filter(Treatment %in% RS_levels) %>%
  filter(Group %in% prey) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = Group, values_from = present, values_fill = list(present = 0))

wg_data %>%
  # mutate names so LAGRHO = Pinfish, etc
  mutate(Species.code = recode(Species.code,
                               LAGRHO = "Pinfish",
                               MICUND = "Croaker",
                               BAICHR = "Silver perch")) %>%
  
  ggplot(aes(x = Length)) +
  geom_histogram(color = "black", fill = "white") +
  facet_wrap(~Species.code + Site) +
  custom_theme()

sup_fig_dir = here::here("res", "figures", "supplementary", "Figures/")


ggsave(file.path(sup_fig_dir, "length_histograms.png"),
       width = 4, height = 5)


#----- Logit MCMC traceplots ---------------------------------------------------

##------ 1. Pinfish --------
ggs_obj <- ggs(fit_pinfish)
ggs_main <- ggs_obj %>%
  filter(
    str_detect(Parameter, "^alpha(\\[|$)") |
      str_detect(Parameter, "^beta_site(\\[|$)") |
      str_detect(Parameter, "^beta_status(\\[|$)") |
      str_detect(Parameter, "^beta_TL(\\[|$)")
  )


main_families = c("alpha", "beta_site", "beta_status", "beta_TL")

# Set number of columns for each parameter family
# Adjust these based on how many parameters you have in each family
ncol_settings = list(
  alpha       = 2,        # For G prey groups
  beta_site   = 3,    # For G x S combinations  
  beta_status = 2,  # For G x K combinations
  beta_TL     = 2       # For G prey groups
)

p = list()
for (fam in main_families) {
  ggs_fam = ggs_main %>%
    filter(str_detect(Parameter, paste0("^", fam, "(\\[|$)")))
  
  if (nrow(ggs_fam) > 0) {
    ncol_fam = ncol_settings[[fam]]
    p[[fam]] = ggs_fam %>%
      ggplot(aes(x = Iteration, y = value, color = as.factor(Chain))) +
      geom_line(alpha = 0.8) +
      facet_wrap(~ Parameter, scales = "free_y", ncol = ncol_fam) +
      labs(
        title = paste("Traceplots for", fam),
        x = "Iteration",
        y = "Parameter Value",
        color = "Chain"
      ) +
      theme_bw() +
      theme(
        strip.text = element_text(size = 8),
        legend.position = "bottom"
      )
    #print(p[[fam]])
    
    ggsave(paste0(sup_fig_dir, "traceplot_", fam, ".png"), p[[fam]], 
            width = 4 * ncol_fam, height = ceiling(length(unique(ggs_fam$Parameter)) / ncol_fam) * 3, 
            dpi = 300)
  }
}

## ----- 2. Croaker --------
ggs_obj = ggs(fit_micund)

ggs_main = ggs_obj %>%
  filter(
    str_detect(Parameter, "^alpha(\\[|$)") |
      str_detect(Parameter, "^beta_site(\\[|$)") |
      str_detect(Parameter, "^beta_status(\\[|$)") |
      str_detect(Parameter, "^beta_TL(\\[|$)")
  )

main_families = c("alpha", "beta_site", "beta_status", "beta_TL")

ncol_settings = list(
  alpha       = 2,        # For G prey groups
  beta_site   = 3,    # For G x S combinations  
  beta_status = 2,  # For G x K combinations
  beta_TL     = 2       # For G prey groups
)

p = list()
for (fam in main_families) {
  ggs_fam = ggs_main %>%
    filter(str_detect(Parameter, paste0("^", fam, "(\\[|$)")))
  
  if (nrow(ggs_fam) > 0) {
    ncol_fam = ncol_settings[[fam]]
    p[[fam]] = ggs_fam %>%
      ggplot(aes(x = Iteration, y = value, color = as.factor(Chain))) +
      geom_line(alpha = 0.8) +
      facet_wrap(~ Parameter, scales = "free_y", ncol = ncol_fam) +
      labs(
        title = paste("Traceplots for", fam),
        x = "Iteration",
        y = "Parameter Value",
        color = "Chain"
      ) +
      theme_bw() +
      theme(
        strip.text = element_text(size = 8),
        legend.position = "bottom"
      )
    #print(p[[fam]])
    
    ggsave(paste0(sup_fig_dir, "traceplot_", fam, ".png"), p[[fam]], 
           width = 4 * ncol_fam, height = ceiling(length(unique(ggs_fam$Parameter)) / ncol_fam) * 3, 
           dpi = 300)
  }
}

## ------ 3. Silver perch --------
ggs_obj = ggs(fit_silverperch)
ggs_main = ggs_obj %>%
  filter(
    str_detect(Parameter, "^alpha(\\[|$)") |
      str_detect(Parameter, "^beta_site(\\[|$)") |
      str_detect(Parameter, "^beta_status(\\[|$)") |
      str_detect(Parameter, "^beta_TL(\\[|$)")
  )

main_families = c("alpha", "beta_site", "beta_status", "beta_TL")

ncol_settings = list(
  alpha       = 2,        # For G prey groups
  beta_site   = 3,    # For G x S combinations  
  beta_status = 2,  # For G x K combinations
  beta_TL     = 2       # For G prey groups
)

p = list()
for (fam in main_families) {
  ggs_fam = ggs_main %>%
    filter(str_detect(Parameter, paste0("^", fam, "(\\[|$)")))
  
  if (nrow(ggs_fam) > 0) {
    ncol_fam = ncol_settings[[fam]]
    p[[fam]] = ggs_fam %>%
      ggplot(aes(x = Iteration, y = value, color = as.factor(Chain))) +
      geom_line(alpha = 0.8) +
      facet_wrap(~ Parameter, scales = "free_y", ncol = ncol_fam) +
      labs(
        title = paste("Traceplots for", fam),
        x = "Iteration",
        y = "Parameter Value",
        color = "Chain"
      ) +
      theme_bw() +
      theme(
        strip.text = element_text(size = 8),
        legend.position = "bottom"
      )
    #print(p[[fam]])
    
    ggsave(paste0(sup_fig_dir, "traceplot_", fam, ".png"), p[[fam]], 
           width = 4 * ncol_fam, height = ceiling(length(unique(ggs_fam$Parameter)) / ncol_fam) * 3, 
           dpi = 300)
  }
}

# Gut fullness mcmc traceplots ---------------
# fit = readRDS(here::here("res", "fit_beta_reg.rds"))
# ggs_obj = ggmcmc::ggs(fit)
# 
# ggs_main = ggs_obj %>%
#   filter(
#     str_detect(Parameter, "^alpha(\\[|$)") |
#       str_detect(Parameter, "^beta_site(\\[|$)") |
#       str_detect(Parameter, "^beta_RS(\\[|$)") |
#       str_detect(Parameter, "^beta_TL(\\[|$)") |
#       str_detect(Parameter, "^beta_FW(\\[|$)") |
#       str_detect(Parameter, "^sigma(\\[|$)")
#   )
# 
# main_families = c("alpha", "beta_site", "beta_RS", "beta_TL", "beta_FW", "sigma")
# 
# ncol_settings = list(
#   alpha = 3,        # S species (3 parameters)
#   beta_site = 3,    # S x (I-1) = 3 x 6 = 18 parameters
#   beta_RS = 3,      # S x (TT-1) = 3 x 1 = 3 parameters
#   beta_TL = 3,      # S species (3 parameters)
#   beta_FW = 3,      # S species (3 parameters)
#   sigma = 3         # S species (3 parameters)
# )
# 
# p = list()
# for (fam in main_families) {
#   ggs_fam = ggs_main %>%
#     filter(str_detect(Parameter, paste0("^", fam, "(\\[|$)")))
#   
#   if (nrow(ggs_fam) > 0) {
#     # Get number of columns for this family
#     ncol_fam = ncol_settings[[fam]]
#     
#     p[[fam]] = ggs_fam %>%
#       ggplot(aes(x = Iteration, y = value, color = as.factor(Chain))) +
#       geom_line(alpha = 0.8) +
#       facet_wrap(~ Parameter, scales = "free_y", ncol = ncol_fam) +
#       labs(
#         title = paste("Traceplots for", fam),
#         x = "Iteration",
#         y = "Parameter Value",
#         color = "Chain"
#       ) +
#       theme_bw() +
#       theme(
#         strip.text = element_text(size = 8),
#         legend.position = "bottom"
#       )
#     
#     print(p[[fam]])
#     
#     n_params = length(unique(ggs_fam$Parameter))
#     plot_width = min(16, 4 * ncol_fam)
#     plot_height = max(6, ceiling(n_params / ncol_fam) * 2.5)
# 
#     ggsave(paste0(sup_fig_dir, "traceplot_", fam, ".png"), p[[fam]],
#            width = plot_width, height = plot_height, dpi = 300)
#   }
# }
# 
# # ---- Hyperprior traceplots (optional) ----
# ggs_hyperpriors = ggs_obj %>%
#   filter(
#     str_detect(Parameter, "^tau_") |
#       str_detect(Parameter, "lp__")  # log posterior
#   )
# 
# if (nrow(ggs_hyperpriors) > 0) {
#   traceplot_hyperpriors = ggs_hyperpriors %>%
#     ggplot(aes(x = Iteration, y = value, color = as.factor(Chain))) +
#     geom_line(alpha = 0.8) +
#     facet_wrap(~ Parameter, scales = "free_y", ncol = 3) +
#     labs(
#       title = "MCMC Traceplots: Hyperpriors & Log Posterior",
#       x = "Iteration",
#       y = "Parameter Value",
#       color = "Chain"
#     ) +
#     theme_bw() +
#     theme(
#       strip.text = element_text(size = 8),
#       legend.position = "bottom"
#     )
#   
#   print(traceplot_hyperpriors)
#   
#   # ggsave("mcmc_traceplots_hyperpriors.png", traceplot_hyperpriors, 
#   #        width = 12, height = 8, dpi = 300)
# }
# 
# # ---- Interaction terms traceplots (optional - can be many parameters) ----
# # ggs_interactions = ggs_obj %>%
# #   filter(str_detect(Parameter, "^beta_site_RS(\\[|$)"))
# # 
# # if (nrow(ggs_interactions) > 0) {
# #   traceplot_interactions = ggs_interactions %>%
# #     ggplot(aes(x = Iteration, y = value, color = as.factor(Chain))) +
# #     geom_line(alpha = 0.8) +
# #     facet_wrap(~ Parameter, scales = "free_y", ncol = 6) +  # More columns for many parameters
# #     labs(
# #       title = "MCMC Traceplots: Site × Status Interactions",
# #       x = "Iteration",
# #       y = "Parameter Value",
# #       color = "Chain"
# #     ) +
# #     theme_bw() +
# #     theme(
# #       strip.text = element_text(size = 6),  # Smaller text for many panels
# #       legend.position = "bottom"
# #     )
# #   
# #   print(traceplot_interactions)
# #   
# #   # ggsave("mcmc_traceplots_interactions.png", traceplot_interactions, 
# #   #        width = 18, height = 12, dpi = 300)
# # }
# 
# # ---- Summary of parameter counts ----
# param_summary = ggs_obj %>%
#   mutate(param_family = str_extract(Parameter, "^[^\\[]+")) %>%
#   count(param_family, name = "n_parameters") %>%
#   arrange(desc(n_parameters))
# 
# print("Parameter family summary:")
# print(param_summary)
# 
# 
# Gut fullness mcmc traceplots ---------------
fit = readRDS(here::here("res", "fit_beta_reg.rds"))
ggs_obj = ggmcmc::ggs(fit)

ggs_main = ggs_obj %>%
  filter(
    str_detect(Parameter, "^alpha(\\[|$)") |
      str_detect(Parameter, "^beta_site(\\[|$)") |
      str_detect(Parameter, "^beta_RS(\\[|$)") |
      str_detect(Parameter, "^beta_TL(\\[|$)") |
      str_detect(Parameter, "^phi(\\[|$)") |
      str_detect(Parameter, "^omega(\\[|$)") |
      str_detect(Parameter, "^omega_TL(\\[|$)") |
      str_detect(Parameter, "^sigma_site(\\[|$)") |
      str_detect(Parameter, "^sigma_RS(\\[|$)") |
      str_detect(Parameter, "^sigma_int(\\[|$)")
  )

main_families = c("alpha", "beta_site", "beta_RS", "beta_TL", "phi", "omega", "omega_TL", "sigma_site", "sigma_RS", "sigma_int")

ncol_settings = list(
  alpha = 3,        # S species (3 parameters)
  beta_site = 3,    # S x (I-1) = 3 x 6 = 18 parameters
  beta_RS = 3,      # S x (TT-1) = 3 x 1 = 3 parameters
  beta_TL = 3,      # S species (3 parameters)
  phi = 3,
  omega = 3,
  omega_TL = 3,
  sigma_site = 3,
  sigma_RS = 3,
  sigma_int = 3
)

p = list()
for (fam in main_families) {
  ggs_fam = ggs_main %>%
    filter(str_detect(Parameter, paste0("^", fam, "(\\[|$)")))
  
  if (nrow(ggs_fam) > 0) {
    ncol_fam = ncol_settings[[fam]]
    
    p[[fam]] = ggs_fam %>%
      ggplot(aes(x = Iteration, y = value, color = as.factor(Chain))) +
      geom_line(alpha = 0.8) +
      facet_wrap(~ Parameter, scales = "free_y", ncol = ncol_fam) +
      labs(
        title = paste("Traceplots for", fam),
        x = "Iteration",
        y = "Parameter Value",
        color = "Chain"
      ) +
      theme_bw() +
      theme(
        strip.text = element_text(size = 8),
        legend.position = "bottom"
      )
    
    print(p[[fam]])
    
    n_params = length(unique(ggs_fam$Parameter))
    plot_width = min(16, 4 * ncol_fam)
    plot_height = max(6, ceiling(n_params / ncol_fam) * 2.5)
    
    ggsave(paste0(sup_fig_dir, "traceplot_", fam, ".png"), p[[fam]],
           width = plot_width, height = plot_height, dpi = 300)
  }
}

# ---- Hyperprior traceplots (optional) ----
ggs_hyperpriors = ggs_obj %>%
  filter(
    str_detect(Parameter, "lp__")  # log posterior
  )

if (nrow(ggs_hyperpriors) > 0) {
  traceplot_hyperpriors = ggs_hyperpriors %>%
    ggplot(aes(x = Iteration, y = value, color = as.factor(Chain))) +
    geom_line(alpha = 0.8) +
    facet_wrap(~ Parameter, scales = "free_y", ncol = 3) +
    labs(
      title = "MCMC Traceplots: Hyperpriors & Log Posterior",
      x = "Iteration",
      y = "Parameter Value",
      color = "Chain"
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(size = 8),
      legend.position = "bottom"
    )
  
  print(traceplot_hyperpriors)
  
  # ggsave("mcmc_traceplots_hyperpriors.png", traceplot_hyperpriors, 
  #        width = 12, height = 8, dpi = 300)
}

# ---- Interaction terms traceplots (optional - can be many parameters) ----
# ggs_interactions = ggs_obj %>%
#   filter(str_detect(Parameter, "^beta_int(\\[|$)"))
# 
# if (nrow(ggs_interactions) > 0) {
#   traceplot_interactions = ggs_interactions %>%
#     ggplot(aes(x = Iteration, y = value, color = as.factor(Chain))) +
#     geom_line(alpha = 0.8) +
#     facet_wrap(~ Parameter, scales = "free_y", ncol = 6) +  # More columns for many parameters
#     labs(
#       title = "MCMC Traceplots: Site × Status Interactions",
#       x = "Iteration",
#       y = "Parameter Value",
#       color = "Chain"
#     ) +
#     theme_bw() +
#     theme(
#       strip.text = element_text(size = 6),  # Smaller text for many panels
#       legend.position = "bottom"
#     )
#   
#   print(traceplot_interactions)
#   
#   # ggsave("mcmc_traceplots_interactions.png", traceplot_interactions, 
#   #        width = 18, height = 12, dpi = 300)
# }

# ---- Summary of parameter counts ----
param_summary = ggs_obj %>%
  mutate(param_family = str_extract(Parameter, "^[^\\[]+")) %>%
  count(param_family, name = "n_parameters") %>%
  arrange(desc(n_parameters))

print("Parameter family summary:")
print(param_summary)