rm()
options(error = NULL)
set.seed(444)
# Packages ---------------------------------------------------------------------
pack_list = c("dplyr", "rstan", "tidyr", "ggplot2", "purrr", "readr",
              "stringr", "here", "reshape2", "PNWColors", "ggpubr", "scales")
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
spps      = sort(c("LAGRHO", "MICUND", "BAICHR"))
sites     = sort(c("AM", "DR", "HWP", "LB", "NEPaP", "SA", "CI"))
RS_levels = sort(c("CT", "LS"))

wg_data =
  raw_data %>%
  filter(Species.code %in% spps) %>%
  filter(Site %in% sites) %>%
  filter(Treatment %in% RS_levels) %>%
  group_by(Fish.ID, Fish_ID_year, Species.code, Site, Treatment, Length, Wet.Weight, Gut.weight) %>%
  summarise(Z = as.numeric(any(Group == "Empty")), .groups = "drop")

n_distinct_fish =
  raw_data %>%
  filter(Species.code %in% spps, Site %in% sites, Treatment %in% RS_levels) %>%
  distinct(Fish.ID, Fish_ID_year, Species.code, Site, Treatment, Length, Wet.Weight, Gut.weight) %>%
  nrow()
stopifnot(n_distinct_fish == nrow(wg_data))

spp  = as.numeric(factor(wg_data$Species.code, levels = spps))
site = as.numeric(factor(wg_data$Site, levels = sites))
RS   = as.numeric(factor(wg_data$Treatment, levels = RS_levels))
GW   = as.numeric(wg_data$Gut.weight)
FW   = as.numeric(wg_data$Wet.Weight)
TL   = as.numeric(wg_data$Length)
Z    = wg_data$Z

# put everything into a dataframe and exclude NA values
df =
  data.frame(
    spp  = spp,
    site = site,
    RS   = RS,
    GW   = GW,
    FW   = FW,
    TL   = TL,
    Z    = Z
  ) %>%
  filter(!is.na(spp) & !is.na(site) & !is.na(RS) & !is.na(FW) & !is.na(TL)) %>%
  filter(FW > 0) %>%
  mutate(GF = GW/FW) %>% # calculate fullness
  filter(!is.na(GF))

global_mean_TL = df %>%
  group_by(spp) %>%
  summarize(mean_TL = mean(TL, na.rm = TRUE)) %>%
  arrange(spp) %>%
  pull(mean_TL)

stan_data = 
  list(
    N       = nrow(df),
    S       = length(spps),         
    I       = length(sites),         
    TT      = length(RS_levels),    
    species = df$spp,
    site    = df$site,              
    status  = df$RS,                
    GF      = df$GF,
    Z       = df$Z,
    TL      = df$TL,
    global_mean_TL = global_mean_TL
  )

# fit --------------------------------------------------------------------------
model_dir = here::here("stan/")
stanc(paste0(model_dir, "GF_beta_reg.stan"))
rstan_options(auto_write = TRUE)

REFIT = FALSE

# load fit object if it exists
if (file.exists(here::here("res", "fit_beta_reg.rds"))) {
  message("Loading existing fit object...")
  fit = readRDS(here::here("res", "fit_beta_reg.rds"))
} else {
  fit = NULL
}

# if fit is NULL or REFIT is TRUE, run the model
if (is.null(fit) || REFIT) {
  message("Running the model...")
  set.seed(413)
  fit = stan(
    file    = paste0(model_dir, "GF_beta_reg.stan"),  
    data    = stan_data,
    chains  = 3,        
    iter    = 10000,      
    warmup  = 1000,        
    cores   = 3,        
    seed    = 427,
    control = list(adapt_delta = 0.95,
                   max_treedepth = 20)
  )
}

# save fit if NULL
if (file.exists(here::here("res", "fit_beta_reg.rds"))) {
  message("Fit already exists, skipping save.")
} else {
  message("Saving fit object...")
  saveRDS(fit, file = here::here("res", "fit_beta_reg.rds"))
}

# plotting ---------------------------------------------------------------------
fit = readRDS(here::here("res", "fit_beta_reg.rds"))

### PROBABILITIES OF LS > CT ###
species_names = c("Pinfish", "Croaker", "Silver perch")
status_names  = c("CT", "LS") 

post = rstan::extract(fit, pars = c("beta_RS", "beta_site", "beta_TL", "beta_int", "alpha", "omega", "omega_TL", "phi", "y_hat_int"))
post = rstan::extract(fit, pars = c("beta_RS", "beta_site", "beta_TL", "beta_int", "alpha", "omega", "omega_TL", "phi", "y_hat_int"))

# calculate overall marginal beta_RS
n_iter = dim(post$beta_RS)[1]
n_species = dim(post$beta_RS)[2]

# create a new array to hold the overall status effect
post$beta_RS_overall = array(NA, dim = c(n_iter, n_species, 1))

for (s in 1:n_species) {
  # average the interaction terms across all sites (including the baseline site, which is 0)
  avg_interaction = rowMeans(cbind(0, post$beta_int[, s, , 1]))
  
  # overall effect = Baseline effect + Average interaction
  post$beta_RS_overall[, s, 1] = post$beta_RS[, s, 1] + avg_interaction
}

post_diffs_plot =
  plot_post_diffs(
    post            = post,
    species_names   = species_names,
    status_names    = status_names,
    prob_text_size  = 3,
    strip_text_size = 10,
    title           = "A",
    absolute_scale  = TRUE # TRUE for %, FALSE for logit
  ) + xlim(-0.1,0.2)

### PREDICTIVE CHECK ###
ppcheck_plot =
  plot_gutweight_ppcheck(post,
                         stan_data,
                         title   = "B",
                         palette = "Bay",
                         upperX = range(stan_data$GF)[2],
                         strip_text_size = 10)

### PLOT COEFFICIENTS ###
coef_plot =
  plot_betas(
    post          = post,
    spps          = c("Pinfish", "Croaker", "Silver perch"),
    site_labels   = c("AM", "DR", "HWP", "LB", "NEPaP", "SA", "CI"),
    status_labels = c("CT", "LS"),
    param_panels  = c("beta_RS_overall", "beta_TL"),
    palette       = "Bay",
    title         = "C",
    ncol          = 3
  )

# arrange 
plot1 =
  ggarrange(
    post_diffs_plot,
    ppcheck_plot,
    nrow    = 2,
    label.y = "Density"
  ); print(plot1)

plot_final = 
  ggarrange(
    plot1,
    coef_plot,
    nrow = 1,
    widths = c(2, 1)
  ); print(plot_final)

fig_dir = here::here("res", "figures")

ggsave(
  filename = file.path(fig_dir, "gut_fullness.pdf"),
  plot = plot_final,
  width = 8.1,
  height = 4.5,
  units = "in",
  device = cairo_pdf
)

##### PLOT GUT WEIGHT ~ FISH WEIGHT #####
spp_labs = c(
  LAGRHO = "Pinfish",
  MICUND = "Croaker",
  BAICHR = "Silver perch"
)

# verify how many empty stomachs have weights
print(
  raw_data %>%
    filter(Species.code %in% c("LAGRHO", "MICUND", "BAICHR")) %>%
    filter(Gut.weight > 0.0005) %>%
    group_by(Group) %>%
    summarize(count = n())
  ,n = 100)

# original scale
raw_data %>%
  filter(Species.code %in% c("LAGRHO", "MICUND", "BAICHR")) %>%
  filter(Gut.weight > 0.005) %>%
  mutate(Species.code = recode(Species.code, !!!spp_labs)) %>%
  mutate(Species.code = factor(Species.code,
                               levels = c("Pinfish", "Croaker", "Silver perch"))) %>%
  group_by(Fish_ID_year, Group, Species.code) %>%
  mutate(Gut_Status = ifelse(Group == "Empty", "Empty", "Non-empty")) %>%
  filter(Gut.weight < 1) %>% # crazy silver perch outlier
  as.data.frame() %>%
  
  ggplot(aes(x = Length, y = Gut.weight, color = Gut_Status, size = Gut_Status)) +
  geom_point(alpha = 0.6, shape = 16) +
  scale_size_manual(values = c("Empty" = 2.0, "Non-empty" = 0.7), guide = "none") +
  facet_wrap(~Species.code, scales = "free_y") +
  xlim(20,125) +
  custom_theme() +
  xlab("Length (mm)") +
  ylab("Gut weight (g)") +
  theme(legend.position = "right", legend.title = element_blank()) +
  scale_color_manual(values = pnw_palette("Bay", 2))


## log-log
raw_data %>%
  filter(Species.code %in% c("LAGRHO", "MICUND", "BAICHR")) %>%
  filter(Gut.weight > 0.005) %>%
  mutate(Species.code = recode(Species.code, !!!spp_labs)) %>%
  mutate(Species.code = factor(Species.code,
                               levels = c("Pinfish", "Croaker", "Silver perch"))) %>%
  group_by(Fish_ID_year, Group, Species.code) %>%
  mutate(Gut_Status = ifelse(Group == "Empty", "Empty", "Non-empty")) %>%
  filter(Gut.weight < 1) %>% 
  as.data.frame() %>%
  
  ggplot(aes(x = log(Length), y = log(Gut.weight), color = Gut_Status, size = Gut_Status)) +
  geom_point(alpha = 0.6, shape = 16) +
  scale_size_manual(values = c("Empty" = 2.0, "Non-empty" = 0.7), guide = "none") +
  facet_wrap(~Species.code, scales = "free_y") +
  custom_theme() +
  xlab("log Length (mm)") +
  ylab("log Gut weight (g)") +
  theme(legend.position = "right", legend.title = element_blank()) +
  scale_color_manual(values = pnw_palette("Bay", 2))

# count the number and proportion of empty stomachs per species
raw_data %>%
  filter(Species.code %in% c("LAGRHO", "MICUND", "BAICHR")) %>%
  mutate(Species.code = recode(Species.code, !!!spp_labs)) %>%
  mutate(Species.code = factor(Species.code,
                               levels = c("Pinfish", "Croaker", "Silver perch"))) %>%
  group_by(Group, Species.code) %>%
  mutate(Gut_Status = ifelse(Group == "Empty", "Empty", "Non-empty")) %>%
  group_by(Species.code) %>%
  summarize(
    total = n(),
    empty = sum(Gut_Status == "Empty"),
    prop_empty = empty / total
  )