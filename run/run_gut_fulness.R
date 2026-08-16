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

# put everything into a dataframe and exclude NA values
df =
data.frame(
  spp  = spp,
  site = site,
  RS   = RS,
  GW   = GW,
  FW   = FW,
  TL   = TL
) %>%
  filter(!is.na(spp) & !is.na(site) & !is.na(RS) & !is.na(GW) & !is.na(FW) & !is.na(TL))

# put into list for stan
stan_data = 
list(
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

# fit --------------------------------------------------------------------------
model_dir = here::here("stan/")
stanc(paste0(model_dir, "GF_centered.stan"))
rstan_options(auto_write = TRUE)

# load fit object if it exists
if (file.exists(here::here("res", "fit_gut_fullness.rds"))) {
  fit = readRDS(here::here("res", "fit_gut_fullness.rds"))
} else {
  fit = NULL
}

# if fit is NULL, run the model
# fit = NULL
if (is.null(fit)) {
  message("Running the model...")
  set.seed(413)
  fit = stan(
    file   = paste0(model_dir, "GF_centered.stan"),  
    data   = stan_data,
    chains = 3,        
    iter   = 30000,      
    warmup = 3000,       
    cores  = 3,        
    seed   = 427,
    control = list(adapt_delta = 0.95,
                   max_treedepth = 20)
  )
}

# save fit if NULL
if (file.exists(here::here("res", "fit_gut_fullness.rds"))) {
  message("Fit already exists, skipping save.")
} else {
  message("Saving fit object...")
  saveRDS(fit, file = here::here("res", "fit_gut_fullness.rds"))
}


# plotting ---------------------------------------------------------------------
fit = readRDS(here::here("res", "fit_gut_fullness.rds"))

### PROBABILITIES OF LS > CT ###
species_names = c("Pinfish", "Croaker", "Silver perch")
status_names  = c("C", "R") # CT is reference (status==1)

post = rstan::extract(fit,
                      pars = c("beta_RS", "beta_TL", "beta_FW"))

post_diffs_plot =
plot_post_diffs(
  beta_RS_array   = post$beta_RS,
  species_names   = species_names,
  status_names    = status_names,
  prob_text_size  = 3,
  strip_text_size = 10,
  title           = "A"
) +
  xlim(-1,2)

### PREDICTIVE CHECK ###
post = rstan::extract(fit, pars = c("alpha", "beta_site", "beta_RS", "beta_TL", "beta_FW"))
ppcheck_plot =
plot_gutweight_ppcheck(post, stan_data, title = "B", palette = "Bay", strip_text_size = 10)

### PLOT COEFFICIENTS ###
coef_plot =
plot_betas(
  post = post,
  spps = c("Pinfish", "Croaker", "Silver perch"),
  site_labels = c("AM", "DR", "HWP", "LB", "NEPaP", "SA", "CI"),
  status_labels = c("CT", "LS"),
  param_panels = c("beta_RS", "beta_TL", "beta_FW"),
  palette = "Bay",
  title = "B",
  ncol = 3
)

# arrange 
plot1 =
ggarrange(
  post_diffs_plot,
  coef_plot,
  nrow = 2,
  label.y = "Density"
  #align = "v"
); print(plot1)


fig_dir = here::here("res", "figures")

ggsave(
  filename = file.path(fig_dir, "gut_fullness.pdf"),
  plot = plot1,
  width = 5,
  height = 4,
  units = "in",
  device = cairo_pdf
)

##### PLOT GUT WEIGHT ~ FISH WEIGHT #####
spp_labs = c(
  LAGRHO = "Pinfish",
  MICUND = "Croaker",
  BAICHR = "Silver perch"
)
colnames(raw_data)

# original scale
raw_data %>%
  filter(Species.code %in% c("LAGRHO", "MICUND", "BAICHR")) %>%
  # filter gut weights < 0.01 g
  filter(Gut.weight > 0.005) %>%
  mutate(Species.code = recode(Species.code, !!!spp_labs)) %>%
  mutate(Species.code = factor(Species.code,
                               levels = c("Pinfish", "Croaker", "Silver perch"))) %>%
  group_by(Fish_ID_year, Group, Species.code) %>%
  mutate(Gut_Status = ifelse(Group == "Empty", "Empty", "Non-empty")) %>%
  filter(Gut.weight < 1) %>%
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
  # filter gut weights < 0.01 g
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
  #xlim(20,125) +
  custom_theme() +
  xlab("log Length (mm)") +
  ylab("log Gut weight (g)") +
  theme(legend.position = "right", legend.title = element_blank()) +
  scale_color_manual(values = pnw_palette("Bay", 2))

# count the number and proportion of empty stomachs per species
raw_data %>%
  filter(Species.code %in% c("LAGRHO", "MICUND", "BAICHR")) %>%
  # filter gut weights < 0.01 g
  filter(Gut.weight > 0.005) %>%
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





