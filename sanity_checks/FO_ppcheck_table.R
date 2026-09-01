library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(PNWColors)

RS_levels = sort(c("CT", "LS"))

species_config = list(
  Pinfish = list(spps = "LAGRHO",
                 sites = sort(c("AM","DR","HWP","LB","NEPaP","SA")),
                 prey = sort(c("Amphipod","Crustacean","Fish","Isopod","Polychaete","SAV","Tanaidacea"))),
  Croaker = list(spps = "MICUND",
                 sites = sort(c("CI","NEPaP","SA")),
                 prey = sort(c("Amphipod","Crustacean","Fish","Isopod","Polychaete","Tanaidacea"))),
  `Silver perch` = list(spps = "BAICHR",
                        sites = sort(c("CI","NEPaP","SA")),
                        prey = sort(c("Amphipod","Crustacean","Fish","Isopod","Polychaete","Tanaidacea")))
)

raw_data = read.csv(here::here("data","rawData.csv"))

# ---- observed FO per site ----
obs_site = imap_dfr(species_config, function(cfg, sp) {
  raw_data %>%
    filter(Species.code %in% cfg$spps, Site %in% cfg$sites, Treatment %in% RS_levels) %>%
    filter(Group %in% cfg$prey) %>%
    mutate(present = 1) %>%
    pivot_wider(
      id_cols = c(Fish.ID, Fish_ID_year, Species.code, Site, Treatment, Length, Wet.Weight, Gut.weight),
      names_from = Group, values_from = present, values_fn = ~1, values_fill = 0
    ) %>%
    pivot_longer(all_of(cfg$prey), names_to = "Prey", values_to = "present") %>%
    group_by(Site, Prey) %>%
    summarise(n = n(), Observed = mean(present), .groups = "drop") %>%
    mutate(Species = sp)
})

# ---- modelled FO per site: p_site, per species fit ----
get_psite = function(fit, cfg, sp) {
  post = rstan::extract(fit, pars = "p_site")$p_site
  expand.grid(Prey = cfg$prey, Site = cfg$sites, stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(g = which(cfg$prey == Prey), s = which(cfg$sites == Site),
           Model_Mean  = mean(post[, g, s]),
           Model_Lower = quantile(post[, g, s], 0.1),
           Model_Upper = quantile(post[, g, s], 0.9)) %>%
    ungroup() %>% select(-g, -s) %>% mutate(Species = sp)
}


fits = list(
  Pinfish        = fit_pinfish,
  Croaker        = fit_micund,
  `Silver perch` = fit_silverperch
)

model_site = imap_dfr(species_config, function(cfg, sp) get_psite(fits[[sp]], cfg, sp))

ppc = obs_site %>%
  left_join(model_site, by = c("Species","Site","Prey")) %>%
  mutate(covered = Observed >= Model_Lower & Observed <= Model_Upper) %>%
  mutate(across(c(Observed, Model_Mean, Model_Lower, Model_Upper), ~round(.x, 3)))

write.csv(ppc, here::here("sanity_checks","model_vs_observed_by_site.csv"), row.names = FALSE)

ppc_plot = ppc %>%
  ggplot(aes(x = Observed, y = Model_Mean)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = Model_Lower, ymax = Model_Upper), alpha = 0.3, width = 0) +
  geom_point(aes(color = Species), size = 1.8, alpha = 0.8) +
  scale_color_manual(values = pnw_palette("Bay", 3)) +
  coord_equal(xlim = c(0,1), ylim = c(0,1)) +
  labs(x = "Observed FO (raw, per site)", y = "Modelled FO (p_site, per site)") +
  theme_bw()
ggsave(here::here("sanity_checks","model_vs_observed_by_site.png"), ppc_plot, width = 6, height = 6, dpi = 300)

# coverage summary
cat("Overall coverage:", round(mean(ppc$covered), 3), "\n")
ppc %>% group_by(Species) %>% summarise(coverage = round(mean(covered), 3))

print(ppc)