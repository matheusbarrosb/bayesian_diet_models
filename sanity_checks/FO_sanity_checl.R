rm(list = ls())
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(PNWColors)

RS_levels = sort(c("CT", "LS"))

species_config = list(
  Pinfish = list(spps = "LAGRHO",
                 sites = sort(c("AM", "DR", "HWP", "LB", "NEPaP", "SA")),
                 prey  = sort(c("Amphipod", "Crustacean", "Fish", "Isopod", "Polychaete", "SAV", "Tanaidacea"))),
  Croaker = list(spps = "MICUND",
                 sites = sort(c("CI", "NEPaP", "SA")),
                 prey  = sort(c("Amphipod", "Crustacean", "Fish", "Isopod", "Polychaete", "Tanaidacea"))),
  `Silver perch` = list(spps = "BAICHR",
                        sites = sort(c("CI", "NEPaP", "SA")),
                        prey  = sort(c("Amphipod", "Crustacean", "Fish", "Isopod", "Polychaete", "Tanaidacea")))
)

 raw_data = read.csv(here::here("data", "rawData.csv"))

 build_presence = function(cfg, sp_name) {
   raw_data %>%
     filter(Species.code %in% cfg$spps, Site %in% cfg$sites, Treatment %in% RS_levels) %>%
     filter(Group %in% cfg$prey) %>%
     mutate(present = 1) %>%
     pivot_wider(
       id_cols = c(Fish.ID, Fish_ID_year, Species.code, Site, Treatment, Length, Wet.Weight, Gut.weight),
       names_from = Group,
       values_from = present,
       values_fn = ~1,
       values_fill = 0
     ) %>%
     mutate(status = ifelse(Treatment == "LS", "Restored", "Control"),
            Species = sp_name)
 }

pres_all = imap(species_config, build_presence)

pres_all$Pinfish %>% count(Treatment)        # expect 116 CT / 208 LS
pres_all$Croaker %>% count(Treatment)        # expect 74 CT / 43 LS
pres_all$`Silver perch` %>% count(Treatment) # expect 50 CT / 169 LS

### 1. FO by status (overall) ###
fo_status = imap_dfr(species_config, function(cfg, sp_name) {
  d = pres_all[[sp_name]]
  map_dfr(cfg$prey, function(p) {
    d %>% group_by(status) %>%
      summarise(n = n(), FO = round(mean(.data[[p]]),3), .groups = "drop") %>%
      mutate(Prey = p, Species = sp_name)
  })
})
write.csv(fo_status, here::here("sanity_checks", "raw_FO_by_status.csv"), row.names = FALSE)

### 2. FO by site and status ###
fo_site_status = imap_dfr(species_config, function(cfg, sp_name) {
  d = pres_all[[sp_name]]
  map_dfr(cfg$prey, function(p) {
    d %>% group_by(Site, status) %>%
      summarise(n = n(), FO = round(mean(.data[[p]]),3), .groups = "drop") %>%
      mutate(Prey = p, Species = sp_name)
  })
})
write.csv(fo_site_status, here::here("sanity_checks", "raw_FO_by_site_status.csv"), row.names = FALSE)

### 3. FO by size class (10 mm bins, pooled across sites & status) ###
fo_size = imap_dfr(species_config, function(cfg, sp_name) {
  d = pres_all[[sp_name]] %>%
    mutate(size_bin = cut(Length, breaks = seq(30, 200, by = 10), right = FALSE))
  map_dfr(cfg$prey, function(p) {
    d %>% filter(!is.na(size_bin)) %>%
      group_by(size_bin) %>%
      summarise(n = n(), FO = round(mean(.data[[p]]),3), .groups = "drop") %>%
      mutate(Prey = p, Species = sp_name)
  })
})
write.csv(fo_size, here::here("sanity_checks", "raw_FO_by_sizeclass.csv"), row.names = FALSE)

sp_levels = c("Pinfish", "Croaker", "Silver perch")

### Figure: FO by status ###
fo_status_plot =
  fo_status %>%
  mutate(Species = factor(Species, levels = sp_levels)) %>%
  ggplot(aes(x = FO, y = Prey, color = status)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  facet_wrap(~ Species, scales = "free_y") +
  scale_color_manual(values = pnw_palette("Bay", 2)) +
  labs(x = "Raw frequency of occurrence", y = "", color = "") +
  theme_bw() +
  theme(legend.position = "top")
ggsave(here::here("sanity_checks", "raw_FO_by_status.png"), fo_status_plot, width = 9, height = 4, dpi = 300)

### Figure: FO by site and status ###
fo_site_plot =
  fo_site_status %>%
  mutate(Species = factor(Species, levels = sp_levels)) %>%
  ggplot(aes(x = FO, y = Prey, color = status)) +
  geom_point(size = 1.8, position = position_dodge(width = 0.5)) +
  facet_grid(Species ~ Site, scales = "free_y") +
  scale_color_manual(values = pnw_palette("Bay", 2)) +
  labs(x = "Raw frequency of occurrence", y = "", color = "") +
  theme_bw() +
  theme(legend.position = "top")
ggsave(here::here("sanity_checks", "raw_FO_by_site_status.png"), fo_site_plot, width = 12, height = 6, dpi = 300)

### Figure: FO by size class ###
fo_size_plot =
  fo_size %>%
  mutate(Species = factor(Species, levels = sp_levels)) %>%
  ggplot(aes(x = size_bin, y = FO, group = Prey, color = Prey)) +
  geom_line() +
  geom_point(aes(size = n), alpha = 0.7) +
  facet_wrap(~ Species, ncol = 1, scales = "free_x") +
  scale_size_continuous(range = c(0.5, 3)) +
  labs(x = "Size class (mm, 10 mm bins)", y = "Raw frequency of occurrence", color = "Prey") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(here::here("sanity_checks", "raw_FO_by_sizeclass.png"), fo_size_plot, width = 9, height = 10, dpi = 300)

print(fo_status)

# check if sample sizes from model matches raw data
##
raw_n = fo_status %>%
  distinct(Species, status, n) %>%
  arrange(Species, status)
print(raw_n)










