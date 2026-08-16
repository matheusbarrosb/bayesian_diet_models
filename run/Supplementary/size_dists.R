#rm(list = ls())
pack_list = c("dplyr", "tidyr", "purrr", "ggplot2", "here")
lapply(pack_list, library, character.only = TRUE)

fun_dir = here::here("R")
fun_files = list.files(fun_dir, pattern = "\\.R$", full.names = TRUE)
lapply(fun_files, source)

data_dir = here::here("data")
raw_data = read.csv(file.path(data_dir, "rawData.csv"))
raw_data = raw_data %>% filter(!(Species.code == "BAICHR" & Site == "CI" & Treatment == "LS"))

RS_levels = sort(c("CT", "LS"))

species_config = list(
  Pinfish = list(spps = sort(c("LAGRHO")),
                 sites = sort(c("AM", "DR", "HWP", "LB", "NEPaP", "SA")),
                 prey  = sort(c("Amphipod", "Crustacean", "Fish", "Isopod", "Polychaete", "SAV", "Tanaidacea"))),
  Croaker = list(spps = sort(c("MICUND")),
                 sites = sort(c("CI", "NEPaP", "SA")),
                 prey  = sort(c("Amphipod", "Crustacean", "Fish", "Isopod", "Polychaete", "Tanaidacea"))),
  `Silver perch` = list(spps = sort(c("BAICHR")),
                        sites = sort(c("CI", "NEPaP", "SA")),
                        prey  = sort(c("Amphipod", "Crustacean", "Fish", "Isopod", "Polychaete", "Tanaidacea")))
)

wg_data_all = purrr::imap_dfr(species_config, function(cfg, sp_name) {
  out =
    raw_data %>%
    filter(Species.code %in% cfg$spps) %>%
    filter(Site %in% cfg$sites) %>%
    filter(Treatment %in% RS_levels) %>%
    filter(Group %in% cfg$prey) %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = Group, values_from = present, values_fill = list(present = 0)) %>%
    select(Site, Treatment, Species.code, Length, Wet.Weight, Gut.weight)
  
  out$Species = sp_name
  out
})

sample_sizes_per_site_spp =
  wg_data_all %>%
  group_by(Site, Species) %>%
  summarise(n = n(), .groups = "drop")

wg_data_all %>%
  ggplot(aes(x = Length, fill = Species)) +
  geom_histogram(binwidth = 1, position = "dodge") +
  facet_grid(Species~Site) +
  labs(x = "Total length (mm)", y = "Count", fill = "Species") +
  custom_theme() +
  geom_text(
    data = sample_sizes_per_site_spp,
    aes(x = 100, y = 10, label = paste0("n = ", n)),
    inherit.aes = FALSE,
    size = 3
  ) +
  theme(legend.position = "none")

ggsave(filename = here::here("res", "figures", "supplementary", "size_distributions.png"), width = 8, height = 6, dpi = 300)