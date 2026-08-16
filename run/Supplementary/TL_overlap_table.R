#rm(list = ls())
pack_list = c("dplyr", "tidyr", "purrr", "here")
lapply(pack_list, library, character.only = TRUE)

data_dir = here::here("data")
raw_data = read.csv(file.path(data_dir, "rawData.csv"))
raw_data = raw_data %>% filter(!(Species.code == "BAICHR" & Site == "CI" & Treatment == "LS"))
RS_levels = sort(c("CT", "LS"))

species_config = list(
  Pinfish = list(spps = sort(c("LAGRHO")),
                 sites = sort(c("AM", "DR", "HWP", "LB", "NEPaP", "SA")),
                 prey  = sort(c("Amphipod","Crustacean", "Fish", "Isopod", "Polychaete", "SAV", "Tanaidacea"))),
  Croaker = list(spps = sort(c("MICUND")),
                 sites = sort(c("CI", "NEPaP", "SA")),
                 prey  = sort(c("Amphipod", "Crustacean", "Fish", "Isopod", "Polychaete", "Tanaidacea"))),
  `Silver perch` = list(spps = sort(c("BAICHR")),
                        sites = sort(c("CI", "NEPaP", "SA")),
                        prey  = sort(c("Amphipod", "Crustacean", "Fish", "Isopod", "Polychaete", "Tanaidacea")))
)

check_TL_range = function(cfg, sp_name) {
  wg_data =
    raw_data %>%
    filter(Species.code %in%cfg$spps) %>%
    filter(Site %in% cfg$sites) %>%
    filter(Treatment %in% RS_levels) %>%
    filter(Group %in% cfg$prey) %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = Group, values_from = present, values_fill = list(present = 0)) 
  
  global_mean = mean(wg_data$Length, na.rm = TRUE)
  
  out =
    wg_data %>%
    group_by(Site, Treatment) %>%
    summarise(
      n = n(),
      min_TL = min(Length, na.rm = TRUE),
      max_TL = max(Length, na.rm = TRUE),
      mean_TL = round(mean(Length, na.rm = TRUE), 1),
      global_mean_TL = round(global_mean, 1),
      mean_within_range = global_mean>= min_TL & global_mean <= max_TL,
      .groups = "drop"
    )
  out$Species =sp_name
  out
}

pinfish_tl = check_TL_range(species_config$Pinfish, "Pinfish")
pinfish_tl

croaker_tl  = check_TL_range(species_config$Croaker, "Croaker")
croaker_tl

silverperch_tl = check_TL_range(species_config$`Silver perch`, "Silver perch")
silverperch_tl