#rm(list = ls())
# Packages ---------------------------------------------------------------------
pack_list = c("dplyr", "tidyr", "purrr", "here")
if (!all(pack_list %in% rownames(installed.packages()))) {
  install.packages(pack_list[!pack_list %in% rownames(installed.packages())])
}
lapply(pack_list, library, character.only = TRUE)

# descriptive FO function -------------------------------------------------------
get_descriptive_FO = function(prey_mat, status, prey_names) {
  purrr::map_dfr(seq_along(prey_names), function(j) {
    x_ct = sum(prey_mat[status == 1, j])
    n_ct = sum(status == 1)
    x_ls = sum(prey_mat[status == 2, j])
    n_ls = sum(status == 2)
    
    ft = fisher.test(matrix(c(x_ls, n_ls - x_ls, x_ct, n_ct - x_ct), nrow = 2))
    
    data.frame(
      Prey = prey_names[j],
      FO_CT = x_ct/n_ct,
      FO_LS = x_ls/n_ls,
      Desc_diff = x_ls/n_ls - x_ct / n_ct,
      Desc_fisher_p = ft$p.value
    )
  })
}

# load data ----------------------------------------------------------------------
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

# compute descriptive FO for each species -----------------------------------------
desc_all = purrr::imap_dfr(species_config, function(cfg, sp_name) {
  wg_data =
    raw_data %>%
    filter(Species.code %in% cfg$spps) %>%
    filter(Site %in% cfg$sites) %>%
    filter(Treatment %in% RS_levels) %>%
    filter(Group %in% cfg$prey) %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = Group, values_from = present,
                values_fill = list(present = 0))
  
  prey_mat = wg_data %>% select(all_of(cfg$prey)) %>% as.matrix()
  status = as.numeric(factor(wg_data$Treatment, levels = RS_levels))
  
  out = get_descriptive_FO(prey_mat, status, cfg$prey)
  out$Species = sp_name
  out
})

model_all = readRDS(here::here("res", "tables", "FO_model_diff_all.rds"))

fo_comparison =
  model_all %>%
  select(Species, Prey, Diff_mean, P_LS_greater, sig) %>%
  rename(Model_diff = Diff_mean, Model_P = P_LS_greater, Model_sig = sig) %>%
  left_join(desc_all, by = c("Species", "Prey")) %>%
  mutate(Desc_fisher_sig = Desc_fisher_p < 0.05) %>%
  mutate(Agreement = case_when(
    Model_sig & Desc_fisher_sig ~ "Both detect difference",
    !Model_sig & !Desc_fisher_sig ~ "Both detect no difference",
    Model_sig & !Desc_fisher_sig ~ "Model only",
    !Model_sig & Desc_fisher_sig ~ "Descriptive only")) %>%
  select(Species, Prey, FO_CT, FO_LS, Desc_diff, Desc_fisher_p, Model_diff, Model_P, Agreement) %>%
  mutate(across(where(is.numeric), ~ round(., 2)))

write.csv(fo_comparison, here::here("res", "tables", "FO_comparison_supplement.csv"), row.names = FALSE)
fo_comparison