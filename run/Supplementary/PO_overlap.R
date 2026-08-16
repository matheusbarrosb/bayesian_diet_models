library(dplyr)
library(tidyr)
library(stringr)
library(ggmcmc)

## Pinfish ##

# ---- Extract posterior draws for p_status ----
df_pstatus = ggmcmc::ggs(fit_pinfish, family = "p_status")

prey_names = sort(c("Amphipod", "Crustacean", "Fish", "Isopod", "Polychaete", "SAV", "Tanaidacea"))

# Parse indices
df_parsed = df_pstatus %>%
  mutate(
    indices = str_extract(Parameter, "\\[(.*?)\\]"),
    indices = str_remove_all(indices, "\\[|\\]"),
    prey_idx   = as.integer(str_split_fixed(indices, ",", 2)[,1]),
    status_idx = as.integer(str_split_fixed(indices, ",", 2)[,2]),
    Prey   = prey_names[prey_idx],
    Status = ifelse(status_idx == 1, "CT", "LS")
  ) %>%
  select(-indices, -prey_idx, -status_idx)

# ---- Compute 80% credible intervals ----
ci_summary = df_parsed %>%
  group_by(Prey, Status) %>%
  summarise(
    lower = quantile(value, 0.20),
    upper = quantile(value, 0.80),
    mean  = mean(value),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Status,
    values_from = c(lower, upper, mean),
    names_glue = "{.value}_{Status}"
  )

# ---- Test interval overlap for each prey ----
overlap_results = ci_summary %>%
  mutate(
    overlap = !(upper_CT < lower_LS | upper_LS < lower_CT),
    overlap_pct = (pmin(upper_CT, upper_LS) - pmax(lower_CT, lower_LS)) /
      (pmax(upper_CT, upper_LS) - pmin(lower_CT, lower_LS)),
    overlap_pct = ifelse(overlap, overlap_pct, 0)
  ) %>%
  mutate(across(where(is.numeric), ~ round(., 2))) %>%
  select(Prey,
         LowerCT = lower_CT, UpperCT = upper_CT, MeanCT = mean_CT,
         LowerLS = lower_LS, UpperLS = upper_LS, MeanLS = mean_LS,
         Overlap = overlap, `Overlap_%` = overlap_pct)

res_dir = file.path(here::here(), "res", "figures", "supplementary", "tables")
write.csv(overlap_results, file.path(res_dir, "pinfish_PO_overlap.csv"), row.names = FALSE)


## Croaker ##

# ---- Extract posterior draws for p_status ----
df_pstatus = ggmcmc::ggs(fit_micund, family = "p_status")

prey_names = sort(c("Amphipod", "Crustacean", "Fish", "Isopod", "Polychaete", "Tanaidacea"))

# Parse indices
df_parsed = df_pstatus %>%
  mutate(
    indices = str_extract(Parameter, "\\[(.*?)\\]"),
    indices = str_remove_all(indices, "\\[|\\]"),
    prey_idx   = as.integer(str_split_fixed(indices, ",", 2)[,1]),
    status_idx = as.integer(str_split_fixed(indices, ",", 2)[,2]),
    Prey   = prey_names[prey_idx],
    Status = ifelse(status_idx == 1, "CT", "LS")
  ) %>%
  select(-indices, -prey_idx, -status_idx)

# ---- Compute 80% credible intervals ----
ci_summary = df_parsed %>%
  group_by(Prey, Status) %>%
  summarise(
    lower = quantile(value, 0.20),
    upper = quantile(value, 0.80),
    mean  = mean(value),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Status,
    values_from = c(lower, upper, mean),
    names_glue = "{.value}_{Status}"
  )

# ---- Test interval overlap for each prey ----
overlap_results = ci_summary %>%
  mutate(
    overlap = !(upper_CT < lower_LS | upper_LS < lower_CT),
    overlap_pct = (pmin(upper_CT, upper_LS) - pmax(lower_CT, lower_LS)) /
      (pmax(upper_CT, upper_LS) - pmin(lower_CT, lower_LS)),
    overlap_pct = ifelse(overlap, overlap_pct, 0)
  ) %>%
  mutate(across(where(is.numeric), ~ round(., 2))) %>%
  select(Prey,
         LowerCT = lower_CT, UpperCT = upper_CT, MeanCT = mean_CT,
         LowerLS = lower_LS, UpperLS = upper_LS, MeanLS = mean_LS,
         Overlap = overlap, `Overlap_%` = overlap_pct)

write.csv(overlap_results, file.path(res_dir, "croaker_PO_overlap.csv"), row.names = FALSE)


## Silver perch ##

# ---- Extract posterior draws for p_status ----
df_pstatus = ggmcmc::ggs(fit_silverperch, family = "p_status")

prey_names = sort(c("Amphipod", "Crustacean", "Fish", "Isopod", "Polychaete", "Tanaidacea"))

# Parse indices
df_parsed = df_pstatus %>%
  mutate(
    indices = str_extract(Parameter, "\\[(.*?)\\]"),
    indices = str_remove_all(indices, "\\[|\\]"),
    prey_idx   = as.integer(str_split_fixed(indices, ",", 2)[,1]),
    status_idx = as.integer(str_split_fixed(indices, ",", 2)[,2]),
    Prey   = prey_names[prey_idx],
    Status = ifelse(status_idx == 1, "CT", "LS")
  ) %>%
  select(-indices, -prey_idx, -status_idx)

# ---- Compute 80% credible intervals ----
ci_summary = df_parsed %>%
  group_by(Prey, Status) %>%
  summarise(
    lower = quantile(value, 0.20),
    upper = quantile(value, 0.80),
    mean  = mean(value),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Status,
    values_from = c(lower, upper, mean),
    names_glue = "{.value}_{Status}"
  )

# ---- Test interval overlap for each prey ----
overlap_results = ci_summary %>%
  mutate(
    overlap = !(upper_CT < lower_LS | upper_LS < lower_CT),
    overlap_pct = (pmin(upper_CT, upper_LS) - pmax(lower_CT, lower_LS)) /
      (pmax(upper_CT, upper_LS) - pmin(lower_CT, lower_LS)),
    overlap_pct = ifelse(overlap, overlap_pct, 0)
  ) %>%
  mutate(across(where(is.numeric), ~ round(., 2))) %>%
  select(Prey,
         LowerCT = lower_CT, UpperCT = upper_CT, MeanCT = mean_CT,
         LowerLS = lower_LS, UpperLS = upper_LS, MeanLS = mean_LS,
         Overlap = overlap, `Overlap_%` = overlap_pct)

write.csv(overlap_results, file.path(res_dir, "silverperch_PO_overlap.csv"), row.names = FALSE)