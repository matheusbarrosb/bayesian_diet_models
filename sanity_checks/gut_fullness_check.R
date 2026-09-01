rm(list = ls())
library(dplyr)
library(ggplot2)
library(PNWColors)

spps       = sort(c("LAGRHO", "MICUND", "BAICHR"))
sites      = sort(c("AM", "DR", "HWP", "LB", "NEPaP", "SA", "CI"))
RS_levels  = sort(c("CT", "LS"))
spp_lookup = c(BAICHR = "Silver perch", LAGRHO = "Pinfish", MICUND = "Croaker")

raw = read.csv(here::here("data", "rawData.csv")) %>%
  filter(Species.code %in% spps, Site %in% sites, Treatment %in% RS_levels) %>%
  filter(!is.na(Gut.weight), !is.na(Wet.Weight), Wet.Weight > 0, !is.na(Length))

fish = raw %>%
  group_by(Fish.ID, Fish_ID_year, Species.code, Site, Treatment, Length, Wet.Weight, Gut.weight) %>%
  summarise(.groups = "drop") %>%
  mutate(species = unname(spp_lookup[Species.code]),
         GF     = Gut.weight / Wet.Weight,
         status = ifelse(Treatment == "LS", "Restored", "Control"))

stopifnot(!anyNA(fish$species))
print(table(fish$species))

nrow(fish)
n_distinct(raw$Fish.ID, raw$Fish_ID_year) 
table(fish$species)

gf_site_status = fish %>%
  group_by(species, Site, status) %>%
  summarise(n = n(),
            mean_GF = round(mean(GF),3),
            sd_GF   = round(sd(GF),4),
            se_GF   = round(sd_GF / sqrt(n),4),
            .groups = "drop")

gf_status = fish %>%
  group_by(species, status) %>%
  summarise(n = n(),
            mean_GF = round(mean(GF),3),
            sd_GF   = round(sd(GF),4),
            se_GF   = round(sd_GF / sqrt(n),4),
            .groups = "drop")

write.csv(gf_site_status, here::here("sanity_checks", "raw_GF_by_site_status.csv"), row.names = FALSE)
write.csv(gf_status,      here::here("sanity_checks", "raw_GF_by_status.csv"),      row.names = FALSE)

gf_status_plot =
  gf_status %>%
  mutate(species = factor(species, levels = c("Pinfish", "Croaker", "Silver perch"))) %>%
  ggplot(aes(x = status, y = mean_GF, color = status)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_GF - 2*se_GF, ymax = mean_GF + 2*se_GF), width = 0.15) +
  facet_wrap(~ species) +
  scale_color_manual(values = pnw_palette("Bay", 2)) +
  labs(x = "", y = "Mean raw gut fullness (\u00b1 SE)") +
  theme_bw() +
  theme(legend.position = "none")

gf_site_plot =
  gf_site_status %>%
  mutate(species = factor(species, levels = c("Pinfish", "Croaker", "Silver perch"))) %>%
  ggplot(aes(x = Site, y = mean_GF, color = status)) +
  geom_point(size = 2, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = mean_GF - se_GF, ymax = mean_GF + se_GF),
                width = 0.2, position = position_dodge(width = 0.4)) +
  facet_wrap(~ species) +
  scale_color_manual(values = pnw_palette("Bay", 2)) +
  labs(x = "", y = "Mean raw gut fullness (\u00b1 SE)", color = "") +
  theme_bw() +
  theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(here::here("sanity_checks", "raw_GF_by_status.png"),      gf_status_plot, width = 7, height = 3.5, dpi = 300)
ggsave(here::here("sanity_checks", "raw_GF_by_site_status.png"), gf_site_plot,   width = 9, height = 4,   dpi = 300)

print(gf_status)

# ---- verify the mechanism: which fish were dropped by adding the length filter? ----

raw_before = read.csv(here::here("data", "rawData.csv")) %>%
  filter(Species.code %in% spps, Site %in% sites, Treatment %in% RS_levels) %>%
  filter(!is.na(Gut.weight), !is.na(Wet.Weight), Wet.Weight > 0)  # no Length filter

fish_before = raw_before %>%
  group_by(Fish.ID, Fish_ID_year, Species.code, Site, Treatment, Length, Wet.Weight, Gut.weight) %>%
  summarise(.groups = "drop") %>%
  mutate(species = unname(spp_lookup[Species.code]),
         status  = ifelse(Treatment == "LS", "Restored", "Control"))

cat("Before Length filter:\n")
print(table(fish_before$species))
cat("\nAfter Length filter (current script):\n")
print(table(fish$species))

dropped = fish_before %>%
  anti_join(fish, by = c("Fish.ID", "Fish_ID_year", "Species.code", "Site", "Treatment", "Wet.Weight", "Gut.weight"))

cat("\nFish dropped by the Length filter:\n")
print(dropped %>% select(Fish.ID, species, Site, Treatment, Length, Wet.Weight, Gut.weight))
cat("\ndropped fish by species x status:\n")
print(dropped %>% count(species, status))
# confirm they all have NA Length (the actual mechanism)
stopifnot(all(is.na(dropped$Length)))
cat("\nConfirmed: all dropped fish have missing lengths.\n")

gf_site_before = fish_before %>%
  group_by(species, Site, status) %>%
  summarise(n_before = n(), .groups = "drop")

gf_site_after = fish %>%
  group_by(species, Site, status) %>%
  summarise(n_after = n(), .groups = "drop")

site_compare = gf_site_before %>%
  full_join(gf_site_after, by = c("species", "Site", "status")) %>%
  mutate(n_before = coalesce(n_before, 0L),
         n_after  = coalesce(n_after, 0L),
         diff = n_before - n_after) %>%
  filter(diff != 0)

cat("Sites where n changed:\n")
print(site_compare)

# cross-check: the sites that changed should exactly match the dropped fish's sites/status
expected = dropped %>% count(species, Site, status = status, name = "expected_diff")
dropped2 = dropped %>% mutate(status = ifelse(Treatment == "LS", "Restored", "Control"))
expected = dropped2 %>% count(species, Site, status, name = "expected_diff")

cat("\nExpected diffs (from dropped fish list):\n")
print(expected)

check = site_compare %>%
  select(species, Site, status, diff) %>%
  full_join(expected, by = c("species","Site","status"))

cat("\nComparison (diff should equal expected_diff, no NAs):\n")
print(check)

stopifnot(all(check$diff == check$expected_diff, na.rm = FALSE))
cat("\nConfirmed: site-level changes exactly match the 3 dropped fish, no unexpected sites affected.\n")