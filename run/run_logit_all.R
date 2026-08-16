# runs all logit models at once, saves figures to directory
rm(list = ls())
for (i in 1:3) {
  run_dir = here::here("run")
  source(list.files(run_dir, pattern = ".R$", full.names = TRUE)[i])
}
