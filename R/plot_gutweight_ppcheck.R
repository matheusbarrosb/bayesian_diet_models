plot_gutweight_ppcheck = function(
    post,            
    stan_data,      
    spps = c("Silver perch", "Pinfish", "Croaker"),
    display_order = c("Silver perch", "Pinfish", "Croaker"),
    palette = "Bay",
    title,
    upperX = 1,
    strip_text_size = 10,
    n_fill = dim(post$alpha)[1]
) {
  library(dplyr)
  library(ggplot2)
  library(PNWColors)
  
  n_iter = dim(post$alpha)[1]
  n_species = length(spps)
  
  rep_TL = mean(stan_data$TL, na.rm = TRUE)
  
  # calculate mean site effect padded with 0 for the baseline
  mean_site_effect = array(NA, dim = c(n_iter, n_species))
  for (s in seq_along(spps)) {
    padded_sites = cbind(0, post$beta_site[, s, ])
    mean_site_effect[, s] = rowMeans(padded_sites)
  }
  
  # calculate mean status effect padded with 0 for the baseline
  mean_status_effect = array(NA, dim = c(n_iter, n_species))
  for (s in seq_along(spps)) {
    padded_status = cbind(0, post$beta_RS[, s, ])
    mean_status_effect[, s] = rowMeans(padded_status)
  }
  
  pred_df = data.frame()
  for (s in seq_along(spps)) {
    
    lp = numeric(n_iter)
    
    for (i in 1:n_fill) {
      lp[i] = post$alpha[i,s] + mean_site_effect[i,s] + mean_status_effect[i, s] + post$beta_TL[i, s]*rep_TL
    }
    
    mu = plogis(lp)
    phi = post$phi[, s]
    
    GF_sim = rbeta(n_iter, mu *phi, (1 - mu) * phi)
    
    pred_df = rbind(pred_df, data.frame(
      species = spps[s],
      draw = 1:n_iter,
      GF = GF_sim
    ))
  }
  
  obs_df = data.frame(GF = stan_data$GF, species = spps[stan_data$species])
  
  pred_df$species = factor(pred_df$species, levels = display_order)
  obs_df$species = factor(obs_df$species,levels = display_order)
  
  p = ggplot() +
    geom_density(data = pred_df, aes(x = GF, fill = "Predicted"), alpha = 0.3, fill = pnw_palette(palette, 2)[1]) +
    geom_density(data = obs_df, aes(x = GF, fill = "Observed"), alpha = 0.3, fill = pnw_palette(palette, 2)[2]) +
    facet_wrap(~ species, scales = "free_y") +
    labs(
      title = title,
      x = "Gut fullness",
      y = "Density",
      fill = ""
    ) +
    custom_theme() +
    theme(legend.position = "top", strip.text = element_blank()) +
    xlim(0, upperX)
  
  return(p)
}