plot_gutweight_ppcheck = function(
    post,           
    stan_data,      
    spps = c("Pinfish", "Croaker", "Silver perch"),
    palette = "Bay",
    title,
    upperX = 1,
    strip_text_size = 10
) {
  library(dplyr)
  library(ggplot2)
  library(PNWColors)
  
  n_iter = dim(post$alpha)[1]
  n_species = length(spps)
  
  rep_TL = mean(stan_data$TL, na.rm = TRUE)
  
  # Calculate mean site effect padded with 0 for the baseline
  mean_site_effect = array(NA, dim = c(n_iter, n_species))
  for (s in seq_along(spps)) {
    padded_sites = cbind(0, post$beta_site[, s, ])
    mean_site_effect[, s] = rowMeans(padded_sites)
  }
  
  # Calculate mean status effect padded with 0 for the baseline
  mean_status_effect = array(NA, dim = c(n_iter, n_species))
  for (s in seq_along(spps)) {
    padded_status = cbind(0, post$beta_RS[, s, ])
    mean_status_effect[, s] = rowMeans(padded_status)
  }
  
  pred_df = data.frame()
  for (s in seq_along(spps)) {
    # Using the marginalized components for an "average" fish
    lp = post$alpha[, s] +
      mean_site_effect[, s] +
      mean_status_effect[, s] +
      post$beta_TL[, s] * rep_TL
    
    mu = plogis(lp)
    phi = post$phi[, s]
    
    theta = plogis(post$omega[, s] + post$omega_TL[, s] * rep_TL)
    is_empty = rbinom(n_iter, 1, theta)

    GF_sim = rbeta(n_iter, mu * phi, (1 - mu) * phi)
    
    pred_df = rbind(pred_df, data.frame(
      species = spps[s],
      draw = 1:n_iter,
      GF = GF_sim
    ))
  }
  
  obs_df = data.frame(
    GF = stan_data$GF,
    species = spps[stan_data$species]
  )
  
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
    theme(legend.position = "top",
          strip.text = element_blank()) +
    xlim(0, upperX)
  
  return(p)
}
