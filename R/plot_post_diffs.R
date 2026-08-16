plot_post_diffs = function(
    post,
    species_names = c("Pinfish", "Croaker", "Silver perch"),
    status_names = c("CT", "LS"),
    palette = "Bay",
    title = NULL,
    strip_text_size = 13,
    prob_text_size = 3,
    absolute_scale = FALSE
) {
  stopifnot(length(status_names) == 2)
  
  library(ggplot2)
  library(dplyr)
  library(PNWColors)
  
  diff_df = data.frame()
  annot_df = data.frame()
  segment_df = data.frame()
  y_max = numeric(length(species_names))
  
  for (s in seq_along(species_names)) {
    
    if (absolute_scale) {
      mu_CT = apply(post$y_hat_int[, s, , 1], 1, mean)
      mu_LS = apply(post$y_hat_int[, s, , 2], 1, mean)
      diff = mu_LS - mu_CT
      x_label = "Absolute difference in % gut fullness"
    } else {
      # In dummy coding, beta_RS[, s, 1] is exactly the logit difference (LS - CT)
      diff = post$beta_RS[, s, 1]
      x_label = "Posterior difference (Logit scale)"
    }
    
    prob_LS_greater = mean(diff > 0)
    
    diff_df = rbind(
      diff_df,
      data.frame(
        species = species_names[s],
        diff = diff
      )
    )
    
    dens = density(diff)
    y_max[s] = max(dens$y)
    yloc = max(dens$y) * 1.25
    x_center = 0 
    
    annot_df = rbind(
      annot_df,
      data.frame(
        species = species_names[s],
        xpos = x_center,
        ypos = yloc,
        label = paste0("P(", status_names[2], " > ", status_names[1], ") = ", round(prob_LS_greater, 3))
      )
    )
    
    segment_df = rbind(
      segment_df,
      data.frame(
        species = species_names[s],
        x = 0,
        xend = 0,
        y = 0,
        yend = yloc * 0.95
      )
    )
  }
  
  diff_df$species = factor(diff_df$species, levels = species_names)
  annot_df$species = factor(annot_df$species, levels = species_names)
  segment_df$species = factor(segment_df$species, levels = species_names)
  colors = pnw_palette(palette, length(species_names))
  global_ymax = max(y_max) * 1.35
  
  p = ggplot(diff_df, aes(x = diff, fill = species)) +
    geom_density(alpha = 0.5, color = "black") +
    facet_wrap(~ species, nrow = 1, strip.position = "top") +
    geom_segment(
      data = segment_df,
      aes(x = x, xend = xend, y = y, yend = yend),
      linetype = "dashed",
      inherit.aes = FALSE
    ) +
    geom_text(
      data = annot_df,
      aes(x = xpos, y = ypos, label = label),
      inherit.aes = FALSE,
      size = prob_text_size,
      fontface = "plain"
    ) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, global_ymax)) +
    labs(
      title = title,
      x = x_label,
      y = "Density"
    ) +
    custom_theme() +
    theme(
      legend.title = element_blank(),
      strip.text = element_text(size = strip_text_size, face = "bold"),
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "none"
    )
  
  return(p)
}
