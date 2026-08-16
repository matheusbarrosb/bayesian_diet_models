plot_betas = function(
    post,
    spps = c("Pinfish", "Croaker", "Silver perch"),
    site_labels = c("AM", "DR", "HWP", "LB", "NEPaP", "SA", "CI"),
    status_labels = c("CT", "LS"),
    param_panels = c("beta_RS_overall", "beta_site", "beta_TL", "beta_int"),
    title = "",
    palette = "Bay",
    ncol
) {
  
  n_species = length(spps)
  species_colors = pnw_palette(palette, n = n_species)
  
  present_sites = list(
    "Pinfish" = site_labels,
    "Croaker" = c("CI", "NEPaP", "SA"),
    "Silver perch" = c("CI", "NEPaP", "SA")
  )
  
  beta_long = data.frame()
  
  for (param in param_panels) {
    if (param == "beta_RS_overall" && !is.null(post$beta_RS_overall)) {
      n_status_params = dim(post$beta_RS_overall)[3]
      n_iter = dim(post$beta_RS_overall)[1]
      for (s in seq_along(spps)) {
        for (tt in seq_len(n_status_params)) {
          df = data.frame(
            draw = seq_len(n_iter),
            value = post$beta_RS_overall[, s, tt],
            species = spps[s],
            index = status_labels[tt + 1], 
            parameter = "beta_RS_overall"
          )
          beta_long = bind_rows(beta_long, df)
        }
      }
    } else if (param == "beta_site" && !is.null(post$beta_site)) {
      n_site_params = dim(post$beta_site)[3]
      n_iter  = dim(post$beta_site)[1]
      for (s in seq_along(spps)) {
        for (site_name in present_sites[[spps[s]]]) {
          if (site_name == site_labels[1]) next 
          i = which(site_labels == site_name) - 1 
          if (length(i) == 1 && i <= n_site_params && i > 0) {
            df = data.frame(
              draw = seq_len(n_iter),
              value = post$beta_site[, s, i],
              species = spps[s],
              index = site_name,
              parameter = "beta_site"
            )
            beta_long = bind_rows(beta_long, df)
          }
        }
      }
    } else if (param == "beta_TL" && !is.null(post$beta_TL)) {
      n_iter = dim(post$beta_TL)[1]
      for (s in seq_along(spps)) {
        df = data.frame(
          draw = seq_len(n_iter),
          value = post$beta_TL[, s],
          species = spps[s],
          index = "",
          parameter = "beta_TL"
        )
        beta_long = bind_rows(beta_long, df)
      }
    } else if (param == "beta_int" && !is.null(post$beta_int)) {
      n_site_params = dim(post$beta_int)[3]
      n_status_params = dim(post$beta_int)[4]
      n_iter = dim(post$beta_int)[1]
      for (s in seq_along(spps)) {
        for (site_name in present_sites[[spps[s]]]) {
          if (site_name == site_labels[1]) next
          i = which(site_labels == site_name) - 1
          if (length(i) == 1 && i <= n_site_params && i > 0) {
            for (tt in seq_len(n_status_params)) {
              df = data.frame(
                draw = seq_len(n_iter),
                value = post$beta_int[, s, i, tt],
                species = spps[s],
                index = paste(site_name, status_labels[tt + 1], sep = " × "),
                parameter = "beta_int"
              )
              beta_long = bind_rows(beta_long, df)
            }
          }
        }
      }
    }
  }
  
  if (nrow(beta_long) == 0) stop("No valid parameters found.")
  
  beta_summary = beta_long %>%
    group_by(parameter, species, index) %>%
    summarise(
      mean = mean(value),
      q025 = quantile(value, 0.1),
      q975 = quantile(value, 0.9),
      .groups = "drop"
    )
  
  pretty_labels = c(
    "beta_RS_overall" = "beta[status]",
    "beta_site" = "beta[site]",
    "beta_TL" = "beta[TL]",
    "beta_int" = "beta[site*','*status]"
  )
  
  beta_summary$parameter = factor(beta_summary$parameter, levels = names(pretty_labels), labels = pretty_labels)
  beta_summary$species = factor(beta_summary$species, levels = c("Silver perch", "Croaker", "Pinfish"))
  
  p = ggplot(beta_summary, aes(x = mean, y = species)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    geom_errorbarh(aes(xmin = q025, xmax = q975), height = 0, position = position_dodge(width = 0.5)) +
    facet_wrap(~ parameter, scales = "free_x", labeller = label_parsed, nrow = 2) +
    labs(
      y = NULL,
      x = "Estimate",
      title = title
    ) +
    custom_theme() +
    theme(
      legend.title = element_blank(),
      legend.position = "top",
      strip.text = element_text(face = "bold", size = 12),
      axis.title.y = element_text(angle = 0, vjust = 0.5),
      axis.text.y = element_text(size = 10),
      panel.spacing = unit(0.5, "lines"),
      strip.placement = "outside"
    ) +
    scale_x_continuous(labels = label_number(accuracy = 0.01))
  
  return(p)
}
