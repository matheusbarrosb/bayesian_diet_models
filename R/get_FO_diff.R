get_FO_diff = function(fit, prey_names, status_levels = c("CT", "LS")) {
  post = rstan::extract(fit, pars = "p_status")$p_status
  ct = which(sort(status_levels) == "CT")
  ls = which(sort(status_levels) == "LS")
  data.frame(
    Prey = prey_names,
    P_LS_greater = sapply(seq_along(prey_names), function(j) mean(post[, j, ls] > post[, j, ct])),
    Diff_mean = sapply(seq_along(prey_names), function(j) mean(post[, j, ls] - post[, j, ct])),
    Diff_lower = sapply(seq_along(prey_names), function(j) quantile(post[, j, ls] - post[, j, ct], 0.1)),
    Diff_upper = sapply(seq_along(prey_names), function(j) quantile(post[, j, ls] - post[, j, ct], 0.9))
  )
}