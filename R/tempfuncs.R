
#' @export
fitted_chi = function(a_func, b_func, thresholds, unique_dist) {
  chi = matrix(0, nrow = length(unique_dist), ncol = length(thresholds))
  count = rep(0, length(thresholds))
  for (y0 in unlist(data$y0)) {
    a = a_func(y0, unique_dist)
    b = b_func(y0, unique_dist)
    for (i in seq_along(thresholds)) {
      if (y0 > thresholds[i]) {
        chi[, i] = chi[, i] + pnorm(thresholds[i], a, b, lower.tail = FALSE)
        count[i] = count[i] + 1
      }
    }
  }
  chi = chi / rep(count, each = length(unique_dist))
  attr(chi, "thresholds") = thresholds
  attr(chi, "unique_dist") = unique_dist
  chi
}

#' @export
empirical_chi = function(data, thresholds, unique_dist) {
  n_big = matrix(0, length(unique_dist), length(thresholds))
  n_all = matrix(0, length(unique_dist), length(thresholds))
  pb = progress_bar(length(data$y))
  for (j in seq_along(data$y)) {
    dist = round(data$dist_to_s0[[j]])
    for (d in seq_along(unique_dist)) {
      d_index = which(unique_dist[d] == dist)
      if (!any(d_index)) next
      for (i in seq_along(thresholds)) {
        y0_index = which(data$y0[[j]] > thresholds[i])
        if (!any(y0_index)) next
        n_big[d, i] = n_big[d, i] +
          sum(data$y[[j]][d_index, y0_index] > thresholds[i], na.rm = TRUE)
        n_all[d, i] = n_all[d, i] +
          sum(!is.na(data$y[[j]][d_index, y0_index]))
      }
    }
    pb$tick()
  }
  pb$terminate()
  chi = n_big / n_all
  attr(chi, "thresholds") = thresholds
  attr(chi, "unique_dist") = unique_dist
  chi
}

plot_chi = function(chi) {
  t = attr(chi, "thresholds")
  unique_dist = attr(chi, "unique_dist")
  data.frame(
    chi = as.numeric(rbind(1, chi)),
    dist = rep(c(0, unique_dist), length(t)),
    t = rep(t, each = length(unique_dist) + 1)) |>
    ggplot() +
    geom_line(aes(x = dist, y = chi, group = t, col = t))
}

#' @export
empirical_moments = function(data, unique_dist, unique_y0) {
  y_sum = matrix(0, length(unique_dist), length(unique_y0))
  y_sq_sum = matrix(0, length(unique_dist), length(unique_y0))
  y_n = matrix(0, length(unique_dist), length(unique_y0))
  pb = progress_bar(length(data$y))
  for (j in seq_along(data$y)) {
    dist = round(data$dist_to_s0[[j]])
    y0 = round(data$y0[[j]], 1)
    for (d in seq_along(unique_dist)) {
      d_index = which(unique_dist[d] == dist)
      if (!any(d_index)) next
      for (i in seq_along(unique_y0)) {
        y0_index = which(y0 == unique_y0[i])
        if (!any(y0_index)) next
        y_sum[d, i] = y_sum[d, i] + sum(data$y[[j]][d_index, y0_index], na.rm = TRUE)
        y_sq_sum[d, i] = y_sq_sum[d, i] + sum(data$y[[j]][d_index, y0_index]^2, na.rm = TRUE)
        y_n[d, i] = y_n[d, i] + sum(!is.na(data$y[[j]][d_index, y0_index]))
      }
    }
    pb$tick()
  }
  pb$terminate()
  y_mean = y_sum / y_n
  y_sq_mean = y_sq_sum / y_n
  y_sd = sqrt(y_sq_mean - y_mean^2)
  res = list(mean = y_mean, sd = y_sd)
  attr(res, "unique_dist") = unique_dist
  attr(res, "unique_y0") = unique_y0
  res
}

plot_moments = function(x) {
  unique_dist = attr(x, "unique_dist")
  unique_y0 = attr(x, "unique_y0")
  sd_plot = data.frame(
    sd = as.numeric(rbind(0, x$sd)),
    dist = rep(c(0, unique_dist), length(unique_y0)),
    y0 = rep(unique_y0, each = length(unique_dist) + 1)) |>
    ggplot() +
    geom_line(aes(x = dist, y = sd, group = y0, col = y0))
  mean_plot = data.frame(
    mean = as.numeric(rbind(unique_y0, x$mean)),
    dist = rep(c(0, unique_dist), length(unique_y0)),
    y0 = rep(unique_y0, each = length(unique_dist) + 1)) |>
    ggplot() +
    geom_line(aes(x = dist, y = mean, group = y0, col = y0))
  list(mean = mean_plot, sd = sd_plot)
}

compare_chi3 = function(theta, data, multimesh_data, p, unique_dist, k) {
  tau = exp(theta[1])
  lambda = exp(theta[2])
  kappa = exp(theta[3])
  sigma = exp(theta[4])
  rho_b = exp(theta[5])
  beta0 = exp(theta[6]) / (1 + exp(theta[6]))
  beta1 = exp(theta[7])
  rho = exp(theta[8])
  chi_emp = empirical_chi(data, p, unique_dist)
  for (i in seq_along(data$y)) {
    dist_from_mesh = multimesh_data$dist[[i]]
    dist = data$dist_to_s0[[i]]
    Q = inla.spde2.precision(multimesh_data$spde[[i]], theta = c(log(rho), log(sigma)))
    z = rnorm_spde(length(data$y0[[i]]), Q)
    for (j in seq_along(data$y0[[i]])) {
      a = data$y0[[i]][j] * exp(-(dist / lambda)^kappa)
      b = sqrt(1 - exp(-2 * dist_from_mesh / rho_b)) *
        data$y0[[i]][j] ^ softplus(beta0 - beta1 * dist_from_mesh, k)
      data$y[[i]][, j] = a +
        as.numeric(multimesh_data$A[[i]] %*% (b * z[, j])) +
        rnorm(length(a), 1 / sqrt(tau))
    }
  }
  chi = empirical_chi(data, p, unique_dist)
  df1 = data.frame(
    chi = as.numeric(rbind(1, chi_emp)),
    dist = rep(c(0, unique_dist), length(p)),
    p = rep(p, each = length(unique_dist) + 1),
    tag = "empirical")
  df2 = data.frame(
    chi = as.numeric(rbind(1, chi)),
    dist = rep(c(0, unique_dist), length(p)),
    p = rep(p, each = length(unique_dist) + 1),
    tag = "model fit")
  rbind(df1, df2) |>
    ggplot() +
    geom_line(aes(x = dist, y = chi, col = p, group = p)) +
    facet_wrap(~tag)
}



compare_chi2 = function(theta, data, p, unique_dist) {
  lambda = exp(theta[1])
  kappa = exp(theta[2])
  sigma = exp(theta[3])
  rho_b = exp(theta[4])
  #beta0 = exp(theta[5]) / (1 + exp(theta[5]))
  beta0 = theta[5]
  beta1 = exp(theta[6])
  chi_emp = empirical_chi(data, p, unique_dist)
  for (i in seq_along(data$y)) {
    dist = data$dist_to_s0[[i]]
    for (j in seq_along(data$y0[[i]])) {
      a = data$y0[[i]][j] * exp(-(dist / lambda)^kappa)
      beta = exp(beta0 - beta1 * dist)
      beta = beta / (1 + beta)
      b = sigma * sqrt(1 - exp(-2 * dist / rho_b)) * data$y0[[i]][j] ^ beta
      data$y[[i]][, j] = rnorm(length(a), a, b)
    }
  }
  chi = empirical_chi(data, p, unique_dist)
  df1 = data.frame(
    chi = as.numeric(rbind(1, chi_emp)),
    dist = rep(c(0, unique_dist), length(p)),
    p = rep(p, each = length(unique_dist) + 1),
    tag = "empirical")
  df2 = data.frame(
    chi = as.numeric(rbind(1, chi)),
    dist = rep(c(0, unique_dist), length(p)),
    p = rep(p, each = length(unique_dist) + 1),
    tag = "model fit")
  rbind(df1, df2) |>
    ggplot() +
    geom_line(aes(x = dist, y = chi, col = p, group = p)) +
    facet_wrap(~tag)
}

compare_chi = function(a_func, b_func, data, p, unique_dist) {
  xx = seq(0, max(unique_dist), length.out = 1000)[-1]
  chi_emp = empirical_chi(data, p, unique_dist)
  thresholds = qlaplace(p)
  chi = matrix(0, nrow = length(xx), ncol = length(p))
  count = rep(0, length(p))
  for (y0 in unlist(data$y0)) {
    a = a_func(y0, xx)
    b = b_func(y0, xx)
    for (i in seq_along(thresholds)) {
      if (y0 > thresholds[i]) {
        chi[, i] = chi[, i] + pnorm(thresholds[i], a, b, lower.tail = FALSE)
        count[i] = count[i] + 1
      }
    }
   }
  chi = chi / rep(count, each = length(xx))
  df1 = data.frame(
    chi = as.numeric(rbind(1, chi_emp)),
    dist = rep(c(0, unique_dist), length(p)),
    p = rep(p, each = length(unique_dist) + 1),
    tag = "empirical")
  df2 = data.frame(
    chi = as.numeric(rbind(1, chi)),
    dist = rep(c(0, xx), length(p)),
    p = rep(p, each = length(xx) + 1),
    tag = "model fit")
  rbind(df1, df2) |>
    ggplot() +
    geom_line(aes(x = dist, y = chi, col = p, group = p)) +
    facet_wrap(~tag)
}

#' @export
get_s0_index = function(coords, delta_s0 = 1) {
  stopifnot(all(as.integer(delta_s0) == delta_s0))
  x_coords = coords[, 1] |> unique() |> sort()
  y_coords = coords[, 2] |> unique() |> sort()
  s0_locs = expand.grid(
    x = x_coords[seq(delta_s0, length(x_coords), by = delta_s0)],
    y = y_coords[seq(delta_s0, length(y_coords), by = delta_s0)]) |>
    as.matrix()
  s0_index = lapply(
    X = 1:nrow(s0_locs),
    FUN = function(i) which(coords[, 1] == s0_locs[i, 1] & coords[, 2] == s0_locs[i, 2]))
  s0_index = s0_index[sapply(s0_index, length) > 0] |>
    unlist() |>
    unname()
  s0_index
}

#' @export
fitted_moments = function(a_func, b_func, unique_dist, unique_y0) {
  mean = list()
  sd = list()
  for (y0 in unique_y0) {
    mean[[length(mean) + 1]] = a_func(y0, unique_dist)
    sd[[length(sd) + 1]] = b_func(y0, unique_dist)
  }
  mean = do.call(cbind, mean)
  sd = do.call(cbind, sd)
  res = list(mean = mean, sd = sd)
  attr(res, "unique_y0") = unique_y0
  attr(res, "unique_dist") = unique_dist
  res
}

#' @export
compare_moments = function(a_func, b_func, data, unique_dist, unique_y0) {
  xx = seq(0, max(unique_dist), length.out = 1000)
  moments1 = empirical_moments(data, unique_dist, unique_y0)
  moments2 = fitted_moments(a_func, b_func, unique_dist, unique_y0)

  
  plot1 = local({
    df1 = lapply(
      unique_y0,
      function(y0) {
        data.frame(
          dist = xx,
          tag = 1,
          y0 = y0,
          sd = b_func(y0, xx))
      }) |>
      do.call(what = rbind)
    df2 = data.frame(
      sd = as.numeric(rbind(0, moments$sd)),
      dist = rep(c(0, unique_dist), length(unique_y0)),
      y0 = rep(unique_y0, each = length(unique_dist) + 1),
      tag = 2)
    rbind(df1, df2) |>
      ggplot() +
      geom_line(aes(x = dist, y = sd, group = y0, col = y0)) +
      facet_wrap(~tag)
  })
  plot2 = local({
    df1 = lapply(
      unique_y0,
      function(y0) {
        data.frame(
          dist = xx,
          tag = 1,
          y0 = y0,
          mean = a_func(y0, xx))
      }) |>
      do.call(what = rbind)
    df2 = data.frame(
      mean = as.numeric(rbind(unique_y0, moments$mean)),
      dist = rep(c(0, unique_dist), length(unique_y0)),
      y0 = rep(unique_y0, each = length(unique_dist) + 1),
      tag = 2)
    rbind(df1, df2) |>
      ggplot() +
      geom_line(aes(x = dist, y = mean, group = y0, col = y0)) +
      facet_wrap(~tag)
  })
  plot3 = local({
    df1 = data.frame(
      dist = xx,
      tag = "sigma",
      value = b_func(1, xx))
    df2 = data.frame(
      dist = xx,
      tag = "beta",
      value = log(b_func(exp(1), xx) / b_func(1, xx)))
    rbind(df1, df2) |>
      ggplot() +
      geom_line(aes(x = dist, y = value)) +
      facet_wrap(~tag, scales = "free")
  })
  list(plot1, plot2, plot3)
}
