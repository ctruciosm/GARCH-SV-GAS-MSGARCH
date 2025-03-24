################################################################################
#####                                Utils                                 #####
################################################################################

loss_mse <- function(h_proxy, h_fore)  mean((h_proxy - h_fore)^2)
loss_qlike <- function(h_proxy, h_fore) mean(log(h_fore) + h_proxy/h_fore)
loss_mse_log <- function(h_proxy, h_fore) mean((log(h_proxy) - log(h_fore))^2)
loss_mse_sd <- function(h_proxy, h_fore) mean((sqrt(h_proxy) - sqrt(h_fore))^2)
loss_mse_prop <- function(h_proxy, h_fore) mean((h_proxy / h_fore - 1)^2)
loss_mae <- function(h_proxy, h_fore)  mean(abs(h_proxy - h_fore))
loss_mae_log <- function(h_proxy, h_fore)  mean(abs(log(h_proxy) - log(h_fore)))
loss_mae_sd <- function(h_proxy, h_fore) mean(abs(sqrt(h_proxy) - sqrt(h_fore)))
loss_mae_prop <- function(h_proxy, h_fore) mean(abs(h_proxy / h_fore - 1))


metrics <- function(h_proxy, h_fore) {
  c(loss_mse(h_proxy, h_fore),
    loss_qlike(h_proxy, h_fore),
    loss_mse_log(h_proxy, h_fore),
    loss_mse_sd(h_proxy, h_fore),
    loss_mse_prop(h_proxy, h_fore),
    loss_mae(h_proxy, h_fore),
    loss_mae_log(h_proxy, h_fore),
    loss_mae_sd(h_proxy, h_fore),
    loss_mae_prop(h_proxy, h_fore))
}


probability_regime_given_time_n <- function(p, q, sigma, r, Pt) {
  numA <- (1 - q) * dnorm(r, 0, sigma[2]) * (1 - Pt)
  numB <- p * dnorm(r, 0, sigma[1]) * Pt
  deno <- dnorm(r, 0, sigma[1]) * Pt + dnorm(r, 0, sigma[2]) * (1 - Pt)
  l <- numA / deno + numB / deno
  return(l)
}


probability_regime_given_time_t <- function(p, q, sigma, r, Pt, nu) {
  numA <- (1 - q) * sqrt(nu / (nu - 2)) / sigma[2] * dt(r * sqrt(nu / (nu - 2)) / sigma[2], nu) * (1 - Pt)
  numB <- p * sqrt(nu / (nu - 2)) / sigma[1] * dt(r * sqrt(nu / (nu - 2)) / sigma[1], nu) * Pt
  deno <- sqrt(nu / (nu - 2)) / sigma[1] * dt(r * sqrt(nu / (nu - 2)) / sigma[1], nu) * Pt +
    sqrt(nu / (nu - 2)) / sigma[2] * dt(r * sqrt(nu / (nu - 2)) / sigma[2], nu) * (1 - Pt)
  l <- numA / deno + numB / deno
  return(l)
}


msgarchfit <- function(spec, data) {
  is_error <- TRUE
  k <- 0
  opt_methods <- c("BFGS", "Nelder-Mead", "CG", "SANN")
  expr <- NULL
  while (is_error == TRUE && k <4) {
    k <- k + 1
    tryCatch(
      expr <- {
        fit_model <- FitML(spec, data, ctr = list(do.se = FALSE, do.plm = FALSE, OptimFUN = function(vPw, f_nll, spec, data, do.plm){
          out <- stats::optim(vPw, f_nll, spec = spec, data = data,
            do.plm = do.plm, method = opt_methods[k])}))
      },
      error = function(e) {
        is_error <- TRUE
      })
    if (!is.null(expr)) {
      is_error <- FALSE
    }
  }
  if (is_error == TRUE) {
    k <- 0
    expr <- NULL
    while (is_error == TRUE && k <4) {
      k <- k + 1
      tryCatch(
        expr <- {
          fit_model <- FitML(spec, data, ctr = list(do.se = FALSE, do.plm = TRUE, OptimFUN = function(vPw, f_nll, spec, data, do.plm){
            out <- stats::optim(vPw, f_nll, spec = spec, data = data,
              do.plm = do.plm, method = opt_methods[k])}))
        },
        error = function(e) {
          is_error <- TRUE
        })
      if (!is.null(expr)) {
        is_error <- FALSE
      }
    }
  }
  return(fit_model)
}





gasfit <- function(spec, data) {
  is_error <- TRUE
  k <- 0
  opt_methods <- c("BFGS", "Nelder-Mead", "CG", "SANN")
  expr <- NULL
  while (is_error == TRUE) {
    k <- k + 1
    tryCatch(
      expr <- {
        fit_model <- UniGASFit(spec, data, Compute.SE = FALSE, fn.optimizer = function(par0, data, GASSpec, FUN) {
          optimizer = optim(par0, FUN, data = data, GASSpec = GASSpec, method = opt_methods[k],
            control = list(trace = 0), hessian = FALSE)
          out = list(pars = optimizer$par,
            value = optimizer$value,
            hessian = optimizer$hessian,
            convergence = optimizer$convergence)
          return(out)
        })
      },
      error = function(e) {
        is_error <- TRUE
      })
    if (!is.null(expr)) {
      is_error <- FALSE
    }
  }
  return(fit_model)
}


