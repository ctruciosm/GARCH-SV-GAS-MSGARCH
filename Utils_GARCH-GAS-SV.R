################################################################################
#####                                Utils                                 #####
################################################################################

loss_mse <- function(h_hat, h)  mean((h_hat - h)^2)
loss_qlike <- function(h_hat, h) mean(log(h) + h_hat/h)
loss_qlike_ult <- function(h_hat, h) mean(h_hat/h - log(h_hat/h) - 1)
loss_mse_log <- function(h_hat, h) mean((log(h_hat) - log(h))^2)
loss_mse_sd <- function(h_hat, h) mean((sqrt(h_hat) - sqrt(h))^2)
loss_mse_prop <- function(h_hat, h) mean((h_hat / h - 1)^2)
loss_mae <- function(h_hat, h)  mean(abs(h_hat - h))
loss_mae_log <- function(h_hat, h)  mean(abs(log(h_hat) - log(h)))
loss_mae_sd <- function(h_hat, h) mean(abs(sqrt(h_hat) - sqrt(h)))
loss_mae_prop <- function(h_hat, h) mean(abs(h_hat / h - 1))


metrics <- function(h_hat, h) {
  c(loss_mse(h_hat, h),
    loss_qlike(h_hat, h),
    loss_qlike_ult(h_hat, h),
    loss_mse_log(h_hat, h),
    loss_mse_sd(h_hat, h),
    loss_mse_prop(h_hat, h),
    loss_mae(h_hat, h),
    loss_mae_log(h_hat, h),
    loss_mae_sd(h_hat, h),
    loss_mae_prop(h_hat, h))
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
  while (is_error == TRUE) {
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
  return(fit_model)
}

