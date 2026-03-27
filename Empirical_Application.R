################################################################################
#####                   Empirical Application                              #####
################################################################################
library(dplyr)
library(readxl)
library(stringr)
library(rugarch)
library(GAS)
library(stochvol)
library(stochvolTMB)
library(MSGARCH)
library(future.apply)
source("DGPs.R")
source("Utils_GARCH-GAS-SV.R")

# Data
data <- read_excel("./Data/capire_daily_returns.xlsx", skip = 3, col_types = c("date", rep("numeric", 30)), na = c("", "-", NA)) |> 
  filter(Data > "2010-01-01" & Data < "2025-01-01") |> 
  filter(!if_all(where(is.numeric), is.na)) |> 
  select(where(~ !any(is.na(.x)))) |> 
  rename_with(~ str_remove(.x, "^(?s).*prov\n"), -Data)

ins <- 1000
oos <- nrow(data) - 2500

# Specs
garch_spec_n <- ugarchspec(variance.model = list(model = 'sGARCH', garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), distribution.model = "norm")
garch_spec_t <- ugarchspec(variance.model = list(model = 'sGARCH', garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), distribution.model = "std")
gas_spec_n <- UniGASSpec(Dist = "norm", ScalingType = "Identity", GASPar = list(locate = FALSE, scale = TRUE, shape = FALSE))
gas_spec_t <- UniGASSpec(Dist = "std", ScalingType = "Identity", GASPar = list(locate = FALSE, scale = TRUE, shape = FALSE))
ms_spec_n <- CreateSpec(variance.spec = list(model = c("sGARCH", "sGARCH")), switch.spec = list(do.mix = FALSE), distribution.spec = list(distribution = c("norm", "norm")))
ms_spec_t <- CreateSpec(variance.spec = list(model = c("sGARCH", "sGARCH")), switch.spec = list(do.mix = FALSE), distribution.spec = list(distribution = c("std", "std")), constraint.spec = list(regime.const = c("nu")))
figarch_spec_n <- ugarchspec(variance.model = list(model = 'fiGARCH', garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0,0), include.mean = FALSE), distribution = 'norm')
figarch_spec_t <- ugarchspec(variance.model = list(model = 'fiGARCH', garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0,0), include.mean = FALSE), distribution = 'std')

garch_n_fore_s <- garch_t_fore_s <- figarch_n_fore_s <- figarch_t_fore_s <- gas_n_fore_s <- gas_t_fore_s <- ms_n_fore_s <- ms_t_fore_s <- sv_n_fore_s <- sv_t_fore_s <- matrix(0, nrow = oos, ncol = ncol(data) - 1, dimnames = list(NULL, colnames(data)[-1]))
garch_n_fore_var1 <- garch_t_fore_var1 <- figarch_n_fore_var1 <- figarch_t_fore_var1 <- gas_n_fore_var1 <- gas_t_fore_var1 <- ms_n_fore_var1 <- ms_t_fore_var1 <- sv_n_fore_var1 <- sv_t_fore_var1 <- matrix(0, nrow = oos, ncol = ncol(data) - 1, dimnames = list(NULL, colnames(data)[-1]))
garch_n_fore_var2 <- garch_t_fore_var2 <- figarch_n_fore_var2 <- figarch_t_fore_var2 <- gas_n_fore_var2 <- gas_t_fore_var2 <- ms_n_fore_var2 <- ms_t_fore_var2 <- sv_n_fore_var2 <- sv_t_fore_var2 <- matrix(0, nrow = oos, ncol = ncol(data) - 1, dimnames = list(NULL, colnames(data)[-1]))
garch_n_fore_es1 <- garch_t_fore_es1 <- figarch_n_fore_es1 <- figarch_t_fore_es1 <- gas_n_fore_es1 <- gas_t_fore_es1 <- ms_n_fore_es1 <- ms_t_fore_es1 <- sv_n_fore_es1 <- sv_t_fore_es1 <- matrix(0, nrow = oos, ncol = ncol(data) - 1, dimnames = list(NULL, colnames(data)[-1]))
garch_n_fore_es2 <- garch_t_fore_es2 <- figarch_n_fore_es2 <- figarch_t_fore_es2 <- gas_n_fore_es2 <- gas_t_fore_es2 <- ms_n_fore_es2 <- ms_t_fore_es2 <- sv_n_fore_es2 <- sv_t_fore_es2 <- matrix(0, nrow = oos, ncol = ncol(data) - 1, dimnames = list(NULL, colnames(data)[-1]))

fhs_garch_n_fore_s2   <- fhs_garch_t_fore_s2    <- fhs_figarch_n_fore_s2    <- fhs_figarch_t_fore_s2    <- fhs_gas_n_fore_s2    <- fhs_gas_t_fore_s2    <- fhs_ms_n_fore_s2   <- fhs_ms_t_fore_s2   <- fhs_sv_n_fore_s2   <- fhs_sv_t_fore_s2   <- matrix(0, nrow = oos, ncol = ncol(data) - 1, dimnames = list(NULL, colnames(data)[-1]))
fhs_garch_n_fore_var1 <- fhs_garch_t_fore_var1  <- fhs_figarch_n_fore_var1  <- fhs_figarch_t_fore_var1  <- fhs_gas_n_fore_var1  <- fhs_gas_t_fore_var1  <- fhs_ms_n_fore_var1 <- fhs_ms_t_fore_var1 <- fhs_sv_n_fore_var1 <- fhs_sv_t_fore_var1 <- matrix(0, nrow = oos, ncol = ncol(data) - 1, dimnames = list(NULL, colnames(data)[-1]))
fhs_garch_n_fore_var2 <- fhs_garch_t_fore_var2  <- fhs_figarch_n_fore_var2  <- fhs_figarch_t_fore_var2  <- fhs_gas_n_fore_var2  <- fhs_gas_t_fore_var2  <- fhs_ms_n_fore_var2 <- fhs_ms_t_fore_var2 <- fhs_sv_n_fore_var2 <- fhs_sv_t_fore_var2 <- matrix(0, nrow = oos, ncol = ncol(data) - 1, dimnames = list(NULL, colnames(data)[-1]))
fhs_garch_n_fore_es1  <- fhs_garch_t_fore_es1   <- fhs_figarch_n_fore_es1   <- fhs_figarch_t_fore_es1   <- fhs_gas_n_fore_es1   <- fhs_gas_t_fore_es1   <- fhs_ms_n_fore_es1  <- fhs_ms_t_fore_es1  <- fhs_sv_n_fore_es1  <- fhs_sv_t_fore_es1  <- matrix(0, nrow = oos, ncol = ncol(data) - 1, dimnames = list(NULL, colnames(data)[-1]))
fhs_garch_n_fore_es2  <- fhs_garch_t_fore_es2   <- fhs_figarch_n_fore_es2   <- fhs_figarch_t_fore_es2   <- fhs_gas_n_fore_es2   <- fhs_gas_t_fore_es2   <- fhs_ms_n_fore_es2  <- fhs_ms_t_fore_es2  <- fhs_sv_n_fore_es2  <- fhs_sv_t_fore_es2  <- matrix(0, nrow = oos, ncol = ncol(data) - 1, dimnames = list(NULL, colnames(data)[-1]))



#plan(multicore, workers = parallel::detectCores() - 4)
plan(sequential)
for (i in 1:oos) {
  print(i)
  returns <- tail(data[i:(i + 2500 - 1), -1], ins) 
  mu <- apply(returns, 2, mean)
  returns_c <- scale(returns, scale = FALSE)

  for (j in 1:ncol(returns_c)) {
    if(abs(acf(returns_c[, j])$acf[2]) > 2/sqrt(nrow(returns_c))) {
      ar_fit <- ar.yw(returns_c[, j], 3, aic = TRUE, se.fit = FALSE)
      mu[j] <- mu[j] + as.numeric(predict(ar_fit, newdata = returns_c[, j], n.ahead = 1, se.fit = FALSE))
      returns_c[, j] <- ar_fit$resid
    }
  }
  
  # Fit Models
  garch_n_fit <- future_apply(returns_c, 2, function(x) ugarchfit(garch_spec_n, na.omit(x), solver = "hybrid"), future.seed = TRUE)
  garch_t_fit <- future_apply(returns_c, 2, function(x) ugarchfit(garch_spec_t, na.omit(x), solver = "hybrid"), future.seed = TRUE)
  figarch_n_fit <- future_apply(returns_c, 2, function(x) figarch_fit(figarch_spec_n, na.omit(x)), future.seed = TRUE)
  figarch_t_fit <- future_apply(returns_c, 2, function(x) figarch_fit(figarch_spec_t, na.omit(x)), future.seed = TRUE)
  gas_n_fit <- future_apply(returns_c, 2, function(x) gas_fit(gas_spec_n, na.omit(x)), future.seed = TRUE)
  gas_t_fit <- future_apply(returns_c, 2, function(x) gas_fit(gas_spec_t, na.omit(x)), future.seed = TRUE)
  ms_n_fit <- future_apply(returns_c, 2, function(x) msgarch_fit(ms_spec_n, na.omit(x)), future.seed = TRUE)
  ms_t_fit <- future_apply(returns_c, 2, function(x) msgarch_fit(ms_spec_t, na.omit(x)), future.seed = TRUE)
  sv_n_fit <- future_apply(returns_c, 2, function(x) estimate_parameters_n(as.numeric(na.omit(x))), future.seed = TRUE)
  sv_t_fit <- future_apply(returns_c, 2, function(x) estimate_parameters_t(as.numeric(na.omit(x))), future.seed = TRUE)

  # Forecasting S2, VaR and ES
  garch_n_fore <- do.call(cbind, future_lapply(garch_n_fit, function(x) {
      alpha <- c(0.01, 0.025)
      e <- x@fit$residuals / x@fit$sigma
      sigma_fore <- ugarchforecast(x, n.ahead = 1)@forecast$sigmaFor[1]
      var <- quantile(e, prob = alpha) * sigma_fore
      es <- sapply(var, function(v) mean(x@fit$residuals[x@fit$residuals < v]))
      var_par <- qdist(distribution = "norm", alpha, mu = 0, sigma = sigma_fore)
      es_par <- -sigma_fore * dnorm(qnorm(alpha)) / alpha
      out <- c(sigma_fore, var, es, var_par, es_par)  
    }, future.seed = TRUE))
  garch_t_fore <- do.call(cbind, future_lapply(garch_t_fit, function(x) {
    alpha <- c(0.01, 0.025)
    e <- x@fit$residuals / x@fit$sigma
    sigma_fore <- ugarchforecast(x, n.ahead = 1)@forecast$sigmaFor[1]
    var <- quantile(e, prob = alpha) * sigma_fore
    es <- sapply(var, function(v) mean(x@fit$residuals[x@fit$residuals < v]))
    nu <- coef(x)["shape"]
    k <- sqrt(nu / (nu - 2))
    var_par <- sigma_fore * qt(alpha, df = nu) / k
    es_par <- -sigma_fore / k * dt(qt(alpha, df = nu), df = nu) / alpha * (nu + qt(alpha, df = nu)^2) / (nu - 1)
    out <- c(sigma_fore, var, es, var_par, es_par)  
  }, future.seed = TRUE))
  figarch_n_fore <- do.call(cbind, future_lapply(figarch_n_fit, function(x) {
    alpha <- c(0.01, 0.025)
    e <- x@fit$residuals / x@fit$sigma
    sigma_fore <- ugarchforecast(x, n.ahead = 1)@forecast$sigmaFor[1]
    var <- quantile(e, prob =  alpha) * sigma_fore
    es <- sapply(var, function(v) mean(x@fit$residuals[x@fit$residuals < v]))
    var_par <- qdist(distribution = "norm", alpha, mu = 0, sigma = sigma_fore)
    es_par <- -sigma_fore * dnorm(qnorm(alpha)) / alpha
    out <- c(sigma_fore, var, es, var_par, es_par)  
  }, future.seed = TRUE))
  figarch_t_fore <- do.call(cbind, future_lapply(figarch_t_fit, function(x) {
    alpha <- c(0.01, 0.025)
    e <- x@fit$residuals / x@fit$sigma
    sigma_fore <- ugarchforecast(x, n.ahead = 1)@forecast$sigmaFor[1]
    var <- quantile(e, prob =  alpha) * sigma_fore
    es <- sapply(var, function(v) mean(x@fit$residuals[x@fit$residuals < v]))
    nu <- coef(x)["shape"]
    k <- sqrt(nu / (nu - 2))
    var_par <- sigma_fore * qt(alpha, df = nu) / k
    es_par <- -sigma_fore / k * dt(qt(alpha, df = nu), df = nu) / alpha * (nu + qt(alpha, df = nu)^2) / (nu - 1)
    out <- c(sigma_fore, var, es, var_par, es_par)
  }, future.seed = TRUE))
  gas_n_fore <- do.call(cbind, future_lapply(gas_n_fit, function(x) {
    alpha <- c(0.01, 0.025)
    sigma <- sqrt(x@GASDyn$mTheta[2, ])
    e <- x@Data$vY / sigma[1 : length(x@Data$vY)]
    sigma_fore <- tail(sigma, 1)
    var <- quantile(e, prob =  alpha) * sigma_fore
    es <- sapply(var, function(v) mean(x@Data$vY[x@Data$vY < v]))
    var_par <- as.numeric(tail(quantile(x, prob =  alpha), 1))
    es_par <- as.numeric(tail(ES(x, prob = alpha), 1))
    out <- c(sigma_fore, var, es, var_par, es_par)  
  }, future.seed = TRUE))
  gas_t_fore <- do.call(cbind, future_lapply(gas_t_fit, function(x) {
    alpha <- c(0.01, 0.025)
    sigma <- sqrt(x@GASDyn$mTheta[2, ] * x@GASDyn$mTheta[3, 1] /(x@GASDyn$mTheta[3, 1] - 2))
    e <- x@Data$vY / sigma[1 : length(x@Data$vY)]
    sigma_fore <- tail(sigma, 1)
    var <- quantile(e, prob = alpha) * sigma_fore
    es <- sapply(var, function(v) mean(x@Data$vY[x@Data$vY < v]))
    var_par <- as.numeric(tail(quantile(x, prob =  alpha), 1))
    es_par <- as.numeric(tail(ES(x, prob = alpha), 1))
    out <- c(sigma_fore, var, es, var_par, es_par)  
  }, future.seed = TRUE))
  ms_n_fore <- do.call(cbind, future_lapply(ms_n_fit, function(x) {
    alpha <- c(0.01, 0.025)
    e <- x$data / Volatility(x)
    sigma_fore <- predict(x, nahead = 1)$vol
    var <- quantile(e, prob = alpha) * sigma_fore
    es <- sapply(var, function(v) mean(x$data[x$data < v]))
    var_es <- Risk(x, nahead = 1, alpha = alpha)
    var_par <- var_es$VaR
    es_par <- var_es$ES
    out <- c(sigma_fore, var, es, var_par, es_par)  
  }, future.seed = TRUE))
  ms_t_fore <- do.call(cbind, future_lapply(ms_t_fit, function(x) {
    alpha <- c(0.01, 0.025)
    e <- x$data / Volatility(x)
    sigma_fore <- predict(x, nahead = 1)$vol
    var <- quantile(e, prob = alpha) * sigma_fore
    es <- sapply(var, function(v) mean(x$data[x$data < v]))
    var_es <- Risk(x, nahead = 1, alpha = alpha)
    var_par <- var_es$VaR
    es_par <- var_es$ES
    out <- c(sigma_fore, var, es, var_par, es_par)  
  }, future.seed = TRUE))
  sv_n_fore <- do.call(cbind, future_lapply(sv_n_fit, function(x) {
    alpha <- c(0.01, 0.025)
    sigma_y <-  as.numeric(summary(x)[1, 2])
    sigma <- exp(0.5 * x$rep$par.random) * sigma_y
    e <- x$data / sigma
    aux <- predict(x, steps = 1)
    sigma_fore <- median(aux$h_exp)
    var <- quantile(e, prob =  alpha) * sigma_fore
    es <- sapply(var, function(v) mean(x$data[x$data < v]))
    one_step_ahead_ret <- aux$y
    var_par <- quantile(one_step_ahead_ret, alpha)
    es_par <- sapply(var_par, function(v) mean(one_step_ahead_ret[one_step_ahead_ret < v]))
    out <- c(sigma_fore, var, es, var_par, es_par)  
  }, future.seed = TRUE))
  sv_t_fore <- do.call(cbind, future_lapply(sv_t_fit, function(x) {
    alpha <- c(0.01, 0.025)
    sigma_y <-  as.numeric(summary(x)[1, 2])
    sigma <- exp(0.5 * x$rep$par.random) * sigma_y
    e <- x$data / sigma
    aux <- predict(x, steps = 1)
    sigma_fore <- median(aux$h_exp)
    var <- quantile(e, prob = alpha) * sigma_fore
    es <- sapply(var, function(v) mean(x$data[x$data < v]))
    one_step_ahead_ret <- aux$y
    var_par <- quantile(one_step_ahead_ret, alpha)
    es_par <- sapply(var_par, function(v) mean(one_step_ahead_ret[one_step_ahead_ret < v]))
    out <- c(sigma_fore, var, es, var_par, es_par)  
  } , future.seed = TRUE))
  
  
  # Saving results
  garch_n_fore_s[i, ]    <-  garch_n_fore[1, ]
  garch_t_fore_s[i, ]    <-  garch_t_fore[1, ]
  figarch_n_fore_s[i, ]  <-  figarch_n_fore[1, ]
  figarch_t_fore_s[i, ]  <-  figarch_t_fore[1, ]
  gas_n_fore_s[i, ]      <-  gas_n_fore[1, ]
  gas_t_fore_s[i, ]      <-  gas_t_fore[1, ]
  sv_n_fore_s[i, ]       <-  sv_n_fore[1, ]
  sv_t_fore_s[i, ]       <-  sv_t_fore[1, ]
  ms_n_fore_s[i, ]       <-  ms_n_fore[1, ]
  ms_t_fore_s[i, ]       <-  ms_t_fore[1, ]
  
  fhs_garch_n_fore_var1[i, ]    <-  garch_n_fore[2, ] + mu
  fhs_garch_t_fore_var1[i, ]    <-  garch_t_fore[2, ] + mu
  fhs_figarch_n_fore_var1[i, ]  <-  figarch_n_fore[2, ] + mu
  fhs_figarch_t_fore_var1[i, ]  <-  figarch_t_fore[2, ] + mu
  fhs_gas_n_fore_var1[i, ]      <-  gas_n_fore[2, ] + mu
  fhs_gas_t_fore_var1[i, ]      <-  gas_t_fore[2, ] + mu
  fhs_sv_n_fore_var1[i, ]       <-  sv_n_fore[2, ] + mu
  fhs_sv_t_fore_var1[i, ]       <-  sv_t_fore[2, ] + mu
  fhs_ms_n_fore_var1[i, ]       <-  ms_n_fore[2, ] + mu
  fhs_ms_t_fore_var1[i, ]       <-  ms_t_fore[2, ] + mu
  
  fhs_garch_n_fore_var2[i, ]    <-  garch_n_fore[3, ] + mu
  fhs_garch_t_fore_var2[i, ]    <-  garch_t_fore[3, ] + mu
  fhs_figarch_n_fore_var2[i, ]  <-  figarch_n_fore[3, ] + mu
  fhs_figarch_t_fore_var2[i, ]  <-  figarch_t_fore[3, ] + mu
  fhs_gas_n_fore_var2[i, ]      <-  gas_n_fore[3, ] + mu
  fhs_gas_t_fore_var2[i, ]      <-  gas_t_fore[3, ] + mu
  fhs_sv_n_fore_var2[i, ]       <-  sv_n_fore[3, ] + mu
  fhs_sv_t_fore_var2[i, ]       <-  sv_t_fore[3, ] + mu
  fhs_ms_n_fore_var2[i, ]       <-  ms_n_fore[3, ] + mu
  fhs_ms_t_fore_var2[i, ]       <-  ms_t_fore[3, ] + mu
  
  fhs_garch_n_fore_es1[i, ]    <-  garch_n_fore[4, ] + mu
  fhs_garch_t_fore_es1[i, ]    <-  garch_t_fore[4, ] + mu
  fhs_figarch_n_fore_es1[i, ]  <-  figarch_n_fore[4, ] + mu
  fhs_figarch_t_fore_es1[i, ]  <-  figarch_t_fore[4, ] + mu
  fhs_gas_n_fore_es1[i, ]      <-  gas_n_fore[4, ] + mu
  fhs_gas_t_fore_es1[i, ]      <-  gas_t_fore[4, ] + mu
  fhs_sv_n_fore_es1[i, ]       <-  sv_n_fore[4, ] + mu
  fhs_sv_t_fore_es1[i, ]       <-  sv_t_fore[4, ] + mu
  fhs_ms_n_fore_es1[i, ]       <-  ms_n_fore[4, ] + mu
  fhs_ms_t_fore_es1[i, ]       <-  ms_t_fore[4, ] + mu
  
  fhs_garch_n_fore_es2[i, ]    <-  garch_n_fore[5, ] + mu
  fhs_garch_t_fore_es2[i, ]    <-  garch_t_fore[5, ] + mu
  fhs_figarch_n_fore_es2[i, ]  <-  figarch_n_fore[5, ] + mu
  fhs_figarch_t_fore_es2[i, ]  <-  figarch_t_fore[5, ] + mu
  fhs_gas_n_fore_es2[i, ]      <-  gas_n_fore[5, ] + mu
  fhs_gas_t_fore_es2[i, ]      <-  gas_t_fore[5, ] + mu
  fhs_sv_n_fore_es2[i, ]       <-  sv_n_fore[5, ] + mu
  fhs_sv_t_fore_es2[i, ]       <-  sv_t_fore[5, ] + mu
  fhs_ms_n_fore_es2[i, ]       <-  ms_n_fore[5, ] + mu
  fhs_ms_t_fore_es2[i, ]       <-  ms_t_fore[5, ] + mu
  
  garch_n_fore_var1[i, ]    <-  garch_n_fore[6, ] + mu
  garch_t_fore_var1[i, ]    <-  garch_t_fore[6, ] + mu
  figarch_n_fore_var1[i, ]  <-  figarch_n_fore[6, ] + mu
  figarch_t_fore_var1[i, ]  <-  figarch_t_fore[6, ] + mu
  gas_n_fore_var1[i, ]      <-  gas_n_fore[6, ] + mu
  gas_t_fore_var1[i, ]      <-  gas_t_fore[6, ] + mu
  sv_n_fore_var1[i, ]       <-  sv_n_fore[6, ] + mu
  sv_t_fore_var1[i, ]       <-  sv_t_fore[6, ] + mu
  ms_n_fore_var1[i, ]       <-  ms_n_fore[6, ] + mu
  ms_t_fore_var1[i, ]       <-  ms_t_fore[6, ] + mu
  
  garch_n_fore_var2[i, ]    <-  garch_n_fore[7, ] + mu
  garch_t_fore_var2[i, ]    <-  garch_t_fore[7, ] + mu
  figarch_n_fore_var2[i, ]  <-  figarch_n_fore[7, ] + mu
  figarch_t_fore_var2[i, ]  <-  figarch_t_fore[7, ] + mu
  gas_n_fore_var2[i, ]      <-  gas_n_fore[7, ] + mu
  gas_t_fore_var2[i, ]      <-  gas_t_fore[7, ] + mu
  sv_n_fore_var2[i, ]       <-  sv_n_fore[7, ] + mu
  sv_t_fore_var2[i, ]       <-  sv_t_fore[7, ] + mu
  ms_n_fore_var2[i, ]       <-  ms_n_fore[7, ] + mu
  ms_t_fore_var2[i, ]       <-  ms_t_fore[7, ] + mu
  
  garch_n_fore_es1[i, ]    <-  garch_n_fore[8, ] + mu
  garch_t_fore_es1[i, ]    <-  garch_t_fore[8, ] + mu
  figarch_n_fore_es1[i, ]  <-  figarch_n_fore[8, ] + mu
  figarch_t_fore_es1[i, ]  <-  figarch_t_fore[8, ] + mu
  gas_n_fore_es1[i, ]      <-  gas_n_fore[8, ] + mu
  gas_t_fore_es1[i, ]      <-  gas_t_fore[8, ] + mu
  sv_n_fore_es1[i, ]       <-  sv_n_fore[8, ] + mu
  sv_t_fore_es1[i, ]       <-  sv_t_fore[8, ] + mu
  ms_n_fore_es1[i, ]       <-  ms_n_fore[8, ] + mu
  ms_t_fore_es1[i, ]       <-  ms_t_fore[8, ] + mu
  
  garch_n_fore_es2[i, ]    <-  garch_n_fore[9, ] + mu
  garch_t_fore_es2[i, ]    <-  garch_t_fore[9, ] + mu
  figarch_n_fore_es2[i, ]  <-  figarch_n_fore[9, ] + mu
  figarch_t_fore_es2[i, ]  <-  figarch_t_fore[9, ] + mu
  gas_n_fore_es2[i, ]      <-  gas_n_fore[9, ] + mu
  gas_t_fore_es2[i, ]      <-  gas_t_fore[9, ] + mu
  sv_n_fore_es2[i, ]       <-  sv_n_fore[9, ] + mu
  sv_t_fore_es2[i, ]       <-  sv_t_fore[9, ] + mu
  ms_n_fore_es2[i, ]       <-  ms_n_fore[9, ] + mu
  ms_t_fore_es2[i, ]       <-  ms_t_fore[9, ] + mu
}


write.csv(garch_n_fore_s, "Empirical_Application/garch_n_fore_s_1000.csv")
write.csv(garch_t_fore_s, "Empirical_Application/garch_t_fore_s_1000.csv")
write.csv(figarch_n_fore_s, "Empirical_Application/figarch_n_fore_s_1000.csv")
write.csv(figarch_t_fore_s, "Empirical_Application/figarch_t_fore_s_1000.csv")
write.csv(gas_n_fore_s, "Empirical_Application/gas_n_fore_s_1000.csv")
write.csv(gas_t_fore_s, "Empirical_Application/gas_t_fore_s_1000.csv")
write.csv(ms_n_fore_s, "Empirical_Application/ms_n_fore_s_1000.csv")
write.csv(ms_t_fore_s, "Empirical_Application/ms_t_fore_s_1000.csv")
write.csv(sv_n_fore_s, "Empirical_Application/sv_n_fore_s_1000.csv")
write.csv(sv_t_fore_s, "Empirical_Application/sv_t_fore_s_1000.csv")

write.csv(garch_n_fore_var1, "Empirical_Application/garch_n_fore_var1_1000.csv")
write.csv(garch_t_fore_var1, "Empirical_Application/garch_t_fore_var1_1000.csv")
write.csv(figarch_n_fore_var1, "Empirical_Application/figarch_n_fore_var1_1000.csv")
write.csv(figarch_t_fore_var1, "Empirical_Application/figarch_t_fore_var1_1000.csv")
write.csv(gas_n_fore_var1, "Empirical_Application/gas_n_fore_var1_1000.csv")
write.csv(gas_t_fore_var1, "Empirical_Application/gas_t_fore_var1_1000.csv")
write.csv(ms_n_fore_var1, "Empirical_Application/ms_n_fore_var1_1000.csv")
write.csv(ms_t_fore_var1, "Empirical_Application/ms_t_fore_var1_1000.csv")
write.csv(sv_n_fore_var1, "Empirical_Application/sv_n_fore_var1_1000.csv")
write.csv(sv_t_fore_var1, "Empirical_Application/sv_t_fore_var1_1000.csv")
  
write.csv(fhs_garch_n_fore_var1, "Empirical_Application/fhs_garch_n_fore_var1_1000.csv")
write.csv(fhs_garch_t_fore_var1, "Empirical_Application/fhs_garch_t_fore_var1_1000.csv")
write.csv(fhs_figarch_n_fore_var1, "Empirical_Application/fhs_figarch_n_fore_var1_1000.csv")
write.csv(fhs_figarch_t_fore_var1, "Empirical_Application/fhs_figarch_t_fore_var1_1000.csv")
write.csv(fhs_gas_n_fore_var1, "Empirical_Application/fhs_gas_n_fore_var1_1000.csv")
write.csv(fhs_gas_t_fore_var1, "Empirical_Application/fhs_gas_t_fore_var1_1000.csv")
write.csv(fhs_ms_n_fore_var1, "Empirical_Application/fhs_ms_n_fore_var1_1000.csv")
write.csv(fhs_ms_t_fore_var1, "Empirical_Application/fhs_ms_t_fore_var1_1000.csv")
write.csv(fhs_sv_n_fore_var1, "Empirical_Application/fhs_sv_n_fore_var1_1000.csv")
write.csv(fhs_sv_t_fore_var1, "Empirical_Application/fhs_sv_t_fore_var1_1000.csv")

write.csv(garch_n_fore_var2, "Empirical_Application/garch_n_fore_var2_1000.csv")
write.csv(garch_t_fore_var2, "Empirical_Application/garch_t_fore_var2_1000.csv")
write.csv(figarch_n_fore_var2, "Empirical_Application/figarch_n_fore_var2_1000.csv")
write.csv(figarch_t_fore_var2, "Empirical_Application/figarch_t_fore_var2_1000.csv")
write.csv(gas_n_fore_var2, "Empirical_Application/gas_n_fore_var2_1000.csv")
write.csv(gas_t_fore_var2, "Empirical_Application/gas_t_fore_var2_1000.csv")
write.csv(ms_n_fore_var2, "Empirical_Application/ms_n_fore_var2_1000.csv")
write.csv(ms_t_fore_var2, "Empirical_Application/ms_t_fore_var2_1000.csv")
write.csv(sv_n_fore_var2, "Empirical_Application/sv_n_fore_var2_1000.csv")
write.csv(sv_t_fore_var2, "Empirical_Application/sv_t_fore_var2_1000.csv")

write.csv(fhs_garch_n_fore_var2, "Empirical_Application/fhs_garch_n_fore_var2_1000.csv")
write.csv(fhs_garch_t_fore_var2, "Empirical_Application/fhs_garch_t_fore_var2_1000.csv")
write.csv(fhs_figarch_n_fore_var2, "Empirical_Application/fhs_figarch_n_fore_var2_1000.csv")
write.csv(fhs_figarch_t_fore_var2, "Empirical_Application/fhs_figarch_t_fore_var2_1000.csv")
write.csv(fhs_gas_n_fore_var2, "Empirical_Application/fhs_gas_n_fore_var2_1000.csv")
write.csv(fhs_gas_t_fore_var2, "Empirical_Application/fhs_gas_t_fore_var2_1000.csv")
write.csv(fhs_ms_n_fore_var2, "Empirical_Application/fhs_ms_n_fore_var2_1000.csv")
write.csv(fhs_ms_t_fore_var2, "Empirical_Application/fhs_ms_t_fore_var2_1000.csv")
write.csv(fhs_sv_n_fore_var2, "Empirical_Application/fhs_sv_n_fore_var2_1000.csv")
write.csv(fhs_sv_t_fore_var2, "Empirical_Application/fhs_sv_t_fore_var2_1000.csv")

write.csv(garch_n_fore_es1, "Empirical_Application/garch_n_fore_es1_1000.csv")
write.csv(garch_t_fore_es1, "Empirical_Application/garch_t_fore_es1_1000.csv")
write.csv(figarch_n_fore_es1, "Empirical_Application/figarch_n_fore_es1_1000.csv")
write.csv(figarch_t_fore_es1, "Empirical_Application/figarch_t_fore_es1_1000.csv")
write.csv(gas_n_fore_es1, "Empirical_Application/gas_n_fore_es1_1000.csv")
write.csv(gas_t_fore_es1, "Empirical_Application/gas_t_fore_es1_1000.csv")
write.csv(ms_n_fore_es1, "Empirical_Application/ms_n_fore_es1_1000.csv")
write.csv(ms_t_fore_es1, "Empirical_Application/ms_t_fore_es1_1000.csv")
write.csv(sv_n_fore_es1, "Empirical_Application/sv_n_fore_es1_1000.csv")
write.csv(sv_t_fore_es1, "Empirical_Application/sv_t_fore_es1_1000.csv")

write.csv(fhs_garch_n_fore_es1, "Empirical_Application/fhs_garch_n_fore_es1_1000.csv")
write.csv(fhs_garch_t_fore_es1, "Empirical_Application/fhs_garch_t_fore_es1_1000.csv")
write.csv(fhs_figarch_n_fore_es1, "Empirical_Application/fhs_figarch_n_fore_es1_1000.csv")
write.csv(fhs_figarch_t_fore_es1, "Empirical_Application/fhs_figarch_t_fore_es1_1000.csv")
write.csv(fhs_gas_n_fore_es1, "Empirical_Application/fhs_gas_n_fore_es1_1000.csv")
write.csv(fhs_gas_t_fore_es1, "Empirical_Application/fhs_gas_t_fore_es1_1000.csv")
write.csv(fhs_ms_n_fore_es1, "Empirical_Application/fhs_ms_n_fore_es1_1000.csv")
write.csv(fhs_ms_t_fore_es1, "Empirical_Application/fhs_ms_t_fore_es1_1000.csv")
write.csv(fhs_sv_n_fore_es1, "Empirical_Application/fhs_sv_n_fore_es1_1000.csv")
write.csv(fhs_sv_t_fore_es1, "Empirical_Application/fhs_sv_t_fore_es1_1000.csv")

write.csv(garch_n_fore_es2, "Empirical_Application/garch_n_fore_es2_1000.csv")
write.csv(garch_t_fore_es2, "Empirical_Application/garch_t_fore_es2_1000.csv")
write.csv(figarch_n_fore_es2, "Empirical_Application/figarch_n_fore_es2_1000.csv")
write.csv(figarch_t_fore_es2, "Empirical_Application/figarch_t_fore_es2_1000.csv")
write.csv(gas_n_fore_es2, "Empirical_Application/gas_n_fore_es2_1000.csv")
write.csv(gas_t_fore_es2, "Empirical_Application/gas_t_fore_es2_1000.csv")
write.csv(ms_n_fore_es2, "Empirical_Application/ms_n_fore_es2_1000.csv")
write.csv(ms_t_fore_es2, "Empirical_Application/ms_t_fore_es2_1000.csv")
write.csv(sv_n_fore_es2, "Empirical_Application/sv_n_fore_es2_1000.csv")
write.csv(sv_t_fore_es2, "Empirical_Application/sv_t_fore_es2_1000.csv")

write.csv(fhs_garch_n_fore_es2, "Empirical_Application/fhs_garch_n_fore_es2_1000.csv")
write.csv(fhs_garch_t_fore_es2, "Empirical_Application/fhs_garch_t_fore_es2_1000.csv")
write.csv(fhs_figarch_n_fore_es2, "Empirical_Application/fhs_figarch_n_fore_es2_1000.csv")
write.csv(fhs_figarch_t_fore_es2, "Empirical_Application/fhs_figarch_t_fore_es2_1000.csv")
write.csv(fhs_gas_n_fore_es2, "Empirical_Application/fhs_gas_n_fore_es2_1000.csv")
write.csv(fhs_gas_t_fore_es2, "Empirical_Application/fhs_gas_t_fore_es2_1000.csv")
write.csv(fhs_ms_n_fore_es2, "Empirical_Application/fhs_ms_n_fore_es2_1000.csv")
write.csv(fhs_ms_t_fore_es2, "Empirical_Application/fhs_ms_t_fore_es2_1000.csv")
write.csv(fhs_sv_n_fore_es2, "Empirical_Application/fhs_sv_n_fore_es2_1000.csv")
write.csv(fhs_sv_t_fore_es2, "Empirical_Application/fhs_sv_t_fore_es2_1000.csv")

#write.csv(data[(1 + 2500):nrow(data), ], "Empirical_Application/data_oos.csv")