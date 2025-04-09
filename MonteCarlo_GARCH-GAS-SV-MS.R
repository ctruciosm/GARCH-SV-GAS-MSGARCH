################################################################################
#####                   Monte Carlo Simulation                             #####
################################################################################
args <- commandArgs(TRUE)
library(rugarch)
library(GAS)
library(stochvol)
library(stochvolTMB)
library(MSGARCH)
source("DGPs.R")
source("Utils_GARCH-GAS-SV.R")


## Setting values
mc <- 1000

# Parse arguments if they exist
if (length(args) > 0) {
  for (i in 1:length(args)) {
    arg_parts <- strsplit(args[[i]], "=")[[1]]
    if (length(arg_parts) == 2) {
      param_name <- arg_parts[1]
      param_value <- arg_parts[2]
      if (param_name == "n") n <- as.integer(param_value)
      if (param_name == "type") type <- param_value
      if (param_name == "outliers") outliers <- param_value
    }
  }
} else {
  n <- 2500
  type <- "BR"
  outliers <- "FALSE"
}


# Estimated a BR time series to obtain parameters
garch_spec_n <- ugarchspec(
  variance.model = list(model = 'sGARCH', garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "norm"
)
garch_spec_t <- ugarchspec(
  variance.model = list(model = 'sGARCH', garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "std"
)
gas_spec_n <- UniGASSpec(
  Dist = "norm",
  ScalingType = "Identity",
  GASPar = list(locate = FALSE, scale = TRUE, shape = FALSE)
)
gas_spec_t <- UniGASSpec(
  Dist = "std",
  ScalingType = "Identity",
  GASPar = list(locate = FALSE, scale = TRUE, shape = FALSE)
)
ms_spec_n <- CreateSpec(
  variance.spec = list(model = c("sGARCH", "sGARCH")),
  switch.spec = list(do.mix = FALSE),
  distribution.spec = list(distribution = c("norm", "norm"))
)
ms_spec_t <- CreateSpec(
  variance.spec = list(model = c("sGARCH", "sGARCH")),
  switch.spec = list(do.mix = FALSE),
  distribution.spec = list(distribution = c("std", "std")),
  constraint.spec = list(regime.const = c("nu"))
)

if (type == "BR") {
  # data <- logret(read.csv("./Data/precos_diarios_ibrx.csv")[, "PETR4"]) * 100
  # data <- logret(read.csv("./Data/BTCUSDT_1d.csv")[, "Close"]) * 100
  garch_params <- c(0.18, 0.09, 0.89)
  gas_params <- c(0.03, 0.22, 0.98)
  sv_params <- c(1.74, 0.97, 0.17)
  ms_params <- c(0.005, 0.025, 0.95, 0.1, 0.25, 0.70)
  P <-  matrix(c(0.75, 0.30, 0.25, 0.70), 2, 2, byrow = TRUE)  #Haas
}
if (type == "US") {
  # data <- readxl::read_xlsx("Data/economatica_nyse_diario.xlsx", skip = 3, col_types = c("text", rep("numeric", 1420)), na = c("-", " ", "NA"))
  # colnames(data) <- stringr::str_replace(colnames(data), "Fechamento\najust p/ prov\nEm moeda orig\n", "")
  # data <- logret(na.omit(data$MCS[5000:8000])) * 100
  garch_params <- c(0.37, 0.14, 0.77)
  gas_params <- c(0.06, 0.34, 0.92)
  sv_params <- c(1.15, 0.90, 0.36)
  ms_params <- c(0.01, 0.16, 0.30, 0.18, 0.46, 0.20)
  P <- matrix(c(0.98, 0.05, 0.02, 0.95), 2, 2, byrow = TRUE) #Haas
}

true_vols_n <- matrix(NA, ncol = 4, nrow = mc)
fore_vols_n <- matrix(NA, ncol = 40, nrow = mc)
true_vols_t <- matrix(NA, ncol = 4, nrow = mc)
fore_vols_t <- matrix(NA, ncol = 40, nrow = mc)
for (i in 1:mc) {
  set.seed(i + 123)
  print(i)

  # Simulate DGPs
  garch_sim_n <- garch_sim(2500 + 1, garch_params, "norm")
  gas_sim_n <- gas_sim(2500 + 1, gas_params, "norm")
  sv_sim_n <- sv_sim(2500 + 1, sv_params, "norm")
  ms_sim_n <- msgarch_sim(2500 + 1, ms_params, "norm", P)

  garch_sim_t <- garch_sim(2500 + 1, c(garch_params, 7), "std")
  gas_sim_t <- gas_sim(2500 + 1, c(gas_params, -2.6625878), "std")
  sv_sim_t <- sv_sim(2500 + 1, c(sv_params, 7), "std")
  ms_sim_t <- msgarch_sim(2500 + 1, c(ms_params, 7), "std", P)

  r_garch_sim_n <- tail(garch_sim_n$returns, n + 1)[1:n]
  r_gas_sim_n <- tail(gas_sim_n$returns, n + 1)[1:n]
  r_sv_sim_n <- tail(sv_sim_n$returns, n + 1)[1:n]
  r_ms_sim_n <- tail(ms_sim_n$returns, n + 1)[1:n]

  r_garch_sim_t <- tail(garch_sim_t$returns, n + 1)[1:n]
  r_gas_sim_t <- tail(gas_sim_t$returns, n + 1)[1:n]
  r_sv_sim_t <- tail(sv_sim_t$returns, n + 1)[1:n]
  r_ms_sim_t <- tail(ms_sim_t$returns, n + 1)[1:n]

  # Contaminate the series with AO
  if (outliers == "TRUE") {
    outlier_position <- floor(runif(1, n - 22, n)) + 1

    r_garch_sim_n[outlier_position] <- r_garch_sim_n[outlier_position] +
      sign(r_garch_sim_n[outlier_position]) * 5 * sd(r_garch_sim_n)
    r_gas_sim_n[outlier_position] <- r_gas_sim_n[outlier_position] +
      sign(r_gas_sim_n[outlier_position]) * 5 * sd(r_gas_sim_n)
    r_sv_sim_n[outlier_position] <- r_sv_sim_n[outlier_position] +
      sign(r_sv_sim_n[outlier_position]) * 5 * sd(r_sv_sim_n)
    r_ms_sim_n[outlier_position] <- r_ms_sim_n[outlier_position] +
      sign(r_ms_sim_n[outlier_position]) * 5 * sd(r_ms_sim_n)

    r_garch_sim_t[outlier_position] <- r_garch_sim_t[outlier_position] +
      sign(r_garch_sim_t[outlier_position]) * 5 * sd(r_garch_sim_t)
    r_gas_sim_t[outlier_position] <- r_gas_sim_t[outlier_position] +
      sign(r_gas_sim_t[outlier_position]) * 5 * sd(r_gas_sim_t)
    r_sv_sim_t[outlier_position] <- r_sv_sim_t[outlier_position] +
      sign(r_sv_sim_t[outlier_position]) * 5 * sd(r_sv_sim_t)
    r_ms_sim_t[outlier_position] <- r_ms_sim_t[outlier_position] +
      sign(r_ms_sim_t[outlier_position]) * 5 * sd(r_ms_sim_t)
  }

  # Fit the models
  garch_n_garch_n <- as.numeric(sigma(ugarchforecast(ugarchfit(garch_spec_n, r_garch_sim_n, solver = "hybrid"), n.ahead = 1)))
  garch_n_gas_n <- sqrt(UniGASFor(gasfit(gas_spec_n, r_garch_sim_n), H = 1)@Forecast$PointForecast[, 2])
  garch_n_sv_n_b <- median(predvola(predict(svsample(r_garch_sim_n, quiet = TRUE), 1)))
  garch_n_sv_n <- median(predict(estimate_parameters(r_garch_sim_n, model = "gaussian", silent = TRUE), steps = 1)$h_exp)
  garch_n_ms_n <- predict(msgarchfit(ms_spec_n, r_garch_sim_n), nahead = 1)$vol
  garch_n_garch_t <- as.numeric(sigma(ugarchforecast(ugarchfit(garch_spec_t, r_garch_sim_n, solver = "hybrid"), n.ahead = 1)))
  aux_gas_t <- gasfit(gas_spec_t, r_garch_sim_n)
    nu <- aux_gas_t@GASDyn$mTheta[3, 1]
  garch_n_gas_t <- sqrt(UniGASFor(aux_gas_t, H = 1)@Forecast$PointForecast[, 2]) * sqrt(nu / (nu - 2))
  garch_n_sv_t_b <- median(predvola(predict(svtsample(r_garch_sim_n, quiet = TRUE), 1)))
  garch_n_sv_t <- median(predict(estimate_parameters_t(r_garch_sim_n), steps = 1)$h_exp)
  garch_n_ms_t <- predict(msgarchfit(ms_spec_t, r_garch_sim_n), nahead = 1)$vol
  
  garch_t_garch_n <- as.numeric(sigma(ugarchforecast(ugarchfit(garch_spec_n, r_garch_sim_t, solver = "hybrid"), n.ahead = 1)))
  garch_t_gas_n <- sqrt(UniGASFor(gasfit(gas_spec_n, r_garch_sim_t), H = 1)@Forecast$PointForecast[, 2])
  garch_t_sv_n_b <- median(predvola(predict(svsample(r_garch_sim_t, quiet = TRUE), 1)))
  garch_t_sv_n <- median(predict(estimate_parameters(r_garch_sim_t, model = "gaussian", silent = TRUE), steps = 1)$h_exp)
  garch_t_ms_n <- predict(msgarchfit(ms_spec_n, r_garch_sim_t), nahead = 1)$vol
  garch_t_garch_t <- as.numeric(sigma(ugarchforecast(ugarchfit(garch_spec_t, r_garch_sim_t, solver = "hybrid"), n.ahead = 1)))
  aux_gas_t <- gasfit(gas_spec_t, r_garch_sim_t)
    nu <- aux_gas_t@GASDyn$mTheta[3, 1]
  garch_t_gas_t <-  sqrt(UniGASFor(aux_gas_t, H = 1)@Forecast$PointForecast[, 2]) * sqrt(nu / (nu - 2))
  garch_t_sv_t_b <- median(predvola(predict(svtsample(r_garch_sim_t, quiet = TRUE), 1)))
  garch_t_sv_t <- median(predict(estimate_parameters_t(r_garch_sim_t), steps = 1)$h_exp)
  garch_t_ms_t <- predict(msgarchfit(ms_spec_t, r_garch_sim_t), nahead = 1)$vol

  gas_n_garch_n <- as.numeric(sigma(ugarchforecast(ugarchfit(garch_spec_n, r_gas_sim_n, solver = "hybrid"), n.ahead = 1)))
  gas_n_gas_n <- sqrt(UniGASFor(gasfit(gas_spec_n, r_gas_sim_n), H = 1)@Forecast$PointForecast[, 2])
  gas_n_sv_n_b <- median(predvola(predict(svsample(r_gas_sim_n, quiet = TRUE), 1)))
  gas_n_sv_n <- median(predict(estimate_parameters(r_gas_sim_n, model = "gaussian", silent = TRUE), steps = 1)$h_exp)
  gas_n_ms_n <- predict(msgarchfit(ms_spec_n, r_gas_sim_n), nahead = 1)$vol
  gas_n_garch_t <- as.numeric(sigma(ugarchforecast(ugarchfit(garch_spec_t, r_gas_sim_n, solver = "hybrid"), n.ahead = 1)))
  aux_gas_t <- gasfit(gas_spec_t, r_gas_sim_n)
    nu <- aux_gas_t@GASDyn$mTheta[3, 1]
  gas_n_gas_t <-  sqrt(UniGASFor(aux_gas_t, H = 1)@Forecast$PointForecast[, 2]) * sqrt(nu / (nu - 2))
  gas_n_sv_t_b <- median(predvola(predict(svtsample(r_gas_sim_n, quiet = TRUE), 1)))
  gas_n_sv_t <- median(predict(estimate_parameters_t(r_gas_sim_n), steps = 1)$h_exp)
  gas_n_ms_t <- predict(msgarchfit(ms_spec_t, r_gas_sim_n), nahead = 1)$vol

  gas_t_garch_n <- as.numeric(sigma(ugarchforecast(ugarchfit(garch_spec_n, r_gas_sim_t, solver = "hybrid"), n.ahead = 1)))
  gas_t_gas_n <- sqrt(UniGASFor(gasfit(gas_spec_n, r_gas_sim_t),H = 1)@Forecast$PointForecast[, 2])
  gas_t_sv_n_b <- median(predvola(predict(svsample(r_gas_sim_t, quiet = TRUE), 1)))
  gas_t_sv_n <-  median(predict(estimate_parameters(r_gas_sim_t, model = "gaussian", silent = TRUE), steps = 1)$h_exp)
  gas_t_ms_n <- predict(msgarchfit(ms_spec_n, r_gas_sim_t), nahead = 1)$vol
  gas_t_garch_t <- as.numeric(sigma(ugarchforecast(ugarchfit(garch_spec_t, r_gas_sim_t, solver = "hybrid"), n.ahead = 1)))
  aux_gas_t <- gasfit(gas_spec_t, r_gas_sim_t)
    nu <- aux_gas_t@GASDyn$mTheta[3, 1]
  gas_t_gas_t <-  sqrt(UniGASFor(aux_gas_t, H = 1)@Forecast$PointForecast[, 2]) * sqrt(nu / (nu - 2))
  gas_t_sv_t_b <- median(predvola(predict(svtsample(r_gas_sim_t, quiet = TRUE), 1)))
  gas_t_sv_t <- median(predict(estimate_parameters_t(r_gas_sim_t), steps = 1)$h_exp)
  gas_t_ms_t <- predict(msgarchfit(ms_spec_t, r_gas_sim_t), nahead = 1)$vol

  sv_n_garch_n <- as.numeric(sigma(ugarchforecast(ugarchfit(garch_spec_n, r_sv_sim_n, solver = "hybrid"), n.ahead = 1)))
  sv_n_gas_n <- sqrt(UniGASFor(gasfit(gas_spec_n, r_sv_sim_n), H = 1)@Forecast$PointForecast[, 2])
  sv_n_sv_n_b <- median(predvola(predict(svsample(r_sv_sim_n, quiet = TRUE), 1)))
  sv_n_sv_n <- median(predict(estimate_parameters(r_sv_sim_n, model = "gaussian", silent = TRUE), steps = 1)$h_exp)
  sv_n_ms_n <- predict(msgarchfit(ms_spec_n, r_sv_sim_n), nahead = 1)$vol
  sv_n_garch_t <- as.numeric(sigma(ugarchforecast(ugarchfit(garch_spec_t, r_sv_sim_n, solver = "hybrid"), n.ahead = 1)))
  aux_gas_t <- gasfit(gas_spec_t, r_sv_sim_n)
    nu <- aux_gas_t@GASDyn$mTheta[3, 1]
  sv_n_gas_t <-  sqrt(UniGASFor(aux_gas_t, H = 1)@Forecast$PointForecast[, 2]) * sqrt(nu / (nu - 2))
  sv_n_sv_t_b <- median(predvola(predict(svtsample(r_sv_sim_n, quiet = TRUE), 1)))
  sv_n_sv_t <- median(predict(estimate_parameters_t(r_sv_sim_n), steps = 1)$h_exp)
  sv_n_ms_t <- predict(msgarchfit(ms_spec_t, r_sv_sim_n), nahead = 1)$vol

  sv_t_garch_n <- as.numeric(sigma(ugarchforecast(ugarchfit(garch_spec_n, r_sv_sim_t, solver = "hybrid"), n.ahead = 1)))
  sv_t_gas_n <- sqrt(UniGASFor(gasfit(gas_spec_n, r_sv_sim_t), H = 1)@Forecast$PointForecast[, 2])
  sv_t_sv_n_b <- median(predvola(predict(svsample(r_sv_sim_t, quiet = TRUE), 1)))
  sv_t_sv_n <- median(predict(estimate_parameters(r_sv_sim_t, model = "gaussian", silent = TRUE), steps = 1)$h_exp)
  sv_t_ms_n <- predict(msgarchfit(ms_spec_n, r_sv_sim_t), nahead = 1)$vol
  sv_t_garch_t <- as.numeric(sigma(ugarchforecast(ugarchfit(garch_spec_t, r_sv_sim_t, solver = "hybrid"), n.ahead = 1)))
  aux_gas_t <- gasfit(gas_spec_t, r_sv_sim_t)
    nu <- aux_gas_t@GASDyn$mTheta[3, 1]
  sv_t_gas_t <-  sqrt(UniGASFor(aux_gas_t, H = 1)@Forecast$PointForecast[, 2]) * sqrt(nu / (nu - 2))
  sv_t_sv_t_b <- median(predvola(predict(svtsample(r_sv_sim_t, quiet = TRUE), 1)))
  sv_t_sv_t <- median(predict(estimate_parameters_t(r_sv_sim_t), steps = 1)$h_exp)
  sv_t_ms_t <- predict(msgarchfit(ms_spec_t, r_sv_sim_t), nahead = 1)$vol

  ms_n_garch_n <- as.numeric(sigma(ugarchforecast(ugarchfit(garch_spec_n, r_ms_sim_n, solver = "hybrid"), n.ahead = 1)))
  ms_n_gas_n <- sqrt(UniGASFor(gasfit(gas_spec_n, r_ms_sim_n), H = 1)@Forecast$PointForecast[, 2])
  ms_n_sv_n_b <- median(predvola(predict(svsample(r_ms_sim_n, quiet = TRUE), 1)))
  ms_n_sv_n <- median(predict(estimate_parameters(r_ms_sim_n, model = "gaussian", silent = TRUE), steps = 1)$h_exp)
  ms_n_ms_n <- predict(msgarchfit(ms_spec_n, r_ms_sim_n), nahead = 1)$vol
  ms_n_garch_t <- as.numeric(sigma(ugarchforecast(ugarchfit(garch_spec_t, r_ms_sim_n, solver = "hybrid"), n.ahead = 1)))
  aux_gas_t <- gasfit(gas_spec_t, r_ms_sim_n)
    nu <- aux_gas_t@GASDyn$mTheta[3, 1]
  ms_n_gas_t <-  sqrt(UniGASFor(aux_gas_t, H = 1)@Forecast$PointForecast[, 2]) * sqrt(nu / (nu - 2))
  ms_n_sv_t_b <- median(predvola(predict(svtsample(r_ms_sim_n, quiet = TRUE), 1)))
  ms_n_sv_t <- median(predict(estimate_parameters_t(r_ms_sim_n), steps = 1)$h_exp)
  ms_n_ms_t <- predict(msgarchfit(ms_spec_t, r_ms_sim_n), nahead = 1)$vol

  ms_t_garch_n <- as.numeric(sigma(ugarchforecast(ugarchfit(garch_spec_n, r_ms_sim_t, solver = "hybrid"), n.ahead = 1)))
  ms_t_gas_n <- sqrt(UniGASFor(gasfit(gas_spec_n, r_ms_sim_t), H = 1)@Forecast$PointForecast[, 2])
  ms_t_sv_n_b <- median(predvola(predict(svsample(r_ms_sim_t, quiet = TRUE), 1)))
  ms_t_sv_n <- median(predict(estimate_parameters(r_ms_sim_t, model = "gaussian", silent = TRUE), steps = 1)$h_exp)
  ms_t_ms_n <- predict(msgarchfit(ms_spec_n, r_ms_sim_t), nahead = 1)$vol
  ms_t_garch_t <- as.numeric(sigma(ugarchforecast(ugarchfit(garch_spec_t, r_ms_sim_t, solver = "hybrid"), n.ahead = 1)))
  aux_gas_t <- gasfit(gas_spec_t, r_ms_sim_t)
    nu <- aux_gas_t@GASDyn$mTheta[3, 1]
  ms_t_gas_t <-  sqrt(UniGASFor(aux_gas_t, H = 1)@Forecast$PointForecast[, 2]) * sqrt(nu / (nu - 2))
  ms_t_sv_t_b <- median(predvola(predict(svtsample(r_ms_sim_t, quiet = TRUE), 1)))
  ms_t_sv_t <- median(predict(estimate_parameters_t(r_ms_sim_t), steps = 1)$h_exp)
  ms_t_ms_t <- predict(msgarchfit(ms_spec_t, r_ms_sim_t), nahead = 1)$vol

  true_vols_n[i, ] <- c(tail(garch_sim_n$volatility, 1), tail(gas_sim_n$volatility, 1), tail(sv_sim_n$volatility, 1), tail(ms_sim_n$volatility[, 3], 1))
  true_vols_t[i, ] <- c(tail(garch_sim_t$volatility, 1), tail(gas_sim_t$volatility, 1), tail(sv_sim_t$volatility, 1), tail(ms_sim_t$volatility[, 3], 1))
  fore_vols_n[i, ] <- c(garch_n_garch_n, garch_n_garch_t, garch_n_gas_n, garch_n_gas_t, garch_n_sv_n, garch_n_sv_t, garch_n_sv_n_b, garch_n_sv_t_b, garch_n_ms_n, garch_n_ms_t,
                        gas_n_garch_n, gas_n_garch_t, gas_n_gas_n, gas_n_gas_t, gas_n_sv_n, gas_n_sv_t, gas_n_sv_n_b, gas_n_sv_t_b, gas_n_ms_n, gas_n_ms_t,
                        sv_n_garch_n, sv_n_garch_t, sv_n_gas_n, sv_n_gas_t, sv_n_sv_n, sv_n_sv_t, sv_n_sv_n_b, sv_n_sv_t_b, sv_n_ms_n, sv_n_ms_t,
                        ms_n_garch_n, ms_n_garch_t, ms_n_gas_n, ms_n_gas_t, ms_n_sv_n, ms_n_sv_t, ms_n_sv_n_b, ms_n_sv_t_b, ms_n_ms_n, ms_n_ms_t)
  fore_vols_t[i, ] <- c(garch_t_garch_n, garch_t_garch_t, garch_t_gas_n, garch_t_gas_t, garch_t_sv_n, garch_t_sv_t, garch_t_sv_n_b, garch_t_sv_t_b, garch_t_ms_n, garch_t_ms_t,
                        gas_t_garch_n, gas_t_garch_t, gas_t_gas_n, gas_t_gas_t, gas_t_sv_n, gas_t_sv_t,  gas_t_sv_n_b, gas_t_sv_t_b, gas_t_ms_n, gas_t_ms_t,
                        sv_t_garch_n, sv_t_garch_t, sv_t_gas_n, sv_t_gas_t, sv_t_sv_n, sv_t_sv_t, sv_t_sv_n_b, sv_t_sv_t_b, sv_t_ms_n, sv_t_ms_t,
                        ms_t_garch_n, ms_t_garch_t, ms_t_gas_n, ms_t_gas_t, ms_t_sv_n, ms_t_sv_t, ms_t_sv_n_b, ms_t_sv_t_b, ms_t_ms_n, ms_t_ms_t)
}

volatilities_n <- cbind(true_vols_n, fore_vols_n)
volatilities_t <- cbind(true_vols_t, fore_vols_t)
colnames(volatilities_n) <- c("garch", "gas", "sv", "ms",
  "garch_n_garch_n", "garch_n_garch_t", "garch_n_gas_n", "garch_n_gas_t", "garch_n_sv_n", "garch_n_sv_t", "garch_n_sv_n_b", "garch_n_sv_t_b", "garch_n_ms_n", "garch_n_ms_t",
  "gas_n_garch_n", "gas_n_garch_t", "gas_n_gas_n", "gas_n_gas_t", "gas_n_sv_n", "gas_n_sv_t", "gas_n_sv_n_b", "gas_n_sv_t_b", "gas_n_ms_n", "gas_n_ms_t",
  "sv_n_garch_n", "sv_n_garch_t", "sv_n_gas_n", "sv_n_gas_t", "sv_n_sv_n", "sv_n_sv_t", "sv_n_sv_n_b", "sv_n_sv_t_b", "sv_n_ms_n", "sv_n_ms_t",
  "ms_n_garch_n", "ms_n_garch_t", "ms_n_gas_n", "ms_n_gas_t", "ms_n_sv_n", "ms_n_sv_t", "ms_n_sv_n_b", "ms_n_sv_t_b", "ms_n_ms_n", "ms_n_ms_t")
colnames(volatilities_t) <- c("garch", "gas", "sv", "ms",
  "garch_t_garch_n", "garch_t_garch_t", "garch_t_gas_n", "garch_t_gas_t", "garch_t_sv_n", "garch_t_sv_t", "garch_t_sv_n_b", "garch_t_sv_t_b", "garch_t_ms_n", "garch_t_ms_t",
  "gas_t_garch_n", "gas_t_garch_t", "gas_t_gas_n", "gas_t_gas_t", "gas_t_sv_n", "gas_t_sv_t", "gas_t_sv_n_b", "gas_t_sv_t_b", "gas_t_ms_n", "gas_t_ms_t",
  "sv_t_garch_n", "sv_t_garch_t", "sv_t_gas_n", "sv_t_gas_t", "sv_t_sv_n", "sv_t_sv_t", "sv_t_sv_n_b", "sv_t_sv_t_b","sv_t_ms_n", "sv_t_ms_t",
  "ms_t_garch_n", "ms_t_garch_t", "ms_t_gas_n", "ms_t_gas_t", "ms_t_sv_n", "ms_t_sv_t", "ms_t_sv_n_b", "ms_t_sv_t_b", "ms_t_ms_n", "ms_t_ms_t")

write.csv(volatilities_n, paste0("volatilities_", n, "_norm_", outliers, "_", type, ".csv"))
write.csv(volatilities_t, paste0("volatilities_", n, "_std_", outliers, "_", type, ".csv"))
