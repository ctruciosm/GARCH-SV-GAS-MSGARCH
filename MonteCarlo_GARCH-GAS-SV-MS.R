################################################################################
#####                   Monte Carlo Simulation                             #####
################################################################################


## APENDIX: COMPARE RUGARCH VS FEGARCH PACKAGES 
## AND GAS VS BETATEGARCH 


args <- commandArgs(TRUE)
library(betategarch)
library(stochvol)
library(stochvolTMB)
library(dplyr)
library(stringr)
library(MSGARCH)
library(Rcpp)
library(fEGarch)
source("./DGPs.R")
source("./Utils_GARCH-GAS-SV.R")
sourceCpp("./utils.cpp")
mc <- 1000
n_dgps <- 5
n_max <- 5000


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
  n <- 5000
  type <- "US"
  outliers <- "FALSE"
}

if (type == "RND"){
  seed <- as.numeric(Sys.Date())
}


# Estimated a BR time series to obtain parameters
ms_spec_n <- CreateSpec(variance.spec = list(model = c("sGARCH", "sGARCH")), switch.spec = list(do.mix = FALSE), distribution.spec = list(distribution = c("norm", "norm")))
ms_spec_t <- CreateSpec(variance.spec = list(model = c("sGARCH", "sGARCH")), switch.spec = list(do.mix = FALSE), distribution.spec = list(distribution = c("std", "std")), constraint.spec = list(regime.const = c("nu")))

if (type == "BR") {
  # data <- logret(read.csv("./Data/precos_diarios_ibrx.csv")[, "PETR4"]) * 100
  # data <- logret(read.csv("./Data/BTCUSDT_1d.csv")[, "Close"]) * 100
  garch_params <- c(0.18, 0.09, 0.89)
  figarch_params <- c(0.72, 0.22, 0.55, 0.44)  #c(0.08, 0.2, 0.5, 0.4)
  dcs_params <- c(0.74, 0.97, 0.05)
  sv_params <- c(1.74, 0.97, 0.17)
  ms_params <- c(0.005, 0.025, 0.95, 0.1, 0.25, 0.70)
  P <-  matrix(c(0.75, 0.30, 0.25, 0.70), 2, 2, byrow = TRUE) 
}
if (type == "US") {
  # data <- readxl::read_xlsx("Data/economatica_nyse_diario.xlsx", skip = 3, col_types = c("text", rep("numeric", 1420)), na = c("-", " ", "NA"))
  # colnames(data) <- stringr::str_replace(colnames(data), "Fechamento\najust p/ prov\nEm moeda orig\n", "")
  # data <- logret(na.omit(data$AAP[3000:10000])) * 100
  garch_params <- c(0.37, 0.14, 0.77)
  figarch_params <- c(0.66, 0.54, 0.70, 0.27) #c(0.02, 0.15, 0.62, 0.48)
  dcs_params <- c(0.42, 0.92, 0.08)
  sv_params <- c(1.15, 0.90, 0.36)
  ms_params <- c(0.01, 0.16, 0.30, 0.18, 0.46, 0.20)
  P <- matrix(c(0.98, 0.05, 0.02, 0.95), 2, 2, byrow = TRUE)
}

if (type == "RND") {
  set.seed(seed)
  if (runif(1)> 0.5) {
    aux <- read.csv("./Data/precos_diarios_ibrx.csv")[, -c(1, 2)]
    ii <- sample(1:ncol(aux), 1)
    data <- stochvol::logret(as.numeric(na.omit(aux[, ii]))) * 100
  } else {
    aux <- readxl::read_xlsx("Data/economatica_nyse_diario.xlsx", skip = 3, col_types = c("date", rep("numeric", 1420)), na = c("-", " ", "NA"))|> 
      filter(Data > "2000-01-01" & Data < "2025-01-01") |> 
      filter(!if_all(where(is.numeric), is.na)) |> 
      select(where(~ mean(!is.na(.x)) >= 0.75)) |> 
      rename_with(~ str_remove(.x, "Fechamento\najust p/ prov\nEm moeda orig\n"), -Data) |> 
      select(!contains("old"))
    ii <- sample(3:ncol(aux), 1)
    data <- stochvol::logret(na.omit(aux[, ii][[1]])) * 100
  }
  if (abs(acf(data, plot = FALSE)$acf[2])> 2/sqrt(length(data))) {
    data <- ar.yw(data, 4, aic = TRUE, se.fit = FALSE)$resid
  }
  
  garch_params <- as.numeric(garch(r_garch_sim_n, cond_dist = "std", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)))[1:3]
  figarch_params <- as.numeric(figarch(data, cond_dist = "std", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE))@pars[1:4])
  dcs_params <- as.numeric(coef(tegarch(data, asym = FALSE, skew = FALSE)))[1:3]
  sv_params <- as.numeric(svtsample(data)$summary$para[c(1, 2, 3), 1])
  aux_ms <- as.numeric(msgarchfit(ms_spec_t, data)$par)
  ms_params <- aux_ms[c(1, 2, 3, 5, 6, 7)]
  P <- matrix(c(aux_ms[9], aux_ms[10], 1 - aux_ms[9], 1- aux_ms[10]), 2, 2, byrow = TRUE) 
}

true_ret_n <- matrix(NA, ncol = n_dgps, nrow = mc)
true_vols_n <- matrix(NA, ncol = n_dgps, nrow = mc)
fore_vols_n <- matrix(NA, ncol = (2 * n_dgps) * n_dgps, nrow = mc)
true_vols_t <- matrix(NA, ncol = n_dgps, nrow = mc)
true_ret_t <- matrix(NA, ncol = n_dgps, nrow = mc)
fore_vols_t <- matrix(NA, ncol = (2 * n_dgps) * n_dgps, nrow = mc)

for (i in 1:mc) {
  set.seed(i + 123)
  print(i)
  
  # Simulate DGPs
  garch_sim_n <- garch_sim(n_max + 1, garch_params, "norm")
  figarch_sim_n <- figarch_sim(pars = list(omega = figarch_params[1], phi = figarch_params[2], beta = figarch_params[3], d = figarch_params[4]), cond_dist = "norm", n = n_max + 1)
  dcs_sim_n <- dcs_sim(n_max + 1, dcs_params, "norm")
  sv_sim_n <- sv_sim(n_max + 1, sv_params, "norm")
  ms_sim_n <- msgarch_sim(n_max + 1, ms_params, "norm", P)
  garch_sim_t <- garch_sim(n_max + 1, c(garch_params, 7), "std")
  figarch_sim_t <- figarch_sim(pars = list(omega = figarch_params[1], phi = figarch_params[2], beta = figarch_params[3], d = figarch_params[4], df = 7), cond_dist = "std", n = n_max + 1)
  dcs_sim_t <- dcs_sim(n_max + 1, c(dcs_params, 7), "std")
  sv_sim_t <- sv_sim(n_max + 1, c(sv_params, 7), "std")
  ms_sim_t <- msgarch_sim(n_max + 1, c(ms_params, 7), "std", P)

  r_garch_sim_n <- tail(garch_sim_n$returns, n + 1)[1:n]
  r_figarch_sim_n <- tail(figarch_sim_n$rt, n + 1)[1:n]
  r_dcs_sim_n <- tail(dcs_sim_n$returns, n + 1)[1:n]
  r_sv_sim_n <- tail(sv_sim_n$returns, n + 1)[1:n]
  r_ms_sim_n <- tail(ms_sim_n$returns, n + 1)[1:n]

  r_garch_sim_t <- tail(garch_sim_t$returns, n + 1)[1:n]
  r_figarch_sim_t <- tail(figarch_sim_t$rt, n + 1)[1:n]
  r_dcs_sim_t <- tail(dcs_sim_t$returns, n + 1)[1:n]
  r_sv_sim_t <- tail(sv_sim_t$returns, n + 1)[1:n]
  r_ms_sim_t <- tail(ms_sim_t$returns, n + 1)[1:n]

  # Contaminate the series with AO
  if (outliers == "TRUE") {
    outlier_position <- floor(runif(1, n - 22, n)) + 1

    r_garch_sim_n[outlier_position] <- r_garch_sim_n[outlier_position] + sign(r_garch_sim_n[outlier_position]) * 5 * sd(r_garch_sim_n)
    r_figarch_sim_n[outlier_position] <- r_figarch_sim_n[outlier_position] + sign(r_figarch_sim_n[outlier_position]) * 5 * sd(r_figarch_sim_n)
    r_dcs_sim_n[outlier_position] <- r_dcs_sim_n[outlier_position] + sign(r_dcs_sim_n[outlier_position]) * 5 * sd(r_dcs_sim_n)
    r_sv_sim_n[outlier_position] <- r_sv_sim_n[outlier_position] + sign(r_sv_sim_n[outlier_position]) * 5 * sd(r_sv_sim_n)
    r_ms_sim_n[outlier_position] <- r_ms_sim_n[outlier_position] + sign(r_ms_sim_n[outlier_position]) * 5 * sd(r_ms_sim_n)

    r_garch_sim_t[outlier_position] <- r_garch_sim_t[outlier_position] + sign(r_garch_sim_t[outlier_position]) * 5 * sd(r_garch_sim_t)
    r_figarch_sim_t[outlier_position] <- r_figarch_sim_t[outlier_position] + sign(r_figarch_sim_t[outlier_position]) * 5 * sd(r_figarch_sim_t)
    r_dcs_sim_t[outlier_position] <- r_dcs_sim_t[outlier_position] + sign(r_dcs_sim_t[outlier_position]) * 5 * sd(r_dcs_sim_t)
    r_sv_sim_t[outlier_position] <- r_sv_sim_t[outlier_position] + sign(r_sv_sim_t[outlier_position]) * 5 * sd(r_sv_sim_t)
    r_ms_sim_t[outlier_position] <- r_ms_sim_t[outlier_position] + sign(r_ms_sim_t[outlier_position]) * 5 * sd(r_ms_sim_t)
  }

  # Fit the models
  garch_n_garch_n <- predict(garch(r_garch_sim_n, cond_dist = "norm", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  garch_n_figarch_n <-  predict(figarch(r_garch_sim_n, cond_dist = "norm", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  garch_n_dcs_n <- tail(vol_dcsm(r_garch_sim_n, dcsn_fit(r_garch_sim_n)), 1)
  garch_n_sv_n <- median(predict(estimate_parameters_n(r_garch_sim_n), steps = 1)$h_exp)
  garch_n_ms_n <- as.numeric(predict(msgarchfit(ms_spec_n, r_garch_sim_n), nahead = 1)$vol)
  garch_n_garch_t <- predict(garch(r_garch_sim_n, cond_dist = "std", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  garch_n_figarch_t <- predict(figarch(r_garch_sim_n, cond_dist = "std", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  garch_n_dcs_t <- as.numeric(predict(tegarch(r_garch_sim_n, asym = FALSE, skew = FALSE, components = 1, hessian = FALSE), n.ahead = 1, verbose=TRUE)$stdev)
  garch_n_sv_t <- median(predict(estimate_parameters_t(r_garch_sim_n), steps = 1)$h_exp)
  garch_n_ms_t <- as.numeric(predict(msgarchfit(ms_spec_t, r_garch_sim_n), nahead = 1)$vol)
  
  garch_t_garch_n <-   predict(garch(r_garch_sim_t, cond_dist = "norm", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  garch_t_figarch_n <- predict(figarch(r_garch_sim_t, cond_dist = "norm", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  garch_t_dcs_n <- tail(vol_dcsm(r_garch_sim_t, dcsn_fit(r_garch_sim_t)), 1)
  garch_t_sv_n <- median(predict(estimate_parameters_n(r_garch_sim_t), steps = 1)$h_exp)
  garch_t_ms_n <- as.numeric(predict(msgarchfit(ms_spec_n, r_garch_sim_t), nahead = 1)$vol)
  garch_t_garch_t <-   predict(garch(r_garch_sim_t, cond_dist = "std", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  garch_t_figarch_t <- predict(figarch(r_garch_sim_t, cond_dist = "std", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  garch_t_dcs_t <- as.numeric(predict(tegarch(r_garch_sim_t, asym = FALSE, skew = FALSE, components = 1, hessian = FALSE), n.ahead = 1, verbose=TRUE)$stdev)
  garch_t_sv_t <- median(predict(estimate_parameters_t(r_garch_sim_t), steps = 1)$h_exp)
  garch_t_ms_t <- as.numeric(predict(msgarchfit(ms_spec_t, r_garch_sim_t), nahead = 1)$vol)
  
  
  figarch_n_garch_n <-   predict(garch(r_figarch_sim_n, cond_dist = "norm", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  figarch_n_figarch_n <- predict(figarch(r_figarch_sim_n, cond_dist = "norm", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  figarch_n_dcs_n <-  tail(vol_dcsm(r_figarch_sim_n, dcsn_fit(r_figarch_sim_n)), 1)
  figarch_n_sv_n <- median(predict(estimate_parameters_n(r_figarch_sim_n), steps = 1)$h_exp)
  figarch_n_ms_n <- as.numeric(predict(msgarchfit(ms_spec_n, r_figarch_sim_n), nahead = 1)$vol)
  figarch_n_garch_t <-  predict(garch(r_figarch_sim_n, cond_dist = "std", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  figarch_n_figarch_t <- predict(figarch(r_figarch_sim_n, cond_dist = "std", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  figarch_n_dcs_t <- as.numeric(predict(tegarch(r_figarch_sim_n, asym = FALSE, skew = FALSE, components = 1, hessian = FALSE), n.ahead = 1, verbose=TRUE)$stdev)
  figarch_n_sv_t <- median(predict(estimate_parameters_t(r_figarch_sim_n), steps = 1)$h_exp)
  figarch_n_ms_t <- as.numeric(predict(msgarchfit(ms_spec_t, r_figarch_sim_n), nahead = 1)$vol)
  
  figarch_t_garch_n <-   predict(garch(r_figarch_sim_t, cond_dist = "norm", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  figarch_t_figarch_n <- predict(figarch(r_figarch_sim_t, cond_dist = "norm", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  figarch_t_dcs_n <- tail(vol_dcsm(r_figarch_sim_t, dcsn_fit(r_figarch_sim_t)), 1)
  figarch_t_sv_n <- median(predict(estimate_parameters_n(r_figarch_sim_t), steps = 1)$h_exp)
  figarch_t_ms_n <- as.numeric(predict(msgarchfit(ms_spec_n, r_figarch_sim_t), nahead = 1)$vol)
  figarch_t_garch_t <- predict(garch(r_figarch_sim_t, cond_dist = "norm", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  figarch_t_figarch_t <- predict(figarch(r_figarch_sim_t, cond_dist = "std", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  figarch_t_dcs_t <- as.numeric(predict(tegarch(r_figarch_sim_t, asym = FALSE, skew = FALSE, components = 1, hessian = FALSE), n.ahead = 1, verbose=TRUE)$stdev)
  figarch_t_sv_t <- median(predict(estimate_parameters_t(r_figarch_sim_t), steps = 1)$h_exp)
  figarch_t_ms_t <- as.numeric(predict(msgarchfit(ms_spec_t, r_figarch_sim_t), nahead = 1)$vol)

  dcs_n_garch_n <- predict(garch(r_dcs_sim_n, cond_dist = "norm", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  dcs_n_figarch_n <- predict(figarch(r_dcs_sim_n, cond_dist = "norm", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  dcs_n_dcs_n <- tail(vol_dcsm(r_dcs_sim_n, dcsn_fit(r_dcs_sim_n)), 1)
  dcs_n_sv_n <- median(predict(estimate_parameters_n(r_dcs_sim_n), steps = 1)$h_exp)
  dcs_n_ms_n <-as.numeric(predict(msgarchfit(ms_spec_n, r_dcs_sim_n), nahead = 1)$vol)
  dcs_n_garch_t <-  predict(garch(r_dcs_sim_n, cond_dist = "std", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  dcs_n_figarch_t <- predict(figarch(r_dcs_sim_n, cond_dist = "std", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  dcs_n_dcs_t <-  as.numeric(predict(tegarch(r_dcs_sim_n, asym = FALSE, skew = FALSE, components = 1, hessian = FALSE), n.ahead = 1, verbose=TRUE)$stdev)
  dcs_n_sv_t <- median(predict(estimate_parameters_t(r_dcs_sim_n), steps = 1)$h_exp)
  dcs_n_ms_t <- as.numeric(predict(msgarchfit(ms_spec_t, r_dcs_sim_n), nahead = 1)$vol)

  dcs_t_garch_n <- predict(garch(r_dcs_sim_t, cond_dist = "norm", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  dcs_t_figarch_n <- predict(figarch(r_dcs_sim_t, cond_dist = "norm", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  dcs_t_dcs_n <- tail(vol_dcsm(r_dcs_sim_t, dcsn_fit(r_dcs_sim_t)), 1)
  dcs_t_sv_n <-  median(predict(estimate_parameters_n(r_dcs_sim_t), steps = 1)$h_exp)
  dcs_t_ms_n <- as.numeric(predict(msgarchfit(ms_spec_n, r_dcs_sim_t), nahead = 1)$vol)
  dcs_t_garch_t <-  predict(garch(r_dcs_sim_t, cond_dist = "std", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  dcs_t_figarch_t <- predict(figarch(r_dcs_sim_t, cond_dist = "std", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  dcs_t_dcs_t <- as.numeric(predict(tegarch(r_dcs_sim_t, asym = FALSE, skew = FALSE, components = 1, hessian = FALSE), n.ahead = 1, verbose=TRUE)$stdev)
  dcs_t_sv_t <- median(predict(estimate_parameters_t(r_dcs_sim_t), steps = 1)$h_exp)
  dcs_t_ms_t <- as.numeric(predict(msgarchfit(ms_spec_t, r_dcs_sim_t), nahead = 1)$vol)

  
  sv_n_garch_n <- predict(garch(r_sv_sim_n, cond_dist = "norm", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  sv_n_figarch_n <- predict(figarch(r_sv_sim_n, cond_dist = "norm", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  sv_n_dcs_n <- tail(vol_dcsm(r_sv_sim_n, dcsn_fit(r_sv_sim_n)), 1)
  sv_n_sv_n <- median(predict(estimate_parameters_n(r_sv_sim_n), steps = 1)$h_exp)
  sv_n_ms_n <- as.numeric(predict(msgarchfit(ms_spec_n, r_sv_sim_n), nahead = 1)$vol)
  sv_n_garch_t <- predict(garch(r_sv_sim_n, cond_dist = "std", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  sv_n_figarch_t <- predict(figarch(r_sv_sim_n, cond_dist = "std", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  sv_n_dcs_t <- as.numeric(predict(tegarch(r_sv_sim_n, asym = FALSE, skew = FALSE, components = 1, hessian = FALSE), n.ahead = 1, verbose=TRUE)$stdev)
  sv_n_sv_t <- median(predict(estimate_parameters_t(r_sv_sim_n), steps = 1)$h_exp)
  sv_n_ms_t <- as.numeric(predict(msgarchfit(ms_spec_t, r_sv_sim_n), nahead = 1)$vol)

  sv_t_garch_n <- predict(garch(r_sv_sim_t, cond_dist = "norm", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  sv_t_figarch_n <- predict(figarch(r_sv_sim_t, cond_dist = "norm", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  sv_t_dcs_n <- tail(vol_dcsm(r_sv_sim_t, dcsn_fit(r_sv_sim_t)), 1)
  sv_t_sv_n <- median(predict(estimate_parameters_n(r_sv_sim_t), steps = 1)$h_exp)
  sv_t_ms_n <- as.numeric(predict(msgarchfit(ms_spec_n, r_sv_sim_t), nahead = 1)$vol)
  sv_t_garch_t <- predict(garch(r_sv_sim_t, cond_dist = "std", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  sv_t_figarch_t <- predict(figarch(r_sv_sim_t, cond_dist = "std", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  sv_t_dcs_t <- as.numeric(predict(tegarch(r_sv_sim_t, asym = FALSE, skew = FALSE, components = 1, hessian = FALSE), n.ahead = 1, verbose=TRUE)$stdev)
  sv_t_sv_t <- median(predict(estimate_parameters_t(r_sv_sim_t), steps = 1)$h_exp)
  sv_t_ms_t <- as.numeric(predict(msgarchfit(ms_spec_t, r_sv_sim_t), nahead = 1)$vol)

  
  ms_n_garch_n <- predict(garch(r_ms_sim_n, cond_dist = "norm", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  ms_n_figarch_n <- predict(figarch(r_ms_sim_n, cond_dist = "norm", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  ms_n_dcs_n <- tail(vol_dcsm(r_ms_sim_n, dcsn_fit(r_ms_sim_n)), 1)
  ms_n_sv_n <- median(predict(estimate_parameters_n(r_ms_sim_n), steps = 1)$h_exp)
  ms_n_ms_n <- as.numeric(predict(msgarchfit(ms_spec_n, r_ms_sim_n), nahead = 1)$vol)
  ms_n_garch_t <- predict(garch(r_ms_sim_n, cond_dist = "std", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  ms_n_figarch_t <- predict(figarch(r_ms_sim_n, cond_dist = "std", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  ms_n_dcs_t <-  as.numeric(predict(tegarch(r_ms_sim_n, asym = FALSE, skew = FALSE, components = 1, hessian = FALSE), n.ahead = 1, verbose=TRUE)$stdev)
  ms_n_sv_t <- median(predict(estimate_parameters_t(r_ms_sim_n), steps = 1)$h_exp)
  ms_n_ms_t <- as.numeric(predict(msgarchfit(ms_spec_t, r_ms_sim_n), nahead = 1)$vol)

  ms_t_garch_n <- predict(garch(r_ms_sim_t, cond_dist = "norm", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  ms_t_figarch_n <-predict(figarch(r_ms_sim_t, cond_dist = "norm", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  ms_t_dcs_n <- tail(vol_dcsm(r_ms_sim_t, dcsn_fit(r_ms_sim_t)), 1)
  ms_t_sv_n <- median(predict(estimate_parameters_n(r_ms_sim_t), steps = 1)$h_exp)
  ms_t_ms_n <- as.numeric(predict(msgarchfit(ms_spec_n, r_ms_sim_t), nahead = 1)$vol)
  ms_t_garch_t <- predict(garch(r_ms_sim_t, cond_dist = "std", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  ms_t_figarch_t <- predict(figarch(r_ms_sim_t, cond_dist = "std", meanspec = mean_spec(orders = c(0, 0), long_memo = FALSE, include_mean = FALSE)), n.ahead = 1)@sigt
  ms_t_dcs_t <-  as.numeric(predict(tegarch(r_ms_sim_t, asym = FALSE, skew = FALSE, components = 1, hessian = FALSE), n.ahead = 1, verbose=TRUE)$stdev)
  ms_t_sv_t <- median(predict(estimate_parameters_t(r_ms_sim_t), steps = 1)$h_exp)
  ms_t_ms_t <- as.numeric(predict(msgarchfit(ms_spec_t, r_ms_sim_t), nahead = 1)$vol)
  

  true_vols_n[i, ] <- c(tail(garch_sim_n$volatility, 1), tail(dcs_sim_n$volatility, 1), tail(sv_sim_n$volatility, 1), tail(ms_sim_n$volatility[, 3], 1), tail(figarch_sim_n$sigt, 1))
  true_vols_t[i, ] <- c(tail(garch_sim_t$volatility, 1), tail(dcs_sim_t$volatility, 1), tail(sv_sim_t$volatility, 1), tail(ms_sim_t$volatility[, 3], 1), tail(figarch_sim_t$sigt, 1))
  fore_vols_n[i, ] <- c(garch_n_garch_n, garch_n_garch_t, garch_n_dcs_n, garch_n_dcs_t, garch_n_sv_n, garch_n_sv_t,  garch_n_ms_n, garch_n_ms_t, garch_n_figarch_n, garch_n_figarch_t,
                        dcs_n_garch_n, dcs_n_garch_t, dcs_n_dcs_n, dcs_n_dcs_t, dcs_n_sv_n, dcs_n_sv_t,  dcs_n_ms_n, dcs_n_ms_t, dcs_n_figarch_n, dcs_n_figarch_t, 
                        sv_n_garch_n, sv_n_garch_t, sv_n_dcs_n, sv_n_dcs_t, sv_n_sv_n, sv_n_sv_t, sv_n_ms_n, sv_n_ms_t, sv_n_figarch_n, sv_n_figarch_t,
                        ms_n_garch_n, ms_n_garch_t, ms_n_dcs_n, ms_n_dcs_t, ms_n_sv_n, ms_n_sv_t, ms_n_ms_n, ms_n_ms_t, ms_n_figarch_n, ms_n_figarch_t,
                        figarch_n_garch_n, figarch_n_garch_t, figarch_n_dcs_n, figarch_n_dcs_t, figarch_n_sv_n, figarch_n_sv_t,  figarch_n_ms_n, figarch_n_ms_t, figarch_n_figarch_n, figarch_n_figarch_t)
  fore_vols_t[i, ] <- c(garch_t_garch_n, garch_t_garch_t, garch_t_dcs_n, garch_t_dcs_t, garch_t_sv_n, garch_t_sv_t, garch_t_ms_n, garch_t_ms_t, garch_t_figarch_n, garch_t_figarch_t,
                        dcs_t_garch_n, dcs_t_garch_t, dcs_t_dcs_n, dcs_t_dcs_t, dcs_t_sv_n, dcs_t_sv_t, dcs_t_ms_n, dcs_t_ms_t, dcs_t_figarch_n, dcs_t_figarch_t,
                        sv_t_garch_n, sv_t_garch_t, sv_t_dcs_n, sv_t_dcs_t, sv_t_sv_n, sv_t_sv_t, sv_t_ms_n, sv_t_ms_t, sv_t_figarch_n, sv_t_figarch_t,
                        ms_t_garch_n, ms_t_garch_t, ms_t_dcs_n, ms_t_dcs_t, ms_t_sv_n, ms_t_sv_t, ms_t_ms_n, ms_t_ms_t, ms_t_figarch_n, ms_t_figarch_t,
                        figarch_t_garch_n, figarch_t_garch_t, figarch_t_dcs_n, figarch_t_dcs_t, figarch_t_sv_n, figarch_t_sv_t, figarch_t_ms_n, figarch_t_ms_t, figarch_t_figarch_n, figarch_t_figarch_t)
  

  true_ret_n[i, ] <- c(tail(garch_sim_n$returns, 1), tail(dcs_sim_n$returns, 1), tail(sv_sim_n$returns, 1), tail(ms_sim_n$returns, 1), tail(figarch_sim_n$rt, 1))
  true_ret_t[i, ] <- c(tail(garch_sim_t$returns, 1), tail(dcs_sim_t$returns, 1), tail(sv_sim_t$returns, 1), tail(ms_sim_t$returns, 1), tail(figarch_sim_t$rt, 1))
  
  
  }

volatilities_n <- cbind(true_vols_n, fore_vols_n)
volatilities_t <- cbind(true_vols_t, fore_vols_t)
colnames(volatilities_n) <- c("garch", "dcs", "sv", "ms", "figarch",
  "garch_n_garch_n", "garch_n_garch_t", "garch_n_dcs_n", "garch_n_dcs_t", "garch_n_sv_n", "garch_n_sv_t", "garch_n_ms_n", "garch_n_ms_t", "garch_n_figarch_n", "garch_n_figarch_t",
  "dcs_n_garch_n", "dcs_n_garch_t", "dcs_n_dcs_n", "dcs_n_dcs_t", "dcs_n_sv_n", "dcs_n_sv_t", "dcs_n_ms_n", "dcs_n_ms_t", "dcs_n_figarch_n", "dcs_n_figarch_t",
  "sv_n_garch_n", "sv_n_garch_t", "sv_n_dcs_n", "sv_n_dcs_t", "sv_n_sv_n", "sv_n_sv_t", "sv_n_ms_n", "sv_n_ms_t", "sv_n_figarch_n", "sv_n_figarch_t",
  "ms_n_garch_n", "ms_n_garch_t", "ms_n_dcs_n", "ms_n_dcs_t", "ms_n_sv_n", "ms_n_sv_t", "ms_n_ms_n", "ms_n_ms_t", "ms_n_figarch_n", "ms_n_figarch_t",
  "figarch_n_garch_n", "figarch_n_garch_t", "figarch_n_dcs_n", "figarch_n_dcs_t", "figarch_n_sv_n", "figarch_n_sv_t", "figarch_n_ms_n", "figarch_n_ms_t", "figarch_n_figarch_n", "figarch_n_figarch_t")
colnames(volatilities_t) <- c("garch", "dcs", "sv", "ms", "figarch",
  "garch_t_garch_n", "garch_t_garch_t", "garch_t_dcs_n", "garch_t_dcs_t", "garch_t_sv_n", "garch_t_sv_t", "garch_t_ms_n", "garch_t_ms_t", "garch_t_figarch_n", "garch_t_figarch_t",
  "dcs_t_garch_n", "dcs_t_garch_t", "dcs_t_dcs_n", "dcs_t_dcs_t", "dcs_t_sv_n", "dcs_t_sv_t", "dcs_t_ms_n", "dcs_t_ms_t", "dcs_t_figarch_n", "dcs_t_figarch_t",
  "sv_t_garch_n", "sv_t_garch_t", "sv_t_dcs_n", "sv_t_dcs_t", "sv_t_sv_n", "sv_t_sv_t", "sv_t_ms_n", "sv_t_ms_t", "sv_t_figarch_n", "sv_t_figarch_t",
  "ms_t_garch_n", "ms_t_garch_t", "ms_t_dcs_n", "ms_t_dcs_t", "ms_t_sv_n", "ms_t_sv_t", "ms_t_ms_n", "ms_t_ms_t", "ms_t_figarch_n", "ms_t_figarch_t", 
  "figarch_t_garch_n", "figarch_t_garch_t", "figarch_t_dcs_n", "figarch_t_dcs_t", "figarch_t_sv_n", "figarch_t_sv_t", "figarch_t_ms_n", "figarch_t_ms_t", "figarch_t_figarch_n", "figarch_t_figarch_t")
colnames(true_ret_n) <- c("garch", "dcs", "sv", "ms", "figarch")
colnames(true_ret_t) <- c("garch", "dcs", "sv", "ms", "figarch")

if (type == "RND") {
  write.csv(volatilities_n, paste0("volatilities_", n, "_norm_", outliers, "_", type, colnames(aux)[ii], ".csv"))
  write.csv(volatilities_t, paste0("volatilities_", n, "_std_", outliers, "_", type, colnames(aux)[ii], ".csv"))
  write.csv(true_ret_n, paste0("true_ret_n_", n, "_norm_", outliers, "_", type, colnames(aux)[ii],  ".csv"))
  write.csv(true_ret_t, paste0("true_ret_t_", n, "_std_", outliers, "_", type, colnames(aux)[ii], ".csv"))
} else {
  write.csv(volatilities_n, paste0("volatilities_", n, "_norm_", outliers, "_", type, ".csv"))
  write.csv(volatilities_t, paste0("volatilities_", n, "_std_", outliers, "_", type, ".csv"))
  write.csv(true_ret_n, paste0("true_ret_n_", n, "_norm_", outliers, "_", type, ".csv"))
  write.csv(true_ret_t, paste0("true_ret_t_", n, "_std_", outliers, "_", type, ".csv"))
}


