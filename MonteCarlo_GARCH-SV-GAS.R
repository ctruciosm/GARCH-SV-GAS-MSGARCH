################################################################################
#####                   Monte Carlo Simulation                             #####
################################################################################
args <- commandArgs(TRUE)
library(rugarch)
library(GAS)
library(stochvol)
library(MSGARCH)
source("DGPs.R")

## Setting values
mc <- 1000


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
        

# Parse arguments if they exist
if(length(args) > 0){
  for(i in 1:length(args)){
    arg_parts <- strsplit(args[[i]], "=")[[1]]
    if(length(arg_parts) == 2){
      param_name <- arg_parts[1]
      param_value <- arg_parts[2]
      if(param_name == "n") n <- as.integer(param_value)
      if(param_name == "type") type <- param_value
      if(param_name == "outliers") outliers <- param_value
    }
  }
} else {
  n <- 1000
  type <- "BR"
  outliers <-  "TRUE"
}


# Estimated a BR time series to obtain parameters 
garch_spec_t <- ugarchspec(variance.model = list(model = 'sGARCH', garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),  distribution.model = "std")
gas_spec_t <- UniGASSpec(Dist = "std", ScalingType = "Identity", GASPar = list(locate = FALSE, scale = TRUE, shape = FALSE))
garch_spec_n <- ugarchspec(variance.model = list(model = 'sGARCH', garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),  distribution.model = "norm")
gas_spec_n <- UniGASSpec(Dist = "norm", ScalingType = "Identity", GASPar = list(locate = FALSE, scale = TRUE, shape = FALSE))
ms_spec_t <- CreateSpec(variance.spec = list(model = c("sGARCH","sGARCH")), switch.spec = list(do.mix = FALSE), distribution.spec = list(distribution = c("std", "std")), constraint.spec = list(regime.const = c("nu")))
ms_spec_n <- CreateSpec(variance.spec = list(model = c("sGARCH","sGARCH")), switch.spec = list(do.mix = FALSE), distribution.spec = list(distribution = c("norm", "norm")))

if (type == "BR") {
  # data <- logret(read.csv("./Data/precos_diarios_ibrx.csv")[, "PETR4"]) * 100
  garch_params <- c(0.18, 0.09, 0.89)
  gas_params <- c(0.03, 0.22, 0.98)
  sv_params <- c(1.74, 0.97, 0.18)
  ms_params <- c(0.01, 0.04, 0.92, 0.57, 0.03, 0.95, 0.73, 0.44)
} 
if (type == "US") {
  # data <- logret(read.csv("./Data/precos_diarios_nyse.csv")[, "Close"]) * 100
  
}


true_vols_n <- matrix(NA, ncol = 4, nrow = mc)
fore_vols_n <- matrix(NA, ncol = 16, nrow = mc)
true_vols_t <- matrix(NA, ncol = 4, nrow = mc)
fore_vols_t <- matrix(NA, ncol = 16, nrow = mc)
for (i in 1:mc) {
  set.seed(i + 123)
  print(i)
  
  garch_sim_n <- garch_sim(2500 + 1, garch_params, "norm")
  gas_sim_n <- gas_sim(2500 + 1, gas_params, "norm")
  sv_sim_n <- sv_sim(2500 + 1, sv_params, "norm")
  ms_sim_n <- simulate(object = ms_spec_n, nsim = 1L, nahead = 2500 + 1, nburn = 500L, par = ms_params)
    
  garch_sim_t <- garch_sim(2500 + 1, c(garch_params, 7), "std")
  gas_sim_t <- gas_sim(2500 + 1, c(gas_params, -2.6625878), "std")
  sv_sim_t <- sv_sim(2500 + 1, c(sv_params, 7), "std")  
  ms_sim_t <- simulate(object = ms_spec_t, nsim = 1L, nahead = 2500 + 1, nburn = 500L, par = c(ms_params[1:3], 7, ms_params[4:6], 7, ms_params[7:8]))
  

  r_garch_sim_n <- tail(garch_sim_n$returns, n + 1)[1:n]
  r_gas_sim_n <- tail(gas_sim_n$returns, n + 1)[1:n]
  r_sv_sim_n <- tail(sv_sim_n$returns, n + 1)[1:n]
  r_ms_sim_n <- tail(ms_sim_n$draw, n + 1)[1:n]
  
  r_garch_sim_t <- tail(garch_sim_t$returns, n + 1)[1:n]
  r_gas_sim_t <- tail(gas_sim_t$returns, n + 1)[1:n]
  r_sv_sim_t <- tail(sv_sim_t$returns, n + 1)[1:n]
  r_ms_sim_t <- tail(ms_sim_t$draw, n + 1)[1:n]
  
  
  if (outliers == "TRUE") {
    outlier_position <- floor(runif(1, n - 22, n)) + 1
    
    r_garch_sim_n[outlier_position] <- r_garch_sim_n[outlier_position] + sign(r_garch_sim_n[outlier_position]) * 5 * sd(r_garch_sim_n)
    r_gas_sim_n[outlier_position] <- r_gas_sim_n[outlier_position] + sign(r_gas_sim_n[outlier_position]) * 5 * sd(r_gas_sim_n)
    r_sv_sim_n[outlier_position] <- r_sv_sim_n[outlier_position] + sign(r_sv_sim_n[outlier_position]) * 5 * sd(r_sv_sim_n)
    r_ms_sim_n[outlier_position] <- r_ms_sim_n[outlier_position] + sign(r_ms_sim_n[outlier_position]) * 5 * sd(r_ms_sim_n)
    
    r_garch_sim_t[outlier_position] <- r_garch_sim_t[outlier_position] + sign(r_garch_sim_t[outlier_position]) * 5 * sd(r_garch_sim_t)
    r_gas_sim_t[outlier_position] <- r_gas_sim_t[outlier_position] + sign(r_gas_sim_t[outlier_position]) * 5 * sd(r_gas_sim_t)
    r_sv_sim_t[outlier_position] <- r_sv_sim_t[outlier_position] + sign(r_sv_sim_t[outlier_position]) * 5 * sd(r_sv_sim_t)
    r_ms_sim_t[outlier_position] <- r_ms_sim_t[outlier_position] + sign(r_ms_sim_t[outlier_position]) * 5 * sd(r_ms_sim_t)
  }
  
  garch_garch_n <- as.numeric(sigma(ugarchforecast(ugarchfit(garch_spec_n, r_garch_sim_n), n.ahead = 1)))
  garch_gas_n <- sqrt(UniGASFor(UniGASFit(gas_spec_n, r_garch_sim_n, Compute.SE = FALSE), H = 1)@Forecast$PointForecast[, 2])
  garch_sv_n <- median(predvola(predict(svsample(r_garch_sim_n), 1)))
  garch_ms_n <- predict(msgarchfit(ms_spec_n, r_garch_sim_n), nahead = 1)$vol

  garch_garch_t <- as.numeric(sigma(ugarchforecast(ugarchfit(garch_spec_t, r_garch_sim_t), n.ahead = 1)))
  garch_gas_t <- sqrt(UniGASFor(UniGASFit(gas_spec_t, r_garch_sim_t, Compute.SE = FALSE), H = 1)@Forecast$PointForecast[, 2])
  garch_sv_t <- median(predvola(predict(svtsample(r_garch_sim_t), 1)))
  garch_ms_t <- predict(msgarchfit(ms_spec_t, r_garch_sim_t), nahead = 1)$vol
  
  gas_garch_n <- as.numeric(sigma(ugarchforecast(ugarchfit(garch_spec_n, r_gas_sim_n), n.ahead = 1)))
  gas_gas_n <- sqrt(UniGASFor(UniGASFit(gas_spec_n, r_gas_sim_n, Compute.SE = FALSE), H = 1)@Forecast$PointForecast[, 2])
  gas_sv_n <- median(predvola(predict(svsample(r_gas_sim_n), 1)))
  gas_ms_n <- predict(msgarchfit(ms_spec_n, r_gas_sim_n), nahead = 1)$vol
  
  gas_garch_t <- as.numeric(sigma(ugarchforecast(ugarchfit(garch_spec_t, r_gas_sim_t), n.ahead = 1)))
  gas_gas_t <- sqrt(UniGASFor(UniGASFit(gas_spec_t, r_gas_sim_t, Compute.SE = FALSE), H = 1)@Forecast$PointForecast[, 2])
  gas_sv_t <- median(predvola(predict(svtsample(r_gas_sim_t), 1)))
  gas_ms_t <- predict(msgarchfit(ms_spec_t, r_gas_sim_t), nahead = 1)$vol
  
  sv_garch_n <- as.numeric(sigma(ugarchforecast(ugarchfit(garch_spec_n, r_sv_sim_n), n.ahead = 1)))
  sv_gas_n <- sqrt(UniGASFor(UniGASFit(gas_spec_n, r_sv_sim_n, Compute.SE = FALSE), H = 1)@Forecast$PointForecast[, 2])
  sv_sv_n <- median(predvola(predict(svsample(r_sv_sim_n), 1)))
  sv_ms_n <- predict(msgarchfit(ms_spec_n, r_sv_sim_n), nahead = 1)$vol
  
  sv_garch_t <- as.numeric(sigma(ugarchforecast(ugarchfit(garch_spec_t, r_sv_sim_t), n.ahead = 1)))
  sv_gas_t <- sqrt(UniGASFor(UniGASFit(gas_spec_t, r_sv_sim_t, Compute.SE = FALSE), H = 1)@Forecast$PointForecast[, 2])
  sv_sv_t <- median(predvola(predict(svtsample(r_sv_sim_t), 1)))
  sv_ms_t <- predict(msgarchfit(ms_spec_t, r_sv_sim_t), nahead = 1)$vol
  
  ms_garch_n <- as.numeric(sigma(ugarchforecast(ugarchfit(garch_spec_n, r_ms_sim_n), n.ahead = 1)))
  ms_gas_n <- sqrt(UniGASFor(UniGASFit(gas_spec_n, r_ms_sim_n, Compute.SE = FALSE), H = 1)@Forecast$PointForecast[, 2])
  ms_sv_n <- median(predvola(predict(svsample(r_ms_sim_n), 1)))
  ms_ms_n <- predict(msgarchfit(ms_spec_n, r_ms_sim_n), nahead = 1)$vol
  
  ms_garch_t <- as.numeric(sigma(ugarchforecast(ugarchfit(garch_spec_t, r_ms_sim_t), n.ahead = 1)))
  ms_gas_t <- sqrt(UniGASFor(UniGASFit(gas_spec_t, r_ms_sim_t, Compute.SE = FALSE), H = 1)@Forecast$PointForecast[, 2])
  ms_sv_t <- median(predvola(predict(svtsample(r_ms_sim_t), 1)))
  ms_ms_t <- predict(msgarchfit(ms_spec_t, r_ms_sim_t), nahead = 1)$vol
  
  true_vols_n[i, ] <- c(tail(garch_sim_n$volatility, 1), tail(gas_sim_n$volatility, 1), tail(sv_sim_n$volatility, 1), tail(ms_sim_n$CondVol[, 1, as.numeric(tail(ms_sim_n$state, 1))], 1))
  true_vols_t[i, ] <- c(tail(garch_sim_t$volatility, 1), tail(gas_sim_t$volatility, 1), tail(sv_sim_t$volatility, 1), tail(ms_sim_t$CondVol[, 1, as.numeric(tail(ms_sim_t$state, 1))], 1))
  fore_vols_n[i, ] <- c(garch_garch_n, garch_gas_n, garch_sv_n, garch_ms_n, gas_garch_n, gas_gas_n, gas_sv_n, gas_ms_n, sv_garch_n, sv_gas_n, sv_sv_n, sv_ms_n, ms_garch_n, ms_gas_n, ms_sv_n, ms_ms_n)
  fore_vols_t[i, ] <- c(garch_garch_t, garch_gas_t, garch_sv_t, garch_ms_t, gas_garch_t, gas_gas_t, gas_sv_t, gas_ms_t, sv_garch_t, sv_gas_t, sv_sv_t, sv_ms_t, ms_garch_t, ms_gas_t, ms_sv_t, ms_ms_t)
}

volatilities_n <- cbind(true_vols_n, fore_vols_n)
volatilities_t <- cbind(true_vols_t, fore_vols_t)
colnames(volatilities_n) <- c("garch", "gas", "sv", "ms", "garch_garch", "garch_gas", "garch_sv", "garch_ms", "gas_garch", "gas_gas", "gas_sv", "gas_ms", "sv_garch", "sv_gas", "sv_sv", "sv_ms", "ms_garch", "ms_gas", "ms_sv", "ms_ms")
colnames(volatilities_t) <- c("garch", "gas", "sv", "ms", "garch_garch", "garch_gas", "garch_sv", "garch_ms", "gas_garch", "gas_gas", "gas_sv", "gas_ms", "sv_garch", "sv_gas", "sv_sv", "sv_ms", "ms_garch", "ms_gas", "ms_sv", "ms_ms")

write.csv(volatilities_n, paste0("volatilities_", n, "_norm_", outliers, "_", type, ".csv"))
write.csv(volatilities_t, paste0("volatilities_", n, "_std_", outliers,"_", type, ".csv"))