################################################################################
#####                   Monte Carlo Simulation                             #####
################################################################################
args <- commandArgs(TRUE)
library(rugarch)
library(GAS)
library(stochvol)
source("DGPs.R")

## Setting values
mc <- 1000

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
  n <- 500
  type <- "BR"
  outliers <-  "FALSE"
}


# Estimated a BR time series to obtain parameters 
garch_spec_t <- ugarchspec(variance.model = list(model = 'sGARCH', garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),  distribution.model = "std")
gas_spec_t <- UniGASSpec(Dist = "std", ScalingType = "Identity", GASPar = list(locate = FALSE, scale = TRUE, shape = FALSE))
garch_spec_n <- ugarchspec(variance.model = list(model = 'sGARCH', garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),  distribution.model = "norm")
gas_spec_n <- UniGASSpec(Dist = "norm", ScalingType = "Identity", GASPar = list(locate = FALSE, scale = TRUE, shape = FALSE))

if (type == "BR") {
  # data <- logret(read.csv("./Data/precos_diarios_ibrx.csv")[, "PETR4"]) * 100
  garch_params <- c(0.18, 0.09, 0.89)
  gas_params <- c(0.03, 0.22, 0.98)
  sv_params <- c(1.74, 0.97, 0.18)
} 
if (type == "BTC") {
  # data <- logret(read.csv("./Data/BTCUSDT_1d.csv")[, "Close"]) * 100
  garch_params <- c(0.14, 0.07, 0.93)
  gas_params <- c(0.03, 0.22, 0.98)
  sv_params <- c(2.47, 0.98, 0.17)
}

if (type == "US") {
  # data <- logret(read.csv("./Data/precos_diarios_nyse.csv")[, "Close"]) * 100
  
}


true_vols_n <- matrix(NA, ncol = 3, nrow = mc)
fore_vols_n <- matrix(NA, ncol = 9, nrow = mc)
true_vols_t <- matrix(NA, ncol = 3, nrow = mc)
fore_vols_t <- matrix(NA, ncol = 9, nrow = mc)
for (i in 1:mc) {
  set.seed(i + 123)
  print(i)
  
  garch_sim_n <- garch_sim(2500 + 1, garch_params, "norm")
  gas_sim_n <- gas_sim(2500 + 1, gas_params, "norm")
  sv_sim_n <- sv_sim(2500 + 1, sv_params, "norm")
    
  garch_sim_t <- garch_sim(2500 + 1, c(garch_params, 7), "std")
  gas_sim_t <- gas_sim(2500 + 1, c(gas_params, -2.6625878), "std")
  sv_sim_t <- sv_sim(2500 + 1, c(sv_params, 7), "std")  

  r_garch_sim_n <- tail(garch_sim_n$returns, n + 1)[1:n]
  r_gas_sim_n <- tail(gas_sim_n$returns, n + 1)[1:n]
  r_sv_sim_n <- tail(sv_sim_n$returns, n + 1)[1:n]
  
  r_garch_sim_t <- tail(garch_sim_t$returns, n + 1)[1:n]
  r_gas_sim_t <- tail(gas_sim_t$returns, n + 1)[1:n]
  r_sv_sim_t <- tail(sv_sim_t$returns, n + 1)[1:n]
  
  
  if (outliers == "TRUE") {
    outlier_position <- floor(runif(1, n - 22, n)) + 1
    
    r_garch_sim_n[outlier_position] <- r_garch_sim_n[outlier_position] + sign(r_garch_sim_n[outlier_position]) * 5 * sd(r_garch_sim_n)
    r_gas_sim_n[outlier_position] <- r_gas_sim_n[outlier_position] + sign(r_gas_sim_n[outlier_position]) * 5 * sd(r_gas_sim_n)
    r_sv_sim_n[outlier_position] <- r_sv_sim_n[outlier_position] + sign(r_sv_sim_n[outlier_position]) * 5 * sd(r_sv_sim_n)
    
    r_garch_sim_t[outlier_position] <- r_garch_sim_t[outlier_position] + sign(r_garch_sim_t[outlier_position]) * 5 * sd(r_garch_sim_t)
    r_gas_sim_t[outlier_position] <- r_gas_sim_t[outlier_position] + sign(r_gas_sim_t[outlier_position]) * 5 * sd(r_gas_sim_t)
    r_sv_sim_t[outlier_position] <- r_sv_sim_t[outlier_position] + sign(r_sv_sim_t[outlier_position]) * 5 * sd(r_sv_sim_t)
  }
  
  garch_garch_n <- as.numeric(sigma(ugarchforecast(ugarchfit(garch_spec_n, r_garch_sim_n), n.ahead = 1)))
  garch_gas_n <- sqrt(UniGASFor(UniGASFit(gas_spec_n, r_garch_sim_n), H = 1)@Forecast$PointForecast[, 2])
  garch_sv_n <- median(predvola(predict(svsample(r_garch_sim_n), 1)))
  garch_garch_t <- as.numeric(sigma(ugarchforecast(ugarchfit(garch_spec_t, r_garch_sim_t), n.ahead = 1)))
  garch_gas_t <- sqrt(UniGASFor(UniGASFit(gas_spec_t, r_garch_sim_t), H = 1)@Forecast$PointForecast[, 2])
  garch_sv_t <- median(predvola(predict(svtsample(r_garch_sim_t), 1)))
  
  gas_garch_n <- as.numeric(sigma(ugarchforecast(ugarchfit(garch_spec_n, r_gas_sim_n), n.ahead = 1)))
  gas_gas_n <- sqrt(UniGASFor(UniGASFit(gas_spec_n, r_gas_sim_n), H = 1)@Forecast$PointForecast[, 2])
  gas_sv_n <- median(predvola(predict(svsample(r_gas_sim_n), 1)))
  gas_garch_t <- as.numeric(sigma(ugarchforecast(ugarchfit(garch_spec_t, r_gas_sim_t), n.ahead = 1)))
  gas_gas_t <- sqrt(UniGASFor(UniGASFit(gas_spec_t, r_gas_sim_t), H = 1)@Forecast$PointForecast[, 2])
  gas_sv_t <- median(predvola(predict(svtsample(r_gas_sim_t), 1)))
  
  sv_garch_n <- as.numeric(sigma(ugarchforecast(ugarchfit(garch_spec_n, r_sv_sim_n), n.ahead = 1)))
  sv_gas_n <- sqrt(UniGASFor(UniGASFit(gas_spec_n, r_sv_sim_n), H = 1)@Forecast$PointForecast[, 2])
  sv_sv_n <- median(predvola(predict(svsample(r_sv_sim_n), 1)))
  sv_garch_t <- as.numeric(sigma(ugarchforecast(ugarchfit(garch_spec_t, r_sv_sim_t), n.ahead = 1)))
  sv_gas_t <- sqrt(UniGASFor(UniGASFit(gas_spec_t, r_sv_sim_t), H = 1)@Forecast$PointForecast[, 2])
  sv_sv_t <- median(predvola(predict(svtsample(r_sv_sim_t), 1)))
  
  true_vols_n[i, ] <- c(tail(garch_sim_n$volatility, 1), tail(gas_sim_n$volatility, 1), tail(sv_sim_n$volatility, 1))
  true_vols_t[i, ] <- c(tail(garch_sim_t$volatility, 1), tail(gas_sim_t$volatility, 1), tail(sv_sim_t$volatility, 1))
  fore_vols_n[i, ] <- c(garch_garch_n, garch_gas_n, garch_sv_n, gas_garch_n, gas_gas_n, gas_sv_n, sv_garch_n, sv_gas_n, sv_sv_n)
  fore_vols_t[i, ] <- c(garch_garch_t, garch_gas_t, garch_sv_t, gas_garch_t, gas_gas_t, gas_sv_t, sv_garch_t, sv_gas_t, sv_sv_t)
}

volatilities_n <- cbind(true_vols_n, fore_vols_n)
volatilities_t <- cbind(true_vols_t, fore_vols_t)
colnames(volatilities_n) <- c("garch", "gas", "sv", "garch_garch", "garch_gas", "garch_sv", "gas_garch", "gas_gas", "gas_sv", "sv_garch", "sv_gas", "sv_sv")
colnames(volatilities_t) <- c("garch", "gas", "sv", "garch_garch", "garch_gas", "garch_sv", "gas_garch", "gas_gas", "gas_sv", "sv_garch", "sv_gas", "sv_sv")

write.csv(volatilities_n, paste0("volatilities_", n, "_norm_", outliers, "_", type, ".csv"))
write.csv(volatilities_t, paste0("volatilities_", n, "_std_", outliers,"_", type, ".csv"))