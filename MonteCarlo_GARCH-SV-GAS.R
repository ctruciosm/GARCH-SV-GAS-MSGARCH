################################################################################
#####                   Monte Carlo Simulation                             #####
################################################################################
args <- commandArgs(TRUE)
library(rugarch)
library(GAS)
library(stochvol)

## Setting values
mc <- 1000
n_burnin <- 1000



# Parse arguments if they exist
if(length(args) > 0){
  for(i in 1:length(args)){
    arg_parts <- strsplit(args[[i]], "=")[[1]]
    if(length(arg_parts) == 2){
      param_name <- arg_parts[1]
      param_value <- arg_parts[2]
      if(param_name == "n") n <- as.integer(param_value)
      if(param_name == "distri") distri <- param_value
    }
  }
} else {
  n <- 1000
  distri <- "std"
}
n
distri


if (distri == "norm") {
  nu <- Inf
  sv_sample <- svsample
} else {
  sv_sample <- svtsample
}


# Estimated a BR time series to obtain parameters 
data <- logret(read.csv("./Data/precos_diarios_ibrx.csv")[, "PETR4"]) * 100
garch_spec <- ugarchspec(variance.model = list(model = 'sGARCH', garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),  distribution.model = distri)
gas_spec <- UniGASSpec(Dist = distri, ScalingType = "Inv", GASPar = list(locate = FALSE, scale = TRUE, shape = FALSE))

garch_fit <- ugarchfit(garch_spec, data)
gas_fit <- UniGASFit(gas_spec, data)
if (distri == "norm") {
  nu <- Inf
  sv_sample <- svsample
  sv_fit <- sv_sample(data)
  sv_params <- sv_fit$summary$para[,1]
} else {
  sv_sample <- svtsample
  sv_fit <- sv_sample(data)
  sv_params <- sv_fit$summary$para[,1]
  nu <- sv_params[4]
}

true_vols <- matrix(NA, ncol = 3, nrow = mc)
fore_vols <- matrix(NA, ncol = 9, nrow = mc)

for (i in 1:mc) {
  set.seed(i + 123)
  print(i)
  garch_sim <- ugarchsim(garch_fit, n.sim = n + 1, n.start = n_burnin)
  gas_sim <- UniGASSim(gas_fit, T.sim = n)
  sv_sim <- svsim(n + 1, mu = sv_params[1], phi = sv_params[2], sigma = sv_params[3], nu)
  
  r_garch_sim <- garch_sim@simulation$seriesSim[1:n]
  r_gas_sim <- gas_sim@GASDyn$vY[1:n, 1]
  r_sv_sim <- sv_sim$y[1:n]
  
  garch_garch <- as.numeric(sigma(ugarchforecast(ugarchfit(garch_spec, r_garch_sim), n.ahead = 1)))
  garch_gas <- sqrt(UniGASFor(UniGASFit(gas_spec, r_garch_sim), H = 1)@Forecast$PointForecast[, 2])
  garch_sv <- median(predvola(predict(sv_sample(r_garch_sim), 1)))
  
  gas_garch <- as.numeric(sigma(ugarchforecast(ugarchfit(garch_spec, r_gas_sim), n.ahead = 1)))
  gas_gas <- sqrt(UniGASFor(UniGASFit(gas_spec, r_gas_sim), H = 1)@Forecast$PointForecast[, 2])
  gas_sv <- median(predvola(predict(sv_sample(r_gas_sim), 1)))
  
  sv_garch <- as.numeric(sigma(ugarchforecast(ugarchfit(garch_spec, r_sv_sim), n.ahead = 1)))
  sv_gas <- sqrt(UniGASFor(UniGASFit(gas_spec, r_sv_sim), H = 1)@Forecast$PointForecast[, 2])
  sv_sv <- median(predvola(predict(sv_sample(r_sv_sim), 1)))
  
  true_vols[i, ] <- c(tail(garch_sim@simulation$sigmaSim, 1), sqrt(tail(gas_sim@GASDyn$mTheta[2, ], 1)), tail(sv_sim$vol, 1))
  fore_vols[i, ] <- c(garch_garch, garch_gas, garch_sv, gas_garch, gas_gas, gas_sv, sv_garch, sv_gas, sv_sv)
}

volatilities <- cbind(true_vols, fore_vols)
colnames(volatilities) <- c("garch", "gas", "sv", "garch_garch", "garch_gas", "garch_sv", "gas_garch", "gas_gas", "gas_sv", "sv_garch", "sv_gas", "sv_sv")

write.csv(volatilities, paste0("volatilities_", n, "_", distri, ".csv"))