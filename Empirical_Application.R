################################################################################
#####                   Monte Carlo Simulation                             #####
################################################################################
library(dplyr)
library(readxl)
library(stringr)
library(rugarch)
library(GAS)
library(stochvol)
library(stochvolTMB)
library(MSGARCH)
source("DGPs.R")
source("Utils_GARCH-GAS-SV.R")


# Data
data <- read_excel("./Data/capire_daily_returns.xlsx", skip = 3, col_types = c("date", rep("numeric", 29)), na = c("", "-", NA)) |> 
  filter(Data > "2010-01-01" & Data < "2025-01-01") |> 
  filter(!if_all(where(is.numeric), is.na)) |> 
  select(where(~ !any(is.na(.x)))) |> 
  rename_with(~ str_remove(.x, "^(?s).*prov\n"), -Data)

ins <- 2500
oos <- nrow(data) - ins

# Specs
garch_spec_n <- ugarchspec(variance.model = list(model = 'sGARCH', garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), distribution.model = "norm")
garch_spec_t <- ugarchspec(variance.model = list(model = 'sGARCH', garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = FALSE), distribution.model = "std")
gas_spec_n <- UniGASSpec(Dist = "norm", ScalingType = "Identity", GASPar = list(locate = FALSE, scale = TRUE, shape = FALSE))
gas_spec_t <- UniGASSpec(Dist = "std", ScalingType = "Identity", GASPar = list(locate = FALSE, scale = TRUE, shape = FALSE))
ms_spec_n <- CreateSpec(variance.spec = list(model = c("sGARCH", "sGARCH")), switch.spec = list(do.mix = FALSE), distribution.spec = list(distribution = c("norm", "norm")))
ms_spec_t <- CreateSpec(variance.spec = list(model = c("sGARCH", "sGARCH")), switch.spec = list(do.mix = FALSE), distribution.spec = list(distribution = c("std", "std")), constraint.spec = list(regime.const = c("nu")))


garch_n_fore <- garch_t_fore <- gas_n_fore <- gas_t_fore <- ms_n_fore <- ms_t_fore <- sv_n_fore <- sv_t_fore <- sv_n_fore_b <- sv_t_fore_b <- matrix(0, nrow = oos, ncol = ncol(data) - 1)

for (i in 1:oos) {
  print(i)
  returns <- data[i:(i + ins - 1), -1]
  mu <- apply(returns, 2, mean)
  returns_c <- scale(returns, scale = FALSE)
  
  garch_n_fore[i, ] <- apply(returns_c, 2, function(x) {
    tryCatch({
      sigma2 <- ugarchforecast(ugarchfit(garch_spec_n, x, solver = "hybrid"), n.ahead = 1)@forecast$sigmaFor[1]^2
      return(sigma2)
    }, error = function(e) { return(NA) }) })
  
  garch_t_fore[i, ] <- apply(returns_c, 2, function(x) {
    tryCatch({
      sigma2 <- ugarchforecast(ugarchfit(garch_spec_t, x, solver = "hybrid"), n.ahead = 1)@forecast$sigmaFor[1]^2
      return(sigma2)
    }, error = function(e) { return(NA) }) })
  
  gas_n_fore[i, ]  <- apply(returns_c, 2, function(x) {
    tryCatch({
      sigma2 <- UniGASFor(gasfit(gas_spec_n, x), H = 1)@Forecast$PointForecast[, 2]
      return(sigma2)
    }, error = function(e) { return(NA) }) })
  
  gas_t_fore[i, ]  <- apply(returns_c, 2, function(x) {
    tryCatch({
      fit <- gasfit(gas_spec_t, x)
      sigma2 <- UniGASFor(fit, H = 1)@Forecast$PointForecast[, 2] * fit@GASDyn$mTheta[3, 1] /(fit@GASDyn$mTheta[3, 1] - 2)
      return(sigma2)
    }, error = function(e) { return(NA) }) })
  
  ms_n_fore[i, ]  <- apply(returns_c, 2, function(x) {
    tryCatch({
      sigma2 <- predict(msgarchfit(ms_spec_n, x), nahead = 1)$vol^2
      return(sigma2)
    }, error = function(e) { return(NA) }) })
  
  ms_t_fore[i, ]  <- apply(returns_c, 2, function(x) {
    tryCatch({
      sigma2 <- predict(msgarchfit(ms_spec_t, x), nahead = 1)$vol^2
      return(sigma2)
    }, error = function(e) { return(NA) }) })
  
  sv_n_fore_b[i, ]  <- apply(returns_c, 2, function(x) {
    tryCatch({
      sigma2 <- median(predvola(predict(svsample(x), 1))^2)
      return(sigma2)
    }, error = function(e) { return(NA) }) })
  
  sv_t_fore_b[i, ]  <- apply(returns_c, 2, function(x) {
    tryCatch({
      sigma2 <- median(predvola(predict(svtsample(x), 1))^2)
      return(sigma2)
    }, error = function(e) { return(NA) }) })
  
  sv_n_fore[i, ]  <- apply(returns_c, 2, function(x) {
    tryCatch({
      sv_n_fit <- estimate_parameters(x, model = "gaussian", silent = TRUE)
      sigma2 <- median(predict(sv_n_fit, steps = 1)$h_exp)
      return(sigma2)
    }, error = function(e) { return(NA) }) })
  
  sv_t_fore[i, ]  <- apply(returns_c, 2, function(x) {
    tryCatch({
      sv_t_fit <- estimate_parameters(x, model = "t", silent = TRUE)
      sigma2 <- median(predict(sv_t_fit, steps = 1)$h_exp)
      return(sigma2)
    }, error = function(e) { return(NA) }) })
}


write.csv(garch_n_fore, "garch_n_fore.csv")
write.csv(garch_t_fore, "garch_t_fore.csv")
write.csv(gas_n_fore, "gas_n_fore.csv")
write.csv(gas_t_fore, "gas_t_fore.csv")
write.csv(ms_n_fore, "ms_n_fore.csv")
write.csv(ms_t_fore, "ms_t_fore.csv")
write.csv(sv_n_fore_b, "sv_n_fore_b.csv")
write.csv(sv_t_fore_b, "sv_t_fore_b.csv")
write.csv(sv_n_fore, "sv_n_fore.csv")
write.csv(sv_t_fore, "sv_t_fore.csv")
write.csv(data[(1 + ins):nrow(data), 1], "datas_oos.csv")
