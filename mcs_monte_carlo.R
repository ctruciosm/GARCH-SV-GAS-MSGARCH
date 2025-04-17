################################################################################
#####                         MCS  Monte Carlo                             #####
################################################################################
library(modelconf)
source("Utils_GARCH-GAS-SV.R")

making_tables_mcs <- function(files_names) {
  data_mc <- list()
  for (name in files_names) data_mc[[name]] <- read.csv(name)
  
  
  MSE <- cbind(loss_mse(data_mc[[1]][, "garch_n_garch_n"], data_mc[[1]][, "garch"]), 
    loss_mse(data_mc[[1]][, "garch_n_garch_t"], data_mc[[1]][, "garch"]), 
    loss_mse(data_mc[[1]][, "garch_n_gas_n"], data_mc[[1]][, "garch"]), 
    loss_mse(data_mc[[1]][, "garch_n_gas_t"], data_mc[[1]][, "garch"]), 
    loss_mse(data_mc[[1]][, "garch_n_sv_n"], data_mc[[1]][, "garch"]), 
    loss_mse(data_mc[[1]][, "garch_n_sv_t"], data_mc[[1]][, "garch"]), 
    loss_mse(data_mc[[1]][, "garch_n_ms_n"], data_mc[[1]][, "garch"]), 
    loss_mse(data_mc[[1]][, "garch_n_ms_n"], data_mc[[1]][, "garch"]))
  colnames(MSE) <- c("GARCH-N", "GARCH-T", "GAS-N", "GAS-T", "SV-N", "SV-T", "MS-N", "MS_T")
  GARCH_MSE_mcs = estMCS.quick(MSE, test = "t.range", B = 10000, l = 1, alpha = 0.05)
  
  MSE <- cbind(loss_mse(data_mc[[1]][, "gas_n_garch_n"], data_mc[[1]][, "gas"]), 
    loss_mse(data_mc[[1]][, "gas_n_garch_t"], data_mc[[1]][, "gas"]), 
    loss_mse(data_mc[[1]][, "gas_n_gas_n"], data_mc[[1]][, "gas"]), 
    loss_mse(data_mc[[1]][, "gas_n_gas_t"], data_mc[[1]][, "gas"]), 
    loss_mse(data_mc[[1]][, "gas_n_sv_n"], data_mc[[1]][, "gas"]), 
    loss_mse(data_mc[[1]][, "gas_n_sv_t"], data_mc[[1]][, "gas"]), 
    loss_mse(data_mc[[1]][, "gas_n_ms_n"], data_mc[[1]][, "gas"]), 
    loss_mse(data_mc[[1]][, "gas_n_ms_n"], data_mc[[1]][, "gas"]))
  colnames(MSE) <- c("GARCH-N", "GARCH-T", "GAS-N", "GAS-T", "SV-N", "SV-T", "MS-N", "MS_T")
  GAS_MSE_mcs = estMCS.quick(MSE, test = "t.range", B = 10000, l = 1, alpha = 0.05)
  
  MSE <- cbind(loss_mse(data_mc[[1]][, "sv_n_garch_n"], data_mc[[1]][, "sv"]), 
    loss_mse(data_mc[[1]][, "sv_n_garch_t"], data_mc[[1]][, "sv"]), 
    loss_mse(data_mc[[1]][, "sv_n_gas_n"], data_mc[[1]][, "sv"]), 
    loss_mse(data_mc[[1]][, "sv_n_gas_t"], data_mc[[1]][, "sv"]), 
    loss_mse(data_mc[[1]][, "sv_n_sv_n"], data_mc[[1]][, "sv"]), 
    loss_mse(data_mc[[1]][, "sv_n_sv_t"], data_mc[[1]][, "sv"]), 
    loss_mse(data_mc[[1]][, "sv_n_ms_n"], data_mc[[1]][, "sv"]), 
    loss_mse(data_mc[[1]][, "sv_n_ms_n"], data_mc[[1]][, "sv"]))
  colnames(MSE) <- c("GARCH-N", "GARCH-T", "GAS-N", "GAS-T", "SV-N", "SV-T", "MS-N", "MS_T")
  SV_MSE_mcs = estMCS.quick(MSE, test = "t.range", B = 10000, l = 1, alpha = 0.05)
  
  MSE <- cbind(loss_mse(data_mc[[1]][, "ms_n_garch_n"], data_mc[[1]][, "ms"]), 
    loss_mse(data_mc[[1]][, "ms_n_garch_t"], data_mc[[1]][, "ms"]), 
    loss_mse(data_mc[[1]][, "ms_n_gas_n"], data_mc[[1]][, "ms"]), 
    loss_mse(data_mc[[1]][, "ms_n_gas_t"], data_mc[[1]][, "ms"]), 
    loss_mse(data_mc[[1]][, "ms_n_sv_n"], data_mc[[1]][, "ms"]), 
    loss_mse(data_mc[[1]][, "ms_n_sv_t"], data_mc[[1]][, "ms"]), 
    loss_mse(data_mc[[1]][, "ms_n_ms_n"], data_mc[[1]][, "ms"]), 
    loss_mse(data_mc[[1]][, "ms_n_ms_n"], data_mc[[1]][, "ms"]))
  colnames(MSE) <- c("GARCH-N", "GARCH-T", "GAS-N", "GAS-T", "SV-N", "SV-T", "MS-N", "MS_T")
  MS_MSE_mcs = estMCS.quick(MSE, test = "t.range", B = 10000, l = 1, alpha = 0.05)
  
  
  
  MSE <- cbind(loss_mse(data_mc[[2]][, "garch_t_garch_n"], data_mc[[2]][, "garch"]), 
    loss_mse(data_mc[[2]][, "garch_t_garch_t"], data_mc[[2]][, "garch"]), 
    loss_mse(data_mc[[2]][, "garch_t_gas_n"], data_mc[[2]][, "garch"]), 
    loss_mse(data_mc[[2]][, "garch_t_gas_t"], data_mc[[2]][, "garch"]), 
    loss_mse(data_mc[[2]][, "garch_t_sv_n"], data_mc[[2]][, "garch"]), 
    loss_mse(data_mc[[2]][, "garch_t_sv_t"], data_mc[[2]][, "garch"]), 
    loss_mse(data_mc[[2]][, "garch_t_ms_n"], data_mc[[2]][, "garch"]), 
    loss_mse(data_mc[[2]][, "garch_t_ms_n"], data_mc[[2]][, "garch"]))
  colnames(MSE) <- c("GARCH-N", "GARCH-T", "GAS-N", "GAS-T", "SV-N", "SV-T", "MS-N", "MS_T")
  GARCH_MSE_mcs = estMCS.quick(MSE, test = "t.range", B = 10000, l = 1, alpha = 0.05)
  
  MSE <- cbind(loss_mse(data_mc[[2]][, "gas_t_garch_n"], data_mc[[2]][, "gas"]), 
    loss_mse(data_mc[[2]][, "gas_t_garch_t"], data_mc[[2]][, "gas"]), 
    loss_mse(data_mc[[2]][, "gas_t_gas_n"], data_mc[[2]][, "gas"]), 
    loss_mse(data_mc[[2]][, "gas_t_gas_t"], data_mc[[2]][, "gas"]), 
    loss_mse(data_mc[[2]][, "gas_t_sv_n"], data_mc[[2]][, "gas"]), 
    loss_mse(data_mc[[2]][, "gas_t_sv_t"], data_mc[[2]][, "gas"]), 
    loss_mse(data_mc[[2]][, "gas_t_ms_n"], data_mc[[2]][, "gas"]), 
    loss_mse(data_mc[[2]][, "gas_t_ms_n"], data_mc[[2]][, "gas"]))
  colnames(MSE) <- c("GARCH-N", "GARCH-T", "GAS-N", "GAS-T", "SV-N", "SV-T", "MS-N", "MS_T")
  GAS_MSE_mcs = estMCS.quick(MSE, test = "t.range", B = 10000, l = 1, alpha = 0.05)
  
  MSE <- cbind(loss_mse(data_mc[[2]][, "sv_t_garch_n"], data_mc[[2]][, "sv"]), 
    loss_mse(data_mc[[2]][, "sv_t_garch_t"], data_mc[[2]][, "sv"]), 
    loss_mse(data_mc[[2]][, "sv_t_gas_n"], data_mc[[2]][, "sv"]), 
    loss_mse(data_mc[[2]][, "sv_t_gas_t"], data_mc[[2]][, "sv"]), 
    loss_mse(data_mc[[2]][, "sv_t_sv_n"], data_mc[[2]][, "sv"]), 
    loss_mse(data_mc[[2]][, "sv_t_sv_t"], data_mc[[2]][, "sv"]), 
    loss_mse(data_mc[[2]][, "sv_t_ms_n"], data_mc[[2]][, "sv"]), 
    loss_mse(data_mc[[2]][, "sv_t_ms_n"], data_mc[[2]][, "sv"]))
  colnames(MSE) <- c("GARCH-N", "GARCH-T", "GAS-N", "GAS-T", "SV-N", "SV-T", "MS-N", "MS_T")
  SV_MSE_mcs = estMCS.quick(MSE, test = "t.range", B = 10000, l = 1, alpha = 0.05)
  
  MSE <- cbind(loss_mse(data_mc[[2]][, "ms_t_garch_n"], data_mc[[2]][, "ms"]), 
    loss_mse(data_mc[[2]][, "ms_t_garch_t"], data_mc[[2]][, "ms"]), 
    loss_mse(data_mc[[2]][, "ms_t_gas_n"], data_mc[[2]][, "ms"]), 
    loss_mse(data_mc[[2]][, "ms_t_gas_t"], data_mc[[2]][, "ms"]), 
    loss_mse(data_mc[[2]][, "ms_t_sv_n"], data_mc[[2]][, "ms"]), 
    loss_mse(data_mc[[2]][, "ms_t_sv_t"], data_mc[[2]][, "ms"]), 
    loss_mse(data_mc[[2]][, "ms_t_ms_n"], data_mc[[2]][, "ms"]), 
    loss_mse(data_mc[[2]][, "ms_t_ms_n"], data_mc[[2]][, "ms"]))
  colnames(MSE) <- c("GARCH-N", "GARCH-T", "GAS-N", "GAS-T", "SV-N", "SV-T", "MS-N", "MS_T")
  MS_MSE_mcs = estMCS.quick(MSE, test = "t.range", B = 10000, l = 1, alpha = 0.05)
  
  
}

# Model Confidence Set
run_mcs <- function(data_index, loss_func, pattern = '*FALSE_BR.csv') {
  files_names <- list.files(path = './MonteCarlo', pattern = pattern, full.names = TRUE)
  data_mc <- lapply(files_names, read.csv)
  names(data_mc) <- files_names
  print(files_names[data_index])
  data <- data_mc[[data_index]]
  target_col <- c("garch", "gas", "sv", "ms")
  for (k in target_col) {
    if (data_index %% 2 == 0) {
      prefix <- paste0(k, "_t")
    } else {
      prefix <- paste0(k, "_n")
    }
    col_names <- paste0(prefix, c("_garch_n", "_garch_t", "_gas_n", "_gas_t", "_sv_n", "_sv_t", "_ms_n", "_ms_n"))
    loss_vals <- sapply(col_names, function(col) loss_func(data[[col]], data[[k]]))
    colnames(loss_vals) <-  c("GARCH-N", "GARCH-T", "GAS-N", "GAS-T", "SV-N", "SV-T", "MS-N", "MS-T")
    mcs <- estMCS.quick(loss_vals, test = "t.range", B = 10000, l = 1, alpha = 0.1)
    print(paste0("DGP_", k))
    print(colnames(loss_vals)[mcs])
  }
}


# Lê arquivos e monta lista nomeada


run_mcs(3, loss_mse, '*FALSE_BR.csv')


