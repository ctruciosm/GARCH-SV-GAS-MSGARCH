################################################################################
#####                       Tables Monte Carlo                             #####
################################################################################
source("Utils_GARCH-GAS-SV.R")


files_names <- list.files(path = './MonteCarlo', pattern = '*.csv', full.names = TRUE)
data_mc <- list()
for (name in files_names) data_mc[[name]] <- read.csv(name)


data_mc[[1]][, "garch"]
data_mc[[1]][, "gas"]
data_mc[[1]][, "sv"]

results <- round(rbind(metrics(data_mc[[1]][, "garch_garch"], data_mc[[1]][, "garch"]),
                 metrics(data_mc[[1]][, "garch_gas"], data_mc[[1]][, "garch"]),
                 metrics(data_mc[[1]][, "garch_sv"], data_mc[[1]][, "garch"]),
                 metrics(data_mc[[1]][, "gas_garch"], data_mc[[1]][, "gas"]),
                 metrics(data_mc[[1]][, "gas_gas"], data_mc[[1]][, "gas"]),
                 metrics(data_mc[[1]][, "gas_sv"], data_mc[[1]][, "gas"]),
                 metrics(data_mc[[1]][, "sv_garch"], data_mc[[1]][, "sv"]),
                 metrics(data_mc[[1]][, "sv_gas"], data_mc[[1]][, "sv"]),
                 metrics(data_mc[[1]][, "sv_sv"], data_mc[[1]][, "sv"])), 4)

colnames(results) <- c("MSE", "QLIKE1", "QLIKE2", "MSE LOG", "MSE SD", "MSE PROP", "MAE", "MAE LOG", "MAE SD", "MAE PROP")
row.names(results) <- c("GARCH - GARCH", "GARCH - GAS", "GARCH - SV",
                        "GAS - GARCH", "GAS - GAS", "GAS - SV",
                        "SV - GARCH", "SV - GAS", "SV - SV")


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