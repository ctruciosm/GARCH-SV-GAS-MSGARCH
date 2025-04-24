################################################################################
#####                  Tables Empirical Application                        #####
################################################################################
library(dplyr)
library(lubridate)
library(modelconf)
source("Utils_GARCH-GAS-SV.R")



datas_oos <- read.csv("./Empirical_Application/datas_oos.csv")[, -1]
garch_n <- read.csv("./Empirical_Application/garch_n_fore.csv")[, -1]
garch_t <- read.csv("./Empirical_Application/garch_t_fore.csv")[, -1]
gas_n <- read.csv("./Empirical_Application/gas_n_fore.csv")[, -1]
gas_t <- read.csv("./Empirical_Application/gas_t_fore.csv")[, -1]
ms_n <- read.csv("./Empirical_Application/ms_n_fore.csv")[, -1]
ms_t <- read.csv("./Empirical_Application/ms_t_fore.csv")[, -1]
sv_n <- read.csv("./Empirical_Application/sv_n_fore.csv")[, -1]
sv_t <- read.csv("./Empirical_Application/sv_t_fore.csv")[, -1]
#sv_n_b <- read.csv("./Empirical_Application/sv_n_fore_b.csv")[, -1]
#sv_t_b <- read.csv("./Empirical_Application/sv_t_fore_b.csv")[, -1]
rv_oos <- read.csv("./Data/capire_bpv_5.csv", sep = ";", dec = ",") |> mutate(Date = lubridate::dmy(Date)) |> filter(Date > "2019-12-08") |> select(-Date)


colnames(garch_n) <- c("MMM", "AMZN", "AXP", "AMGN", "AAPL", "BA", "CAT", "CVX", "CSCO", "KO", "HD", "HON", "INTC", "IBM", "JNJ", "JPM", "MCD", "MRK", "MSFT", "NKE", "PG", "GS", "TRV", "UNH", "VZ", "V", "WMT", "DIS", "CRM")
colnames(garch_t) <- c("MMM", "AMZN", "AXP", "AMGN", "AAPL", "BA", "CAT", "CVX", "CSCO", "KO", "HD", "HON", "INTC", "IBM", "JNJ", "JPM", "MCD", "MRK", "MSFT", "NKE", "PG", "GS", "TRV", "UNH", "VZ", "V", "WMT", "DIS", "CRM")
colnames(gas_n) <- c("MMM", "AMZN", "AXP", "AMGN", "AAPL", "BA", "CAT", "CVX", "CSCO", "KO", "HD", "HON", "INTC", "IBM", "JNJ", "JPM", "MCD", "MRK", "MSFT", "NKE", "PG", "GS", "TRV", "UNH", "VZ", "V", "WMT", "DIS", "CRM")
colnames(gas_t) <- c("MMM", "AMZN", "AXP", "AMGN", "AAPL", "BA", "CAT", "CVX", "CSCO", "KO", "HD", "HON", "INTC", "IBM", "JNJ", "JPM", "MCD", "MRK", "MSFT", "NKE", "PG", "GS", "TRV", "UNH", "VZ", "V", "WMT", "DIS", "CRM")
colnames(ms_n) <- c("MMM", "AMZN", "AXP", "AMGN", "AAPL", "BA", "CAT", "CVX", "CSCO", "KO", "HD", "HON", "INTC", "IBM", "JNJ", "JPM", "MCD", "MRK", "MSFT", "NKE", "PG", "GS", "TRV", "UNH", "VZ", "V", "WMT", "DIS", "CRM")
colnames(ms_t) <- c("MMM", "AMZN", "AXP", "AMGN", "AAPL", "BA", "CAT", "CVX", "CSCO", "KO", "HD", "HON", "INTC", "IBM", "JNJ", "JPM", "MCD", "MRK", "MSFT", "NKE", "PG", "GS", "TRV", "UNH", "VZ", "V", "WMT", "DIS", "CRM")
colnames(sv_n) <- c("MMM", "AMZN", "AXP", "AMGN", "AAPL", "BA", "CAT", "CVX", "CSCO", "KO", "HD", "HON", "INTC", "IBM", "JNJ", "JPM", "MCD", "MRK", "MSFT", "NKE", "PG", "GS", "TRV", "UNH", "VZ", "V", "WMT", "DIS", "CRM")
colnames(sv_t) <- c("MMM", "AMZN", "AXP", "AMGN", "AAPL", "BA", "CAT", "CVX", "CSCO", "KO", "HD", "HON", "INTC", "IBM", "JNJ", "JPM", "MCD", "MRK", "MSFT", "NKE", "PG", "GS", "TRV", "UNH", "VZ", "V", "WMT", "DIS", "CRM")
#colnames(sv_n_b) <- c("MMM", "AMZN", "AXP", "AMGN", "AAPL", "BA", "CAT", "CVX", "CSCO", "KO", "HD", "HON", "INTC", "IBM", "JNJ", "JPM", "MCD", "MRK", "MSFT", "NKE", "PG", "GS", "TRV", "UNH", "VZ", "V", "WMT", "DIS", "CRM")
#colnames(sv_t_b) <- c("MMM", "AMZN", "AXP", "AMGN", "AAPL", "BA", "CAT", "CVX", "CSCO", "KO", "HD", "HON", "INTC", "IBM", "JNJ", "JPM", "MCD", "MRK", "MSFT", "NKE", "PG", "GS", "TRV", "UNH", "VZ", "V", "WMT", "DIS", "CRM")


m_mse <- matrix(0, ncol = 8, nrow = ncol(garch_n))
m_qlike <- matrix(0, ncol = 8, nrow = ncol(garch_n))

for (i in 1:ncol(garch_n)) {
  m_mse[i, ] <- c(mean(loss_mse(rv_oos[, colnames(garch_n)[i]], garch_n[, colnames(garch_n)[i]])),
                  mean(loss_mse(rv_oos[, colnames(garch_n)[i]], garch_t[, colnames(garch_n)[i]])),
                  mean(loss_mse(rv_oos[, colnames(garch_n)[i]], gas_n[, colnames(garch_n)[i]])),
                  mean(loss_mse(rv_oos[, colnames(garch_n)[i]], gas_t[, colnames(garch_n)[i]])),
                  mean(loss_mse(rv_oos[, colnames(garch_n)[i]], ms_n[, colnames(garch_n)[i]])),
                  mean(loss_mse(rv_oos[, colnames(garch_n)[i]], ms_t[, colnames(garch_n)[i]])),
                  mean(loss_mse(rv_oos[, colnames(garch_n)[i]], sv_n[, colnames(garch_n)[i]])),
                  mean(loss_mse(rv_oos[, colnames(garch_n)[i]], sv_t[, colnames(garch_n)[i]])))
  m_qlike[i, ] <- c(mean(loss_qlike(rv_oos[, colnames(garch_n)[i]], garch_n[, colnames(garch_n)[i]])),
                    mean(loss_qlike(rv_oos[, colnames(garch_n)[i]], garch_t[, colnames(garch_n)[i]])),
                    mean(loss_qlike(rv_oos[, colnames(garch_n)[i]], gas_n[, colnames(garch_n)[i]])),
                    mean(loss_qlike(rv_oos[, colnames(garch_n)[i]], gas_t[, colnames(garch_n)[i]])),
                    mean(loss_qlike(rv_oos[, colnames(garch_n)[i]], ms_n[, colnames(garch_n)[i]])),
                    mean(loss_qlike(rv_oos[, colnames(garch_n)[i]], ms_t[, colnames(garch_n)[i]])),
                    mean(loss_qlike(rv_oos[, colnames(garch_n)[i]], sv_n[, colnames(garch_n)[i]])),
                    mean(loss_qlike(rv_oos[, colnames(garch_n)[i]], sv_t[, colnames(garch_n)[i]])))
  
}


col_names <- c("GARCH-N", "GARCH-T", "GAS-N", "GAS-T", "MS-N", "MS-T", "SV-N", "SV-T")
colnames(m_mse) <- col_names
colnames(m_qlike) <- col_names
row.names(m_mse) <- colnames(garch_n)
row.names(m_qlike) <- colnames(garch_n)
xtable::xtable(m_mse, digits = 3)
xtable::xtable(m_qlike, digits = 3)


################################################################################
#####                              MSC                                     #####
################################################################################

stock_names <- colnames(garch_n)
for (i in 1:ncol(garch_n)) {
  m_mse <- cbind(loss_mse(rv_oos[, stock_names[i]], garch_n[, stock_names[i]]),
    loss_mse(rv_oos[, stock_names[i]], garch_t[, stock_names[i]]),
    loss_mse(rv_oos[, stock_names[i]], gas_n[, stock_names[i]]),
    loss_mse(rv_oos[, stock_names[i]], gas_t[, stock_names[i]]),
    loss_mse(rv_oos[, stock_names[i]], ms_n[, stock_names[i]]),
    loss_mse(rv_oos[, stock_names[i]], ms_t[, stock_names[i]]),
    loss_mse(rv_oos[, stock_names[i]], sv_n[, stock_names[i]]),
    loss_mse(rv_oos[, stock_names[i]], sv_t[, stock_names[i]]))
  
  m_qlike <- cbind(loss_qlike(rv_oos[, stock_names[i]], garch_n[, stock_names[i]]),
    loss_qlike(rv_oos[, stock_names[i]], garch_t[, stock_names[i]]),
    loss_qlike(rv_oos[, stock_names[i]], gas_n[, stock_names[i]]),
    loss_qlike(rv_oos[, stock_names[i]], gas_t[, stock_names[i]]),
    loss_qlike(rv_oos[, stock_names[i]], ms_n[, stock_names[i]]),
    loss_qlike(rv_oos[, stock_names[i]], ms_t[, stock_names[i]]),
    loss_qlike(rv_oos[, stock_names[i]], sv_n[, stock_names[i]]),
    loss_qlike(rv_oos[, stock_names[i]], sv_t[, stock_names[i]]))
  
  colnames(m_mse) <- col_names[1:8]
  colnames(m_qlike) <- col_names[1:8]
  
  mcs_mse <- estMCS.quick(m_mse, test = "t.range", B = 10000, l = 21, alpha = 0.25)
  mcs_qlike <- estMCS.quick(m_qlike, test = "t.range", B = 10000, l = 21, alpha = 0.25)

  print(stock_names[i])
  print("MSE")
  print(colnames(m_mse)[mcs_mse])
  print("QLIKE")
  print(colnames(m_qlike)[mcs_qlike])
}
