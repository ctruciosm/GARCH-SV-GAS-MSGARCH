################################################################################
#####                       Tables Monte Carlo                             #####
################################################################################
library(dplyr)
library(ggplot2)
library(tidyr)
library(modelconf)
library(wrMisc)
source("Utils_GARCH-GAS-SV.R")


making_tables <- function(files_names) {
  data_mc <- list()
  for (name in files_names) data_mc[[name]] <- read.csv(name)
  
  
  simul <- c(rep(data_mc[[1]][, "garch"], 8), rep(data_mc[[1]][, "gas"], 8), rep(data_mc[[1]][, "sv"], 8), rep(data_mc[[1]][, "ms"], 8))
  dgp <- c(rep("GARCH", 8000), rep("GAS", 8000), rep("SV", 8000), rep("MS", 8000))
  fore <- c(data_mc[[1]][, "garch_n_garch_n"], data_mc[[1]][, "garch_n_garch_t"], data_mc[[1]][, "garch_n_gas_n"], data_mc[[1]][, "garch_n_gas_t"],
    data_mc[[1]][, "garch_n_sv_n"], data_mc[[1]][, "garch_n_sv_t"], data_mc[[1]][, "garch_n_ms_n"],  data_mc[[1]][, "garch_n_ms_t"],
    data_mc[[1]][, "gas_n_garch_n"], data_mc[[1]][, "gas_n_garch_t"], data_mc[[1]][, "gas_n_gas_n"], data_mc[[1]][, "gas_n_gas_t"],
    data_mc[[1]][, "gas_n_sv_n"], data_mc[[1]][, "gas_n_sv_t"], data_mc[[1]][, "gas_n_ms_n"], data_mc[[1]][, "gas_n_ms_t"],
    data_mc[[1]][, "sv_n_garch_n"], data_mc[[1]][, "sv_n_garch_t"], data_mc[[1]][, "sv_n_gas_n"], data_mc[[1]][, "sv_n_gas_t"],
    data_mc[[1]][, "sv_n_sv_n"], data_mc[[1]][, "sv_n_sv_t"], data_mc[[1]][, "sv_n_ms_n"], data_mc[[1]][, "sv_n_ms_t"],
    data_mc[[1]][, "ms_n_garch_n"], data_mc[[1]][, "ms_n_garch_t"], data_mc[[1]][, "ms_n_gas_n"], data_mc[[1]][, "ms_n_gas_t"],
    data_mc[[1]][, "ms_n_sv_n"], data_mc[[1]][, "ms_n_sv_t"], data_mc[[1]][, "ms_n_ms_n"], data_mc[[1]][, "ms_n_ms_t"])
  model <- rep(c(rep("GARCH-N", 1000), rep("GARCH-T", 1000), rep("GAS-N", 1000), rep("GAS-T", 1000), rep("SV-N", 1000), rep("SV-T", 1000), rep("MS-N", 1000), rep("MS-T", 1000)), 4)
  n <- rep(1000, 32000)
  distri <- rep("N", 32000)
  mc_1 <- data.frame(simul, dgp, fore, model, n, distri)
  
  simul <- c(rep(data_mc[[2]][, "garch"], 8), rep(data_mc[[2]][, "gas"], 8), rep(data_mc[[2]][, "sv"], 8), rep(data_mc[[2]][, "ms"], 8))
  dgp <- c(rep("GARCH", 8000), rep("GAS", 8000), rep("SV", 8000), rep("MS", 8000))
  fore <- c(data_mc[[2]][, "garch_t_garch_n"], data_mc[[2]][, "garch_t_garch_t"], data_mc[[2]][, "garch_t_gas_n"], data_mc[[2]][, "garch_t_gas_t"],
    data_mc[[2]][, "garch_t_sv_n"], data_mc[[2]][, "garch_t_sv_t"], data_mc[[2]][, "garch_t_ms_n"],  data_mc[[2]][, "garch_t_ms_t"],
    data_mc[[2]][, "gas_t_garch_n"], data_mc[[2]][, "gas_t_garch_t"], data_mc[[2]][, "gas_t_gas_n"], data_mc[[2]][, "gas_t_gas_t"],
    data_mc[[2]][, "gas_t_sv_n"], data_mc[[2]][, "gas_t_sv_t"], data_mc[[2]][, "gas_t_ms_n"], data_mc[[2]][, "gas_t_ms_t"],
    data_mc[[2]][, "sv_t_garch_n"], data_mc[[2]][, "sv_t_garch_t"], data_mc[[2]][, "sv_t_gas_n"], data_mc[[2]][, "sv_t_gas_t"],
    data_mc[[2]][, "sv_t_sv_n"], data_mc[[2]][, "sv_t_sv_t"], data_mc[[2]][, "sv_t_ms_n"], data_mc[[2]][, "sv_t_ms_t"],
    data_mc[[2]][, "ms_t_garch_n"], data_mc[[2]][, "ms_t_garch_t"], data_mc[[2]][, "ms_t_gas_n"], data_mc[[2]][, "ms_t_gas_t"],
    data_mc[[2]][, "ms_t_sv_n"], data_mc[[2]][, "ms_t_sv_t"], data_mc[[2]][, "ms_t_ms_n"], data_mc[[2]][, "ms_t_ms_t"])
  model <- rep(c(rep("GARCH-N", 1000), rep("GARCH-T", 1000), rep("GAS-N", 1000), rep("GAS-T", 1000), rep("SV-N", 1000), rep("SV-T", 1000), rep("MS-N", 1000), rep("MS-T", 1000)), 4)
  n <- rep(1000, 32000)
  distri <- rep("T", 32000)
  mc_2 <- data.frame(simul, dgp, fore, model, n, distri)
  
  simul <- c(rep(data_mc[[3]][, "garch"], 8), rep(data_mc[[3]][, "gas"], 8), rep(data_mc[[3]][, "sv"], 8), rep(data_mc[[3]][, "ms"], 8))
  dgp <- c(rep("GARCH", 8000), rep("GAS", 8000), rep("SV", 8000), rep("MS", 8000))
  fore <- c(data_mc[[3]][, "garch_n_garch_n"], data_mc[[3]][, "garch_n_garch_t"], data_mc[[3]][, "garch_n_gas_n"], data_mc[[3]][, "garch_n_gas_t"],
    data_mc[[3]][, "garch_n_sv_n"], data_mc[[3]][, "garch_n_sv_t"], data_mc[[3]][, "garch_n_ms_n"],  data_mc[[3]][, "garch_n_ms_t"],
    data_mc[[3]][, "gas_n_garch_n"], data_mc[[3]][, "gas_n_garch_t"], data_mc[[3]][, "gas_n_gas_n"], data_mc[[3]][, "gas_n_gas_t"],
    data_mc[[3]][, "gas_n_sv_n"], data_mc[[3]][, "gas_n_sv_t"], data_mc[[3]][, "gas_n_ms_n"], data_mc[[3]][, "gas_n_ms_t"],
    data_mc[[3]][, "sv_n_garch_n"], data_mc[[3]][, "sv_n_garch_t"], data_mc[[3]][, "sv_n_gas_n"], data_mc[[3]][, "sv_n_gas_t"],
    data_mc[[3]][, "sv_n_sv_n"], data_mc[[3]][, "sv_n_sv_t"], data_mc[[3]][, "sv_n_ms_n"], data_mc[[3]][, "sv_n_ms_t"],
    data_mc[[3]][, "ms_n_garch_n"], data_mc[[3]][, "ms_n_garch_t"], data_mc[[3]][, "ms_n_gas_n"], data_mc[[3]][, "ms_n_gas_t"],
    data_mc[[3]][, "ms_n_sv_n"], data_mc[[3]][, "ms_n_sv_t"], data_mc[[3]][, "ms_n_ms_n"], data_mc[[3]][, "ms_n_ms_t"])
  model <- rep(c(rep("GARCH-N", 1000), rep("GARCH-T", 1000), rep("GAS-N", 1000), rep("GAS-T", 1000), rep("SV-N", 1000), rep("SV-T", 1000), rep("MS-N", 1000), rep("MS-T", 1000)), 4)
  n <- rep(2500, 32000)
  distri <- rep("N", 32000)
  mc_3 <- data.frame(simul, dgp, fore, model, n, distri)
  
  simul <- c(rep(data_mc[[4]][, "garch"], 8), rep(data_mc[[4]][, "gas"], 8), rep(data_mc[[4]][, "sv"], 8), rep(data_mc[[4]][, "ms"], 8))
  dgp <- c(rep("GARCH", 8000), rep("GAS", 8000), rep("SV", 8000), rep("MS", 8000))
  fore <- c(data_mc[[4]][, "garch_t_garch_n"], data_mc[[4]][, "garch_t_garch_t"], data_mc[[4]][, "garch_t_gas_n"], data_mc[[4]][, "garch_t_gas_t"],
    data_mc[[4]][, "garch_t_sv_n"], data_mc[[4]][, "garch_t_sv_t"], data_mc[[4]][, "garch_t_ms_n"],  data_mc[[4]][, "garch_t_ms_t"],
    data_mc[[4]][, "gas_t_garch_n"], data_mc[[4]][, "gas_t_garch_t"], data_mc[[4]][, "gas_t_gas_n"], data_mc[[4]][, "gas_t_gas_t"],
    data_mc[[4]][, "gas_t_sv_n"], data_mc[[4]][, "gas_t_sv_t"], data_mc[[4]][, "gas_t_ms_n"], data_mc[[4]][, "gas_t_ms_t"],
    data_mc[[4]][, "sv_t_garch_n"], data_mc[[4]][, "sv_t_garch_t"], data_mc[[4]][, "sv_t_gas_n"], data_mc[[4]][, "sv_t_gas_t"],
    data_mc[[4]][, "sv_t_sv_n"], data_mc[[4]][, "sv_t_sv_t"], data_mc[[4]][, "sv_t_ms_n"], data_mc[[4]][, "sv_t_ms_t"],
    data_mc[[4]][, "ms_t_garch_n"], data_mc[[4]][, "ms_t_garch_t"], data_mc[[4]][, "ms_t_gas_n"], data_mc[[4]][, "ms_t_gas_t"],
    data_mc[[4]][, "ms_t_sv_n"], data_mc[[4]][, "ms_t_sv_t"], data_mc[[4]][, "ms_t_ms_n"], data_mc[[4]][, "ms_t_ms_t"])
  model <- rep(c(rep("GARCH-N", 1000), rep("GARCH-T", 1000), rep("GAS-N", 1000), rep("GAS-T", 1000), rep("SV-N", 1000), rep("SV-T", 1000), rep("MS-N", 1000), rep("MS-T", 1000)), 4)
  n <- rep(2500, 32000)
  distri <- rep("T", 32000)
  mc_4 <- data.frame(simul, dgp, fore, model, n, distri)
  
  simul <- c(rep(data_mc[[5]][, "garch"], 8), rep(data_mc[[5]][, "gas"], 8), rep(data_mc[[5]][, "sv"], 8), rep(data_mc[[5]][, "ms"], 8))
  dgp <- c(rep("GARCH", 8000), rep("GAS", 8000), rep("SV", 8000), rep("MS", 8000))
  fore <- c(data_mc[[5]][, "garch_n_garch_n"], data_mc[[5]][, "garch_n_garch_t"], data_mc[[5]][, "garch_n_gas_n"], data_mc[[5]][, "garch_n_gas_t"],
    data_mc[[5]][, "garch_n_sv_n"], data_mc[[5]][, "garch_n_sv_t"], data_mc[[5]][, "garch_n_ms_n"],  data_mc[[5]][, "garch_n_ms_t"],
    data_mc[[5]][, "gas_n_garch_n"], data_mc[[5]][, "gas_n_garch_t"], data_mc[[5]][, "gas_n_gas_n"], data_mc[[5]][, "gas_n_gas_t"],
    data_mc[[5]][, "gas_n_sv_n"], data_mc[[5]][, "gas_n_sv_t"], data_mc[[5]][, "gas_n_ms_n"], data_mc[[5]][, "gas_n_ms_t"],
    data_mc[[5]][, "sv_n_garch_n"], data_mc[[5]][, "sv_n_garch_t"], data_mc[[5]][, "sv_n_gas_n"], data_mc[[5]][, "sv_n_gas_t"],
    data_mc[[5]][, "sv_n_sv_n"], data_mc[[5]][, "sv_n_sv_t"], data_mc[[5]][, "sv_n_ms_n"], data_mc[[5]][, "sv_n_ms_t"],
    data_mc[[5]][, "ms_n_garch_n"], data_mc[[5]][, "ms_n_garch_t"], data_mc[[5]][, "ms_n_gas_n"], data_mc[[5]][, "ms_n_gas_t"],
    data_mc[[5]][, "ms_n_sv_n"], data_mc[[5]][, "ms_n_sv_t"], data_mc[[5]][, "ms_n_ms_n"], data_mc[[5]][, "ms_n_ms_t"])
  model <- rep(c(rep("GARCH-N", 1000), rep("GARCH-T", 1000), rep("GAS-N", 1000), rep("GAS-T", 1000), rep("SV-N", 1000), rep("SV-T", 1000), rep("MS-N", 1000), rep("MS-T", 1000)), 4)
  n <- rep(500, 32000)
  distri <- rep("N", 32000)
  mc_5 <- data.frame(simul, dgp, fore, model, n, distri)
  
  simul <- c(rep(data_mc[[6]][, "garch"], 8), rep(data_mc[[6]][, "gas"], 8), rep(data_mc[[6]][, "sv"], 8), rep(data_mc[[6]][, "ms"], 8))
  dgp <- c(rep("GARCH", 8000), rep("GAS", 8000), rep("SV", 8000), rep("MS", 8000))
  fore <- c(data_mc[[6]][, "garch_t_garch_n"], data_mc[[6]][, "garch_t_garch_t"], data_mc[[6]][, "garch_t_gas_n"], data_mc[[6]][, "garch_t_gas_t"],
    data_mc[[6]][, "garch_t_sv_n"], data_mc[[6]][, "garch_t_sv_t"], data_mc[[6]][, "garch_t_ms_n"],  data_mc[[6]][, "garch_t_ms_t"],
    data_mc[[6]][, "gas_t_garch_n"], data_mc[[6]][, "gas_t_garch_t"], data_mc[[6]][, "gas_t_gas_n"], data_mc[[6]][, "gas_t_gas_t"],
    data_mc[[6]][, "gas_t_sv_n"], data_mc[[6]][, "gas_t_sv_t"], data_mc[[6]][, "gas_t_ms_n"], data_mc[[6]][, "gas_t_ms_t"],
    data_mc[[6]][, "sv_t_garch_n"], data_mc[[6]][, "sv_t_garch_t"], data_mc[[6]][, "sv_t_gas_n"], data_mc[[6]][, "sv_t_gas_t"],
    data_mc[[6]][, "sv_t_sv_n"], data_mc[[6]][, "sv_t_sv_t"], data_mc[[6]][, "sv_t_ms_n"], data_mc[[6]][, "sv_t_ms_t"],
    data_mc[[6]][, "ms_t_garch_n"], data_mc[[6]][, "ms_t_garch_t"], data_mc[[6]][, "ms_t_gas_n"], data_mc[[6]][, "ms_t_gas_t"],
    data_mc[[6]][, "ms_t_sv_n"], data_mc[[6]][, "ms_t_sv_t"], data_mc[[6]][, "ms_t_ms_n"], data_mc[[6]][, "ms_t_ms_t"])
  model <- rep(c(rep("GARCH-N", 1000), rep("GARCH-T", 1000), rep("GAS-N", 1000), rep("GAS-T", 1000), rep("SV-N", 1000), rep("SV-T", 1000), rep("MS-N", 1000), rep("MS-T", 1000)), 4)
  n <- rep(500, 32000)
  distri <- rep("T", 32000)
  mc_6 <- data.frame(simul, dgp, fore, model, n, distri)
  
  results_1 <- round(
    rbind(
      metrics(data_mc[[1]][, "garch_n_garch_n"], data_mc[[1]][, "garch"]),
      metrics(data_mc[[1]][, "garch_n_garch_t"], data_mc[[1]][, "garch"]),
      metrics(data_mc[[1]][, "garch_n_gas_n"], data_mc[[1]][, "garch"]),
      metrics(data_mc[[1]][, "garch_n_gas_t"], data_mc[[1]][, "garch"]),
      metrics(data_mc[[1]][, "garch_n_sv_n"], data_mc[[1]][, "garch"]),
      metrics(data_mc[[1]][, "garch_n_sv_t"], data_mc[[1]][, "garch"]),
      metrics(data_mc[[1]][, "garch_n_ms_n"], data_mc[[1]][, "garch"]),
      metrics(data_mc[[1]][, "garch_n_ms_t"], data_mc[[1]][, "garch"]),
      metrics(data_mc[[1]][, "gas_n_garch_n"], data_mc[[1]][, "gas"]),
      metrics(data_mc[[1]][, "gas_n_garch_t"], data_mc[[1]][, "gas"]),
      metrics(data_mc[[1]][, "gas_n_gas_n"], data_mc[[1]][, "gas"]),
      metrics(data_mc[[1]][, "gas_n_gas_t"], data_mc[[1]][, "gas"]),
      metrics(data_mc[[1]][, "gas_n_sv_n"], data_mc[[1]][, "gas"]),
      metrics(data_mc[[1]][, "gas_n_sv_t"], data_mc[[1]][, "gas"]),
      metrics(data_mc[[1]][, "gas_n_ms_n"], data_mc[[1]][, "gas"]),
      metrics(data_mc[[1]][, "gas_n_ms_t"], data_mc[[1]][, "gas"]),
      metrics(data_mc[[1]][, "sv_n_garch_n"], data_mc[[1]][, "sv"]),
      metrics(data_mc[[1]][, "sv_n_garch_t"], data_mc[[1]][, "sv"]),
      metrics(data_mc[[1]][, "sv_n_gas_n"], data_mc[[1]][, "sv"]),
      metrics(data_mc[[1]][, "sv_n_gas_t"], data_mc[[1]][, "sv"]),
      metrics(data_mc[[1]][, "sv_n_sv_n"], data_mc[[1]][, "sv"]),
      metrics(data_mc[[1]][, "sv_n_sv_t"], data_mc[[1]][, "sv"]),
      metrics(data_mc[[1]][, "sv_n_ms_n"], data_mc[[1]][, "sv"]),
      metrics(data_mc[[1]][, "sv_n_ms_t"], data_mc[[1]][, "sv"]),
      metrics(data_mc[[1]][, "ms_n_garch_n"], data_mc[[1]][, "ms"]),
      metrics(data_mc[[1]][, "ms_n_garch_t"], data_mc[[1]][, "ms"]),
      metrics(data_mc[[1]][, "ms_n_gas_n"], data_mc[[1]][, "ms"]),
      metrics(data_mc[[1]][, "ms_n_gas_t"], data_mc[[1]][, "ms"]),
      metrics(data_mc[[1]][, "ms_n_sv_n"], data_mc[[1]][, "ms"]),
      metrics(data_mc[[1]][, "ms_n_sv_t"], data_mc[[1]][, "ms"]),
      metrics(data_mc[[1]][, "ms_n_ms_n"], data_mc[[1]][, "ms"]),
      metrics(data_mc[[1]][, "ms_n_ms_t"], data_mc[[1]][, "ms"])
    ),
    4
  ) |>
    data.frame() |>
    mutate(
      n = 1000,
      distri = "N",
      DGP = c(rep("GARCH", 8), rep("GAS", 8), rep("SV", 8), rep("MS", 8)),
      Estim = rep(
        c("GARCH-N", "GARCH-T", "GAS-N", "GAS-T", "SV-N", "SV-T", "MS-N", "MS-T"),
        4
      )
    )
  
  results_2 <- round(
    rbind(
      metrics(data_mc[[2]][, "garch_t_garch_n"], data_mc[[2]][, "garch"]),
      metrics(data_mc[[2]][, "garch_t_garch_t"], data_mc[[2]][, "garch"]),
      metrics(data_mc[[2]][, "garch_t_gas_n"], data_mc[[2]][, "garch"]),
      metrics(data_mc[[2]][, "garch_t_gas_t"], data_mc[[2]][, "garch"]),
      metrics(data_mc[[2]][, "garch_t_sv_n"], data_mc[[2]][, "garch"]),
      metrics(data_mc[[2]][, "garch_t_sv_t"], data_mc[[2]][, "garch"]),
      metrics(data_mc[[2]][, "garch_t_ms_n"], data_mc[[2]][, "garch"]),
      metrics(data_mc[[2]][, "garch_t_ms_t"], data_mc[[2]][, "garch"]),
      metrics(data_mc[[2]][, "gas_t_garch_n"], data_mc[[2]][, "gas"]),
      metrics(data_mc[[2]][, "gas_t_garch_t"], data_mc[[2]][, "gas"]),
      metrics(data_mc[[2]][, "gas_t_gas_n"], data_mc[[2]][, "gas"]),
      metrics(data_mc[[2]][, "gas_t_gas_t"], data_mc[[2]][, "gas"]),
      metrics(data_mc[[2]][, "gas_t_sv_n"], data_mc[[2]][, "gas"]),
      metrics(data_mc[[2]][, "gas_t_sv_t"], data_mc[[2]][, "gas"]),
      metrics(data_mc[[2]][, "gas_t_ms_n"], data_mc[[2]][, "gas"]),
      metrics(data_mc[[2]][, "gas_t_ms_t"], data_mc[[2]][, "gas"]),
      metrics(data_mc[[2]][, "sv_t_garch_n"], data_mc[[2]][, "sv"]),
      metrics(data_mc[[2]][, "sv_t_garch_t"], data_mc[[2]][, "sv"]),
      metrics(data_mc[[2]][, "sv_t_gas_n"], data_mc[[2]][, "sv"]),
      metrics(data_mc[[2]][, "sv_t_gas_t"], data_mc[[2]][, "sv"]),
      metrics(data_mc[[2]][, "sv_t_sv_n"], data_mc[[2]][, "sv"]),
      metrics(data_mc[[2]][, "sv_t_sv_t"], data_mc[[2]][, "sv"]),
      metrics(data_mc[[2]][, "sv_t_ms_n"], data_mc[[2]][, "sv"]),
      metrics(data_mc[[2]][, "sv_t_ms_t"], data_mc[[2]][, "sv"]),
      metrics(data_mc[[2]][, "ms_t_garch_n"], data_mc[[2]][, "ms"]),
      metrics(data_mc[[2]][, "ms_t_garch_t"], data_mc[[2]][, "ms"]),
      metrics(data_mc[[2]][, "ms_t_gas_n"], data_mc[[2]][, "ms"]),
      metrics(data_mc[[2]][, "ms_t_gas_t"], data_mc[[2]][, "ms"]),
      metrics(data_mc[[2]][, "ms_t_sv_n"], data_mc[[2]][, "ms"]),
      metrics(data_mc[[2]][, "ms_t_sv_t"], data_mc[[2]][, "ms"]),
      metrics(data_mc[[2]][, "ms_t_ms_n"], data_mc[[2]][, "ms"]),
      metrics(data_mc[[2]][, "ms_t_ms_t"], data_mc[[2]][, "ms"])
    ),
    4
  ) |>
    data.frame() |>
    mutate(
      n = 1000,
      distri = "T",
      DGP = c(rep("GARCH", 8), rep("GAS", 8), rep("SV", 8), rep("MS", 8)),
      Estim = rep(
        c("GARCH-N", "GARCH-T", "GAS-N", "GAS-T", "SV-N", "SV-T", "MS-N", "MS-T"),
        4
      )
    )
  
  
  results_3 <- round(
    rbind(
      metrics(data_mc[[3]][, "garch_n_garch_n"], data_mc[[3]][, "garch"]),
      metrics(data_mc[[3]][, "garch_n_garch_t"], data_mc[[3]][, "garch"]),
      metrics(data_mc[[3]][, "garch_n_gas_n"], data_mc[[3]][, "garch"]),
      metrics(data_mc[[3]][, "garch_n_gas_t"], data_mc[[3]][, "garch"]),
      metrics(data_mc[[3]][, "garch_n_sv_n"], data_mc[[3]][, "garch"]),
      metrics(data_mc[[3]][, "garch_n_sv_t"], data_mc[[3]][, "garch"]),
      metrics(data_mc[[3]][, "garch_n_ms_n"], data_mc[[3]][, "garch"]),
      metrics(data_mc[[3]][, "garch_n_ms_t"], data_mc[[3]][, "garch"]),
      metrics(data_mc[[3]][, "gas_n_garch_n"], data_mc[[3]][, "gas"]),
      metrics(data_mc[[3]][, "gas_n_garch_t"], data_mc[[3]][, "gas"]),
      metrics(data_mc[[3]][, "gas_n_gas_n"], data_mc[[3]][, "gas"]),
      metrics(data_mc[[3]][, "gas_n_gas_t"], data_mc[[3]][, "gas"]),
      metrics(data_mc[[3]][, "gas_n_sv_n"], data_mc[[3]][, "gas"]),
      metrics(data_mc[[3]][, "gas_n_sv_t"], data_mc[[3]][, "gas"]),
      metrics(data_mc[[3]][, "gas_n_ms_n"], data_mc[[3]][, "gas"]),
      metrics(data_mc[[3]][, "gas_n_ms_t"], data_mc[[3]][, "gas"]),
      metrics(data_mc[[3]][, "sv_n_garch_n"], data_mc[[3]][, "sv"]),
      metrics(data_mc[[3]][, "sv_n_garch_t"], data_mc[[3]][, "sv"]),
      metrics(data_mc[[3]][, "sv_n_gas_n"], data_mc[[3]][, "sv"]),
      metrics(data_mc[[3]][, "sv_n_gas_t"], data_mc[[3]][, "sv"]),
      metrics(data_mc[[3]][, "sv_n_sv_n"], data_mc[[3]][, "sv"]),
      metrics(data_mc[[3]][, "sv_n_sv_t"], data_mc[[3]][, "sv"]),
      metrics(data_mc[[3]][, "sv_n_ms_n"], data_mc[[3]][, "sv"]),
      metrics(data_mc[[3]][, "sv_n_ms_t"], data_mc[[3]][, "sv"]),
      metrics(data_mc[[3]][, "ms_n_garch_n"], data_mc[[3]][, "ms"]),
      metrics(data_mc[[3]][, "ms_n_garch_t"], data_mc[[3]][, "ms"]),
      metrics(data_mc[[3]][, "ms_n_gas_n"], data_mc[[3]][, "ms"]),
      metrics(data_mc[[3]][, "ms_n_gas_t"], data_mc[[3]][, "ms"]),
      metrics(data_mc[[3]][, "ms_n_sv_n"], data_mc[[3]][, "ms"]),
      metrics(data_mc[[3]][, "ms_n_sv_t"], data_mc[[3]][, "ms"]),
      metrics(data_mc[[3]][, "ms_n_ms_n"], data_mc[[3]][, "ms"]),
      metrics(data_mc[[3]][, "ms_n_ms_t"], data_mc[[3]][, "ms"])
    ),
    4
  ) |>
    data.frame() |>
    mutate(
      n = 2500,
      distri = "N",
      DGP = c(rep("GARCH", 8), rep("GAS", 8), rep("SV", 8), rep("MS", 8)),
      Estim = rep(
        c("GARCH-N", "GARCH-T", "GAS-N", "GAS-T", "SV-N", "SV-T", "MS-N", "MS-T"),
        4
      )
    )
  
  
  results_4 <- round(
    rbind(
      metrics(data_mc[[4]][, "garch_t_garch_n"], data_mc[[4]][, "garch"]),
      metrics(data_mc[[4]][, "garch_t_garch_t"], data_mc[[4]][, "garch"]),
      metrics(data_mc[[4]][, "garch_t_gas_n"], data_mc[[4]][, "garch"]),
      metrics(data_mc[[4]][, "garch_t_gas_t"], data_mc[[4]][, "garch"]),
      metrics(data_mc[[4]][, "garch_t_sv_n"], data_mc[[4]][, "garch"]),
      metrics(data_mc[[4]][, "garch_t_sv_t"], data_mc[[4]][, "garch"]),
      metrics(data_mc[[4]][, "garch_t_ms_n"], data_mc[[4]][, "garch"]),
      metrics(data_mc[[4]][, "garch_t_ms_t"], data_mc[[4]][, "garch"]),
      metrics(data_mc[[4]][, "gas_t_garch_n"], data_mc[[4]][, "gas"]),
      metrics(data_mc[[4]][, "gas_t_garch_t"], data_mc[[4]][, "gas"]),
      metrics(data_mc[[4]][, "gas_t_gas_n"], data_mc[[4]][, "gas"]),
      metrics(data_mc[[4]][, "gas_t_gas_t"], data_mc[[4]][, "gas"]),
      metrics(data_mc[[4]][, "gas_t_sv_n"], data_mc[[4]][, "gas"]),
      metrics(data_mc[[4]][, "gas_t_sv_t"], data_mc[[4]][, "gas"]),
      metrics(data_mc[[4]][, "gas_t_ms_n"], data_mc[[4]][, "gas"]),
      metrics(data_mc[[4]][, "gas_t_ms_t"], data_mc[[4]][, "gas"]),
      metrics(data_mc[[4]][, "sv_t_garch_n"], data_mc[[4]][, "sv"]),
      metrics(data_mc[[4]][, "sv_t_garch_t"], data_mc[[4]][, "sv"]),
      metrics(data_mc[[4]][, "sv_t_gas_n"], data_mc[[4]][, "sv"]),
      metrics(data_mc[[4]][, "sv_t_gas_t"], data_mc[[4]][, "sv"]),
      metrics(data_mc[[4]][, "sv_t_sv_n"], data_mc[[4]][, "sv"]),
      metrics(data_mc[[4]][, "sv_t_sv_t"], data_mc[[4]][, "sv"]),
      metrics(data_mc[[4]][, "sv_t_ms_n"], data_mc[[4]][, "sv"]),
      metrics(data_mc[[4]][, "sv_t_ms_t"], data_mc[[4]][, "sv"]),
      metrics(data_mc[[4]][, "ms_t_garch_n"], data_mc[[4]][, "ms"]),
      metrics(data_mc[[4]][, "ms_t_garch_t"], data_mc[[4]][, "ms"]),
      metrics(data_mc[[4]][, "ms_t_gas_n"], data_mc[[4]][, "ms"]),
      metrics(data_mc[[4]][, "ms_t_gas_t"], data_mc[[4]][, "ms"]),
      metrics(data_mc[[4]][, "ms_t_sv_n"], data_mc[[4]][, "ms"]),
      metrics(data_mc[[4]][, "ms_t_sv_t"], data_mc[[4]][, "ms"]),
      metrics(data_mc[[4]][, "ms_t_ms_n"], data_mc[[4]][, "ms"]),
      metrics(data_mc[[4]][, "ms_t_ms_t"], data_mc[[4]][, "ms"])
    ),
    4
  ) |>
    data.frame() |>
    mutate(
      n = 2500,
      distri = "T",
      DGP = c(rep("GARCH", 8), rep("GAS", 8), rep("SV", 8), rep("MS", 8)),
      Estim = rep(
        c("GARCH-N", "GARCH-T", "GAS-N", "GAS-T", "SV-N", "SV-T", "MS-N", "MS-T"),
        4
      )
    )
  
  
  results_5 <- round(
    rbind(
      metrics(data_mc[[5]][, "garch_n_garch_n"], data_mc[[5]][, "garch"]),
      metrics(data_mc[[5]][, "garch_n_garch_t"], data_mc[[5]][, "garch"]),
      metrics(data_mc[[5]][, "garch_n_gas_n"], data_mc[[5]][, "garch"]),
      metrics(data_mc[[5]][, "garch_n_gas_t"], data_mc[[5]][, "garch"]),
      metrics(data_mc[[5]][, "garch_n_sv_n"], data_mc[[5]][, "garch"]),
      metrics(data_mc[[5]][, "garch_n_sv_t"], data_mc[[5]][, "garch"]),
      metrics(data_mc[[5]][, "garch_n_ms_n"], data_mc[[5]][, "garch"]),
      metrics(data_mc[[5]][, "garch_n_ms_t"], data_mc[[5]][, "garch"]),
      metrics(data_mc[[5]][, "gas_n_garch_n"], data_mc[[5]][, "gas"]),
      metrics(data_mc[[5]][, "gas_n_garch_t"], data_mc[[5]][, "gas"]),
      metrics(data_mc[[5]][, "gas_n_gas_n"], data_mc[[5]][, "gas"]),
      metrics(data_mc[[5]][, "gas_n_gas_t"], data_mc[[5]][, "gas"]),
      metrics(data_mc[[5]][, "gas_n_sv_n"], data_mc[[5]][, "gas"]),
      metrics(data_mc[[5]][, "gas_n_sv_t"], data_mc[[5]][, "gas"]),
      metrics(data_mc[[5]][, "gas_n_ms_n"], data_mc[[5]][, "gas"]),
      metrics(data_mc[[5]][, "gas_n_ms_t"], data_mc[[5]][, "gas"]),
      metrics(data_mc[[5]][, "sv_n_garch_n"], data_mc[[5]][, "sv"]),
      metrics(data_mc[[5]][, "sv_n_garch_t"], data_mc[[5]][, "sv"]),
      metrics(data_mc[[5]][, "sv_n_gas_n"], data_mc[[5]][, "sv"]),
      metrics(data_mc[[5]][, "sv_n_gas_t"], data_mc[[5]][, "sv"]),
      metrics(data_mc[[5]][, "sv_n_sv_n"], data_mc[[5]][, "sv"]),
      metrics(data_mc[[5]][, "sv_n_sv_t"], data_mc[[5]][, "sv"]),
      metrics(data_mc[[5]][, "sv_n_ms_n"], data_mc[[5]][, "sv"]),
      metrics(data_mc[[5]][, "sv_n_ms_t"], data_mc[[5]][, "sv"]),
      metrics(data_mc[[5]][, "ms_n_garch_n"], data_mc[[5]][, "ms"]),
      metrics(data_mc[[5]][, "ms_n_garch_t"], data_mc[[5]][, "ms"]),
      metrics(data_mc[[5]][, "ms_n_gas_n"], data_mc[[5]][, "ms"]),
      metrics(data_mc[[5]][, "ms_n_gas_t"], data_mc[[5]][, "ms"]),
      metrics(data_mc[[5]][, "ms_n_sv_n"], data_mc[[5]][, "ms"]),
      metrics(data_mc[[5]][, "ms_n_sv_t"], data_mc[[5]][, "ms"]),
      metrics(data_mc[[5]][, "ms_n_ms_n"], data_mc[[5]][, "ms"]),
      metrics(data_mc[[5]][, "ms_n_ms_t"], data_mc[[5]][, "ms"])
    ),
    4
  ) |>
    data.frame() |>
    mutate(
      n = 500,
      distri = "N",
      DGP = c(rep("GARCH", 8), rep("GAS", 8), rep("SV", 8), rep("MS", 8)),
      Estim = rep(
        c("GARCH-N", "GARCH-T", "GAS-N", "GAS-T", "SV-N", "SV-T", "MS-N", "MS-T"),
        4
      )
    )
  
  
  results_6 <- round(
    rbind(
      metrics(data_mc[[6]][, "garch_t_garch_n"], data_mc[[6]][, "garch"]),
      metrics(data_mc[[6]][, "garch_t_garch_t"], data_mc[[6]][, "garch"]),
      metrics(data_mc[[6]][, "garch_t_gas_n"], data_mc[[6]][, "garch"]),
      metrics(data_mc[[6]][, "garch_t_gas_t"], data_mc[[6]][, "garch"]),
      metrics(data_mc[[6]][, "garch_t_sv_n"], data_mc[[6]][, "garch"]),
      metrics(data_mc[[6]][, "garch_t_sv_t"], data_mc[[6]][, "garch"]),
      metrics(data_mc[[6]][, "garch_t_ms_n"], data_mc[[6]][, "garch"]),
      metrics(data_mc[[6]][, "garch_t_ms_t"], data_mc[[6]][, "garch"]),
      metrics(data_mc[[6]][, "gas_t_garch_n"], data_mc[[6]][, "gas"]),
      metrics(data_mc[[6]][, "gas_t_garch_t"], data_mc[[6]][, "gas"]),
      metrics(data_mc[[6]][, "gas_t_gas_n"], data_mc[[6]][, "gas"]),
      metrics(data_mc[[6]][, "gas_t_gas_t"], data_mc[[6]][, "gas"]),
      metrics(data_mc[[6]][, "gas_t_sv_n"], data_mc[[6]][, "gas"]),
      metrics(data_mc[[6]][, "gas_t_sv_t"], data_mc[[6]][, "gas"]),
      metrics(data_mc[[6]][, "gas_t_ms_n"], data_mc[[6]][, "gas"]),
      metrics(data_mc[[6]][, "gas_t_ms_t"], data_mc[[6]][, "gas"]),
      metrics(data_mc[[6]][, "sv_t_garch_n"], data_mc[[6]][, "sv"]),
      metrics(data_mc[[6]][, "sv_t_garch_t"], data_mc[[6]][, "sv"]),
      metrics(data_mc[[6]][, "sv_t_gas_n"], data_mc[[6]][, "sv"]),
      metrics(data_mc[[6]][, "sv_t_gas_t"], data_mc[[6]][, "sv"]),
      metrics(data_mc[[6]][, "sv_t_sv_n"], data_mc[[6]][, "sv"]),
      metrics(data_mc[[6]][, "sv_t_sv_t"], data_mc[[6]][, "sv"]),
      metrics(data_mc[[6]][, "sv_t_ms_n"], data_mc[[6]][, "sv"]),
      metrics(data_mc[[6]][, "sv_t_ms_t"], data_mc[[6]][, "sv"]),
      metrics(data_mc[[6]][, "ms_t_garch_n"], data_mc[[6]][, "ms"]),
      metrics(data_mc[[6]][, "ms_t_garch_t"], data_mc[[6]][, "ms"]),
      metrics(data_mc[[6]][, "ms_t_gas_n"], data_mc[[6]][, "ms"]),
      metrics(data_mc[[6]][, "ms_t_gas_t"], data_mc[[6]][, "ms"]),
      metrics(data_mc[[6]][, "ms_t_sv_n"], data_mc[[6]][, "ms"]),
      metrics(data_mc[[6]][, "ms_t_sv_t"], data_mc[[6]][, "ms"]),
      metrics(data_mc[[6]][, "ms_t_ms_n"], data_mc[[6]][, "ms"]),
      metrics(data_mc[[6]][, "ms_t_ms_t"], data_mc[[6]][, "ms"])
    ),
    4
  ) |>
    data.frame() |>
    mutate(
      n = 500,
      distri = "T",
      DGP = c(rep("GARCH", 8), rep("GAS", 8), rep("SV", 8), rep("MS", 8)),
      Estim = rep(
        c("GARCH-N", "GARCH-T", "GAS-N", "GAS-T", "SV-N", "SV-T", "MS-N", "MS-T"),
        4
      )
    )
  
  results <- rbind(results_5, results_6, results_1, results_2, results_3, results_4) |> select(DGP, n, distri, Estim, everything())
  colnames(results) <- c("DGP", "N", "DIST", "Estim", "MSE", "QLIKE", "MSE LOG", "MSE SD", "MSE PROP", "MAE", "MAE LOG", "MAE SD", "MAE PROP")
  return(results)
}


making_plots <- function(files_names_F, files_names_T) {
  
  # Outliers FALSE
  data_mc <- list()
  for (name in files_names_F) data_mc[[name]] <- read.csv(name)
  
  simul <- c(rep(data_mc[[1]][, "garch"], 8), rep(data_mc[[1]][, "gas"], 8), rep(data_mc[[1]][, "sv"], 8), rep(data_mc[[1]][, "ms"], 8))
  dgp <- c(rep("GARCH", 8000), rep("GAS", 8000), rep("SV", 8000), rep("MS", 8000))
  fore <- c(data_mc[[1]][, "garch_n_garch_n"], data_mc[[1]][, "garch_n_garch_t"], data_mc[[1]][, "garch_n_gas_n"], data_mc[[1]][, "garch_n_gas_t"],
    data_mc[[1]][, "garch_n_sv_n"], data_mc[[1]][, "garch_n_sv_t"], data_mc[[1]][, "garch_n_ms_n"],  data_mc[[1]][, "garch_n_ms_t"],
    data_mc[[1]][, "gas_n_garch_n"], data_mc[[1]][, "gas_n_garch_t"], data_mc[[1]][, "gas_n_gas_n"], data_mc[[1]][, "gas_n_gas_t"],
    data_mc[[1]][, "gas_n_sv_n"], data_mc[[1]][, "gas_n_sv_t"], data_mc[[1]][, "gas_n_ms_n"], data_mc[[1]][, "gas_n_ms_t"],
    data_mc[[1]][, "sv_n_garch_n"], data_mc[[1]][, "sv_n_garch_t"], data_mc[[1]][, "sv_n_gas_n"], data_mc[[1]][, "sv_n_gas_t"],
    data_mc[[1]][, "sv_n_sv_n"], data_mc[[1]][, "sv_n_sv_t"], data_mc[[1]][, "sv_n_ms_n"], data_mc[[1]][, "sv_n_ms_t"],
    data_mc[[1]][, "ms_n_garch_n"], data_mc[[1]][, "ms_n_garch_t"], data_mc[[1]][, "ms_n_gas_n"], data_mc[[1]][, "ms_n_gas_t"],
    data_mc[[1]][, "ms_n_sv_n"], data_mc[[1]][, "ms_n_sv_t"], data_mc[[1]][, "ms_n_ms_n"], data_mc[[1]][, "ms_n_ms_t"])
  model <- rep(c(rep("GARCH-N", 1000), rep("GARCH-T", 1000), rep("GAS-N", 1000), rep("GAS-T", 1000), rep("SV-N", 1000), rep("SV-T", 1000), rep("MS-N", 1000), rep("MS-T", 1000)), 4)
  n <- rep(1000, 32000)
  distri <- rep("N", 32000)
  mc_1 <- data.frame(simul, dgp, fore, model, n, distri)
  
  simul <- c(rep(data_mc[[2]][, "garch"], 8), rep(data_mc[[2]][, "gas"], 8), rep(data_mc[[2]][, "sv"], 8), rep(data_mc[[2]][, "ms"], 8))
  dgp <- c(rep("GARCH", 8000), rep("GAS", 8000), rep("SV", 8000), rep("MS", 8000))
  fore <- c(data_mc[[2]][, "garch_t_garch_n"], data_mc[[2]][, "garch_t_garch_t"], data_mc[[2]][, "garch_t_gas_n"], data_mc[[2]][, "garch_t_gas_t"],
    data_mc[[2]][, "garch_t_sv_n"], data_mc[[2]][, "garch_t_sv_t"], data_mc[[2]][, "garch_t_ms_n"],  data_mc[[2]][, "garch_t_ms_t"],
    data_mc[[2]][, "gas_t_garch_n"], data_mc[[2]][, "gas_t_garch_t"], data_mc[[2]][, "gas_t_gas_n"], data_mc[[2]][, "gas_t_gas_t"],
    data_mc[[2]][, "gas_t_sv_n"], data_mc[[2]][, "gas_t_sv_t"], data_mc[[2]][, "gas_t_ms_n"], data_mc[[2]][, "gas_t_ms_t"],
    data_mc[[2]][, "sv_t_garch_n"], data_mc[[2]][, "sv_t_garch_t"], data_mc[[2]][, "sv_t_gas_n"], data_mc[[2]][, "sv_t_gas_t"],
    data_mc[[2]][, "sv_t_sv_n"], data_mc[[2]][, "sv_t_sv_t"], data_mc[[2]][, "sv_t_ms_n"], data_mc[[2]][, "sv_t_ms_t"],
    data_mc[[2]][, "ms_t_garch_n"], data_mc[[2]][, "ms_t_garch_t"], data_mc[[2]][, "ms_t_gas_n"], data_mc[[2]][, "ms_t_gas_t"],
    data_mc[[2]][, "ms_t_sv_n"], data_mc[[2]][, "ms_t_sv_t"], data_mc[[2]][, "ms_t_ms_n"], data_mc[[2]][, "ms_t_ms_t"])
  model <- rep(c(rep("GARCH-N", 1000), rep("GARCH-T", 1000), rep("GAS-N", 1000), rep("GAS-T", 1000), rep("SV-N", 1000), rep("SV-T", 1000), rep("MS-N", 1000), rep("MS-T", 1000)), 4)
  n <- rep(1000, 32000)
  distri <- rep("T", 32000)
  mc_2 <- data.frame(simul, dgp, fore, model, n, distri)
  
  simul <- c(rep(data_mc[[3]][, "garch"], 8), rep(data_mc[[3]][, "gas"], 8), rep(data_mc[[3]][, "sv"], 8), rep(data_mc[[3]][, "ms"], 8))
  dgp <- c(rep("GARCH", 8000), rep("GAS", 8000), rep("SV", 8000), rep("MS", 8000))
  fore <- c(data_mc[[3]][, "garch_n_garch_n"], data_mc[[3]][, "garch_n_garch_t"], data_mc[[3]][, "garch_n_gas_n"], data_mc[[3]][, "garch_n_gas_t"],
    data_mc[[3]][, "garch_n_sv_n"], data_mc[[3]][, "garch_n_sv_t"], data_mc[[3]][, "garch_n_ms_n"],  data_mc[[3]][, "garch_n_ms_t"],
    data_mc[[3]][, "gas_n_garch_n"], data_mc[[3]][, "gas_n_garch_t"], data_mc[[3]][, "gas_n_gas_n"], data_mc[[3]][, "gas_n_gas_t"],
    data_mc[[3]][, "gas_n_sv_n"], data_mc[[3]][, "gas_n_sv_t"], data_mc[[3]][, "gas_n_ms_n"], data_mc[[3]][, "gas_n_ms_t"],
    data_mc[[3]][, "sv_n_garch_n"], data_mc[[3]][, "sv_n_garch_t"], data_mc[[3]][, "sv_n_gas_n"], data_mc[[3]][, "sv_n_gas_t"],
    data_mc[[3]][, "sv_n_sv_n"], data_mc[[3]][, "sv_n_sv_t"], data_mc[[3]][, "sv_n_ms_n"], data_mc[[3]][, "sv_n_ms_t"],
    data_mc[[3]][, "ms_n_garch_n"], data_mc[[3]][, "ms_n_garch_t"], data_mc[[3]][, "ms_n_gas_n"], data_mc[[3]][, "ms_n_gas_t"],
    data_mc[[3]][, "ms_n_sv_n"], data_mc[[3]][, "ms_n_sv_t"], data_mc[[3]][, "ms_n_ms_n"], data_mc[[3]][, "ms_n_ms_t"])
  model <- rep(c(rep("GARCH-N", 1000), rep("GARCH-T", 1000), rep("GAS-N", 1000), rep("GAS-T", 1000), rep("SV-N", 1000), rep("SV-T", 1000), rep("MS-N", 1000), rep("MS-T", 1000)), 4)
  n <- rep(2500, 32000)
  distri <- rep("N", 32000)
  mc_3 <- data.frame(simul, dgp, fore, model, n, distri)
  
  simul <- c(rep(data_mc[[4]][, "garch"], 8), rep(data_mc[[4]][, "gas"], 8), rep(data_mc[[4]][, "sv"], 8), rep(data_mc[[4]][, "ms"], 8))
  dgp <- c(rep("GARCH", 8000), rep("GAS", 8000), rep("SV", 8000), rep("MS", 8000))
  fore <- c(data_mc[[4]][, "garch_t_garch_n"], data_mc[[4]][, "garch_t_garch_t"], data_mc[[4]][, "garch_t_gas_n"], data_mc[[4]][, "garch_t_gas_t"],
    data_mc[[4]][, "garch_t_sv_n"], data_mc[[4]][, "garch_t_sv_t"], data_mc[[4]][, "garch_t_ms_n"],  data_mc[[4]][, "garch_t_ms_t"],
    data_mc[[4]][, "gas_t_garch_n"], data_mc[[4]][, "gas_t_garch_t"], data_mc[[4]][, "gas_t_gas_n"], data_mc[[4]][, "gas_t_gas_t"],
    data_mc[[4]][, "gas_t_sv_n"], data_mc[[4]][, "gas_t_sv_t"], data_mc[[4]][, "gas_t_ms_n"], data_mc[[4]][, "gas_t_ms_t"],
    data_mc[[4]][, "sv_t_garch_n"], data_mc[[4]][, "sv_t_garch_t"], data_mc[[4]][, "sv_t_gas_n"], data_mc[[4]][, "sv_t_gas_t"],
    data_mc[[4]][, "sv_t_sv_n"], data_mc[[4]][, "sv_t_sv_t"], data_mc[[4]][, "sv_t_ms_n"], data_mc[[4]][, "sv_t_ms_t"],
    data_mc[[4]][, "ms_t_garch_n"], data_mc[[4]][, "ms_t_garch_t"], data_mc[[4]][, "ms_t_gas_n"], data_mc[[4]][, "ms_t_gas_t"],
    data_mc[[4]][, "ms_t_sv_n"], data_mc[[4]][, "ms_t_sv_t"], data_mc[[4]][, "ms_t_ms_n"], data_mc[[4]][, "ms_t_ms_t"])
  model <- rep(c(rep("GARCH-N", 1000), rep("GARCH-T", 1000), rep("GAS-N", 1000), rep("GAS-T", 1000), rep("SV-N", 1000), rep("SV-T", 1000), rep("MS-N", 1000), rep("MS-T", 1000)), 4)
  n <- rep(2500, 32000)
  distri <- rep("T", 32000)
  mc_4 <- data.frame(simul, dgp, fore, model, n, distri)
  
  simul <- c(rep(data_mc[[5]][, "garch"], 8), rep(data_mc[[5]][, "gas"], 8), rep(data_mc[[5]][, "sv"], 8), rep(data_mc[[5]][, "ms"], 8))
  dgp <- c(rep("GARCH", 8000), rep("GAS", 8000), rep("SV", 8000), rep("MS", 8000))
  fore <- c(data_mc[[5]][, "garch_n_garch_n"], data_mc[[5]][, "garch_n_garch_t"], data_mc[[5]][, "garch_n_gas_n"], data_mc[[5]][, "garch_n_gas_t"],
    data_mc[[5]][, "garch_n_sv_n"], data_mc[[5]][, "garch_n_sv_t"], data_mc[[5]][, "garch_n_ms_n"],  data_mc[[5]][, "garch_n_ms_t"],
    data_mc[[5]][, "gas_n_garch_n"], data_mc[[5]][, "gas_n_garch_t"], data_mc[[5]][, "gas_n_gas_n"], data_mc[[5]][, "gas_n_gas_t"],
    data_mc[[5]][, "gas_n_sv_n"], data_mc[[5]][, "gas_n_sv_t"], data_mc[[5]][, "gas_n_ms_n"], data_mc[[5]][, "gas_n_ms_t"],
    data_mc[[5]][, "sv_n_garch_n"], data_mc[[5]][, "sv_n_garch_t"], data_mc[[5]][, "sv_n_gas_n"], data_mc[[5]][, "sv_n_gas_t"],
    data_mc[[5]][, "sv_n_sv_n"], data_mc[[5]][, "sv_n_sv_t"], data_mc[[5]][, "sv_n_ms_n"], data_mc[[5]][, "sv_n_ms_t"],
    data_mc[[5]][, "ms_n_garch_n"], data_mc[[5]][, "ms_n_garch_t"], data_mc[[5]][, "ms_n_gas_n"], data_mc[[5]][, "ms_n_gas_t"],
    data_mc[[5]][, "ms_n_sv_n"], data_mc[[5]][, "ms_n_sv_t"], data_mc[[5]][, "ms_n_ms_n"], data_mc[[5]][, "ms_n_ms_t"])
  model <- rep(c(rep("GARCH-N", 1000), rep("GARCH-T", 1000), rep("GAS-N", 1000), rep("GAS-T", 1000), rep("SV-N", 1000), rep("SV-T", 1000), rep("MS-N", 1000), rep("MS-T", 1000)), 4)
  n <- rep(500, 32000)
  distri <- rep("N", 32000)
  mc_5 <- data.frame(simul, dgp, fore, model, n, distri)
  
  simul <- c(rep(data_mc[[6]][, "garch"], 8), rep(data_mc[[6]][, "gas"], 8), rep(data_mc[[6]][, "sv"], 8), rep(data_mc[[6]][, "ms"], 8))
  dgp <- c(rep("GARCH", 8000), rep("GAS", 8000), rep("SV", 8000), rep("MS", 8000))
  fore <- c(data_mc[[6]][, "garch_t_garch_n"], data_mc[[6]][, "garch_t_garch_t"], data_mc[[6]][, "garch_t_gas_n"], data_mc[[6]][, "garch_t_gas_t"],
    data_mc[[6]][, "garch_t_sv_n"], data_mc[[6]][, "garch_t_sv_t"], data_mc[[6]][, "garch_t_ms_n"],  data_mc[[6]][, "garch_t_ms_t"],
    data_mc[[6]][, "gas_t_garch_n"], data_mc[[6]][, "gas_t_garch_t"], data_mc[[6]][, "gas_t_gas_n"], data_mc[[6]][, "gas_t_gas_t"],
    data_mc[[6]][, "gas_t_sv_n"], data_mc[[6]][, "gas_t_sv_t"], data_mc[[6]][, "gas_t_ms_n"], data_mc[[6]][, "gas_t_ms_t"],
    data_mc[[6]][, "sv_t_garch_n"], data_mc[[6]][, "sv_t_garch_t"], data_mc[[6]][, "sv_t_gas_n"], data_mc[[6]][, "sv_t_gas_t"],
    data_mc[[6]][, "sv_t_sv_n"], data_mc[[6]][, "sv_t_sv_t"], data_mc[[6]][, "sv_t_ms_n"], data_mc[[6]][, "sv_t_ms_t"],
    data_mc[[6]][, "ms_t_garch_n"], data_mc[[6]][, "ms_t_garch_t"], data_mc[[6]][, "ms_t_gas_n"], data_mc[[6]][, "ms_t_gas_t"],
    data_mc[[6]][, "ms_t_sv_n"], data_mc[[6]][, "ms_t_sv_t"], data_mc[[6]][, "ms_t_ms_n"], data_mc[[6]][, "ms_t_ms_t"])
  model <- rep(c(rep("GARCH-N", 1000), rep("GARCH-T", 1000), rep("GAS-N", 1000), rep("GAS-T", 1000), rep("SV-N", 1000), rep("SV-T", 1000), rep("MS-N", 1000), rep("MS-T", 1000)), 4)
  n <- rep(500, 32000)
  distri <- rep("T", 32000)
  mc_6 <- data.frame(simul, dgp, fore, model, n, distri)
  
  mc_FALSE <- rbind(mc_5, mc_6, mc_1, mc_2, mc_3, mc_4) |>  mutate(ratio = simul/fore, outliers = FALSE)
  
  
  # Outliers TRUE
  data_mc <- list()
  for (name in files_names_T) data_mc[[name]] <- read.csv(name)
  
  simul <- c(rep(data_mc[[1]][, "garch"], 8), rep(data_mc[[1]][, "gas"], 8), rep(data_mc[[1]][, "sv"], 8), rep(data_mc[[1]][, "ms"], 8))
  dgp <- c(rep("GARCH", 8000), rep("GAS", 8000), rep("SV", 8000), rep("MS", 8000))
  fore <- c(data_mc[[1]][, "garch_n_garch_n"], data_mc[[1]][, "garch_n_garch_t"], data_mc[[1]][, "garch_n_gas_n"], data_mc[[1]][, "garch_n_gas_t"],
    data_mc[[1]][, "garch_n_sv_n"], data_mc[[1]][, "garch_n_sv_t"], data_mc[[1]][, "garch_n_ms_n"],  data_mc[[1]][, "garch_n_ms_t"],
    data_mc[[1]][, "gas_n_garch_n"], data_mc[[1]][, "gas_n_garch_t"], data_mc[[1]][, "gas_n_gas_n"], data_mc[[1]][, "gas_n_gas_t"],
    data_mc[[1]][, "gas_n_sv_n"], data_mc[[1]][, "gas_n_sv_t"], data_mc[[1]][, "gas_n_ms_n"], data_mc[[1]][, "gas_n_ms_t"],
    data_mc[[1]][, "sv_n_garch_n"], data_mc[[1]][, "sv_n_garch_t"], data_mc[[1]][, "sv_n_gas_n"], data_mc[[1]][, "sv_n_gas_t"],
    data_mc[[1]][, "sv_n_sv_n"], data_mc[[1]][, "sv_n_sv_t"], data_mc[[1]][, "sv_n_ms_n"], data_mc[[1]][, "sv_n_ms_t"],
    data_mc[[1]][, "ms_n_garch_n"], data_mc[[1]][, "ms_n_garch_t"], data_mc[[1]][, "ms_n_gas_n"], data_mc[[1]][, "ms_n_gas_t"],
    data_mc[[1]][, "ms_n_sv_n"], data_mc[[1]][, "ms_n_sv_t"], data_mc[[1]][, "ms_n_ms_n"], data_mc[[1]][, "ms_n_ms_t"])
  model <- rep(c(rep("GARCH-N", 1000), rep("GARCH-T", 1000), rep("GAS-N", 1000), rep("GAS-T", 1000), rep("SV-N", 1000), rep("SV-T", 1000), rep("MS-N", 1000), rep("MS-T", 1000)), 4)
  n <- rep(1000, 32000)
  distri <- rep("N", 32000)
  mc_1 <- data.frame(simul, dgp, fore, model, n, distri)
  
  simul <- c(rep(data_mc[[2]][, "garch"], 8), rep(data_mc[[2]][, "gas"], 8), rep(data_mc[[2]][, "sv"], 8), rep(data_mc[[2]][, "ms"], 8))
  dgp <- c(rep("GARCH", 8000), rep("GAS", 8000), rep("SV", 8000), rep("MS", 8000))
  fore <- c(data_mc[[2]][, "garch_t_garch_n"], data_mc[[2]][, "garch_t_garch_t"], data_mc[[2]][, "garch_t_gas_n"], data_mc[[2]][, "garch_t_gas_t"],
    data_mc[[2]][, "garch_t_sv_n"], data_mc[[2]][, "garch_t_sv_t"], data_mc[[2]][, "garch_t_ms_n"],  data_mc[[2]][, "garch_t_ms_t"],
    data_mc[[2]][, "gas_t_garch_n"], data_mc[[2]][, "gas_t_garch_t"], data_mc[[2]][, "gas_t_gas_n"], data_mc[[2]][, "gas_t_gas_t"],
    data_mc[[2]][, "gas_t_sv_n"], data_mc[[2]][, "gas_t_sv_t"], data_mc[[2]][, "gas_t_ms_n"], data_mc[[2]][, "gas_t_ms_t"],
    data_mc[[2]][, "sv_t_garch_n"], data_mc[[2]][, "sv_t_garch_t"], data_mc[[2]][, "sv_t_gas_n"], data_mc[[2]][, "sv_t_gas_t"],
    data_mc[[2]][, "sv_t_sv_n"], data_mc[[2]][, "sv_t_sv_t"], data_mc[[2]][, "sv_t_ms_n"], data_mc[[2]][, "sv_t_ms_t"],
    data_mc[[2]][, "ms_t_garch_n"], data_mc[[2]][, "ms_t_garch_t"], data_mc[[2]][, "ms_t_gas_n"], data_mc[[2]][, "ms_t_gas_t"],
    data_mc[[2]][, "ms_t_sv_n"], data_mc[[2]][, "ms_t_sv_t"], data_mc[[2]][, "ms_t_ms_n"], data_mc[[2]][, "ms_t_ms_t"])
  model <- rep(c(rep("GARCH-N", 1000), rep("GARCH-T", 1000), rep("GAS-N", 1000), rep("GAS-T", 1000), rep("SV-N", 1000), rep("SV-T", 1000), rep("MS-N", 1000), rep("MS-T", 1000)), 4)
  n <- rep(1000, 32000)
  distri <- rep("T", 32000)
  mc_2 <- data.frame(simul, dgp, fore, model, n, distri)
  
  simul <- c(rep(data_mc[[3]][, "garch"], 8), rep(data_mc[[3]][, "gas"], 8), rep(data_mc[[3]][, "sv"], 8), rep(data_mc[[3]][, "ms"], 8))
  dgp <- c(rep("GARCH", 8000), rep("GAS", 8000), rep("SV", 8000), rep("MS", 8000))
  fore <- c(data_mc[[3]][, "garch_n_garch_n"], data_mc[[3]][, "garch_n_garch_t"], data_mc[[3]][, "garch_n_gas_n"], data_mc[[3]][, "garch_n_gas_t"],
    data_mc[[3]][, "garch_n_sv_n"], data_mc[[3]][, "garch_n_sv_t"], data_mc[[3]][, "garch_n_ms_n"],  data_mc[[3]][, "garch_n_ms_t"],
    data_mc[[3]][, "gas_n_garch_n"], data_mc[[3]][, "gas_n_garch_t"], data_mc[[3]][, "gas_n_gas_n"], data_mc[[3]][, "gas_n_gas_t"],
    data_mc[[3]][, "gas_n_sv_n"], data_mc[[3]][, "gas_n_sv_t"], data_mc[[3]][, "gas_n_ms_n"], data_mc[[3]][, "gas_n_ms_t"],
    data_mc[[3]][, "sv_n_garch_n"], data_mc[[3]][, "sv_n_garch_t"], data_mc[[3]][, "sv_n_gas_n"], data_mc[[3]][, "sv_n_gas_t"],
    data_mc[[3]][, "sv_n_sv_n"], data_mc[[3]][, "sv_n_sv_t"], data_mc[[3]][, "sv_n_ms_n"], data_mc[[3]][, "sv_n_ms_t"],
    data_mc[[3]][, "ms_n_garch_n"], data_mc[[3]][, "ms_n_garch_t"], data_mc[[3]][, "ms_n_gas_n"], data_mc[[3]][, "ms_n_gas_t"],
    data_mc[[3]][, "ms_n_sv_n"], data_mc[[3]][, "ms_n_sv_t"], data_mc[[3]][, "ms_n_ms_n"], data_mc[[3]][, "ms_n_ms_t"])
  model <- rep(c(rep("GARCH-N", 1000), rep("GARCH-T", 1000), rep("GAS-N", 1000), rep("GAS-T", 1000), rep("SV-N", 1000), rep("SV-T", 1000), rep("MS-N", 1000), rep("MS-T", 1000)), 4)
  n <- rep(2500, 32000)
  distri <- rep("N", 32000)
  mc_3 <- data.frame(simul, dgp, fore, model, n, distri)
  
  simul <- c(rep(data_mc[[4]][, "garch"], 8), rep(data_mc[[4]][, "gas"], 8), rep(data_mc[[4]][, "sv"], 8), rep(data_mc[[4]][, "ms"], 8))
  dgp <- c(rep("GARCH", 8000), rep("GAS", 8000), rep("SV", 8000), rep("MS", 8000))
  fore <- c(data_mc[[4]][, "garch_t_garch_n"], data_mc[[4]][, "garch_t_garch_t"], data_mc[[4]][, "garch_t_gas_n"], data_mc[[4]][, "garch_t_gas_t"],
    data_mc[[4]][, "garch_t_sv_n"], data_mc[[4]][, "garch_t_sv_t"], data_mc[[4]][, "garch_t_ms_n"],  data_mc[[4]][, "garch_t_ms_t"],
    data_mc[[4]][, "gas_t_garch_n"], data_mc[[4]][, "gas_t_garch_t"], data_mc[[4]][, "gas_t_gas_n"], data_mc[[4]][, "gas_t_gas_t"],
    data_mc[[4]][, "gas_t_sv_n"], data_mc[[4]][, "gas_t_sv_t"], data_mc[[4]][, "gas_t_ms_n"], data_mc[[4]][, "gas_t_ms_t"],
    data_mc[[4]][, "sv_t_garch_n"], data_mc[[4]][, "sv_t_garch_t"], data_mc[[4]][, "sv_t_gas_n"], data_mc[[4]][, "sv_t_gas_t"],
    data_mc[[4]][, "sv_t_sv_n"], data_mc[[4]][, "sv_t_sv_t"], data_mc[[4]][, "sv_t_ms_n"], data_mc[[4]][, "sv_t_ms_t"],
    data_mc[[4]][, "ms_t_garch_n"], data_mc[[4]][, "ms_t_garch_t"], data_mc[[4]][, "ms_t_gas_n"], data_mc[[4]][, "ms_t_gas_t"],
    data_mc[[4]][, "ms_t_sv_n"], data_mc[[4]][, "ms_t_sv_t"], data_mc[[4]][, "ms_t_ms_n"], data_mc[[4]][, "ms_t_ms_t"])
  model <- rep(c(rep("GARCH-N", 1000), rep("GARCH-T", 1000), rep("GAS-N", 1000), rep("GAS-T", 1000), rep("SV-N", 1000), rep("SV-T", 1000), rep("MS-N", 1000), rep("MS-T", 1000)), 4)
  n <- rep(2500, 32000)
  distri <- rep("T", 32000)
  mc_4 <- data.frame(simul, dgp, fore, model, n, distri)
  
  simul <- c(rep(data_mc[[5]][, "garch"], 8), rep(data_mc[[5]][, "gas"], 8), rep(data_mc[[5]][, "sv"], 8), rep(data_mc[[5]][, "ms"], 8))
  dgp <- c(rep("GARCH", 8000), rep("GAS", 8000), rep("SV", 8000), rep("MS", 8000))
  fore <- c(data_mc[[5]][, "garch_n_garch_n"], data_mc[[5]][, "garch_n_garch_t"], data_mc[[5]][, "garch_n_gas_n"], data_mc[[5]][, "garch_n_gas_t"],
    data_mc[[5]][, "garch_n_sv_n"], data_mc[[5]][, "garch_n_sv_t"], data_mc[[5]][, "garch_n_ms_n"],  data_mc[[5]][, "garch_n_ms_t"],
    data_mc[[5]][, "gas_n_garch_n"], data_mc[[5]][, "gas_n_garch_t"], data_mc[[5]][, "gas_n_gas_n"], data_mc[[5]][, "gas_n_gas_t"],
    data_mc[[5]][, "gas_n_sv_n"], data_mc[[5]][, "gas_n_sv_t"], data_mc[[5]][, "gas_n_ms_n"], data_mc[[5]][, "gas_n_ms_t"],
    data_mc[[5]][, "sv_n_garch_n"], data_mc[[5]][, "sv_n_garch_t"], data_mc[[5]][, "sv_n_gas_n"], data_mc[[5]][, "sv_n_gas_t"],
    data_mc[[5]][, "sv_n_sv_n"], data_mc[[5]][, "sv_n_sv_t"], data_mc[[5]][, "sv_n_ms_n"], data_mc[[5]][, "sv_n_ms_t"],
    data_mc[[5]][, "ms_n_garch_n"], data_mc[[5]][, "ms_n_garch_t"], data_mc[[5]][, "ms_n_gas_n"], data_mc[[5]][, "ms_n_gas_t"],
    data_mc[[5]][, "ms_n_sv_n"], data_mc[[5]][, "ms_n_sv_t"], data_mc[[5]][, "ms_n_ms_n"], data_mc[[5]][, "ms_n_ms_t"])
  model <- rep(c(rep("GARCH-N", 1000), rep("GARCH-T", 1000), rep("GAS-N", 1000), rep("GAS-T", 1000), rep("SV-N", 1000), rep("SV-T", 1000), rep("MS-N", 1000), rep("MS-T", 1000)), 4)
  n <- rep(500, 32000)
  distri <- rep("N", 32000)
  mc_5 <- data.frame(simul, dgp, fore, model, n, distri)
  
  simul <- c(rep(data_mc[[6]][, "garch"], 8), rep(data_mc[[6]][, "gas"], 8), rep(data_mc[[6]][, "sv"], 8), rep(data_mc[[6]][, "ms"], 8))
  dgp <- c(rep("GARCH", 8000), rep("GAS", 8000), rep("SV", 8000), rep("MS", 8000))
  fore <- c(data_mc[[6]][, "garch_t_garch_n"], data_mc[[6]][, "garch_t_garch_t"], data_mc[[6]][, "garch_t_gas_n"], data_mc[[6]][, "garch_t_gas_t"],
    data_mc[[6]][, "garch_t_sv_n"], data_mc[[6]][, "garch_t_sv_t"], data_mc[[6]][, "garch_t_ms_n"],  data_mc[[6]][, "garch_t_ms_t"],
    data_mc[[6]][, "gas_t_garch_n"], data_mc[[6]][, "gas_t_garch_t"], data_mc[[6]][, "gas_t_gas_n"], data_mc[[6]][, "gas_t_gas_t"],
    data_mc[[6]][, "gas_t_sv_n"], data_mc[[6]][, "gas_t_sv_t"], data_mc[[6]][, "gas_t_ms_n"], data_mc[[6]][, "gas_t_ms_t"],
    data_mc[[6]][, "sv_t_garch_n"], data_mc[[6]][, "sv_t_garch_t"], data_mc[[6]][, "sv_t_gas_n"], data_mc[[6]][, "sv_t_gas_t"],
    data_mc[[6]][, "sv_t_sv_n"], data_mc[[6]][, "sv_t_sv_t"], data_mc[[6]][, "sv_t_ms_n"], data_mc[[6]][, "sv_t_ms_t"],
    data_mc[[6]][, "ms_t_garch_n"], data_mc[[6]][, "ms_t_garch_t"], data_mc[[6]][, "ms_t_gas_n"], data_mc[[6]][, "ms_t_gas_t"],
    data_mc[[6]][, "ms_t_sv_n"], data_mc[[6]][, "ms_t_sv_t"], data_mc[[6]][, "ms_t_ms_n"], data_mc[[6]][, "ms_t_ms_t"])
  model <- rep(c(rep("GARCH-N", 1000), rep("GARCH-T", 1000), rep("GAS-N", 1000), rep("GAS-T", 1000), rep("SV-N", 1000), rep("SV-T", 1000), rep("MS-N", 1000), rep("MS-T", 1000)), 4)
  n <- rep(500, 32000)
  distri <- rep("T", 32000)
  mc_6 <- data.frame(simul, dgp, fore, model, n, distri)
  
  mc_TRUE <- rbind(mc_5, mc_6, mc_1, mc_2, mc_3, mc_4) |>  mutate(ratio = simul/fore, outliers = TRUE)
  
  mc <- rbind(mc_FALSE, mc_TRUE) 
  
  
  mc_full <- mc |> filter(n > 500) |>
    mutate(outliers = as.factor(outliers), n = factor(n), model = factor(model))
  
  ggplot(mc_full) +
    geom_boxplot(aes(
      y = ratio,
      x = factor(n),
      colour = model,
      fill = outliers
    ), position = position_dodge(width = 0.75)) +
    facet_grid(distri ~ dgp, scales = "free_y") +
    labs(
      x = "Sample Size",
      y = expression(sigma[T + 1] / hat(sigma)[T + 1]),
      colour = "Estimated Model: ",
      fill = "Outliers"
    ) +
    scale_colour_brewer(palette = "Set2") +
    theme(legend.position = "bottom")
  
  
}



files_false <- list.files(path = './MonteCarlo', pattern = '*FALSE_BR.csv', full.names = TRUE)
results_FALSE <- making_tables(files_false) |> mutate(Outliers = FALSE)

files_true <- list.files(path = './MonteCarlo', pattern = '*TRUE_BR.csv', full.names = TRUE)
results_TRUE <- making_tables(files_true)  |> mutate(Outliers = TRUE)


results <- rbind(results_FALSE, results_TRUE)



options(pillar.sigfig = 5)
results |> 
  mutate(N = factor(N, levels = c(500, 1000, 2500))) |> 
  filter(DIST == "T") |> 
  group_by(DGP, N, Estim, Outliers) |> 
  reframe(MSE, QLIKE, `MSE LOG`, `MSE SD`, `MSE PROP`, MAE, `MAE LOG`, `MAE SD`, `MAE PROP`) |> 
  ungroup() |> 
  pivot_wider(
    names_from = Outliers, 
    values_from = c(MSE, QLIKE, `MSE LOG`, `MSE SD`, `MSE PROP`, MAE, `MAE LOG`, `MAE SD`, `MAE PROP`)) |> 
  select(DGP, N, Estim, ends_with("_FALSE"), ends_with("_TRUE")) |> 
  xtable::xtable(digits = 4) |> 
  print(file = "results_MC_T_BR.tex", , include.rownames = FALSE)



