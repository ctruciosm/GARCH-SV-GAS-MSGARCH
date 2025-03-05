################################################################################
#####                       Tables Monte Carlo                             #####
################################################################################
library(dplyr)
library(ggplot2)
source("Utils_GARCH-GAS-SV.R")


files_names <- list.files(path = './MonteCarlo', pattern = '*TRUE_BR.csv', full.names = TRUE)
data_mc <- list()
for (name in files_names) data_mc[[name]] <- read.csv(name)


simul <- c(rep(data_mc[[1]][, "garch"], 8), rep(data_mc[[1]][, "gas"], 8), rep(data_mc[[1]][, "sv"], 8), rep(data_mc[[1]][, "ms"], 8))
dgp <- c(rep("GARCH", 8000), rep("GAS", 8000), rep("SV", 8000), rep("MS", 8000))
fore <- c(data_mc[[1]][, "garch_n_garch_n"], data_mc[[1]][, "garch_n_garch_t"], data_mc[[1]][, "garch_n_gas_n"], data_mc[[1]][, "garch_n_gas_t"],
          data_mc[[1]][, "garch_n_sv_n"], data_mc[[1]][, "garch_n_sv_t"], data_mc[[1]][, "garch_n_ms_n"], data_mc[[1]][, "garch_n_ms_t"],
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
  data_mc[[2]][, "garch_t_sv_n"], data_mc[[2]][, "garch_t_sv_t"], data_mc[[2]][, "garch_t_ms_n"], data_mc[[2]][, "garch_t_ms_t"],
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
  data_mc[[3]][, "garch_n_sv_n"], data_mc[[3]][, "garch_n_sv_t"], data_mc[[3]][, "garch_n_ms_n"], data_mc[[3]][, "garch_n_ms_t"],
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
  data_mc[[4]][, "garch_t_sv_n"], data_mc[[4]][, "garch_t_sv_t"], data_mc[[4]][, "garch_t_ms_n"], data_mc[[4]][, "garch_t_ms_t"],
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
  data_mc[[5]][, "garch_n_sv_n"], data_mc[[5]][, "garch_n_sv_t"], data_mc[[5]][, "garch_n_ms_n"], data_mc[[5]][, "garch_n_ms_t"],
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
  data_mc[[6]][, "garch_t_sv_n"], data_mc[[6]][, "garch_t_sv_t"], data_mc[[6]][, "garch_t_ms_n"], data_mc[[6]][, "garch_t_ms_t"],
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

mc <- rbind(mc_5, mc_6, mc_1, mc_2, mc_3, mc_4) |> 
  mutate(ratio = fore/simul)

mc |> 
ggplot() +
  geom_boxplot(aes(y = ratio, x = factor(n), colour = model)) +
  facet_grid(distri~dgp, scales = "free_y") +
  labs(x = "Sample Size", 
    y = expression(hat(sigma)[T+1] / sigma[T+1]),
    colour = "Estimated Model: ") +
  ylim(c(0, 3.5)) + 
  scale_colour_brewer(palette = "Set2") +
  theme(legend.position = "bottom")


results_1 <- round(rbind(metrics(data_mc[[1]][, "garch_garch"], data_mc[[1]][, "garch"]),
                 metrics(data_mc[[1]][, "garch_gas"], data_mc[[1]][, "garch"]),
                 metrics(data_mc[[1]][, "garch_sv"], data_mc[[1]][, "garch"]),
                 metrics(data_mc[[1]][, "gas_garch"], data_mc[[1]][, "gas"]),
                 metrics(data_mc[[1]][, "gas_gas"], data_mc[[1]][, "gas"]),
                 metrics(data_mc[[1]][, "gas_sv"], data_mc[[1]][, "gas"]),
                 metrics(data_mc[[1]][, "sv_garch"], data_mc[[1]][, "sv"]),
                 metrics(data_mc[[1]][, "sv_gas"], data_mc[[1]][, "sv"]),
                 metrics(data_mc[[1]][, "sv_sv"], data_mc[[1]][, "sv"])), 4) |> data.frame() |> 
  mutate(n = 1000, distri = "N", 
    DGP = c(rep("GARCH", 3), rep("GAS", 3), rep("SV", 3)),
    Estim = rep(c("GARCH", "GAS", "SV"), 3))

results_2 <- round(rbind(metrics(data_mc[[2]][, "garch_garch"], data_mc[[2]][, "garch"]),
  metrics(data_mc[[2]][, "garch_gas"], data_mc[[2]][, "garch"]),
  metrics(data_mc[[2]][, "garch_sv"], data_mc[[2]][, "garch"]),
  metrics(data_mc[[2]][, "gas_garch"], data_mc[[2]][, "gas"]),
  metrics(data_mc[[2]][, "gas_gas"], data_mc[[2]][, "gas"]),
  metrics(data_mc[[2]][, "gas_sv"], data_mc[[2]][, "gas"]),
  metrics(data_mc[[2]][, "sv_garch"], data_mc[[2]][, "sv"]),
  metrics(data_mc[[2]][, "sv_gas"], data_mc[[2]][, "sv"]),
  metrics(data_mc[[2]][, "sv_sv"], data_mc[[2]][, "sv"])), 4) |> data.frame() |> 
  mutate(n = 1000, distri = "T",
    DGP = c(rep("GARCH", 3), rep("GAS", 3), rep("SV", 3)),
    Estim = rep(c("GARCH", "GAS", "SV"), 3))


results_3 <- round(rbind(metrics(data_mc[[3]][, "garch_garch"], data_mc[[3]][, "garch"]),
  metrics(data_mc[[3]][, "garch_gas"], data_mc[[3]][, "garch"]),
  metrics(data_mc[[3]][, "garch_sv"], data_mc[[3]][, "garch"]),
  metrics(data_mc[[3]][, "gas_garch"], data_mc[[3]][, "gas"]),
  metrics(data_mc[[3]][, "gas_gas"], data_mc[[3]][, "gas"]),
  metrics(data_mc[[3]][, "gas_sv"], data_mc[[3]][, "gas"]),
  metrics(data_mc[[3]][, "sv_garch"], data_mc[[3]][, "sv"]),
  metrics(data_mc[[3]][, "sv_gas"], data_mc[[3]][, "sv"]),
  metrics(data_mc[[3]][, "sv_sv"], data_mc[[3]][, "sv"])), 4)|> data.frame() |> 
  mutate(n = 2500, distri = "N", 
    DGP = c(rep("GARCH", 3), rep("GAS", 3), rep("SV", 3)),
    Estim = rep(c("GARCH", "GAS", "SV"), 3))

results_4 <- round(rbind(metrics(data_mc[[4]][, "garch_garch"], data_mc[[4]][, "garch"]),
  metrics(data_mc[[4]][, "garch_gas"], data_mc[[4]][, "garch"]),
  metrics(data_mc[[4]][, "garch_sv"], data_mc[[4]][, "garch"]),
  metrics(data_mc[[4]][, "gas_garch"], data_mc[[4]][, "gas"]),
  metrics(data_mc[[4]][, "gas_gas"], data_mc[[4]][, "gas"]),
  metrics(data_mc[[4]][, "gas_sv"], data_mc[[4]][, "gas"]),
  metrics(data_mc[[4]][, "sv_garch"], data_mc[[4]][, "sv"]),
  metrics(data_mc[[4]][, "sv_gas"], data_mc[[4]][, "sv"]),
  metrics(data_mc[[4]][, "sv_sv"], data_mc[[4]][, "sv"])), 4)|> data.frame() |> 
  mutate(n = 2500, distri = "T", 
    DGP = c(rep("GARCH", 3), rep("GAS", 3), rep("SV", 3)),
    Estim = rep(c("GARCH", "GAS", "SV"), 3))


results_5 <- round(rbind(metrics(data_mc[[5]][, "garch_garch"], data_mc[[5]][, "garch"]),
  metrics(data_mc[[5]][, "garch_gas"], data_mc[[5]][, "garch"]),
  metrics(data_mc[[5]][, "garch_sv"], data_mc[[5]][, "garch"]),
  metrics(data_mc[[5]][, "gas_garch"], data_mc[[5]][, "gas"]),
  metrics(data_mc[[5]][, "gas_gas"], data_mc[[5]][, "gas"]),
  metrics(data_mc[[5]][, "gas_sv"], data_mc[[5]][, "gas"]),
  metrics(data_mc[[5]][, "sv_garch"], data_mc[[5]][, "sv"]),
  metrics(data_mc[[5]][, "sv_gas"], data_mc[[5]][, "sv"]),
  metrics(data_mc[[5]][, "sv_sv"], data_mc[[5]][, "sv"])), 4)|> data.frame() |> 
  mutate(n = 500, distri = "N", 
    DGP = c(rep("GARCH", 3), rep("GAS", 3), rep("SV", 3)),
    Estim = rep(c("GARCH", "GAS", "SV"), 3))


results_6 <- round(rbind(metrics(data_mc[[6]][, "garch_garch"], data_mc[[6]][, "garch"]),
  metrics(data_mc[[6]][, "garch_gas"], data_mc[[6]][, "garch"]),
  metrics(data_mc[[6]][, "garch_sv"], data_mc[[6]][, "garch"]),
  metrics(data_mc[[6]][, "gas_garch"], data_mc[[6]][, "gas"]),
  metrics(data_mc[[6]][, "gas_gas"], data_mc[[6]][, "gas"]),
  metrics(data_mc[[6]][, "gas_sv"], data_mc[[6]][, "gas"]),
  metrics(data_mc[[6]][, "sv_garch"], data_mc[[6]][, "sv"]),
  metrics(data_mc[[6]][, "sv_gas"], data_mc[[6]][, "sv"]),
  metrics(data_mc[[6]][, "sv_sv"], data_mc[[6]][, "sv"])), 4)|> data.frame() |> 
  mutate(n = 500, distri = "T", 
    DGP = c(rep("GARCH", 3), rep("GAS", 3), rep("SV", 3)),
    Estim = rep(c("GARCH", "GAS", "SV"), 3))

results <- rbind(results_1, results_2, results_3, results_4, results_5, results_6)

colnames(results) <- c("MSE", "QLIKE1", "QLIKE2", "MSE LOG", "MSE SD", "MSE PROP", "MAE", "MAE LOG", "MAE SD", "MAE PROP", "N", "DIST", "DGP", "Estim")


