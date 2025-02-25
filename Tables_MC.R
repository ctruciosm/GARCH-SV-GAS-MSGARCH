################################################################################
#####                       Tables Monte Carlo                             #####
################################################################################
library(dplyr)
library(ggplot2)
source("Utils_GARCH-GAS-SV.R")


files_names <- list.files(path = './MonteCarlo', pattern = '*FALSE_BR.csv', full.names = TRUE)
data_mc <- list()
for (name in files_names) data_mc[[name]] <- read.csv(name)


simul <- c(rep(data_mc[[1]][, "garch"], 4), rep(data_mc[[1]][, "gas"], 4), rep(data_mc[[1]][, "sv"], 4), rep(data_mc[[1]][, "ms"], 4))
dgp <- c(rep("GARCH", 4000), rep("GAS", 4000), rep("SV", 4000), rep("MS", 4000))
fore <- c(data_mc[[1]][, "garch_garch"], data_mc[[1]][, "garch_gas"], data_mc[[1]][, "garch_sv"], data_mc[[1]][, "garch_ms"],
          data_mc[[1]][, "gas_garch"], data_mc[[1]][, "gas_gas"], data_mc[[1]][, "gas_sv"], data_mc[[1]][, "gas_ms"],
          data_mc[[1]][, "sv_garch"], data_mc[[1]][, "sv_gas"], data_mc[[1]][, "sv_sv"], data_mc[[1]][, "sv_ms"],
          data_mc[[1]][, "ms_garch"], data_mc[[1]][, "ms_gas"], data_mc[[1]][, "ms_sv"], data_mc[[1]][, "ms_ms"])
model <- rep(c(rep("GARCH", 1000), rep("GAS", 1000), rep("SV", 1000), rep("MS", 1000)), 4)
n <- rep(1000, 16000)
distri <- rep("N", 16000)
mc_1 <- data.frame(simul, dgp, fore, model, n, distri)

simul <- c(rep(data_mc[[2]][, "garch"], 4), rep(data_mc[[2]][, "gas"], 4), rep(data_mc[[2]][, "sv"], 4), rep(data_mc[[2]][, "ms"], 4))
dgp <- c(rep("GARCH", 4000), rep("GAS", 4000), rep("SV", 4000), rep("MS", 4000))
fore <- c(data_mc[[2]][, "garch_garch"], data_mc[[2]][, "garch_gas"], data_mc[[2]][, "garch_sv"], data_mc[[2]][, "garch_ms"],
  data_mc[[2]][, "gas_garch"], data_mc[[2]][, "gas_gas"], data_mc[[2]][, "gas_sv"], data_mc[[2]][, "gas_ms"],
  data_mc[[2]][, "sv_garch"], data_mc[[2]][, "sv_gas"], data_mc[[2]][, "sv_sv"], data_mc[[2]][, "sv_ms"],
  data_mc[[2]][, "ms_garch"], data_mc[[2]][, "ms_gas"], data_mc[[2]][, "ms_sv"], data_mc[[2]][, "ms_ms"])
model <- rep(c(rep("GARCH", 1000), rep("GAS", 1000), rep("SV", 1000), rep("MS", 1000)), 4)
n <- rep(1000, 16000)
distri <- rep("T", 16000)
mc_2 <- data.frame(simul, dgp, fore, model, n, distri)

simul <- c(rep(data_mc[[3]][, "garch"], 4), rep(data_mc[[3]][, "gas"], 4), rep(data_mc[[3]][, "sv"], 4), rep(data_mc[[3]][, "ms"], 4))
dgp <- c(rep("GARCH", 4000), rep("GAS", 4000), rep("SV", 4000), rep("MS", 4000))
fore <- c(data_mc[[3]][, "garch_garch"], data_mc[[3]][, "garch_gas"], data_mc[[3]][, "garch_sv"], data_mc[[3]][, "garch_ms"],
  data_mc[[3]][, "gas_garch"], data_mc[[3]][, "gas_gas"], data_mc[[3]][, "gas_sv"], data_mc[[3]][, "gas_ms"],
  data_mc[[3]][, "sv_garch"], data_mc[[3]][, "sv_gas"], data_mc[[3]][, "sv_sv"], data_mc[[3]][, "sv_ms"],
  data_mc[[3]][, "ms_garch"], data_mc[[3]][, "ms_gas"], data_mc[[3]][, "ms_sv"], data_mc[[3]][, "ms_ms"])
model <- rep(c(rep("GARCH", 1000), rep("GAS", 1000), rep("SV", 1000), rep("MS", 1000)), 4)
n <- rep(2500, 16000)
distri <- rep("N", 16000)
mc_3 <- data.frame(simul, dgp, fore, model, n, distri)

simul <- c(rep(data_mc[[4]][, "garch"], 4), rep(data_mc[[4]][, "gas"], 4), rep(data_mc[[4]][, "sv"], 4), rep(data_mc[[4]][, "ms"], 4))
dgp <- c(rep("GARCH", 4000), rep("GAS", 4000), rep("SV", 4000), rep("MS", 4000))
fore <- c(data_mc[[4]][, "garch_garch"], data_mc[[4]][, "garch_gas"], data_mc[[4]][, "garch_sv"], data_mc[[4]][, "garch_ms"],
  data_mc[[4]][, "gas_garch"], data_mc[[4]][, "gas_gas"], data_mc[[4]][, "gas_sv"], data_mc[[4]][, "gas_ms"],
  data_mc[[4]][, "sv_garch"], data_mc[[4]][, "sv_gas"], data_mc[[4]][, "sv_sv"], data_mc[[4]][, "sv_ms"],
  data_mc[[4]][, "ms_garch"], data_mc[[4]][, "ms_gas"], data_mc[[4]][, "ms_sv"], data_mc[[4]][, "ms_ms"])
model <- rep(c(rep("GARCH", 1000), rep("GAS", 1000), rep("SV", 1000), rep("MS", 1000)), 4)
n <- rep(2500, 16000)
distri <- rep("T", 16000)
mc_4 <- data.frame(simul, dgp, fore, model, n, distri)

simul <- c(rep(data_mc[[5]][, "garch"], 4), rep(data_mc[[5]][, "gas"], 4), rep(data_mc[[5]][, "sv"], 4), rep(data_mc[[5]][, "ms"], 4))
dgp <- c(rep("GARCH", 4000), rep("GAS", 4000), rep("SV", 4000), rep("MS", 4000))
fore <- c(data_mc[[5]][, "garch_garch"], data_mc[[5]][, "garch_gas"], data_mc[[5]][, "garch_sv"], data_mc[[5]][, "garch_ms"],
  data_mc[[5]][, "gas_garch"], data_mc[[5]][, "gas_gas"], data_mc[[5]][, "gas_sv"], data_mc[[5]][, "gas_ms"],
  data_mc[[5]][, "sv_garch"], data_mc[[5]][, "sv_gas"], data_mc[[5]][, "sv_sv"], data_mc[[5]][, "sv_ms"],
  data_mc[[5]][, "ms_garch"], data_mc[[5]][, "ms_gas"], data_mc[[5]][, "ms_sv"], data_mc[[5]][, "ms_ms"])
model <- rep(c(rep("GARCH", 1000), rep("GAS", 1000), rep("SV", 1000), rep("MS", 1000)), 4)
n <- rep(500, 16000)
distri <- rep("N", 16000)
mc_5 <- data.frame(simul, dgp, fore, model, n, distri)

simul <- c(rep(data_mc[[6]][, "garch"], 4), rep(data_mc[[6]][, "gas"], 4), rep(data_mc[[6]][, "sv"], 4), rep(data_mc[[6]][, "ms"], 4))
dgp <- c(rep("GARCH", 4000), rep("GAS", 4000), rep("SV", 4000), rep("MS", 4000))
fore <- c(data_mc[[6]][, "garch_garch"], data_mc[[6]][, "garch_gas"], data_mc[[6]][, "garch_sv"], data_mc[[6]][, "garch_ms"],
  data_mc[[6]][, "gas_garch"], data_mc[[6]][, "gas_gas"], data_mc[[6]][, "gas_sv"], data_mc[[6]][, "gas_ms"],
  data_mc[[6]][, "sv_garch"], data_mc[[6]][, "sv_gas"], data_mc[[6]][, "sv_sv"], data_mc[[6]][, "sv_ms"],
  data_mc[[6]][, "ms_garch"], data_mc[[6]][, "ms_gas"], data_mc[[6]][, "ms_sv"], data_mc[[6]][, "ms_ms"])
model <- rep(c(rep("GARCH", 1000), rep("GAS", 1000), rep("SV", 1000), rep("MS", 1000)), 4)
n <- rep(500, 16000)
distri <- rep("T", 16000)
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


