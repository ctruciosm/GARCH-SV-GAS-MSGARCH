################################################################################
#####                    Descriptive Statistics                           #####
################################################################################
library(dplyr)
library(readxl)
library(stringr)
library(moments)
library(xtable)


descriptive_statistics <- function(x) {
  out <- c(min(x), quantile(x, 0.25), median(x), mean(x), quantile(x, 0.75),
  max(x), skewness(x), kurtosis(x), sd(x), acf(x, plot = FALSE)$acf[2])
  names(out) <- c("Min", "Q1", "Med", "Mean", "Q3", "Max", "Skew", "Kurt", "Sd", "ACF")
  return(out)
}


# Data
data <- read_excel("./Data/capire_daily_returns.xlsx", skip = 3, col_types = c("date", rep("numeric", 30)), na = c("", "-", NA)) |> 
  filter(Data > "2010-01-01" & Data < "2025-01-01") |> 
  filter(!if_all(where(is.numeric), is.na)) |> 
  select(where(~ !any(is.na(.x)))) |> 
  rename_with(~ str_remove(.x, "^(?s).*prov\n"), -Data) |> 
  select(-Data)



apply(data, 2, descriptive_statistics) |> t() |> round(4) |> xtable(digits = 4)