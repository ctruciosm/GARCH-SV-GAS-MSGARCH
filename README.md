# Replication Codes

This repository contains the codes used to replicate the results in Marques and Trucíos (2025).

### Empirical Application

- `Empirical_Application.R` computes one-step-ahead conditional variances using GARCH, SV, GAS, and MSGARCH models, each with both standardized Normal and Student-t innovation distributions.
- `Tables_App.R` generates the tables and performs the Model Confidence Set procedure for the empirical application results.

> Five-minute realized variances are freely available from the [CaPiRe](https://capire.stat.unipd.it/) database. Daily returns were obtained from Economatica.

### Monte Carlo Simulation

- `MonteCarlo_GARCH-GAS-SV-MS.R` runs the one-step-ahead forecasting experiment. To use the code, modify the parameters accordingly, or execute it in batch mode using the following command:  
  `R CMD BATCH "--args n=2500 type=BR outliers=FALSE" MonteCarlo_GARCH-GAS-SV-MS.R MonteCarlo_2500_BR_FALSE.txt &`  
  (You can change `BR` to `US`, `FALSE` to `TRUE`, or adjust the sample size as desired.)
- `Tables_MC.R` reproduces the results shown in Tables 3 to 6 of the paper.
- `Model_Confidence_Set_MC.R` performs the Model Confidence Set procedure for the simulation study.

### Auxiliary Functions

- `DGPs.R` defines the data-generating processes used in the simulations.
- `Utils_GARCH-GAS-SV.R` contains additional functions for model estimation and forecasting.
- `Descriptive_Statistics` displays the descriptive statistics in Table 7.


## References

Marques, F. and Trucíos C. (2025). *"GARCH, GAS, SV, and MSGARCH models: Do we really need all of them for forecasting daily volatility?".* Submitted
