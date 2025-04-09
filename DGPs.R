################################################################################
#####                                DGPs                                 #####
################################################################################


garch_sim <- function(n, params, distri) { 
  n_burnin <- 500
  n_tot <- n_burnin + n
  sigma2 <- rep(NA, n_tot)
  ret <- rep(NA, n_tot)
  
  if (distri == "std") {
    epsilon <- rt(n_tot, df = params[4]) * sqrt((params[4] - 2) / params[4])
  } else {
    epsilon <-  rnorm(n_tot)
  }

  sigma2[1] <- params[1] / (1 - params[2] - params[3])
  ret[1] <- sqrt(sigma2[1])*epsilon[1]
  for (i in 2:n_tot) {
    sigma2[i] <- params[1] + params[2] * ret[i - 1]^2 + params[3] *  sigma2[i - 1]
    ret[i] <- sqrt(sigma2[i])*epsilon[i]
  }
  return(list(returns = ret[-c(1:n_burnin)], volatility = sqrt(sigma2[-c(1:n_burnin)]), e = epsilon[-c(1:n_burnin)]))
}
#dados <- garch_sim(10000, c(0.01, 0.1, 0.86, 7), "std")

sv_sim <- function(n, params, distri) {
  n_burnin <- 500
  n_tot <- n_burnin + n
  h <- rep(NA, n_tot)
  ret <- rep(NA, n_tot)
  
  eta <- rnorm(n_tot)
  if (distri == "std") {
    epsilon <- rt(n_tot, df = params[4]) * sqrt((params[4] - 2) / params[4])
  } else {
    epsilon <-  rnorm(n_tot)
  }
  
  h[1] <- rnorm(params[1], params[3]^2 / (1 - params[2]^2))
  ret[1] <- exp(h[1]/2)*epsilon[1]
  for (i in 2:n_tot) {
    h[i] <- params[1] + params[2] * (h[i - 1] - params[1]) + params[3] * eta[i - 1]
    ret[i] <- exp(h[i]/2)*epsilon[i]
  }
  return(list(returns = ret[-c(1:n_burnin)], volatility = exp(h[-c(1:n_burnin)]/2), e = epsilon[-c(1:n_burnin)]))
}
#dados <- sv_sim(10000, c(1.68, 0.95, 0.23, 7), "std")

sv_sim_package <- function(n, params, distri) {
  if (distri == "std") {
    aux <- svsim(n, params[1], params[2], params[3], params[4])
    ret <- aux$y
    vol <- aux$vol
    e <- ret/vol
  } else {
    aux <- svsim(n, params[1], params[2], params[3])
    ret <- aux$y
    vol <- aux$vol
    e <- ret/vol
  }
  return(list(returns = ret, volatility = vol, e = e))
}

gas_sim <- function(n, params, distri) {
  n_burnin <- 500
  n_tot <- n_burnin + n
 
  if (distri == "norm") {
    A <- diag(c(0.0, params[2]))
    B <- diag(c(0.0, params[3]))
    k <- c(0, params[1])
    Sim <- UniGASSim(T.sim = n_tot, kappa = k, A = A, B = B, Dist = "norm", ScalingType = "Identity")
    out <- list(returns = Sim@GASDyn$vY[-c(1:n_burnin)], volatility = head(sqrt(Sim@GASDyn$mTheta[2, -c(1:n_burnin)]), n), e = Sim@GASDyn$mInnovations[1,-c(1:n_burnin)])
  } else {
    A <- diag(c(0.0, params[2], 0.0))
    B <- diag(c(0.0, params[3], 0.0))
    k <- c(0, params[1], params[4])
    Sim <- UniGASSim(T.sim = n_tot, kappa = k, A = A, B = B, Dist = "std", ScalingType = "Identity")
    out <- list(returns = Sim@GASDyn$vY[-c(1:n_burnin)], volatility = head(sqrt(Sim@GASDyn$mTheta[2, -c(1:n_burnin)]) * sqrt(7/5), n), e = Sim@GASDyn$mInnovations[1,-c(1:n_burnin)]/sqrt(7/5))
  }
  return(out)
}
#dados <- gas_sim(100000, c(0.04, 0.22, 0.96, -2.6625878), "std"). # -2.6625878 is equivalento to 7 d.f


msgarch_sim <- function(n, params, distri, P) {
  k <- 2
  n_burnin <- 500
  n_tot <- n_burnin + n 
  h <- matrix(NA, n_tot, k + 1)
  ret <- numeric(n_tot)
  Pt <- numeric(n_tot)
  M <- matrix(NA, 4, 4)
  I4 <- diag(1, 4)
  
  omega <- params[c(1, 4)]
  alpha <- params[c(2, 5)]
  beta  <- params[c(3, 6)]
  if (distri == "std") {
    nu    <- params[7]
    epsilon <- rt(n_tot, df = nu) * sqrt((nu - 2) / nu)
  } else {
    epsilon <-  rnorm(n_tot)
  }
  
  p <- P[1, 1]
  q <- P[2, 2]
  
  M[1, 1] <- P[1, 1] * (alpha[1] + beta[1])
  M[1, 2] <- 0
  M[1, 3] <- P[1, 2] * (alpha[1] + beta[1])
  M[1, 4] <- 0
  M[2, 1] <- P[1, 1] * alpha[2]
  M[2, 2] <- P[1, 1] * beta[2]
  M[2, 3] <- P[1, 2] * alpha[2]
  M[2, 4] <- P[1, 2] * beta[2]
  M[3, 1] <- P[2, 1] * beta[1]
  M[3, 2] <- P[2, 1] * alpha[1]
  M[3, 3] <- P[2, 2] * beta[1]
  M[3, 4] <- P[2, 2] * alpha[1]
  M[4, 1] <- 0;
  M[4, 2] <- P[2, 1] * (alpha[2] + beta[2])
  M[4, 3] <- 0
  M[4, 4] <- P[2, 2] * (alpha[2] + beta[2])
    
  Pt[1] <- (1 - q) / (2 - p - q)       
  pi_inf <- c(Pt[1], 1 - Pt[1])     
    
  h[1, 1:k] <- matrix(c(1, 0, 1, 0, 0, 1, 0, 1), 2, 4, byrow = TRUE) %*% solve(I4 - M) %*% kronecker(pi_inf, omega)
  h[1, k + 1] <- Pt[1] * h[1, 1] + (1 - Pt[1]) * h[1, 2]
    
  s <- numeric(n_tot)
  s[1] <- 1
  ret[1] <- epsilon[1] * sqrt(h[1, s[1]])
    
  if (distri == "norm") {
    for (i in 2:n_tot) {
      h[i, 1:k] <- omega + alpha * ret[i - 1]^2 + beta * h[i - 1, 1:k]
      Pt[i] <- probability_regime_given_time_n(p, q, sqrt(h[i - 1, 1:k]), ret[i - 1], Pt[i - 1])
      h[i, k + 1] <- Pt[i] * h[i, 1] + (1 - Pt[i]) * h[i, 2]
      s[i] <- sample(1:2, 1, prob = P[, s[i - 1]])
      ret[i] <- epsilon[i] * sqrt(h[i, s[i]])
    }
  } else {
    for (i in 2:n_tot) {
      h[i, 1:k] <- omega + alpha * ret[i - 1]^2 + beta * h[i - 1, 1:k]
      Pt[i] <- probability_regime_given_time_t(p, q, sqrt(h[i - 1, 1:k]), ret[i - 1], Pt[i - 1], 7)
      h[i, k + 1] <- Pt[i] * h[i, 1] + (1 - Pt[i]) * h[i, 2]
      s[i] <- sample(1:2, 1, prob = P[, s[i - 1]])
      ret[i] <- epsilon[i] * sqrt(h[i, s[i]])
    }
  }
  return(list(returns = ret[-c(1:n_burnin)], volatility = sqrt(h[-c(1:n_burnin), ]), e = epsilon[-c(1:n_burnin)], s = s[-c(1:n_burnin)]))
}
  