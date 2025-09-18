#include<RcppArmadillo.h>
#include<Rmath.h>

using namespace Rcpp;

// [[Rcpp::export]]
SEXP dcsn_like(NumericVector params, NumericVector r){
  int n = r.size();
  NumericVector lambda(n), lambda_star(n), log_lik(n - 1);
  lambda_star[0] = 0.0;
  lambda[0] = params[0];
  for(int i = 1; i < n; i++){
    lambda_star[i] = params[1] * lambda_star[i - 1] + params[2] * (pow(r[i - 1], 2) / exp(2 * lambda[i - 1]) - 1);
    lambda[i] = params[0] + lambda_star[i];
    log_lik[i - 1] = 0.5* pow(r[i], 2) / exp(2 * lambda[i]) + lambda[i];
  }
  double out = sum(log_lik);
  return Rcpp::wrap(out);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
SEXP grid_dcsn(NumericVector y){
  NumericVector coeff(3), vi(3);
  double omega, phi, kapa;
  double omega_min = -0.01, omega_max= 0.05, phi_min = 0.25,  phi_max = 0.99, kapa_min = 0.05,  kapa_max = 0.25, n_omega = 5, n_phi = 5,  n_kapa = 5;
  double ml = 100000000, nml;
  double lm_omega = (omega_max - omega_min) / n_omega;
  double lm_phi = (phi_max - phi_min) / n_phi;
  double lm_kapa = (kapa_max - kapa_min) / n_kapa;
  Rcpp::Function dcsn_like("dcsn_like");
  
  for(int no = 0; no < n_omega; no++){
    for(int np = 0; np < n_phi; np++){
      for(int nk =0; nk < n_kapa; nk++){
        omega = omega_min + no * lm_omega;
        phi = phi_min + np * lm_phi; 
        kapa = kapa_min +nk * lm_kapa;
        coeff[0] = omega;
        coeff[1] = phi;
        coeff[2] = kapa;
        nml = Rcpp::as<double>(dcsn_like(coeff, y));
        if (nml < ml){
          vi[0] = coeff[0];
          vi[1] = coeff[1];
          vi[2] = coeff[2];
          ml=nml;
        }
      }
    }
  }
  return(vi);
}
