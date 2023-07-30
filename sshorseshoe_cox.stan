// generated with brms 2.16.3
functions {
  /* Efficient computation of the horseshoe prior
   * see Appendix C.1 in https://projecteuclid.org/euclid.ejs/1513306866
   * Args:
   *   z: standardized population-level coefficients
   *   lambda: local shrinkage parameters
   *   tau: global shrinkage parameter
   *   c2: slap regularization parameter
   * Returns:
   *   population-level coefficients following the horseshoe prior
   */
  vector horseshoe(vector z, vector lambda, real tau, real c2) {
    int K = rows(z);
    vector[K] lambda2 = square(lambda);
    vector[K] lambda_tilde = sqrt(c2 * lambda2 ./ (c2 + tau^2 * lambda2));
    return z .* lambda_tilde * tau;
  }
  vector lambdatilde(vector lambda, real c2, real tau) {
    int Kk = rows(lambda);
    vector[Kk] lambda2 = square(lambda);
    vector[Kk] lambda_tilde = sqrt(c2 * lambda2 ./ (c2 + tau^2 * lambda2));
    return lambda_tilde;
  }
  real zlambdatilde(real lambda, real c2, real tau) {
    real lambda2 = square(lambda);
    real lambda_tilde = sqrt(c2 * lambda2 ./ (c2 + tau^2 * lambda2));
    return lambda_tilde;
  }
  /* distribution functions of the Cox proportional hazards model
   * parameterize hazard(t) = baseline(t) * mu 
   * so that higher values of 'mu' imply lower survival times
   * Args:
   *   y: the response value; currently ignored as the relevant 
   *     information is passed via 'bhaz' and 'cbhaz'
   *   mu: positive location parameter
   *   bhaz: baseline hazard
   *   cbhaz: cumulative baseline hazard
   */
  real cox_lhaz(real y, real mu, real bhaz, real cbhaz) {
    return log(bhaz) + log(mu);
  }
  real cox_lccdf(real y, real mu, real bhaz, real cbhaz) {
    // equivalent to the log survival function
    return - cbhaz * mu;
  }
  real cox_lcdf(real y, real mu, real bhaz, real cbhaz) {
    return log1m_exp(cox_lccdf(y | mu, bhaz, cbhaz));
  }
  real cox_lpdf(real y, real mu, real bhaz, real cbhaz) { 
    return cox_lhaz(y, mu, bhaz, cbhaz) + cox_lccdf(y | mu, bhaz, cbhaz);
  }
  // Distribution functions of the Cox model in log parameterization
  real cox_log_lhaz(real y, real log_mu, real bhaz, real cbhaz) {
    return log(bhaz) + log_mu;
  }
  real cox_log_lccdf(real y, real log_mu, real bhaz, real cbhaz) {
    return - cbhaz * exp(log_mu);
  }
  real cox_log_lcdf(real y, real log_mu, real bhaz, real cbhaz) {
    return log1m_exp(cox_log_lccdf(y | log_mu, bhaz, cbhaz));
  }
  real cox_log_lpdf(real y, real log_mu, real bhaz, real cbhaz) { 
    return cox_log_lhaz(y, log_mu, bhaz, cbhaz) + 
           cox_log_lccdf(y | log_mu, bhaz, cbhaz);
  }
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=-1,upper=2> cens[N];  // indicates censoring
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  // data for flexible baseline functions
  int Kbhaz;  // number of basis functions
  // design matrix of the baseline function
  matrix[N, Kbhaz] Zbhaz;
  // design matrix of the cumulative baseline function
  matrix[N, Kbhaz] Zcbhaz;
  // a-priori concentration vector of baseline coefficients
  vector<lower=0>[Kbhaz] con_sbhaz;
  // data for the horseshoe prior
  real<lower=0> hs_df;  // local degrees of freedom
  real<lower=0> hs_df_global;  // global degrees of freedom
  real<lower=0> hs_df_slab;  // slab degrees of freedom
  real<lower=0> hs_scale_global;  // global prior scale
  real<lower=0> hs_scale_slab;  // slab prior scale
  int prior_only;  // should the likelihood be ignored?
  int nsnp;
  int nvar;
  int nz;
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  // local parameters for horseshoe prior
  vector[Kc] zb;
  vector<lower=0>[Kc] hs_local;
  real Intercept;  // temporary intercept for centered predictors
  simplex[Kbhaz] sbhaz;  // baseline coefficients
  // horseshoe shrinkage parameters
  real<lower=0> hs_global;  // global shrinkage parameters
  real<lower=0> hs_global_inter;
  real<lower=0> hs_global_z;
  real<lower=0> hs_slab;  // slab regularization parameter
}
transformed parameters {
  vector[Kc] b;  // population-level effects
  vector[nsnp] hlambdatilde=lambdatilde(hs_local[1:nsnp], hs_scale_slab^2 * hs_slab, hs_global) * zlambdatilde(hs_local[nz], hs_scale_slab^2 * hs_slab, hs_global);
  // compute actual regression coefficients
  b[1:nsnp] = horseshoe(zb[1:nsnp], hs_local[1:nsnp], hs_global, hs_scale_slab^2 * hs_slab);
  b[(nsnp+1):nvar] = horseshoe(zb[(nsnp+1):nvar], hs_local[(nsnp+1):nvar], hs_global_z, hs_scale_slab^2 * hs_slab);
  b[(nvar+1):Kc] = horseshoe(hlambdatilde[1:nsnp] .* zb[(nvar+1):Kc], hs_local[(nvar+1):Kc], hs_global_inter, hs_scale_slab^2 * hs_slab);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // compute values of baseline function
    vector[N] bhaz = Zbhaz * sbhaz;
    // compute values of cumulative baseline function
    vector[N] cbhaz = Zcbhaz * sbhaz;
    // initialize linear predictor term
    vector[N] mu = Intercept + Xc * b;
    for (n in 1:N) {
    // special treatment of censored data
      if (cens[n] == 0) {
        target += cox_log_lpdf(Y[n] | mu[n], bhaz[n], cbhaz[n]);
      } else if (cens[n] == 1) {
        target += cox_log_lccdf(Y[n] | mu[n], bhaz[n], cbhaz[n]);
      } else if (cens[n] == -1) {
        target += cox_log_lcdf(Y[n] | mu[n], bhaz[n], cbhaz[n]);
      }
    }
  }
  // priors including constants
  target += std_normal_lpdf(zb);
  target += student_t_lpdf(hs_local | hs_df, 0, 1)
    - rows(hs_local) * log(0.5);
  target += student_t_lpdf(Intercept | 3, -1.4, 2.5);
  target += dirichlet_lpdf(sbhaz | con_sbhaz);
  target += student_t_lpdf(hs_global | hs_df_global, 0, hs_scale_global) - 1 * log(0.5);
  target += student_t_lpdf(hs_global_z | hs_df_global, 0, hs_scale_global) - 1 * log(0.5);
  target += student_t_lpdf(hs_global_inter | hs_df_global, 0, hs_scale_global) - 1 * log(0.5);
  target += inv_gamma_lpdf(hs_slab | 0.5 * hs_df_slab, 0.5 * hs_df_slab);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
}
