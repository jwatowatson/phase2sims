require(rstan)

lin_model =  '
data {
    int<lower=1> Nobs;            //the number of observations
    int<lower=1> Ncens;           //the number of observations
    int<lower=1> J;               //the number of individuals
    real CT_threshold;

    int<lower=1,upper=J> id_obs[Nobs];   //vector of individual indices
    vector[Nobs] t_obs;                  //time since start of treatment
    vector[Nobs] Trt_obs;                //treatment assigned (0: control; 1:intervention)
    vector[Nobs] y_obs;                  //the raw ct values

    int<lower=1,upper=J> id_cens[Ncens];   //vector of individual indices
    vector[Ncens] t_cens;                  //time since start of treatment
    vector[Ncens] Trt_cens;                //treatment assigned (0: control; 1:intervention)
    vector[Ncens] y_cens;                  //the raw ct values

    real<lower=0> prior_Trt_sd;
    real<lower=0> prior_CT_intercept;
    real<lower=0> prior_sigma_sd;
    real prior_mean_slope;
    real<lower=0> prior_tau_intercept_sd;
    real<lower=0> prior_tau_slope_sd;
    real<lower=0> prior_tau_intercept_mean;
    real prior_tau_slope_mean;
  }
parameters {
  real beta_0;                // population-level regression intercept
  real beta_t;                // population-level regression time coef
  real beta_Trt;              // population-level regression treatment coef
  real<lower=0> tau_0;        //the standard deviation of the regression coefficients
  real<lower=0> tau_t;        //the standard deviation of the regression coefficients
  vector[J] beta_0_id;
  vector[J] beta_t_id;
  real<lower=0> sigma;        //standard deviation of the CT values
}
transformed parameters {
  vector[Nobs] mu_obs;     //linear predictor
  vector[Ncens] mu_cens;     //linear predictor
  for(n in 1:Nobs){
    mu_obs[n] = (beta_0 + beta_0_id[id_obs[n]]) + (beta_t*(1+beta_Trt*Trt_obs[n]) + beta_t_id[id_obs[n]])*t_obs[n];
  }
  for(n in 1:Ncens){
    mu_cens[n] = (beta_0 + beta_0_id[id_cens[n]]) + (beta_t + beta_t_id[id_cens[n]])*(1+beta_Trt*Trt_cens[n])*t_cens[n];
  }
}
model {
  //priors
  tau_0 ~ normal(prior_tau_intercept_mean, prior_tau_intercept_sd);
  tau_t ~ normal(prior_tau_slope_mean, prior_tau_slope_sd);
  sigma ~ normal(1, prior_sigma_sd);
  beta_0 ~ normal(prior_CT_intercept, 3);
  beta_t ~ normal(prior_mean_slope, 1);
  beta_Trt ~ normal(0, prior_Trt_sd);
  beta_0_id ~ normal(0, tau_0);
  beta_t_id ~ normal(0, tau_t);

  //likelihood
  // Observed ct values
  y_obs ~ student_t(7, mu_obs, sigma);
  // censored CT values
  for(n in 1:Ncens){
    target += student_t_lccdf(CT_threshold | 7, mu_cens[n], sigma);
  }

}

'

stan_lin_model = stan_model(model_code = lin_model)
