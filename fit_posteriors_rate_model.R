library(rstan)
options(mc.cores = 4)

load(file = 'input_stan_model.RData')
### change with the rate model
ct_model = stan_model("fit_posteriors_symptomatic_rate.stan")
ct_data=list(
  N=nrow(indiv_data),
  n_id=length(unique(indiv_data$id_clean)),
  lod=as.list(global_pars)$lod,
  id=indiv_data$id_clean,
  symp=as.list(prior_pars)$symp,
  t=indiv_data$t,
  y=indiv_data$y,
  tpsd=as.list(prior_pars)$tpsd,
  dpmean_prior=as.list(prior_pars)$dpmean_prior,
  dpsd_prior=as.list(prior_pars)$dpsd_prior,
  wpmax=as.list(prior_pars)$wpmax,
  wpmean_prior=as.list(prior_pars)$wpmean_prior,
  wpsd_prior=as.list(prior_pars)$wpsd_prior,
  alphamean_prior = -1.5,
  alphasd_prior = .5,
  sigma_max=as.list(prior_pars)$sigma_max,
  sigma_prior_scale=as.list(prior_pars)$sigma_prior_scale,
  lambda=as.list(prior_pars)$lambda,
  fpmean=as.list(prior_pars)$fpmean,
  epsilon=(indiv_data$adjusted)*(as.list(global_pars)$adjusted_sd))


ct_fit_rate = sampling(ct_model,
                       data = ct_data,
                       iter=5000, chains=4,
                       control = list(adapt_delta=0.99))

thetas = extract(ct_fit_rate)
save(thetas, file = 'stan_out_rate.RData')

pdf('Best_profiles.pdf', width = 12, height = 12)
par(mfrow=c(5,5), las=1)
for(id in unique(ct_data$id)){
  ind = which(ct_data$id == id & ct_data$t>= -7)
  if(length(ind)>4){
    plot(ct_data$t[ind], ct_data$y[ind], ylim = c(40,15),pch=16,
         xlab = 'Time from peak observed Ct', ylab='Ct')
    lines(ct_data$t[ind], 40-colMeans(thetas$mu)[ind])
    }
}
dev.off()
