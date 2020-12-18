### Generate synthetic data from the posterior distribution
my_mufun = function(t, alpha, dp){
  out = rep(0, length(t))
  ind2 = t>=0
  out[ind2] = dp - alpha*t[ind2]
  return(out)
}

# simulate Ct data given a set of parameters
simulate_data = function(N, effect, thetas, xs){
  sim_dat = array(dim = c(N, length(xs)))
  if(length(effect)==1) effect = rep(effect, N)
  if(length(effect) != N) stop('Effect needs to be of length N')

  for(i in 1:N){
    ind1 = sample(1:length(thetas$dpmeanS),1)
    sim_dat[i,] = round(my_mufun(t = xs,
                                 alpha = thetas$alphameanS[ind1]*effect[i],
                                 dp = thetas$dpmeanS[ind1]) +
                          rnorm(length(xs), sd=thetas$sigma[ind1]))
  }
  sim_dat[ sim_dat > 40 ] = 40
  return(sim_dat)
}


# estimate rate of clearance from Ct data
compute_rate = function(sim_data,xs){
  sim_rate = array(dim = nrow(sim_data))
  for(k in 1:nrow(sim_data)){
    ind = sim_data[k,] < 40
    sim_rate[k] = coef(lm(sim_data[k,ind] ~ xs[ind]))[2]
  }
  return(sim_rate)
}


# estimate clearance time from simulated Ct data - defined as first day to reach 40
compute_clearance_time = function(sim_data,xs){
  obs_clearance = xs[apply(sim_data, 1, function(x) which(x==40)[1])]
  cens = as.numeric(!is.na(obs_clearance))
  obs_clearance[cens==0] = max(xs)
  out = data.frame(t=obs_clearance,s=cens)
  return(out)
}

# wrapper function for computing test statistics for a set of sample sizes, effect sizes and follow-times
compute_pvalues = function(Effect_sizes, Sample_sizes, Follow_up_days,
                           Nsim, summary_function, thetas, endpoint){

  power_matrix = array(dim = c(length(Effect_sizes),length(Sample_sizes),Nsim))
  tic()
  for(i in 1:length(Effect_sizes)){
    effect_size = Effect_sizes[i]
    writeLines(sprintf('Doing for effect size %s',effect_size))

    for(j in 1:length(Sample_sizes)){
      Nsample = Sample_sizes[j]
      ## Loop over simulated trials
      power_matrix[i,j,] = foreach(n = 1:Nsim, .combine = c) %dopar% {
        # simulate data
        sim_dat_cont = simulate_data(N = Nsample, effect = 1, thetas = thetas, xs = Follow_up_days)
        sim_dat_drug = simulate_data(N = Nsample, effect = effect_size, thetas = thetas, xs = Follow_up_days)

        # compute summaries
        sum_control = summary_function(sim_dat_cont, Follow_up_days)
        sum_drug = summary_function(sim_dat_drug, Follow_up_days)

        if(endpoint == 'rate'){
          # compute p-value
          pval = t.test(sum_control, sum_drug)$p.value
        }
        if(endpoint == 'timetoevent'){
          ct_data = rbind(sum_control,sum_drug)
          ct_data$intervention = c(rep(0,nrow(sum_control)),rep(1,nrow(sum_drug)))
          sdf = survdiff(Surv(t, s) ~ intervention, data = ct_data)
          pval = 1 - pchisq(sdf$chisq, length(sdf$n)-1)
        }
        pval
      }
    }
  }
  toc()
  return(power_matrix)
}

# run a single trial with adaptive randomisation
run_trial = function(Trt_max, Nmax, Nbatch,
                     effects, max_follow, init_probs,
                     success_stop_prob = 0.95,
                     futility_stop_prob = 0.2,
                     minimum_effect = 0.01){

  if(length(init_probs) != Trt_max) stop('init_probs needs to be of length = Trt_max')

  # formula for stan_lmer regression model
  reg_form = as.formula('y ~ 1 + t*Trt + t - Trt + (1 + t | id)')

  # Load model parameters used to simulate data
  load('Rout/stan_out_rate.RData')

  # Set simulation parameters
  niter = 1000
  nchains = 1
  Nbatch_1 = 50 # burn-in
  p_id = 1
  current_probs = init_probs

  # simulate first batch
  out = simulate_batch(Trt_max = Trt_max, Nbatch = Nbatch_1, effects = effects,
                       max_follow = max_follow, probs = current_probs,
                       thetas = thetas)
  sim_dat = data.frame(t = rep(0:max_follow, Nbatch_1),
                       Trt = as.factor(as.vector(sapply(out$rand_ind, rep, max_follow+1))),
                       id = as.vector(sapply(p_id:(p_id+Nbatch_1-1), rep, max_follow+1 )),
                       y = out$ct_vals)
  sim_dat = sim_dat[sim_dat$y<40, ]
  p_id = p_id + Nbatch_1

  # Fit model to data accrued so far
  interim_fit = stan_lmer(reg_form, chains = nchains, iter = niter,
                          data = sim_dat, algorithm = 'sampling')
  # get parameter estimates for the trt effects
  xx=as.matrix(interim_fit, pars = names(interim_fit$coefficients)[grep('t:Trt',names(interim_fit$coefficients))])

  # compute randomisation probabilities
  current_probs[2:Trt_max] = compute_probs(xx, minimum_effect, Trt_max)
  reached_stopping = any(c(compute_succes(xx, success_stop_prob, minimum_effect),
                           compute_futility(xx, futility_stop_prob, minimum_effect)))

  # keep running with interim analyses every Nbatch
  while(p_id < Nmax & !reached_stopping){

    out = simulate_batch(Trt_max = Trt_max, Nbatch = Nbatch, effects = effects,
                         max_follow = max_follow, probs = current_probs,
                         thetas = thetas)
    sim_dat_new = data.frame(t = rep(0:max_follow, Nbatch),
                             Trt = as.factor(as.vector(sapply(out$rand_ind, rep, max_follow+1))),
                             id = as.vector(sapply(p_id:(p_id+Nbatch-1), rep, max_follow+1 )),
                             y = out$ct_vals)
    sim_dat = rbind(sim_dat, sim_dat_new)
    sim_dat = sim_dat[sim_dat$y<40, ] # avoid having to truncate

    p_id = p_id + Nbatch

    interim_fit = stan_lmer(reg_form, chains = nchains,
                            data = sim_dat, iter = niter,
                            algorithm = 'sampling')
    xx=as.matrix(interim_fit, pars = names(interim_fit$coefficients)[grep('t:Trt',names(interim_fit$coefficients))])
    current_probs[2:Trt_max] = compute_probs(xx, minimum_effect, Trt_max)
    reached_stopping = any(c(compute_succes(xx, success_stop_prob, minimum_effect),
                             compute_futility(xx, futility_stop_prob, minimum_effect)))
  }

  print(sum(is.na(sim_dat)))
  if(sum(is.na(sim_dat))>0){
    print(sim_dat)
  }
  final_out = list(interim_fit=interim_fit, sim_dat=sim_dat)
  return(final_out)
}

# simulate a single patient batch
simulate_batch = function(Trt_max, Nbatch, effects, max_follow, probs, thetas){
  if(!(length(effects)==Trt_max)) stop('effects has to be of length Trt_max')
  rand_ind = sample(x = 1:Trt_max, size = Nbatch, replace = T, prob = probs)
  ct_vals = as.vector(t(simulate_data(N = Nbatch, effect = effects[rand_ind], thetas = thetas, xs = 0:max_follow)))
  out = list(rand_ind=rand_ind, ct_vals=ct_vals)
  return(out)
}

# adaptive randomization probabilities
compute_probs = function(xx, minimum_effect, Trt_max){
  pps = apply(xx, 2, function(x) mean(x > minimum_effect))
  pps = (1-1/Trt_max) * pps/sum(pps)
  return(pps)
}

# declare success
compute_succes = function(xx, success_stop_prob, minimum_effect){
  # compute probability that each arm has effect greater than the minimum specified effect
  pps = apply(xx, 2, function(x) mean(x > minimum_effect))
  return(any(pps>success_stop_prob))
}

# declare futility
compute_futility = function(xx, futility_stop_prob, minimum_effect){
  # compute probability that each arm has effect less than the minimum specified effect
  pps = apply(xx, 2, function(x) mean(x <= minimum_effect))
  return(all(pps < futility_stop_prob))
}

# extract trial summaries
summary_trial = function(trial_out, Trt_max, success_stop_prob, futility_stop_prob, minimum_effect){

  sim_data = trial_out$sim_dat
  N_patients = length(unique((sim_data$id)))
  trt_numbers = table(sim_data$Trt[!duplicated(sim_data$id)])
  xx=as.matrix(trial_out$interim_fit, pars = names(trial_out$interim_fit$coefficients)[grep('t:Trt',names(trial_out$interim_fit$coefficients))])

  winner = 0
  if(compute_succes(xx, success_stop_prob, minimum_effect)) {
    winner = which.max(apply(xx, 2, function(x) mean(x>minimum_effect)))
  }
  if(compute_futility(xx, futility_stop_prob, minimum_effect)) {
    winner = -1
  }

  prob_superior = apply(xx, 2, function(x) mean(x > minimum_effect))
  prob_futile = apply(xx, 2, function(x) mean(x < minimum_effect))

  res = c(N_patients, trt_numbers, winner, prob_superior, prob_futile)

  return(res)
}

