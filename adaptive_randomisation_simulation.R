library(doParallel)
library(tictoc)
library(rstanarm)
registerDoParallel(cores = 8)
FORCE_RERUN = T
source('functions.R')
# simulate an adaptive trial with 1 drug that works, effect size of 10%, others same as control
Trt_max = 5
effects = c(rep(1,Trt_max-1), 1.1)
max_follow = 10  # follow-up in days
Ntrials = 200  # number of trials to simulate
Nbatch = 20 # number of patients per batch
Nmax = 400 # max number of patients
init_probs = rep(1/Trt_max, Trt_max)
success_stop_prob = 0.99
futility_stop_prob = 0.1
minimum_effect = 0.01

sim_path = 'adaptive_trial_simulation1.RData'

if(FORCE_RERUN | !file.exists(sim_path)){
  tic()
  sims = foreach(i = 1:Ntrials, .combine = rbind) %dopar% {
    out = run_trial(Trt_max = Trt_max,
                    Nmax = Nmax,
                    Nbatch = Nbatch,
                    effects = effects,
                    max_follow = max_follow,
                    init_probs = init_probs,
                    success_stop_prob = success_stop_prob,
                    futility_stop_prob = futility_stop_prob,
                    minimum_effect = minimum_effect)
    my_summary = summary_trial(out, Trt_max, success_stop_prob, futility_stop_prob, minimum_effect)
    my_summary
  }
  toc()
  save(sims, file = sim_path)
}


######################################################
# simulate an adaptive trial with 0 drugs that work
effects = rep(1,Trt_max)
sim_path = 'adaptive_trial_simulation2.RData'

Ntrials = 20  # number of trials to simulate
Nbatch = 40 # number of patients per batch
Nmax = 200 # max number of patients

if(FORCE_RERUN | !file.exists(sim_path)){
  tic()
  sims = foreach(i = 1:Ntrials, .combine = rbind) %dopar% {
    out = run_trial(Trt_max = Trt_max,
                    Nmax = 200,
                    Nbatch = Nbatch,
                    effects = effects,
                    max_follow = max_follow,
                    init_probs = init_probs,
                    success_stop_prob = success_stop_prob,
                    futility_stop_prob = futility_stop_prob,
                    minimum_effect = minimum_effect)
    my_summary = summary_trial(out, Trt_max,success_stop_prob, futility_stop_prob, minimum_effect)
    my_summary
  }
  toc()
  save(sims, file = sim_path)
}


#### Plot results
load('adaptive_trial_simulation1.RData')
sims = as.data.frame(sims)
100*table(sims$winner)/nrow(sims)

pdf('Adaptive_trial1_simulation.pdf', height = 7, width = 10)
par(las =1 , mfrow=c(1,2), family='serif')
hist(sims$Ntotal, breaks = seq(50,410,by=20), freq = T, xlab = 'Total trial size',
     ylab = 'Number of trials', xaxt='n', main = '')
axis(1, at = seq(50,400, by = 50))

hist(100*sims$t5/sims$Ntotal, main = '', xlab = "Randomised to active drug (%)",
     ylab = 'Number of trials', breaks = seq(0,70, by=5))
abline(v = 100/5, col='red',lwd=2)

# hist(sims$prob_superior)
dev.off()



sims = as.data.frame(sims)
100*table(sims$winner)/nrow(sims)

pdf('Adaptive_trial1_simulation.pdf', height = 7, width = 10)
par(las =1 , mfrow=c(1,2), family='serif')
hist(sims$Ntotal, breaks = seq(50,410,by=20), freq = T, xlab = 'Total trial size',
     ylab = 'Number of trials', xaxt='n', main = '')
axis(1, at = seq(50,400, by = 50))

hist(100*sims$t5/sims$Ntotal, main = '', xlab = "Randomised to active drug (%)",
     ylab = 'Number of trials')
abline(v = 100/5, col='red',lwd=2)

# hist(sims$prob_superior)
