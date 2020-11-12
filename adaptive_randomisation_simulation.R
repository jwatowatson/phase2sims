library(doParallel)
library(tictoc)
library(rstanarm)
registerDoParallel(cores = 8)
source('functions.R')
# simulate an adaptive trial with 1 drug that works, effect size of 10%, others same as control
Tmax = 5
effects = c(rep(1,Tmax-1), 1.1)
max_follow = 10
Ntrials = 500
Nbatch = 20
Nmax = 400
init_probs = rep(1/Tmax, Tmax)
my_success_threshold = 0.9975

tic()
sims = foreach(i = 1:Ntrials, .combine = rbind) %dopar% {
  out = run_trial(Tmax = Tmax,
                  Nmax = Nmax,
                  Nbatch = Nbatch,
                  effects = effects,
                  max_follow = max_follow,
                  init_probs = init_probs,
                  success_threshold = my_success_threshold)
  my_summary = summary_trial(out, Tmax)
  my_summary
}
toc()

save(sims, file = 'adaptive_trial_simulation1.RData')

######################################################
# simulate an adaptive trial with 0 drugs that work
effects = rep(1,Tmax)

tic()
sims = foreach(i = 1:Ntrials, .combine = rbind) %dopar% {
  out = run_trial(Tmax = Tmax,
                  Nmax = Nmax,
                  Nbatch = Nbatch,
                  effects = effects,
                  max_follow = max_follow,
                  init_probs = init_probs,
                  success_threshold = my_success_threshold)
  my_summary = summary_trial(out, Tmax)
  my_summary
}
toc()

save(sims, file = 'adaptive_trial_simulation2.RData')

#### Plot results
load('adaptive_trial_simulation1.RData')
sims = as.data.frame(sims)
table(sims$`t:Trt4`)
pdf('Adaptive_trial1_simulation.pdf', height = 7, width = 10)
par(las =1 , mfrow=c(1,2), family='serif')
hist(sims$Ntotal, breaks = seq(0,400,by=10), freq = T, xlab = 'Total trial size',
     ylab = 'Number of trials', xaxt='n', main = '')
axis(1, at = seq(40,400, by = 80))

hist(100*sims[,5]/sims$Ntotal, main = '', xlab = "Randomised to best intervention (%)",ylab = 'Number of trials')
dev.off()
