library(doParallel)
library(tictoc)
library(rstanarm)
registerDoParallel(cores = 8)
FORCE_RERUN = F
source('functions.R')
# simulate an adaptive trial with 1 drug that works, effect size of 10%, others same as control
Trt_max = 5
effects = c(rep(1,Trt_max-1), 1.1)
max_follow = 10  # follow-up in days
Ntrials = 1000  # number of trials to simulate
Nbatch = 30 # number of patients per batch
Nmax = 400 # max number of patients
init_probs = rep(1/Trt_max, Trt_max)
success_stop_prob = 0.99
futility_stop_prob = 0.1
minimum_effect = 0.01

sim_path = 'Rout/adaptive_trial_simulation1.txt'

if(FORCE_RERUN | !file.exists(sim_path)){
  tic()
  my_cols = c("Ntotal", paste('t', 1:Trt_max, sep = ''), 'winner',
              paste('prob_superior', 2:Trt_max, sep = ''),
              paste('prob_futile', 2:Trt_max, sep = ''))
  write.table(t(my_cols), sim_path, row.names = F, col.names = F, quote = F)
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
    # write results to file
    write.table(t(my_summary), sim_path, append = T, row.names = F, col.names = F, quote = F)
    my_summary
  }
  toc()
}


######################################################
# simulate an adaptive trial with 0 drugs that work
effects = rep(1,Trt_max)
sim_path = 'Rout/adaptive_trial_simulation2.txt'

# Make this smaller - otherwise takes aages
Ntrials = 100  # number of trials to simulate
Append_results = F
if(FORCE_RERUN | !file.exists(sim_path) | Append_results){
  tic()
  if(!Append_results){
    my_cols = c("Ntotal", paste('t', 1:Trt_max, sep = ''), 'winner',
                paste('prob_superior', 2:Trt_max, sep = ''),
                paste('prob_futile', 2:Trt_max, sep = ''))
    write.table(t(my_cols), sim_path, row.names = F, col.names = F, quote = F)
  }

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
    # write results to file
    write.table(t(my_summary), sim_path, append = T, row.names = F, col.names = F, quote = F)
    my_summary
  }
  toc()
}


#### Plot results
sims1 = read.table('Rout/adaptive_trial_simulation1.txt', header = F,skip = 1,
                   col.names = c("Ntotal", paste('t', 1:Trt_max, sep = ''), 'winner',
                                 paste('prob_superior', 2:Trt_max, sep = ''),
                                 paste('prob_futile', 2:Trt_max, sep = '')), fill = T)
sims2 = read.table('Rout/adaptive_trial_simulation2.txt', header = F,skip = 1,
                   col.names = c("Ntotal", paste('t', 1:Trt_max, sep = ''), 'winner',
                                 paste('prob_superior', 2:Trt_max, sep = ''),
                                 paste('prob_futile', 2:Trt_max, sep = '')), fill = T)


# sims1 give us an estimate of type 2 error for this effect size
type_2_error = 100*sum(sims1$winner != 4)/nrow(sims1)
type_1_error = 100*sum(sims2$winner > 0)/nrow(sims2)

round(100*table(sims1$winner)/nrow(sims1))
round(100*table(sims2$winner)/nrow(sims2))

apply(sims1[, paste('t', 1:Trt_max, sep = '')], 2, median)
apply(sims2[, paste('t', 1:Trt_max, sep = '')], 2, median)

pdf('Plots/Adaptive_trial_simulation.pdf', height = 7, width = 10)
par(las = 1 , mfrow=c(2,2), family='serif')
hist(sims1$Ntotal, breaks = seq(50,410,by=30), freq = T, xlab = 'Total trial size',
     ylab = 'Number of trials', xaxt='n', main = '')
axis(1, at = seq(50,400, by = 50))

hist(100*sims1$t5/sims1$Ntotal, main = '', xlab = "Randomised to active drug (%)",
     ylab = 'Number of trials', breaks = seq(0,70, by=5))
abline(v = 100/5, col='red',lwd=2)

hist(sims2$Ntotal, breaks = seq(50,200,by=10), freq = T, xlab = 'Total trial size',
     ylab = 'Number of trials', xaxt='n', main = '')
axis(1, at = seq(50,300, by = 50))

barplot(table(sims2$winner), main = '', xlab = "",
        ylab = 'Number of trials',
        names.arg = c('Futile','None',paste('T', 1:(Trt_max-1), sep = '')))

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
