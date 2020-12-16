#*** This script does the power simulations
#* outputs Figures 1 - 3 in manuscript

library(doParallel)
library(tictoc)
library(rstan)
library(survival)
doParallel::registerDoParallel(cores = 8)
FORCE_RE_RUN = F

source('functions.R')
# get posterior samples from the pharmacodynamic model
load('temp_stan_rate.out')
thetas = extract(ct_fit_rate)

# Parameters for the simulations
Sample_sizes = seq(20,100,by=10)
Effect_sizes = c(1.05,1.075,1.1,1.15,1.2,1.25)
Max_follows = c(7,10,14) # Follow-up durations
Nsim = 20000 # number of trials
my_alpha = 0.05 # significance level to compute type 2 error

#####*********************************************
#####* Compare power for a bunch of sample sizes, effect sizes and differing follow-up durations
#####*

if(FORCE_RE_RUN){
  power_summaries = list(); j = 1
  for(i in 1:length(Max_follows)){
    pp1 = compute_pvalues(Effect_sizes = Effect_sizes,
                          Sample_sizes = Sample_sizes,
                          Follow_up_days = 0:Max_follows[i],
                          Nsim = Nsim,
                          summary_function = compute_clearance_time,
                          thetas=thetas, endpoint = 'timetoevent')
    pp2 = compute_pvalues(Effect_sizes = Effect_sizes,
                          Sample_sizes = Sample_sizes,
                          Follow_up_days = 0:Max_follows[i],
                          Nsim = Nsim,
                          summary_function = compute_rate,
                          thetas=thetas, endpoint = 'rate')
    power_summaries[[j]] = pp1; j = j+1;
    power_summaries[[j]] = pp2; j = j+1;
  }
  save(power_summaries, file = 'power_summaries.RData')
} else {
  load(file = 'power_summaries.RData')
}

# compute average 1 - type 2 error
power_list = lapply(power_summaries, function(ll) apply(ll, 1:2, function(x) mean(x<my_alpha)))

cols=RColorBrewer::brewer.pal(7,'Dark2')[c(4,1)]
my_cols = rep(cols, length(Max_follows))
my_ltys = foreach(x = 1:length(Max_follows), .combine = c) %do% rep(x,2)
pdf('Plots/Rate_better_timetoevent.pdf', width = 10, height = 8)
par(mfrow=c(2,3),las=1,lwd=1.5, cex.lab=1.3,cex.axis=1.3,family='serif',mar=c(5,5,3,2),bty='n')
for(i in 1:length(Effect_sizes)) {
  plot(Sample_sizes, power_list[[1]][i,], ylim = c(0,1),type='l',
       xlab='Sample size per arm',ylab='Power',col=my_cols[1], lty=my_ltys[1])
  for(j in 2:length(power_list)){
    lines(Sample_sizes, power_list[[j]][i,], col=my_cols[j], lty=my_ltys[j])
  }
  title(paste(100*(Effect_sizes[i]-1),'%',sep=''))
  if(i==1){
    legend('topleft',col=c(cols, rep('black',length(Max_follows))),bty='n',
           lty=c(1,1,1:length(Max_follows)),inset=0.03,title = 'Endpoint',
           legend = c('Time to clearance','Rate of clearance'))
    legend('left',col='black',title = 'Duration of follow-up (days)',
           lty=1:length(Max_follows),inset=0.03,bty='n',
           legend = Max_follows)
  }
}

dev.off()


#####*********************************************
#####* Twice daily versus once daily for rate endpoints

Effect_sizes = 1.075
if(FORCE_RE_RUN){
  power_Rates_10_once = compute_pvalues(Effect_sizes = Effect_sizes,
                                        Sample_sizes = Sample_sizes,
                                        Follow_up_days = 0:10,
                                        Nsim = Nsim,
                                        summary_function = compute_rate,
                                        thetas=thetas, endpoint = 'rate')[1,,]
  power_Rates_10_twice = compute_pvalues(Effect_sizes = Effect_sizes,
                                         Sample_sizes = Sample_sizes,
                                         Follow_up_days = seq(0,10,by=0.5),
                                         Nsim = Nsim,
                                         summary_function = compute_rate,
                                         thetas=thetas, endpoint = 'rate')[1,,]
  save(power_Rates_10_twice, power_Rates_10_once, file = 'once_vs_twice.RData')
} else {
  load(file = 'once_vs_twice.RData')
}

pdf('Plots/Once_versus_twice.pdf', width = 8, height = 8)
par(mfrow=c(1,1),las=1,lwd=1.5, cex.lab=1.3,cex.axis=1.3,family='serif',mar=c(5,5,3,2),bty='n')
plot(Sample_sizes, apply(power_Rates_10_once,1,function(x) mean(x<my_alpha)),type='l',
     ylab='Power', xlab= 'Sample size per arm', ylim=c(0,1))
lines(Sample_sizes, apply(power_Rates_10_twice,1,function(x) mean(x<my_alpha)),lty=2)
legend('bottomright',legend = c('Once daily','Twice daily'),lty=1:2,inset=0.03,title = 'qPCR sampling')
dev.off()

#####*********************************************
#####* Fix effect and sample size and look at how duration influences power

Effect_sizes = 1.075
Sample_sizes = 50
Max_follow = seq(3, 21, by=1)
if(FORCE_RE_RUN){
  power_TE_follow = power_rates_follow = array(dim = c(length(Max_follow), Nsim))
  for(mm in 1:length(Max_follow)){
    power_rates_follow[mm,]=compute_pvalues(Effect_sizes = Effect_sizes,
                                            Sample_sizes = Sample_sizes,
                                            Follow_up_days = 0:Max_follow[mm],
                                            Nsim = Nsim,
                                            summary_function = compute_rate,
                                            thetas=thetas, endpoint = 'rate')[1,1,,drop=T]
  }
  for(mm in 5:length(Max_follow)){
    power_TE_follow[mm,]=compute_pvalues(Effect_sizes = Effect_sizes,
                                         Sample_sizes = Sample_sizes,
                                         Follow_up_days = 0:Max_follow[mm],
                                         Nsim = Nsim,
                                         summary_function = compute_clearance_time,
                                         thetas=thetas, endpoint = 'timetoevent')[1,1,,drop=T]
  }
  save(power_rates_follow, power_TE_follow, file = 'duration.RData')
} else {
  load(file = 'duration.RData')
}


pdf('Plots/Duration_follow_up.pdf', width = 8, height = 8)
par(mfrow=c(1,1),las=1,lwd=1.5, cex.lab=1.3,cex.axis=1.3,family='serif',mar=c(5,5,3,2),bty='n')
plot(Max_follow, apply(power_rates_follow, 1, function(x) mean(x<my_alpha)),
     type='l', ylim = c(0,1), ylab='Power',panel.first=grid(),
     xlab= 'Duration of follow-up',col=cols[2])
lines(Max_follow,apply(power_TE_follow, 1, function(x) mean(x<my_alpha)),col=cols[1])
abline(h=c(0.6,.44),v=c(10,12),lty=2)

legend('topright',legend = c('Time-to-clearance', 'Rate-of-clearance'),
       col=cols,lwd=2, title = 'Endpoint', bty='n')
dev.off()




### Figure to demonstrate what the model looks like
sim_dat_null = simulate_data(N = 10^4, effect = 1, thetas = thetas, xs = 0:21)
sim_dat_5 = simulate_data(N = 10^4, effect = 1.05, thetas = thetas, xs = 0:21)
sim_dat_75 = simulate_data(N = 10^4, effect = 1.075, thetas = thetas, xs = 0:21)
sim_dat_15 = simulate_data(N = 10^4, effect = 1.15, thetas = thetas, xs = 0:21)
sim_dat_25 = simulate_data(N = 10^4, effect = 1.25, thetas = thetas, xs = 0:21)

pdf('Plots/Example.pdf', width = 10, height = 10)
par(las = 1, mfrow=c(2,2), family='serif', bty='n',cex.axis=1.3,cex.lab=1.3)
set.seed(485)
plot(0:21, sim_dat_null[1,], col = adjustcolor('grey',.3),type='l',ylim = c(40,13),
     ylab='Ct',xlab = 'Days since peak viral load')
points(0:21, sim_dat_null[1,], pch='.')
cols = RColorBrewer::brewer.pal(n = 11, name = 'Set3'); j = 1
for(i in sample(0:1000, 9)){
  lines(0:21, sim_dat_null[i,], col = adjustcolor(cols[j%%11 + 1],.4))
  points(0:21, sim_dat_null[i,], pch='.')
  j = j+1
}


cols = rev(RColorBrewer::brewer.pal(n = 5, name = 'RdYlBu'))
plot(0:21, colMeans(sim_dat_null), type='l', ylab='Ct', col=cols[1],
     xlab = 'Days since peak viral load', ylim = c(40,13),lwd=2)
lines(0:21, colMeans(sim_dat_5),col=cols[2],lwd=2)
lines(0:21, colMeans(sim_dat_75),col=cols[3],lwd=2)
lines(0:21, colMeans(sim_dat_15),col=cols[4],lwd=2)
lines(0:21, colMeans(sim_dat_25),col=cols[5],lwd=2)
legend('topright', legend = c(0,5,7.5,15,25),col=cols, bty='n',
       inset=0.03,lwd=3, title = 'Effect size (%)')


rates_15 = compute_rate(sim_data = sim_dat_15, xs = 0:21)
rates_0 = compute_rate(sim_data = sim_dat_null, xs = 0:21)

plot(density(rates_0), main = '', ylab = '', col = cols[1],lwd=3,
     yaxt='n', xaxt='n', xlab='Clearance half-life', xlim = c(3,1))
lines(density(rates_15), col=cols[4], lwd=3)
axis(1, at = seq(1,3,by=.5), labels = 24/seq(1,3,by=.5))
legend('topright', legend = c(0,15),col=cols[c(1,4)], bty='n',
       inset=0.03,lwd=3, title = 'Effect size (%)')

TC_0 = compute_clearance_time(sim_dat_null,xs = 0:21)[,1]
TC_15 = compute_clearance_time(sim_dat_15,xs = 0:21)[,1]
plot(density(TC_15[complete.cases(TC_15)], bw = .6),
     main = '', ylab = '', col = cols[4],lwd=3,
     yaxt='n',  xlab='Clearance time')
lines(density(TC_0[complete.cases(TC_0)],bw = 0.6), col=cols[1], lwd=3)
legend('topright', legend = c(0,15),col=cols[c(1,4)], bty='n',
       inset=0.03,lwd=3, title = 'Effect size (%)')
dev.off()
