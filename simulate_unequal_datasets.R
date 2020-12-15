source('functions.R')
# Load model parameters used to simulate data
load('stan_out_rate.RData')
Nsims = 100000
N = 250
time_pvals = time_ct_diff = list()
K = 15
ts = rep(0:K, N)
pvals = mydiff = array(dim = Nsims)

eps=.3
i_dat = 1
for(i in 1:Nsims){
  dat = simulate_data(N = N, effect = 1,xs = 0:K,thetas = thetas)
  start = truncdist::rtrunc(N, 'pois',lambda = 1.5, a = 0, b = 3)
  for(j in 1:N){
    dat[j, 1:(K+1-start[j]+1)] = dat[j, start[j]:(K+1)]
    if(start[j]>1) dat[j, (K-start[j]+2):(K+1)] = NA
  }
  dat = dat[, 1:(11)]
  trt = sample(1:2, size = N, replace = T)
  mydiff[i] = mean(dat[trt==1,1]) - mean(dat[trt==2,1])
  pvals[i] = t.test(dat[trt==1,1] , dat[trt==2,1])$p.value

  if(abs(mydiff[i]) > eps & pvals[i] > 0.05 ){
    if(sign(mydiff[i]) == 1){
      time_ct_diff[[i_dat]] = colMeans(dat[trt==1, ]) - colMeans(dat[trt==2, ])
    } else {
      time_ct_diff[[i_dat]] = colMeans(dat[trt==2, ]) - colMeans(dat[trt==1, ])
    }
    time_pvals[[i_dat]] = apply(dat, 2, function(x, trt) t.test(x[trt==1] , x[trt==2])$p.value, trt=trt)
    i_dat = i_dat+1
  }
}
hist(mydiff); abline(v=c(-1,1) * eps,lwd=3, col='red')


xx = sapply(time_ct_diff, rbind)
yy = (sapply(time_pvals, rbind))
plot(1:11, rowMeans(xx), ylim = range(xx), type='l', lwd=3)
for(i in 1:ncol(xx)){
  lines(1:11, xx[,i])
}

par(las=1)
plot(1:11, 100*apply(yy, 1, function(x) mean(x<0.05)), type='l', lwd=3,
     xlab = 'Time from enrollment',
     ylab = 'Probability of seeing statistically significant difference (%)')



