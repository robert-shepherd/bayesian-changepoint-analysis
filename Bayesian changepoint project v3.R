#Set seed so results can be reproduced
set.seed(123)

#Reading in data
library(boot)
data(coal)
x <- 1851:1962
y <- tabulate(floor(coal[[1]]))
y <- y[x]


#Set up
n <- length(y)
m <- 10000
lambda1 <- lambda2 <- T1 <- numeric(m)
L <- numeric(n)
T1[1] <- sample(1:n,1)
lambda1[1] <- 1
lambda2[1] <- 1
b1 <- 1
b2 <- 1

#Gibbs sampler
for (i in 2:m) {
  T1t <- T1[i-1]
  
  #generate lambda1
  r <- 2 + sum(y[1:T1t])
  lambda1[i] <- rgamma(1,shape = r, rate = T1t + b1)
  
  #generate lambda2
  if (T1t + 1 > n) r <- 2 + sum(y) else
    r <- 2 + sum(y[(T1t+1):n])
  lambda2[i] <- rgamma(1,shape = r, rate = n-T1t+b2)
  
  #generate b1 and b2
  b1 <- rgamma(1,shape = 2, rate = lambda1[i]+1)
  b2 <- rgamma(1,shape = 2, rate = lambda2[i]+1)
  
  for (j in 1:n) {
    L[j] <- exp((lambda2[i] - lambda1[i])* j) *
      (lambda1[i]/lambda2[i])^sum(y[1:j])
  }
  L <- L / sum(L)
  
  #generate T from discrete distribution L on 1:n
  T1[i] <- sample(1:n, prob=L,size = 1)
}

#Checking convergence
par(mfrow=c(3,1))
plot(T1, type="l")
plot(lambda1,type="l")
plot(lambda2,type="l")

require(coda)
summary(mcmc(T1))
summary(mcmc(lambda1))
summary(mcmc(lambda2))

#Creating a matrix
gibbs_sum <- matrix(c(T1,lambda1,lambda2),nrow=10000,ncol=3)
#head(gibbs_sum)

#Removing burn in 
b <- 201
gibbs_sum_no_burn <- gibbs_sum[b:m,]
#dim(gibbs_sum_no_burn)

#Check short term correlation
par(mfrow=c(1,1))
acf(T1,lag.max = 20)
acf(lambda1,lag.max = 20)
acf(lambda2,lag.max = 20)

#Taking every second draw
gibbs_sample <- gibbs_sum_no_burn[seq(0,9800,2),]
#dim(gibbs_sample)
#head(gibbs_sample)

print(mean(gibbs_sample[,1]))
#39.83
print(mean(gibbs_sample[,2]))
#3.144551
print(mean(gibbs_sample[,3]))
#0.9396872

#Plotting histograms
par(mfrow=c(3,1))
hist(gibbs_sample[,1],main="Posterior distribution of first changepoint",xlab="T1")
hist(gibbs_sample[,2],main="Posterior distribution of the rate of diasters prior to the change",xlab="lambda 1")
hist(gibbs_sample[,3],main="Posterior distribution of the rate of diasters after the change",xlab="lambda 2")

#Plotting predictions
pred <- matrix(0,nrow=4900,ncol=112)
for (sample in 1:4900) {
  for (year in 1:112) {
    if (year < gibbs_sample[sample,1]) {
      pred[sample, year] <- rpois(1,gibbs_sample[sample,2]) 
    } else {
      pred[sample,year] <- rpois(1,gibbs_sample[sample,3])
    }
  }
}

cp <- t(apply(pred,1,cumsum))
dim(cp)

median <- rep(0,112) 
lower <- rep(0,112)
upper <- rep(0,112)


for (year in 1:112) {
  median[year] <- quantile(cp[,year],probs <- c(0.5))
  lower[year] <- quantile(cp[,year],probs <- c(0.025))
  upper[year] <- quantile(cp[,year],probs <- c(0.975))
}

#Plot predictions
par(mfrow=c(1,1))
plot(1851:1962,cumsum(y),type="l"
     ,main="One changepoint predictions vs. observed"
     ,ylim=c(0,250)
     ,xlab="year")
lines(1851:1962,median,col="red")
lines(1851:1962,lower,col="grey")
lines(1851:1962,upper,col="grey")
abline(v=1890)

#Changepoint year
changepoint <- round(mean(gibbs_sample[,1],0))

changepoint
#40

1850 + changepoint
#1890

#Checking RMSE
pred_median <- rep(0,112) 

for (year in 1:112) {
  pred_median[year] <- quantile(pred[,year],probs <- c(0.5))
}

rmse_pre <- (sum((y[1:(changepoint-1)]-pred_median[1:(changepoint-1)])^2)/length(y[1:(changepoint-1)]))^0.5
rmse_post <- (sum((y[changepoint:length(y)]-pred_median[changepoint:length(y)])^2)/length(y[changepoint:length(y)]))^0.5
rmse_total <- (sum((y-pred_median)^2)/n)^0.5

rmse_pre
#1.577079

rmse_post
#1.013606

rmse_total
#1.239239

#########
#Using fixed constants to generate multiple touchpoints
#########

#Set up

#Pre changepoint

y_pre <- y[1:(changepoint-1)]
n_pre <- length(y_pre)
m <- 10000
lambda1_pre <- lambda2_pre <- T2_pre <- numeric(m)
L_pre <- numeric(n_pre)
T2_pre[1] <- sample(1:(changepoint-1),1)
lambda1_pre[1] <- 1
lambda2_pre[1] <- 1
b1_pre <- 1
b2_pre <- 1

#Gibbs sampler for pre
for (i in 2:m) {
  T2t_pre <- if (T2_pre[i-1] == n_pre) n_pre-1 else T2_pre[i-1]
  
  #generate lambda1
  r <- 2 + sum(y_pre[1:T2t_pre])
  lambda1_pre[i] <- rgamma(1,shape = r, rate = T2t_pre + b1_pre)
  
  #generate lambda2
  if (T2t_pre + 1 > n_pre) r <- 2 + sum(y_pre) else
    r <- 2 + sum(y[(T2t_pre+1):n_pre])
  lambda2_pre[i] <- rgamma(1,shape = r, rate = n_pre-T2t_pre+b2_pre)
  
  #generate b1 and b2
  b1_pre <- rgamma(1,shape = 2, rate = lambda1_pre[i]+1)
  b2_pre <- rgamma(1,shape = 2, rate = lambda2_pre[i]+1)
  
  for (j in 1:n_pre) {
    L_pre[j] <- exp((lambda2_pre[i] - lambda1_pre[i])* j) *
      (lambda1_pre[i]/lambda2_pre[i])^sum(y_pre[1:j])
  }
  L_pre <- L_pre / sum(L_pre)
  
  #generate T from discrete distribution L on 1:n
  T2_pre[i] <- sample(1:n_pre, prob=L_pre,size = 1)
}

#Checking convergence
plot(T2_pre, type="l")
plot(lambda1_pre,type="l")
plot(lambda2_pre,type="l")

require(coda)
summary(mcmc(T2_pre))
summary(mcmc(lambda1_pre))
summary(mcmc(lambda2_pre))

#Creating a matrix
gibbs_sum_pre <- matrix(c(T2_pre,lambda1_pre,lambda2_pre),nrow=10000,ncol=3)
#head(gibbs_sum_pre)

#Removing burn in 
b <- 201
gibbs_sum_no_burn_pre <- gibbs_sum_pre[b:m,]
#dim(gibbs_sum_no_burn_pre)

#Check short term correlation
acf(T2_pre,lag.max = 20)
acf(lambda1_pre,lag.max = 20)
acf(lambda2_pre,lag.max = 20)

#Taking every fourth draw
gibbs_sample_pre <- gibbs_sum_no_burn_pre[seq(0,9800,2),]
#dim(gibbs_sample)
#head(gibbs_sample)

print(mean(gibbs_sample_pre[,1]))
#21.35
print(mean(gibbs_sample_pre[,2]))
#3.394308
print(mean(gibbs_sample_pre[,3]))
#2.982563

length(gibbs_sample_pre) 

#Plotting predictions
pred_pre <- matrix(0,nrow=4900,ncol=n_pre)
for (sample in 1:4900) {
  for (year in 1:n_pre) {
    if (year < gibbs_sample_pre[sample,1]) {
      pred_pre[sample, year] <- rpois(1,gibbs_sample_pre[sample,2]) 
    } else {
      pred_pre[sample,year] <- rpois(1,gibbs_sample_pre[sample,3])
    }
  }
}

cp_pre <- t(apply(pred_pre,1,cumsum))
dim(cp_pre)

median_pre <- rep(0,n_pre) 
lower_pre <- rep(0,n_pre)
upper_pre <- rep(0,n_pre)




for (year in 1:n_pre) {
  median_pre[year] <- quantile(cp_pre[,year],probs <- c(0.5))
  lower_pre[year] <- quantile(cp_pre[,year],probs <- c(0.025))
  upper_pre[year] <- quantile(cp_pre[,year],probs <- c(0.975))
}


#Plot predictions
plot(1851:1889,cumsum(y_pre),type="l"
     ,main="One changepoint predictions vs. actual"
     ,ylim=c(0,250)
     ,xlab="year")
lines(1851:1889,median_pre,col="red")
lines(1851:1889,lower_pre,col="grey")
lines(1851:1889,upper_pre,col="grey")

#Changepoint year
changepoint2 <- round(mean(gibbs_sample_pre[,1],0))

changepoint2
#21

1850 + changepoint2
#1871

#Checking RMSE
pred_median_pre <- rep(0,n_pre) 

for (year in 1:n_pre) {
  pred_median_pre[year] <- quantile(pred_pre[,year],probs <- c(0.5))
}

rmse_pre_pre <- (sum((y_pre[1:(changepoint2-1)]-pred_median_pre[1:(changepoint2-1)])^2)/length(y_pre[1:(changepoint2-1)]))^0.5
rmse_post_pre <- (sum((y_pre[changepoint2:n_pre]-pred_median_pre[changepoint2:n_pre])^2)/length(y_pre[changepoint2:n_pre]))^0.5
rmse_total_pre <- (sum((y_pre-pred_median_pre)^2)/n_pre)^0.5


rmse_pre_pre
#1.788854

rmse_post_pre
#1.376494

rmse_total_pre
#1.601282

#Post changepoint
y_post <- y[changepoint:n]
n_post <- length(y_post)
m <- 10000
lambda1_post <- lambda2_post <- T2_post <- numeric(m)
L_post <- numeric(n_post)
T2_post[1] <- sample(1:n_post,1)
lambda1_post[1] <- 1
lambda2_post[1] <- 1
b1_post <- 1
b2_post <- 1

#Gibbs sampler for post
for (i in 2:m) {
  T2t_post <- if (T2_post[i-1] == n_post) n_post-1 else T2_post[i-1]
  
  #generate lambda1
  r <- 2 + sum(y_post[1:T2t_post])
  lambda1_post[i] <- rgamma(1,shape = r, rate = T2t_post + b1_post)
  
  #generate lambda2
  if (T2t_post + 1 > n_post) r <- 2 + sum(y_post) else
    r <- 2 + sum(y[(T2t_post+1):n_post])
  lambda2_post[i] <- rgamma(1,shape = r, rate = n_post-T2t_post+b2_post)
  
  #generate b1 and b2
  b1_post <- rgamma(1,shape = 2, rate = lambda1_post[i]+1)
  b2_post <- rgamma(1,shape = 2, rate = lambda2_post[i]+1)
  
  for (j in 1:n_post) {
    L_post[j] <- exp((lambda2_post[i] - lambda1_post[i])* j) *
      (lambda1_post[i]/lambda2_post[i])^sum(y_post[1:j])
  }
  L_post <- L_post / sum(L_post)
  
  #generate T from discrete distribution L on 1:n
  T2_post[i] <- sample(1:n_post, prob=L_post,size = 1)
}

#Checking convergence
plot(T2_post, type="l")
plot(lambda1_post,type="l")
plot(lambda2_post,type="l")

require(coda)
summary(mcmc(T2_post))
summary(mcmc(lambda1_post))
summary(mcmc(lambda2_post))

#Creating a matrix
gibbs_sum_post <- matrix(c(T2_post,lambda1_post,lambda2_post),nrow=10000,ncol=3)
#head(gibbs_sum_post)

#Removing burn in 
b <- 201
gibbs_sum_no_burn_post <- gibbs_sum_post[b:m,]
#dim(gibbs_sum_no_burn_post)

#Check short term correlation
acf(T2_post,lag.max = 20)
acf(lambda1_post,lag.max = 20)
acf(lambda2_post,lag.max = 20)

#Taking every second draw
gibbs_sample_post <- gibbs_sum_no_burn_post[seq(0,9800,2),]
#dim(gibbs_sample)
#head(gibbs_sample)

print(mean(gibbs_sample_post[,1]))
#54
print(mean(gibbs_sample_post[,2]))
#1.065339
print(mean(gibbs_sample_post[,3]))
#1.439215

length(gibbs_sample_post) 

#Plotting predictions
pred_post <- matrix(0,nrow=4900,ncol=n_post)
for (sample in 1:4900) {
  for (year in 1:n_post) {
    if (year < gibbs_sample_post[sample,1]) {
      pred_post[sample, year] <- rpois(1,gibbs_sample_post[sample,2]) 
    } else {
      pred_post[sample,year] <- rpois(1,gibbs_sample_post[sample,3])
    }
  }
}

cp_post <- t(apply(pred_post,1,cumsum))
dim(cp_post)

median_post <- rep(0,n_post) 
lower_post <- rep(0,n_post)
upper_post <- rep(0,n_post)



for (year in 1:n_post) {
  median_post[year] <- quantile(cp_post[,year],probs <- c(0.5))
  lower_post[year] <- quantile(cp_post[,year],probs <- c(0.025))
  upper_post[year] <- quantile(cp_post[,year],probs <- c(0.975))
}

#Plot predictions
plot(1890:1962,cumsum(y_post),type="l"
     ,main="One changepoint predictions vs. actual"
     ,ylim=c(0,250)
     ,xlab="year")
lines(1890:1962,median_post,col="red")
lines(1890:1962,lower_post,col="grey")
lines(1890:1962,upper_post,col="grey")

#Changepoint year
changepoint3 <- round(mean(gibbs_sample_post[,1],0))

changepoint3
#54

1850 + changepoint + changepoint3
#1944

#Checking RMSE
pred_median_post <- rep(0,n_post) 

for (year in 1:n_post) {
  pred_median_post[year] <- quantile(pred_post[,year],probs <- c(0.5))
}

rmse_post_pre <- (sum((y_post[1:(changepoint2-1)]-pred_median_post[1:(changepoint2-1)])^2)/length(y_post[1:(changepoint2-1)]))^0.5
rmse_post_post <- (sum((y_post[changepoint2:n_post]-pred_median_post[changepoint2:n_post])^2)/length(y_post[changepoint2:n_post]))^0.5
rmse_total_post <- (sum((y_post-pred_median_post)^2)/n_post)^0.5

rmse_post_pre
#1.024695

rmse_post_post
#1.00939

rmse_total_post
#1.013606

#Combining predictions
combined_median <- c(median_pre,median_post+123)
combined_lower <- c(lower_pre,lower_post+94)
combined_upper <- c(upper_pre,upper_post+156)

#Plot predictions
par(mfrow=c(1,1))
plot(1851:1962,cumsum(y),type="l"
     ,main="Three changepoint predictions vs. observed"
     ,ylim=c(0,250)
     ,xlab="year")
lines(1851:1962,combined_median,col="red")
lines(1851:1962,combined_lower,col="grey")
lines(1851:1962,combined_upper,col="grey")
abline(v=1890)
abline(v=1871)
abline(v=1944)

