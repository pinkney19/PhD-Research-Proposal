setwd("C:/Users/pinkney/Downloads/STOR603/Neuronal Synch in High Dimensions/Week 9/Asymptotic-Distributions")
I_4_1000 = read.table("Corr_old_beta_1000_4.txt")[,1]
I_4_10000 = read.table("Corr_old_beta_10000_4.txt")[,1]

I_10_1000 = read.table("Corr_old_beta_1000_10.txt")[,1]
I_10_10000 = read.table("Corr_old_beta_10000_10.txt")[,1]

R = 0 #set  true coherence = 0 since they are independent processes
n = 4 #number of tapers
p = 2

library(hypergeo)

Goodman_eq = function(n,p,R,R_hat){
  t1 = gamma(n) / (gamma(p-1)*gamma(n-p+1))
  t2 = (1-R)
  t3 = (R_hat)^(p-2)
  t4 = (1-R_hat)^(n-p)
  t5 = hypergeo(n,n,(p-1), (R*R_hat)) #??
  return(t1*t2*t3*t4*t5)
}

plot_both_corr =function(chosen_corr, n,p,R){
  #par(mfrow=c(1,1))
  x=seq(0,1,0.001) #points at which to evaluate Goodman
  t=Goodman_eq(n,p,R,x)
  hist(chosen_corr,freq=F, xlim=c(0,1), main="", breaks=20, col="tomato", xlab = expression(hat(R)[12](omega)),cex.lab=2,cex.axis=2)
  lines(x,t,type='l', main="Theoretical Distribution", col='black')
}
par(mfrow=c(2,2))
plot_both_corr(I_4_1000, n , p, R)
#title(main="K=4, T=1000")
plot_both_corr(I_4_10000, n , p, R)
#title(main="K=4, T=10,000")

#k=10
plot_both_corr(I_10_1000, n=10 , p, R)
#title(main="K=10, T=1000")
plot_both_corr =function(chosen_corr, n,p,R){
  #par(mfrow=c(1,1))
  x=seq(0,1,0.001) #points at which to evaluate Goodman
  t=Goodman_eq(n,p,R,x)
  hist(chosen_corr,freq=F, xlim=c(0,1), main="", breaks=10, col="tomato", xlab = expression(hat(R)[12](omega)), cex.lab=2, cex.axis=2)
  lines(x,t,type='l', main="Theoretical Distribution", col='black')
}
plot_both_corr(I_10_10000, n=10 ,p, R)
#title(main="K=10, T=10,000")

#Now with inflated betas
I_4_1000_2 = read.table("Corr_new_beta_1000_4.txt")[,1]
I_4_10000_2 = read.table("Corr_new_beta_10000_4.txt")[,1]

I_10_1000_2 = read.table("Corr_new_beta_1000_10.txt")[,1]
I_10_10000_2 = read.table("Corr_new_beta_10000_10.txt")[,1]

R = 0 #set  true coherence = 0 since they are independent processes
n = 4 #number of tapers
p = 2

library(hypergeo)

Goodman_eq = function(n,p,R,R_hat){
  t1 = gamma(n) / (gamma(p-1)*gamma(n-p+1))
  t2 = (1-R)
  t3 = (R_hat)^(p-2)
  t4 = (1-R_hat)^(n-p)
  t5 = hypergeo(n,n,(p-1), (R*R_hat)) #??
  return(t1*t2*t3*t4*t5)
}

plot_both_corr =function(chosen_corr, n,p,R){
  #par(mfrow=c(1,1))
  x=seq(0,1,0.001) #points at which to evaluate Goodman
  t=Goodman_eq(n,p,R,x)
  hist(chosen_corr,freq=F, xlim=c(0,1), main="", breaks=20, col="tomato")
  lines(x,t,type='l', main="Theoretical Distribution", col='black')
}
par(mfrow=c(2,2))
plot_both_corr(I_4_1000_2, n , p, R)
title(main="K=4, T=1000")
plot_both_corr(I_4_10000_2, n , p, R)
title(main="K=4, T=10,000")

#k=10
plot_both_corr =function(chosen_corr, n,p,R){
  #par(mfrow=c(1,1))
  x=seq(0,1,0.001) #points at which to evaluate Goodman
  t=Goodman_eq(n,p,R,x)
  hist(chosen_corr,freq=F, xlim=c(0,1), main="", breaks=10, col="tomato")
  lines(x,t,type='l', main="Theoretical Distribution", col='black')
}
plot_both_corr(I_10_1000_2, n=10 , p, R)
title(main="K=10, T=1000")

plot_both_corr(I_10_10000_2, n=10 ,p, R)
title(main="K=10, T=10,000")


#Dependent Plots 
setwd("C:/Users/pinkney/Downloads/STOR603/Neuronal Synch in High Dimensions/Week 9/Asymptotic-Distributions/Dependent")
C_1000_10_dep = read.table("Corr_1000_10.txt")[,1]
C_1000_4_dep = read.table("Corr_1000_4.txt")[,1]
C_10000_10_dep = read.table("Corr_10000_10.txt")[,1]
C_10000_4_dep = read.table("Corr_10000_4.txt")[,1]

plot_both_corr =function(chosen_corr, n,p,R){
  #par(mfrow=c(1,1))
  x=seq(0,1,0.001) #points at which to evaluate Goodman
  t=Goodman_eq(n,p,R,x)
  hist(chosen_corr,freq=F, xlim=c(0,1), main="", breaks=20, col="tomato", xlab= expression(hat(R)[12](omega)), cex.lab=1.7, cex.axis=1.7)
  #lines(x,t,type='l', main="Theoretical Distribution", col='black')
}

par(mfrow=c(1,2))
plot_both_corr(C_1000_4_dep, n=4 , p, R)
plot_both_corr(C_1000_10_dep, n=10 , p, R)

plot_both_corr(C_10000_4_dep, n=4 , p, R)
plot_both_corr(C_1000_10_dep, n=10 , p, R)
plot_both_corr(C_10000_10_dep, n=10 , p, R)
C_10000_10_dep
