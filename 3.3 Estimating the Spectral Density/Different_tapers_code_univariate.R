#library(hawkesbow)
library(phonTools) # for sinc function
library(hawkes)

alpha=0.5#amplitude of spike
lambda0=1 #baseline intensity
beta=0.55 #speed of exponential decay 
k.sections = 10 #the number of sections you want to split your data into in order to average over
set.seed(1)
#T_s=hawkes_ogata(10000, nu, alpha, beta)$p  #simulate spike times via ogata's method
T_s = simulateHawkes(lambda0,alpha,beta,10000)
T_s = unlist(T_s)




#function to obtain  periodogram estimate
periodogram_function = function(group, k.sections){
  J=NULL;I=NULL;omega=NULL;x=NULL; H_omega=NULL;
  Big_T = 100#ceiling(max(T_s)/k.sections) 
  estim_lambda = length(group)/Big_T #N(T)/T
  img=sqrt(as.complex(-1)) #imaginary number i=-1
  
  for(k in 1:256){ #unsure about this / arbitrary 
    omega[k] = 2*pi*k/Big_T
    H_omega[k] = (Big_T/sqrt(Big_T)) * (sinc((omega[k]*Big_T)/(2))) * (exp(-img*omega[k]*Big_T/2)) #FT of the taper - pi or no pi?
    J[k] = (1/sqrt(Big_T))*(sum(exp(-img*group*omega[k]))) - (estim_lambda*H_omega[k]) #equation 4.4
    I[k] = J[k]*Conj(J[k])
    
  }
  
  return(I)
}

#function to obtain averaged estimate and plot vs theoretical spectrum

periodogram = function(T_s, nu, alpha, beta,k.sections){
  
  
  #split data into k sections of equal time length
  
  n=100#ceiling(max(T_s)/k.sections) #if k = 10, then n = 1000
  groups = list()
  
  groups[[1]] = T_s[T_s<=n]
  for (i in 2:k.sections){
    groups[[i]] = T_s[T_s<i*n & T_s>(i-1)*n]-((i-1)*n) #greater than those in previous group but less than i*n
  }
  
  #get periodogram estimates for each group
  
  I = list()
  for (i in 1:k.sections){
    I[[i]] = periodogram_function(groups[[i]], k.sections)
  }
  
  #get smoothed estimate - i.e. average over the individual estimates.
  
  combined_Is = sapply(1:k.sections, function(i) rep(I[[i]]))
  smoothed_I=apply(combined_Is, 1, mean)
  
  
  k_seq=seq(1,256, by=1)
  omega = 2*pi*k_seq/n #n=Big_t here
  
  
  #theoretical spectrum 
  t1 =  nu*beta/((2*pi)*(beta-alpha))
  t2 = alpha*((2*beta)-alpha)
  t3 = ((beta-alpha)^2) + (omega)^2
  f_omega = t1* (1+(t2/t3))
  
  #plot
  plot(omega,smoothed_I/(2*pi),type='l',col='red', xlab=expression(omega), ylab=expression(hat(I)(omega)), ylim=c(0,10), cex.lab=3, cex.axis=2)
  lines(omega,f_omega, col="blue")
  legend("topright", col=c("red","blue"), lty=c(1,1), legend = c("Estimated", "Theoretical"), cex=2)
  
  
}
par(mfrow=c(1,3))
par(mar=c(8,8,1,1))

#try with one taper
periodogram(T_s, nu=lambda0, alpha, beta,1)

#try with 10 tapers
periodogram(T_s, nu=lambda0, alpha, beta,10)

#try with 1000 tapers
periodogram(T_s, nu=lambda0, alpha, beta,1000)

#uncorrected periodogram 

