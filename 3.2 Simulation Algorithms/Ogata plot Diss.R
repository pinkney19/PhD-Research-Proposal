#Simulating Univariate Hawkes process by thinning algorithm
par(mfrow=c(1,1))
#Set current time and event counter
t=0
i=1

#time to simulate until
big_T=10;

#Hawkes Process initial parameters
lambda = 1; alpha = 0.5; beta=2;
#lambda = 0.5; alpha = 4; beta=4.5;
#Current maximum of the process
M = lambda;

#function for intensity rate
lambda_t = function(t, t_i, lambda, alpha, beta){lambda + sum(alpha*exp(-beta*(t-t_i[t_i<=t])))}

#checks
#t=4; t_i=c(1,2,3,4,5); alpha=1; beta=1.2;
#1+exp(-1.2*(4-1))+exp(-1.2*(4-2))+exp(-1.2*(4-3))
#lambda_t(t, t_i, lambda, alpha, beta)

P = NULL; # empty vector for accepted times
M_out=NULL; t_out=NULL; E_out=NULL; U_out=NULL; lambdas=NULL;#set up empty vectors for outputs
e=0
#e = 10^(-10) # a tiny value > 0

set.seed(1) #for checking with Hawkes_ogata package
while (t<big_T){
  
  #find new upper bound
  M = lambda_t(t+e, P, lambda, alpha, beta)
  #save upper bounds for plots
  M_out = append(M_out, M)
  
  #generate next candidate point
  E = rexp(1, M) #i.e. sample inter-arrival time 
  t = t + E #update current time
  
  #save values for plots
  E_out = append(E_out, E) #all arrival times
  t_out = append(t_out, t) #all current times
  
  #accept with some probabilty 
  U = runif(1,0,M)
  U_out = append(U_out, U)
  
  
  if (U<= lambda_t(t, P, lambda, alpha, beta) && t<big_T){
    P = append(P, t) #accepted arrivals
    i = i+1 #update counter
  }
  
}
print(P)

# #check 
# library(hawkesbow)
# set.seed(1)
# hawkes_ogata(10, lambda=1, alpha=0.5, beta=2)$p
# set.seed(1)
# #par(mfrow=c(1,2))
# plot(hawkes_ogata(10, lambda=1, alpha=0.5, beta=2))

#Plots 

times = seq(0,big_T, by = 0.01)
lambdas = numeric();
for (i in 1:length(times)){
  lambdas[i] = lambda_t(times[i], P, lambda, alpha, beta)
}


plot(times, lambdas, type='l', lwd=1.5,xlab='t', ylab=(expression(D~bar(Lambda))), ylim=c(0, max(U_out)+0.7), main ="", cex.lab=1.2, cex.axis=1.2)

#which are accepted?
v1 = match(P,t_out)
points(P, U_out[v1], col='black', pch=1)

#show time points on x axist
points(P, rep(0, length(U_out[v1])), col="black", pch=15)

#add in dashed vertical lines
for(i in 1:length(P)){
  segments(P[i],0,P[i], U_out[v1][i], col="green", lty=3)
}

#which points are not accepted?
not_accepted = setdiff(t_out,P)
#get indices of those not accepted
v2 = match(not_accepted, t_out)
points(not_accepted, U_out[v2], col='red', pch=17, lwd=3)

#add in upper bounds
segments(0,M_out[1],t_out[1],M_out[1], col = 'blue',lwd=2)
for (i in 1:length(M_out)){
  segments(t_out[i], M_out[i+1],t_out[i+1],M_out[i+1], col = 'blue',lwd=3)
}
legend("topright", legend=c(expression(paste(bar(Lambda))),expression(paste(Lambda(t))), "Accepted", "Rejected"),  pch = c(NA,NA,1,17),lwd=c(3,1,NA,NA), col=c("blue", "black", "black", "red"))


