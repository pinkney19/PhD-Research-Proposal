#Ogata method 
lambda_t= function(t, t_i, lambda, alpha, beta){lambda + sum(alpha*exp(-beta*(t-t_i[t_i<t])))}

Hawkes_sim_ogata = function(big_T, lambda,alpha, beta){
  
  #Simulating Univariate Hawkes process by thinning algorithm
  #big_T = time to simulate until
  #alpha,beta and lambda are the initial parameters of the process
  
  #Set current time and event counter
  t=0
  i=0
  
  
  #Current maximum of the process
  M = lambda;
  
  #function for intensity rate
  #lambda_t = function(t, t_i, lambda, alpha, beta){lambda + sum(alpha*exp(-beta*(t-t_i[t_i<t])))}
  
  P = NULL; # empty vector for accepted times
  M_out=NULL; t_out=NULL; E_out=NULL; U_out=NULL; lambdas=NULL;#set up empty vectors for outputs
  e = 10^(-10) # a tiny value > 0
  
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
  
  
  #plotting N(t) and E(N(t))
  Nt=seq(1,length(P), by=1)
  times=seq(0,max(P), by=0.01)
  
  m=alpha/beta
  
  T1 = (lambda*times)/(1-m)
  T2 = m/(1-m)^2
  T3 = lambda/beta
  T4 = 1-exp(-(beta)*(1-m)*times)
  
  EN = T1 - (T2*T3*T4)
  EN2 = lambda/(1-m)*times
  
  plot(P,Nt ,xlab = "t", ylab='Count', type='n',ylim=c(0,max(EN)), main = "",cex.lab=1.5, cex.axis=1.5) #s for stair steps
  lines(times, EN, lty=2, col="blue")
  lines(P,Nt, type='s')
  
  
  return(P)
}
par(mfrow=c(1,2))
set.seed(25)


P5 = Hawkes_sim_ogata(1000, 1,0.5, 2) 
legend("topleft", c("N(t)", "E[N(t)]"), bty = "n", lty=c(1,2), col=c("black","blue"), cex=1.5)


#####################################################################################
####################################################################################
####################################DASSIOS METHOD###################################


#Intensity Function 
lambda_t_func = function(a, lambda0, alpha, beta, t, t_k){
  l_t = a + (lambda0-a)*(exp(-beta*t)) + sum(alpha*exp(-beta*(t-t_k[t_k<t])))
  return(l_t)
}


Exact_sim = function(N,lambda, alpha, beta){
  #Initial Conditions
  a=lambda #background intensity
  #set up empty vectors
  #T_plus = NULL; T_minus = NULL; D=NULL; s1=NULL; s2=NULL; L_plus=NULL; L_minus=NULL; u=NULL;
  #s=NULL; lambda0 = 1; Time=NULL;
  #T_plus[1] = 0; T_minus[1]=0; L_plus[1]=lambda0; L_minus[1]=lambda0; Time[1]=0;
  
  L_plus=NULL; L_minus=NULL;D=NULL;s1=NULL; s2=NULL; s=NULL;lambda0 = 1; Time=NULL;Time[1]=0;
  L_plus[1]=lambda0; L_minus[1]=lambda0; 
  #Set current time 
  
  count = 0
  t=0
  
  #main loop
  #set.seed(1)
  while(t<N){
    
    u = runif(1,0,1)
    s2 = rexp(1,a) #sample from the exponential with rate a
    D = 1 + (beta*log(u))/(tail(L_plus, n=1)-a)
    
    if(D>0){
      s1 = (-1/beta)*log(D) #sample via the inverse CDF
      s = min(s1,s2) #get arrival time
    }
    else{
      s = s2
    }
    
    #update current time
    t = t+s
    #Record the ith jump time
    Time = append(Time, t)
    #print(t)
    
    #update intensities
    L_plus_i_1 = tail(L_plus, n=1)
    term_to_append = ((L_plus_i_1-a)*(exp(-beta*s))) + a 
    L_minus = append(L_minus, term_to_append)
    L_plus = append(L_plus, L_minus+alpha)
    
    #L_minus[i] = (L_plus[i-1]-a)*(exp(-beta*(Time[i]-Time[i-1]))) + a
    # L_plus[i] = L_minus[i] + alpha
    
    #update counter
    count=count+1
  }
  
  #return(Time)
  
  #plotting N(t) and E(N(t))
  x=seq(1,length(Time), by=1)
  times=seq(0,max(Time), by=0.01)
  
  m=alpha/beta
  
  T1 = (lambda*times)/(1-m)
  T2 = m/(1-m)^2
  T3 = lambda/beta
  T4 = 1-exp(-(beta)*(1-m)*times)
  
  EN = T1 - (T2*T3*T4)
  
  plot(Time,x ,xlab = "t", ylab='Count', type='n',ylim=c(0,max(EN)), main="", cex.lab=1.5, cex.axis=1.5) #s for stair steps
  lines(times, EN, lty=2, col="blue")
  lines(Time,x, type='s')
  legend("topleft", c("N(t)", "E[N(t)]"), col=c("black", "blue"), lty=c(1,2), bty="n", cex=1.5)
  
}

set.seed(25)
Exact_sim(1000, 1,0.5, 2)



#############compare run times##################
# set.seed(25)
# system.time(Hawkes_sim_ogata(1000, 1,0.5, 2))
# set.seed(25)
# system.time(Exact_sim(1000, 1,0.5, 2))
# 
# 
# set.seed(25)
# system.time(Hawkes_sim_ogata(10000, 1,0.5, 2))
# set.seed(25)
# system.time(Exact_sim(10000, 1,0.5, 2))


