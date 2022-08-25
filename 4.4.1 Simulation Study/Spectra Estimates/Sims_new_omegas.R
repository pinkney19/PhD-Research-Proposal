library(QZ)
library(phonTools)
#Multivariate Hawkes process
library(hawkes)
lambda0<-c(0.5,0.5)
#alpha<-matrix(c(0.1,0.5,0.7,0.2),byrow=TRUE,nrow=2)
alpha<-matrix(c(4,0,0,3),byrow=TRUE,nrow=2)
beta<-c(5,4)
set.seed(1)

T_end = 1000 #time length 
k.sections = 4 #number of tapers/blocks we average over -smoother here

T_s<-simulateHawkes(lambda0,alpha,beta,T_end)
#function for groups for smoothing
groups = function(k.sections){
  
  
  #split data into k sections of equal time length
  
  n=list(); 
  
  for(i in 1:2){
    
    n[[i]]=ceiling(max(T_s[[i]])/k.sections) 
  }
  
  groups = list(vector(mode = "list", length = k.sections), vector(mode = "list", length = k.sections))
  
  groups[[1]][[1]] = T_s[[1]][T_s[[1]]<=n[[1]]]
  for (j in 2:k.sections){
    groups[[1]][[j]] = T_s[[1]][T_s[[1]]<j*n[[1]] & T_s[[1]]>(j-1)*n[[1]]]   #-((j-1)*n[[1]])
  }
  groups[[2]][[1]] = T_s[[2]][T_s[[2]]<=n[[2]]]
  for(j in 2:k.sections){
    groups[[2]][[j]] = T_s[[2]][T_s[[2]]<j*n[[2]] & T_s[[2]]>(j-1)*n[[2]]]#-((j-1)*n[[2]]) #greater than those in previous group but less than j*n[[i]]
  }
  
  return(groups)
}


groups_new = groups(k.sections)
group1=groups_new[[1]] #groups for process 1
group2=groups_new[[2]] #groups for process 2


#function to obtain  periodogram estimate
periodogram_function = function(group1, group2, k.sections, T_s){ 
  
  J_1=vector(mode = "list", length=length(group1));I=list(); H_omega=vector(mode = "list", length = length(group2)); Big_T = NULL;estim_lambda_1=NULL;estim_lambda_2=NULL;
  omega=vector(mode = "list", length(group2)); J_2=vector(mode = "list", length = length(group2));
  
  
  
  img=sqrt(as.complex(-1)) #imaginary number i^2=-1
  Big_T = T_end/k.sections #1000 #ceiling(max(T_s[[i]])/k.sections) 
  
  
  for(i in 1:length(group1)){
    estim_lambda_1[i] = length(group1[[i]])/Big_T #N(T)/T
    estim_lambda_2[i] = length(group2[[i]])/Big_T #N(T)/T
    
    # points = seq(0,1000,length.out = 100) #adapting omegas for Alex
    
    for(k in 1:256){
      omega[[i]][k] = 2*pi*k/100
      H_omega[[i]][k] = (Big_T/sqrt(Big_T)) * (sinc((omega[[i]][k]*Big_T)/(2))) * (exp(-img*omega[[i]][k]*Big_T/2)) #FT of the taper 
      J_1[[i]][k] = (1/sqrt(Big_T))*(sum(exp(-img*group1[[i]]*omega[[i]][k]))) - (estim_lambda_1[i]*H_omega[[i]][k]) #equation 4.4
      J_2[[i]][k] = (1/sqrt(Big_T))*(sum(exp(-img*group2[[i]]*omega[[i]][k]))) - (estim_lambda_2[i]*H_omega[[i]][k]) #equation 4.4
    }
  }
  
  
  h = vector(mode = "list", length = 4)
  I = rep(list(h),k.sections)
  
  for(i in 1:length(group1)){
    I[[i]][[1]] = (J_1[[i]] * H(J_1[[i]]))[1,] #just de-lists it
    I[[i]][[2]] = (J_1[[i]] * H(J_2[[i]]))[1,]
    I[[i]][[3]] = (J_1[[i]] * H(J_2[[i]]))[1,]
    I[[i]][[4]] = (J_2[[i]] * H(J_2[[i]]))[1,]
  }
  
  
  
  #populate matrices 
  m1 = matrix(NA, ncol=length(group1), nrow=length(J_1[[1]]));m2 = matrix(NA, ncol=length(group1), nrow=length(J_1[[1]]));
  m3 = matrix(NA, ncol=length(group1), nrow=length(J_1[[1]]));m4 = matrix(NA, ncol=length(group1), nrow=length(J_1[[1]]));
  for(i in 1:ncol(m1)){
    m1[,i] =  I[[i]][[1]]
    m2[,i] =  I[[i]][[2]]
    m3[,i] =  I[[i]][[3]]
    m4[,i] =  I[[i]][[4]]
  }
  
  #averages
  smoothed_I1=apply(m1, 1, mean) #I_11 - periodogram
  smoothed_I2=apply(m2, 1, mean) #I_12
  smoothed_I3=apply(m3, 1, mean) #I_21
  smoothed_I4=apply(m4, 1, mean) #I_22 - periodogram 
  
  # par(mfrow=c(1,2))
  # plot(omega[[1]],smoothed_I1/(2*pi),type='l',xlab=expression(omega), ylab=expression(hat(I[11](omega))))
  # plot(omega[[1]],smoothed_I4/(2*pi),type='l',xlab=expression(omega), ylab=expression(hat(I[22](omega))))
  
  
  
  co_1_2 =  Mod(smoothed_I2)^2/(smoothed_I1*smoothed_I4)
  
  corrected_co =  Mod(smoothed_I2/(2*pi))^2/(smoothed_I1/(2*pi)*smoothed_I4/(2*pi))
  
  return(list(smoothed_I1, smoothed_I4, corrected_co, co_1_2))
  
  
}

p=periodogram_function(group1, group2, k.sections, T_s)

k=seq(1,256,1)
Big_T = 100
omega=2*pi*k/Big_T
#write.csv(omega, "omegas_used.csv")


MC_periodograms = matrix(NA, nrow = 256, ncol=1000)
for(i in 1:1000){
  T_s<-simulateHawkes(lambda0,alpha,beta,T_end)
  groups_new = groups(k.sections)
  group1=groups_new[[1]]
  group2=groups_new[[2]]
  MC_periodograms[,i] = periodogram_function(group1,group2, k.sections,T_s)[[1]] #gets I_11 estimate in each column
  print(i) #checker
}

#do the same for second process
MC_periodograms2 = matrix(NA, nrow = 256, ncol=1000)
for(i in 1:1000){
  T_s<-simulateHawkes(lambda0,alpha,beta,T_end)
  groups_new = groups(k.sections)
  group1=groups_new[[1]]
  group2=groups_new[[2]]
  MC_periodograms2[,i] = periodogram_function(group1,group2, k.sections,T_s)[[2]]
  print(i)
}

MC_coherence = matrix(NA, nrow = 256, ncol=1000)
for(i in 1:1000){
  T_s<-simulateHawkes(lambda0,alpha,beta,T_end)
  groups_new = groups(k.sections)
  group1=groups_new[[1]]
  group2=groups_new[[2]]
  MC_coherence[,i] = periodogram_function(group1,group2, k.sections,T_s)[[3]]
  print(i)
}

chosen_corr = as.numeric(MC_coherence[1,])

#Saving Data 
write.table(MC_periodograms, file="MC_periodograms_p1_1000_4.txt")
write.table(MC_periodograms2, file="MC_periodograms_p2_1000_4.txt")
write.table(chosen_corr, file="Corr_1000_4.txt")
