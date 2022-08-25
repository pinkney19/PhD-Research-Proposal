# Script to process Schultz data.
# Original data from the paper
# J. Tang, S. C. Ardila Jimenez, S. Chakraborty, and S. R.
# Schultz. Visual receptive field properties of neurons in
# the mouse lateral geniculate nucleus. Plos One, 11(1):
#  1-34, 2015

# Note: I have only extracted two neurons 
# We can maybe use this as an example of data for now, but the data likely
# exhibits non-stationary behaviour because it is stimulus driven

# I used the same data in Figure 4, of this paper:
# https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8645237

# We can treat the periodogram obtained across each of the 200 trials as a tapered estimate
# simply calculate periodogram for each trial, and then average across trials

# calculate the range of omega by looking at max T across all trials.
# use same range of omega for each trial
par(mar=c(7,7,1,1))
setwd("C:/Users/pinkney/Downloads/STOR603/Neuronal Synch in High Dimensions/Week 7/Real Data")
library(phonTools)
library(QZ)
library(R.matlab)
data <- readMat("schultz_127_130.mat")

E <- data$E

Ntrial = 200
P = 2

E1 = list()
E2 = list()

Tmax = 0

for (i  in 1:P){
  for (trial in 1:Ntrial){
  idx = (trial-1)*P + i
  # Data is interlaced, i.e. i=1, then i=2
    if(i ==1 ){
      E1[[trial]] = E[[idx]][[1]]
    }else if( i==2){
      E2[[trial]] = E[[idx]][[1]]
    }
    Tmax = max(Tmax, E[[idx]][[1]])
  }
}


#function to obtain  periodogram estimate
periodogram_function_real_data = function(E1, E2, Ntrial){ 
  
  J_1=vector(mode = "list", length=length(E1));I=list(); Big_T = NULL;estim_lambda_1=NULL;estim_lambda_2=NULL;
  omega=vector(mode = "list", length(E1)); J_2=vector(mode = "list", length = length(E1));H_omega=vector(mode = "list", length=length(E1))
  
  
  
  img=sqrt(as.complex(-1)) #imaginary number i^2=-1
  Big_T = Tmax
  
  for(i in 1:length(E1)){
    estim_lambda_1[i] = length(E1[[i]][1,])/Big_T #N(T)/T
    estim_lambda_2[i] = length(E2[[i]][1,])/Big_T #N(T)/T
    
    # points = seq(0,1000,length.out = 100) #adapting omegas for Alex
    
    for(k in 1:40){
      omega[[i]][k] = 2*pi*k/Big_T
      H_omega[[i]][k] = (Big_T/sqrt(Big_T)) * (sinc((omega[[i]][k]*Big_T)/(2))) * (exp(-img*omega[[i]][k]*Big_T/2)) #FT of the taper 
      J_1[[i]][k] = (1/sqrt(Big_T))*(sum(exp(-img*E1[[i]][1,]*omega[[i]][k]))) - (estim_lambda_1[i]*H_omega[[i]][k]) #equation 4.4
      J_2[[i]][k] = (1/sqrt(Big_T))*(sum(exp(-img*E2[[i]][1,]*omega[[i]][k]))) - (estim_lambda_2[i]*H_omega[[i]][k]) #equation 4.4
    }
  }
  
  
  h = vector(mode = "list", length = 4)
  I = rep(list(h),Ntrial)
  
  for(i in 1:length(E1)){
    I[[i]][[1]] = (J_1[[i]] * H(J_1[[i]]))[1,] #just de-lists it
    I[[i]][[2]] = (J_1[[i]] * H(J_2[[i]]))[1,]
    I[[i]][[3]] = (J_1[[i]] * H(J_2[[i]]))[1,]
    I[[i]][[4]] = (J_2[[i]] * H(J_2[[i]]))[1,]
  }
  
  
  
  #populate matrices 
  m1 = matrix(NA, ncol=length(E1), nrow=length(J_1[[1]]));m2 = matrix(NA, ncol=length(E1), nrow=length(J_1[[1]]));
  m3 = matrix(NA, ncol=length(E1), nrow=length(J_1[[1]]));m4 = matrix(NA, ncol=length(E1), nrow=length(J_1[[1]]));
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
  
  
  
  #co_1_2 =  Mod(smoothed_I2)^2/(smoothed_I1*smoothed_I4)
  
  corrected_co =  Mod(smoothed_I2/(2*pi))^2/(smoothed_I1/(2*pi)*smoothed_I4/(2*pi))
  
  return(list(smoothed_I1, smoothed_I4, corrected_co))
  
  
}

p=periodogram_function_real_data(E1, E2, Ntrial)
k=seq(1,40,1)
omega = 2*pi*k/Tmax
par(mfrow=c(1,3))
plot(k, (p[[1]]/(2*pi)), type='l', xlab="", ylab="", main = "", xaxt="n", yaxt="n", cex.axis=2)
axis(2,cex.axis=2, cex=2)
axis(1,cex.axis=2, cex=2)
mtext(expression(I[11](f)), side=2, line=2.2, cex=2)
mtext("Frequency (Hz)", side=1, line=2.2, cex=2)


plot(k, (p[[2]]/(2*pi)), type='l', xlab="", ylab="", main = "", xaxt="n", yaxt="n", cex.axis=2)
axis(2,cex.axis=3)
axis(1,cex.axis=3)
mtext(expression(I[22](f)), side=2, line=2.2, cex=2)
mtext("Frequency (Hz)", side=1, line=2.2, cex=2)
plot(k, p[[3]], type='l', ylim=c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", cex.axis=2)
axis(2,cex.axis=2)
axis(1,cex.axis=2)
mtext(expression(hat(R)[12](omega)), side=2, line=2.2, cex=2)
mtext("Frequency (Hz)", side=1, line=2.2, cex=2)

par(mfrow=c(1,2))
x=seq(0,2,length.out=30)
y=seq(0, 200, length.out=30)
plot(x,y, pch="", bty="n", xlab="Time (ms)", ylab="Trial", main = "", cex.lab=2, cex.axis=2)
#h=rep(1, length(E1[[1]]))
#points(E1[[1]],h, pch=16)
#points(E1[[2]],rep(2, length(E1[[2]])) )
#Raster plot
x=list();
for(i in 1:length(E1)){
  x[[i]] = rep(i,length(E1[[i]]))
  points(E1[[i]], x[[i]], pch=16)
  
}

x=seq(0,2,length.out=30)
y=seq(0, 200, length.out=30)
plot(x,y, pch="", bty="n",  xlab="Time (ms)", ylab="Trial", main="", cex.lab=2, cex.axis=2)
#h=rep(1, length(E1[[1]]))
#points(E1[[1]],h, pch=16)
#points(E1[[2]],rep(2, length(E1[[2]])) )
#Raster plot
x=list();
for(i in 1:length(E2)){
  x[[i]] = rep(i,length(E2[[i]]))
  points(E2[[i]], x[[i]], pch=16)
  
}







# 
# New_E1 = list();
# New_E1[[1]] = E1[[1]][1,]
# for(i in 2:length(E1)){
#   New_E1[[i]] = E1[[i]][1,]+(i-1)
# }
# #plot the process
# New_E1 = unlist(New_E1)
# x=seq(1,length(New_E1), by=1)
# plot(sort(New_E1),x,type='s')
# 
# plot(k, (p[[1]]/(2*pi)), type='l', xlab="", ylab="", main = "")
# mtext(expression(I[11](f)), side=2, line=2.2, cex=1.5)
# mtext("Frequency (Hz)", side=1, line=2.2, cex=1.5)
