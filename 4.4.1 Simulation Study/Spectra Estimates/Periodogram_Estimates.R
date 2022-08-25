setwd("C:/Users/pinkney/Downloads/STOR603/Neuronal Synch in High Dimensions/Week 8/New Omega Simulations/Sims_new_omegas")

#Plots for Dissertation 

#Required Parameters and plotting Functions
library(QZ)
library(phonTools)
#Multivariate Hawkes process
library(hawkes)
lambda0<-c(0.5,0.5)
#alpha<-matrix(c(0.1,0.5,0.7,0.2),byrow=TRUE,nrow=2)
alpha<-matrix(c(4,0,0,3),byrow=TRUE,nrow=2)
beta<-c(5,4)
set.seed(1)

plotting_point_wise_averaged_periodograms = function(I_matrix, omega, nu, beta, alpha,f,j){
  p = apply(I_matrix,1,mean) #point wise average
  p2 = rowQuantiles(Re(as.matrix(I_matrix)), probs=0.975)
  p3 = rowQuantiles(Re(as.matrix(I_matrix)), probs=0.025)
  
  #theoretical spectrum
  t1 = nu*beta/((2*pi)*(beta-alpha))
  t2 = alpha *((2*beta)-alpha)
  t3 = ((beta-alpha)^2) + (omega^2)
  f_omega = t1*(1+(t2/t3))
  plot(omega, f_omega, col='blue',type='l', xlab=expression(omega), ylab=expression(hat(I)["11"](omega)), ylim=c(0, max(p2/(2*pi))), cex.axis=2, cex.lab=2)
  #segments(x0=omega, y0=p3/(2*pi), x1=omega, y1=p2/(2*pi) ,code=3,angle = 90, length = 0.1)
  lines(omega, (p/(2*pi)), type='l')
  lines(omega, p2/(2*pi), type='l', col="red")
  lines(omega, p3/(2*pi), type='l', col="red")
  legend("topright", c("Periodogram", "Theoretical", "95% CI"), lty=c(1,1,1), col=c("black", "blue", "red"), cex=1.7, lwd=c(1,1,1))
}


#Scenario 1 - Independent Processes T=1,000 k=10

I1_1000_10 = read.table("MC_periodograms_p1_1000_10.txt")
I2_1000_10 = read.table("MC_periodograms_p2_1000_10.txt")

# 
# p = apply(I1_1000_10,1,mean) #point wise average
# library(matrixStats)
# p2 = colQuantiles(Re(as.matrix(I1_1000_10)), probs=0.95)
# p3 = colQuantiles(Re(as.matrix(I1_1000_10)), probs=0.05)


#Scenario 2 - T=1000, k=4

I1_1000_4 = read.table("MC_periodograms_p1_1000_4.txt")
I2_1000_4 = read.table("MC_periodograms_p2_1000_4.txt")


I1_10000_10 = read.table("MC_periodograms_p1_10000_10.txt")
I2_10000_10 = read.table("MC_periodograms_p2_10000_10.txt")


I1_10000_4 = read.table("MC_periodograms_p1_10000_4.txt")
I2_10000_4 = read.table("MC_periodograms_p2_10000_4.txt")


k=seq(1,256,1)
Big_T = 100
omega = 2*pi*k/Big_T

#plots for I_{11}(\omega)
par(mfrow=c(2,2))
plotting_point_wise_averaged_periodograms(I1_1000_4, omega, lambda0[1], beta[1], alpha[1,1])
#title("Process 1: T=1000 and k=4")
plotting_point_wise_averaged_periodograms(I1_10000_4, omega, lambda0[1], beta[1], alpha[1,1])
#title("Process 1: T=10000 and k=4")
plotting_point_wise_averaged_periodograms(I1_1000_10, omega, lambda0[1], beta[1], alpha[1,1])
#title("Process 1: T=1000 and k=10")
plotting_point_wise_averaged_periodograms(I1_10000_10, omega, lambda0[1], beta[1], alpha[1,1])
#title("Process 1: T=10000 and k=10")