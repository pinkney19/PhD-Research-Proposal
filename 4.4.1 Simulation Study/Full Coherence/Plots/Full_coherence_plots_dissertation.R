#full coherence plots for T=1000
par(mar=c(7,7,1,1))
setwd("C:/Users/pinkney/Downloads/STOR603/Neuronal Synch in High Dimensions/Week 9/Full-Coherence")
C_10_1000 = read.table("Full_Corr_1000_10.txt")
C_4_1000 = read.table("Full_Corr_1000_4.txt")

C_10_1000_dep = read.table("Full_Corr_dep_1000_10.txt")
C_4_1000_dep = read.table("Full_Corr_dep_1000_4.txt")

k=seq(1,256,1)
Big_T = 100
omega=2*pi*k/Big_T


library(matrixStats)
coherence_plot=function(MC_coherence,omega){
c=apply(MC_coherence, 1, mean) #point wise average
c2 = rowQuantiles(Re(as.matrix(MC_coherence)), probs=0.975)
c3 = rowQuantiles(Re(as.matrix(MC_coherence)), probs=0.025)
plot(omega, c, type='l', ylim=c(0,1), xlab=expression(omega), ylab=expression(hat(R)["12"](omega)), cex.lab=2, cex.axis=2)
lines(omega, c2, col='red')
lines(omega, c3, col='red')
legend("topright", c("Estimated Coherence", "95% Confidence Interval"), col=c("black", "red"), lty=c(1,1), cex=1.7)
}
par(mfrow=c(2,2))
coherence_plot(C_4_1000, omega)
coherence_plot(C_10_1000, omega)
coherence_plot(C_4_1000_dep, omega)
coherence_plot(C_10_1000_dep, omega)
