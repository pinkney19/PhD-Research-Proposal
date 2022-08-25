k_seq=seq(1,256, by=1)
n=500
omega = 2*pi*k_seq/n #n=Big_t here

nu=0.5
alpha=4
beta=5


#theoretical spectrum 
t1 =  nu*beta/((2*pi)*(beta-alpha))
t2 = alpha*((2*beta)-alpha)
t3 = ((beta-alpha)^2) + (omega)^2
f_omega = t1* (1+(t2/t3))

#plot
par(mar=c(5,5,1,1))
plot(omega,f_omega,type='l',col='red', xlab=expression(omega), ylab=expression(f+(omega)), ylim=c(0,20), cex.lab=1.5, cex.axis=1.5)

#different example 

nu=0.5
alpha=4
beta=4.8


#theoretical spectrum 
t1 =  nu*beta/((2*pi)*(beta-alpha))
t2 = alpha*((2*beta)-alpha)
t3 = ((beta-alpha)^2) + (omega)^2
f_omega = t1* (1+(t2/t3))
lines(omega,f_omega,type='l',col='blue')

#different example 
nu=0.5
alpha=4
beta=9
t1 =  nu*beta/((2*pi)*(beta-alpha))
t2 = alpha*((2*beta)-alpha)
t3 = ((beta-alpha)^2) + (omega)^2
f_omega = t1* (1+(t2/t3))
lines(omega,f_omega,type='l',col='green')

legend("topright", c(expression(beta~ "= 4.8 "), expression(beta~"= 5 "), expression(beta~"= 9 ")), lty=c(1,1,1), col=c("blue", "red", "green"), cex=1.5) 
