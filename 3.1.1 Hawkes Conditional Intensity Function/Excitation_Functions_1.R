#Excitation Functions
alpha = 0.5
beta = 2
u=seq(0,1,0.01)
exp_k = alpha*exp(-beta*u)
plot(u, exp_k, type='l', xlab = "Time", cex.lab=2, cex.axis=2, ylab="Excitation Function")
p=2.5
c=1
k=0.5
pow_k = (k)/(u+c)^(p) 
lines(u, pow_k, type='l',col="red")
legend("topright", c("Exponential", "Power Law"), col=c("black", "red"),  lty=c(1,1), cex=2)


