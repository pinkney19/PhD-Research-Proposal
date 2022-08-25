par(mfrow=c(1,1))
set.seed(12)
times = rexp(5,1)
y = seq(1,5,1)
times=sort(times)
plot(times, y, pch=16, ylab=expression(N(t)), xlab=expression(t), xaxt = "n", axes=F, cex.lab=2, cex.axis=2)
for (j in 1:length(y)){
  segments(times[j], y[j],times[j+1],y[j], col = 'red',lwd=3)
}

segments(times[5],y[5],times[5]+1,y[5], col = 'red',lwd=2)

for (i in 1:length(y)){
  segments(times[i], y[i],times[i],0, col = 'black',lty=2)
}
axis(1, at = c(0,times), labels = c(NA,expression(t[1]),expression(t[2]), expression(t[3]), expression(t[4]), expression(t[5])), cex.axis=2, cex=2)
axis(2, at=c(0,y), cex=2, cex.axis=2)
