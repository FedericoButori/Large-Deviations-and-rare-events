#The following code reproduces the image 1 and 2. The first part produces the flow
#lines of the deterministic ODE, the second sample the process X_\varepsilon, record and plot only the rare
#realizations of interest. The tolerance for the event `F(X_\varepsilon(T)) = 0.286' is 10−3. 
#The second part of the code reproduce image 2.

Tsim=2
plot(x = 1,                 
            xlab = "Prey", 
            ylab = "Predators",
            xlim = c(0, 0.4), 
            ylim = c(0, 0.6),
            main = "Lotka-Volterra",
            type = "n")
  for(i in seq(0, 0.4, by=0.02)){
  for(j in seq(0, 0.6, by=0.02)){
  x_0=c(i,j)
  timeS=seq(0, Tsim, by=dt)
  path=matrix(NA, 2, length(timeS))
  path[1,]=timeS
  path[2,]=timeS
  path[,1]=x_0
  for (k in 2:length(timeS)) {
    path[,k]=path[,k-1] + B(path[,k-1], params)*dt
  #points(path[1,k], path[2,k], col="red")
  }
 points(path[1,], path[2,], type="l", col="gray")
  }
}
abline(v=0.286)


rare=0
pos=matrix(c(NA, NA), 2, 1)
for(ntimes in 1:10^6){ 
  samp.path=simulate(barx)
  if(!(T %in% is.nan(samp.path))){
  if(max(abs(samp.path[1, L] - 0.286)) < 10^-3){
    points(samp.path[1,], samp.path[2,], type="l", 
            col=rgb(red = 0, green = 0, blue = 1, alpha = 0.1))
    pos=cbind(pos, samp.path)
    rare=rare + 1
    }
  }
}

rare
points(phi[1,], phi[2,], type="l", col="orange", lwd="2")

############################################################

u=pos
df=data.frame(x=u[1,-1], y=u[2,-1])
dfphi = data.frame(z=phi[1,], t=phi[2,])
require(ggplot2)
p <- ggplot(df, aes(x, y))
p <- (p + stat_bin2d(bins = 150) + scale_fill_gradient(trans="log10") 
         + geom_path(dfphi, mapping=aes(x=z,y=t), color="orange") 
          + coord_fixed(ratio=0.7))
p
