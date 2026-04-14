#The following code compute the solutions to the Hamiltonian system for the action minimizing path. 
#By simply substituting at the definition of the variable val an arrow of values, 
#it can be used for a line search for the value of lambda for a desired output of F(\phi_T).

Tmax=5 # 0.0901596 -> lambda = 0.286
dt=0.01
eps=0.004
maxiter=25000
time=seq(0, Tmax, by=dt)
L=length(time)

params=c(1,5,1,0.1)
barx=c(0.1414214,0.3414214)
val= 0.0901596 #seq(0.09015971, 0.09015973, by=0.000000005)
z_lambda=rep(NA, length(val))
c=1

for(t in val){
lambda=t  
phi=matrix(rep(barx, L), 2, L)
theta=matrix(rep(0, L), 2, L)
phi[,1]=barx
theta[, L]=c(lambda,0)

i=0
C=F
N=F

while(i < maxiter && C==F && N==F){
  oldtheta=theta
  oldphi=phi
  for (k in 1:(L-1)) {
    theta[,L-k] = theta[,(L-k+1)] +
    gradphiH(phi[,(L-k+1)], theta[,(L-k+1)], params)*dt
  }
  
  for (k in 2:L) {
    phi[,k] = phi[,k-1] + gradthetaH(phi[,k-1], theta[,k-1], params)*dt
  }
  
  if(T %in% is.nan(phi)){
    z_lambda[c]=0
    N=T} 
     
  else if(max(abs(oldphi-phi)) < 10^-13){
    z_lambda[c]=phi[1,L]
    C=T
  }
  i=i+1
}
c=c+1

}
plot(val, z_lambda)
