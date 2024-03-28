##########################################################
####
####  Fictitious BMI data 
####
##########################################################
x = c ( 26, 32, 19, 23, 22, 25, 28, 19, 23, 27, 26, 22, 18, 26, 24 )     ## Data

####  Prior Parameters (assuming a NSinv prior)
mu0 = 29
sigma0 = 8
phi0 = sigma0^2
kappa0 = 2
nu0 = 1

####  Get ready to plot the prior distributions
muVec = seq( 20 , 38 , 0.1 )
phiVec = seq(0.5, 200, by=0.25 )
muPrior = (1/sqrt(phi0/kappa0)) * dt( (muVec-mu0)/(sqrt(phi0/kappa0)) , nu0 )

logConst = (nu0/2)*log(nu0/2) + (nu0/2)*log(phi0) - lgamma(nu0/2)
phiPrior = exp( logConst - ((nu0/2)+1)*log(phiVec) - nu0*phi0/(2*phiVec) )

# windows( 9 , 9 )
par( mfrow=c(2,1) )
plot( muVec , muPrior , type="l" , col="red" , ylim=c(0,0.08) , lwd=2 )
abline( h=0 , col="gray")
plot( phiVec , phiPrior , type="l" , col="red" , ylim=c(0,0.01) , lwd=2 ) 
abline( h=0 , col="gray" )


xbar = mean(x)
s2 = var(x)
n = length(x)

####  Bayesian updating
mu1    = (kappa0/(kappa0+n))*mu0 + (n/(kappa0+n))*xbar
kappa1 = kappa0 + n
nu1    = nu0 + n
phi1   = (nu0/(nu0+n))*phi0 + ((n-1)/(nu0+n))*s2 +
         (kappa0*n/((nu0+n)*(kappa0+n)))
print( c(mu1,kappa1,phi1,nu1) )

####  Function for normal-inverse-chi-square PDF
NormalInvChiSq = function(mu,phi,mu0,kappa0,phi0,nu0) 
{
  z = sqrt(phi)^(-1) * phi^(-(1+nu0/2)) * 
      exp( - (nu0*phi0 + kappa0*(mu0-mu)^2 ) / (2*phi) )
  return( z )
}
 

muVec  = seq(20, 36, by=0.2 )
phiVec = seq(0, 80, by=0.5 )
xy = outer( muVec , phiVec )
z0 = matrix( 0 , nrow=length(muVec) , ncol=length(phiVec) )
for (i in 1:length(muVec))
{
  for (j in 1:length(phiVec))
  {
    z0[i,j] = NormalInvChiSq(muVec[i],phiVec[j],mu0,kappa0,phi0,nu0)
  }
}
# windows( 9 , 9 )
filled.contour( muVec , phiVec , z0 , nlevels=20 ) 

muVec  = seq(20, 36, by=0.2 )
phiVec = seq(0, 80, by=0.5 )
xy = outer( muVec , phiVec )
z1 = matrix( 0 , nrow=length(muVec) , ncol=length(phiVec) )
for (i in 1:length(muVec))
{
  for (j in 1:length(phiVec))
  {
    z1[i,j] = NormalInvChiSq(muVec[i],phiVec[j],mu1,kappa1,phi1,nu1)
  }
}
# windows( 9 , 9 )
filled.contour( muVec , phiVec , z1 , nlevels=30 )

####
####  Prior and Posterior for mu and sigma 
####

muVec = seq( 20 , 38 , 0.1 )
phiVec = seq(0.5, 80, by=0.25 )
muPrior = (1/sqrt(phi0/kappa0)) * 
                  dt( (muVec-mu0)/(sqrt(phi0/kappa0)) , nu0 )
muPost  = (1/sqrt(phi1/kappa1)) * 
                  dt( (muVec-mu1)/(sqrt(phi1/kappa1)) , nu1 )

logConst = (nu0/2)*log(nu0/2) + (nu0/2)*log(phi0) - lgamma(nu0/2)
phiPrior = exp( logConst - ((nu0/2)+1)*log(phiVec) - nu0*phi0/(2*phiVec) )

logConst = (nu1/2)*log(nu1/2) + (nu1/2)*log(s2) - lgamma(nu1/2)
phiPost = exp( logConst + (-(nu1/2+1))*log(phiVec) - nu1*phi1/(2*phiVec) )

# windows( 9 , 9 )
par( mfrow=c(2,1) )
plot( muVec , muPrior , type="l" , col="red" , ylim=c(0,0.4) , lwd=2 )
lines( muVec , muPost , col="blue" , lwd=2 )
abline( h=0 , col="gray")
plot( phiVec , phiPrior , type="l" , col="red" , ylim=c(0,0.03) , lwd=2 ) 
lines( phiVec , phiPost , type="l" , col="blue" , lwd=2 )
abline( h=0 , col="gray" )

