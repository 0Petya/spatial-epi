##########  FUNCTIONS ###############################################
pBeta = function(p,a,b) 
{ 
  z = 0
  if (p >= 0 && p <= 1) {  z = p^(a-1) * (1-p)^(b-1) }
  return(z)
}

pLik = function(x,n,p) { p^x * (1-p)^(n-x) }

pProp = function(p)  
{ 
  rnorm(n=1,mean=p,sd=0.05) 
}

####################################################################

n = 6
x = 1
a = 2 
b = 2
set.seed(123)
nsim = 10000
pvec = rep(0,nsim)
pvec[1] = 0.50
acc = 0
for (i in 2:nsim)
{
  pp = pProp( pvec[i-1] )
  pc = pvec[i-1]
  num = pBeta(pp,a,b)*pLik(x,n,pp)
  den = pBeta(pc,a,b)*pLik(x,n,pc)
  alpha = num/den
  probAcc = min(1,alpha)
  u = runif(1)
  if ( u < probAcc ) {  pvec[i] = pp; acc = acc + 1  }
  else pvec[i] = pc
}

# windows( 12 , 12 )
par( mfrow=c(4,1) , mar=c(3,3,1,1) )
plot( pvec , type="l" )
plot( pvec , type="l" , xlim=c(0,100))
hist( pvec , xlim=c(0,1) , breaks=seq(0,1,0.02) , col="red" )
plot( (0:100)/100 , pBeta( (0:100)/100 , a+x , b+n-x ) , type="l" )
