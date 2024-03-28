a = 4
b = 16

xVec = seq( 0 , 1 , 0.001 )
yVec = dbeta( xVec , a , b )
plot( xVec , yVec , type="l" , col="red" , lwd=2 )

plot( xVec , gamma(a+b)/(gamma(a)*gamma(b)) * xVec^(a-1) * (1-xVec)^(b-1), 
      type="l" , col="blue" , lwd=2 )

plot( xVec ,  xVec^(a-1) * (1-xVec)^(b-1) , type="l" , col="blue" , lwd=2 )


####  Prior distribution
a0 = 5
b0 = 5

####  Data
n = 10
x = 7 

####  Posterior distribution
a1 = a0 + x
b1 = b0 + n - x

####  Plot prior and posterior
# windows( 12 , 12 )
par( mfrow=c(2,1) )

priorVec = dbeta( xVec , a0 , b0 )
plot( xVec , priorVec , type="l" , col="red" , lwd=2 )

postVec = dbeta( xVec , a1 , b1 )
plot( xVec , postVec , type="l" , col="green3" , lwd=2 )



