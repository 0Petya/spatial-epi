x = seq( 0 , 10 , 0.1 )
shape = 4
rate = 1.5
y = dgamma( x , shape , rate )
f = function(x) return( dgamma(x,shape,rate) - a )

# windows( 12 , 12 )
par( mfrow=c(2,1) )

L1 = qgamma( 0.025 , shape , rate )
U1 = qgamma( 0.975 , shape , rate )
plot( x , y , type="l" )
abline( h=0 )
xx = seq( L1 , U1 , 0.1 )
polygon( c( L1 , xx              , rev(xx) ,            L1) , 
         c( 0  , dgamma(xx,shape,rate) , rep(0,length(xx)) , 0 ) , col="gray" ) 

plot( x , y , type="l" )
a = 0.2       ####  Try a = 0.044
abline( h=a )
abline( h=0 )
L = uniroot(  f , lower = 0.1, upper = 2.2 )$root
U = uniroot(  f , lower = 2.2, upper = 10 )$root
containment = pgamma( U , shape , rate ) - pgamma( L , shape , rate )
text( 8 , 0.3 , paste0("Containment = ",round(containment,3)) , 
      adj=c(0,0) )
lines( c(L,L) , c(0,dgamma(L,shape,rate)) )
lines( c(U,U) , c(0,dgamma(L,shape,rate)) )
polygon( c(L,L,U,U,L) , c(0,dgamma(L,shape,rate),dgamma(U,shape,rate),0,0) , 
         col="gray" )


