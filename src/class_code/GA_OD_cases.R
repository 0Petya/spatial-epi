library( tmap )
library( tmaptools )
library( sf )
library( dplyr )
library( ggplot2 )
library( sp )
library( reshape )   
library( grid )
library( SpatialEpi )
library( spdep )
library( nimble )
library( coda )
library( tigris )

GA = counties( "Georgia" )
class( GA )

# setwd("C:/Users/serig/Dropbox/MyCourses/BST_5610_Fall2022")
OD = read.delim( "./data/drugOD.csv", sep = "\t" )
Pop = read.csv( "./data/georgia_population_22.csv" )

OD$Cases = as.numeric( gsub( "," , "" , OD$Cases ) )
Pop$Pop = as.numeric( gsub( "," , "" , Pop$Pop ) )

OD$County = toupper( OD$County )
GA$County  = toupper( GA$NAME )

GA1 = GA %>% 
  left_join( OD , by="County" ) %>% 
  left_join( Pop , by="County" )

GA1$ODRate = GA1$Cases / GA1$Pop

tm_shape( GA1 ) +
  tm_fill( col="ODRate" ) +
  tm_polygons()

GA_sp = as( GA1 , "Spatial" )
class( GA_sp )
GA.nb = poly2nb( GA_sp )


#GA.net = nb2lines( PA.nb , coords=coordinates(PA_sp) )
#class( PA.net )

#PA.nb2 = poly2nb( PA_sp , queen=FALSE )
#PA.net2 = nb2lines( PA.nb2 , coords=coordinates(PA_sp) )

k = nrow( GA1 ) 
num = rep(0,k)
for (i in 1:k) num[i] = length( GA.nb[[i]] )
adj = c()
for (i in 1:k) adj = c(adj,GA.nb[[i]] )
L = length(adj)

################################################################################
####
####  Model 1:  v_i  only
####
################################################################################

spatialCode = nimbleCode({
  for ( i in 1:k )
  {
    y[i] ~ dpois( pop[i]*theta[i] )
  }
  for ( i in 1:k )
  {
    log(theta[i])  <-  mu + v[i]
  }
  mu ~ dnorm( -3 , 0.001 )
  tauv ~ dgamma( 1 , 0.00001 )
  for (i in 1:L)
  {
    weights[i]  <-  1
  }
  for ( i in 1:k )
  {
    v[i] ~ dnorm( 0 , tauv)
  }
})

constants = list( k=k , L=L , num=num , adj=adj , pop=GA1$Pop  )
data = list( y=GA1$Cases )
inits = list( mu=-2 , tauv=0.1 , v=rep(0,k) )

GA.spatial.model = nimbleModel( code=spatialCode , 
                                constants=constants , 
                                data=data ,
                                inits=inits )

compile.GA.spatial.model = compileNimble( GA.spatial.model )

GA.spatial.model.Conf = configureMCMC( GA.spatial.model , print = TRUE , 
                                       enableWAIC = TRUE , thin=10 )

GA.spatial.model.Conf$addMonitors(c("mu","tauv","theta"))

GA.spatial.model.MCMC = buildMCMC( GA.spatial.model.Conf )

compile.GA.spatial.MCMC = compileNimble( GA.spatial.model.MCMC )

niter = 210000
nburn =  10000
set.seed(1)

start.time = proc.time()
samples.spatial = runMCMC( compile.GA.spatial.MCMC, niter = niter, 
                           nburnin = nburn, inits = inits, nchains = 1, 
                           samplesAsCodaMCMC = TRUE , WAIC=TRUE)
stop.time = proc.time()
time.elapsed = stop.time - start.time
print( time.elapsed )

waic1 = samples.spatial$WAIC


samples1 = samples.spatial[[1]]
plot( samples1[ , "mu" ] )
plot( samples1[ , "tauv" ] )

plot( samples1[ , "theta[1]"])

# windows( 9 , 12 )
par( mfrow=c(3,1) )
ts.plot( samples1[ , 'mu'], xlab = 'iteration', col="red" , lwd=1.5 ,
         ylab = expression(mu), main = expression(beta0) )
ts.plot( samples1[ , 'tauv'], xlab = 'iteration', col="blue" , lwd=1.5 ,
         ylab = expression(tauv), main = expression(tauv) )
ts.plot( as.numeric( samples1[ , 'theta[1]' ] ) , main = 'theta[1]' , 
         ylab='theta[1]' , ylim=c(0,0.030) , col="brown" )

head( samples1 )
GA1$UH.model = apply( samples1 , 2 , mean )[3:161]

p1 = tm_shape( GA1 ) + 
  tm_fill( col="ODRate" ) +
  tm_borders()

p2 = tm_shape( GA1 ) +
  tm_fill( col="UH.model" , breaks=seq(0,0.005,0.001) ) + 
  tm_borders()

grid.newpage()
pushViewport(viewport(layout=grid.layout(1,2)))
print(p1, vp=viewport(layout.pos.col = 1, height = 5))
print(p2, vp=viewport(layout.pos.col = 2, height = 5))

plot( GA1$ODRate , GA1$UH.model )
abline( a=0 , b=1 )

################################################################################
####
####  Model 3:  v_i  and  s_i
####
################################################################################

spatialCode = nimbleCode({
  for ( i in 1:k )
  {
    y[i] ~ dpois( pop[i]*theta[i] )
  }
  for ( i in 1:k )
  {
    log(theta[i])  <-  mu + s[i] + v[i]
    ypred[i] ~ dpois( pop[i]*theta[i] ) 
  }
  mu ~ dnorm( -4 , 0.01 )
  taus ~ dgamma( 1 , 0.0001 )
  tauv ~ dgamma( 1 , 0.0001 )
  for (i in 1:L)
  {
    weights[i]  <-  1
  }
  s[1:k] ~ dcar_normal(adj[1:L],weights[1:L],num[1:k],taus,zero_mean=1)
  for ( i in 1:k )
  {
    v[i] ~ dnorm( 0 , tauv)
  }
})

constants = list( k=k , L=L , num=num , adj=adj , pop=GA1$Pop  )
data = list( y=GA1$Cases )
inits = list( mu=-2 , taus=0.1 , tauv=0.1 , s=rep(0,k) , v=rep(0,k) , ypred=GA1$Cases )

GA.spatial.model = nimbleModel( code=spatialCode , 
                                constants=constants , 
                                data=data ,
                                inits=inits )

compile.GA.spatial.model = compileNimble( GA.spatial.model )

GA.spatial.model.Conf = configureMCMC( GA.spatial.model , print = TRUE , 
                                       enableWAIC = TRUE , thin=100 )

GA.spatial.model.Conf$addMonitors( c("mu","taus","tauv","theta","ypred") )

GA.spatial.model.MCMC = buildMCMC( GA.spatial.model.Conf )

compile.GA.spatial.MCMC = compileNimble( GA.spatial.model.MCMC )

niter =  100000
nburn =   10000
set.seed(1)

start.time = proc.time()
samples.spatial = runMCMC( compile.GA.spatial.MCMC, niter = niter, 
                           nburnin = nburn, inits = inits, nchains = 1, 
                           samplesAsCodaMCMC = TRUE , WAIC=TRUE)
stop.time = proc.time()
time.elapsed = stop.time - start.time
print( time.elapsed )

waic3 = samples.spatial$WAIC

samples1 = samples.spatial[[1]]
# write.csv( samples1 , "./data/GA_OD_Modele3.csv" , row.names = FALSE )

samples.model3 = samples1

plot( samples1[ , "mu" ] )
plot( samples1[ , "tauv" ] )
plot( samples1[ , "taus" ] )
plot( samples1[ , "theta[1]"])

theta.matrix = samples1[ , 4:162 ]
ypred.matrix = samples1[ , 163:321 ]

GA1$CH.UH.model = apply( theta.matrix , 2 , mean )

# windows( 12 , 12 )
par( mfrow=c(3,2) )
ts.plot( samples1[ , 'mu'], xlab = 'iteration', col="red" , lwd=1.5 ,
         ylab = expression(mu), main = expression(beta0) )
ts.plot( samples1[ , 'taus'], xlab = 'iteration', col="blue" , lwd=1.5 ,
         ylab = expression(taus), main = expression(taus) )
ts.plot( samples1[ , 'tauv'], xlab = 'iteration', col="blue" , lwd=1.5 ,
         ylab = expression(tauv), main = expression(tauv) )
ts.plot( theta.matrix[,1] , main = 'theta for Appling Co.' , 
         ylab='theta[1]' , ylim=c(0,0.01) , col="brown" )
hist( ypred.matrix[,1] , main = "Predicted for Co 1" )
  abline( v=GA1$Cases[1] )
hist( ypred.matrix[,2] , main = "Predicted for Co 2" )
  abline( v=GA1$Cases[2] )
  
p1 = tm_shape( GA1 ) + 
  tm_fill( col="ODRate" ) +
  tm_borders()

p2 = tm_shape( GA1 ) +
  tm_fill( col="UH.model" , breaks=seq(0,0.005,0.001) ) + 
  tm_borders()

p3 = tm_shape( GA1 ) +
  tm_fill( col="CH.UH.model" , breaks=seq(0,0.005,0.001) ) + 
  tm_borders()

# windows( 12 , 9 )
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,3)))
print(p1, vp=viewport(layout.pos.col = 1, height = 5))
print(p2, vp=viewport(layout.pos.col = 2, height = 5))
print(p3, vp=viewport(layout.pos.col = 3, height = 5))

plot( GA1$ODRate , GA1$CH.UH.model , xlim=c(0,0.005) , ylim=c(0,0.005) , cex=0.005*sqrt(GA1$Pop) )
abline( a=0 , b=1 )

#########################################################################
####
####  Which Dominates?  Uncorrelated or Correlated Heterogeneity
####
#########################################################################

# windows( 9 , 12 )
par( mfrow=c(2,1) )
taus = as.numeric( samples1[,"taus"] )
tauv = as.numeric( samples1[,"tauv"] )
taushat = mean( taus )
tauvhat = mean( tauv )
sigmashat = mean( 1/sqrt(taus) )
sigmavhat = mean( 1/sqrt(tauv) )

plot( 1/sqrt(taus) , 1/sqrt(tauv) )
plot( 1/sqrt(taus) , 1/sqrt(tauv) , log="xy" )

print( c( sigmashat , sigmavhat ) )

#########################################################################
####
####  Outlier Detection via Residuals
####
#########################################################################
  
GA1$CH.UH.PredCases = apply( ypred.matrix , 2 , mean )
GA1$CH.UH.SD.PredCases = apply( ypred.matrix , 2 , sd )
GA1$residuals = ( GA1$Cases - GA1$CH.UH.PredCases ) / GA1$CH.UH.SD.PredCases

plot( GA1$Cases , GA1$CH.UH.PredCases )
  abline( a=0 , b=1 )

plot( GA1$residuals )

tm_shape( GA1 ) +
  tm_borders() +
  tm_polygons( col="residuals" )


#########################################################################
####
####  Excedence Probability
####
#########################################################################

threshold = sum(GA1$Cases) / sum(GA1$Pop)
GA1$exceed.mean = apply( theta.matrix > threshold , 2 , mean )

p.exceed1 = tm_shape( GA1 ) +
  tm_borders() +
  tm_fill( col="exceed.mean" ) +
  tm_layout( title="Exceedence Probability for Threshold = Avg for State" ,
             title.position = c('left', 'bottom'))

threshold = 0.004    #### Arbitrary threshold
GA1$exceed.004 = apply( theta.matrix > threshold , 2 , mean )

p.exceed2 = tm_shape( GA1 ) +
  tm_borders() +
  tm_fill( col="exceed.004" ) +
  tm_layout( title="Exceedence Probability for Threshold = 0.04" ,
             title.position = c('left', 'bottom'))

# windows( 12 , 9 )
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,2)))
print(p.exceed1, vp=viewport(layout.pos.col = 1, height = 5))
print(p.exceed2, vp=viewport(layout.pos.col = 2, height = 5))

