# setwd("C:/Users/serig/Dropbox/MyCourses/BST_5610_Fall2022")

library( nimble )
library( spdep )
library( dplyr )
library( tmap )
library( MASS )
library( rgdal )

data( nc.sids ) 
example( nc.sids )

nc.sids = readOGR(system.file("shapes/sids.shp", package="spData")[1])
proj4string(nc.sids) = CRS("+proj=longlat +ellps=clrk66")
row.names(nc.sids) = as.character(nc.sids$FIPS)

# windows( 14 , 8 )
tm_shape( nc.sids ) +
  tm_polygons( col="gray" ) +
  tm_text( "NAME" , "lon" , "lat" , size=0.8 , col="black")

 
###################################################
####  Model 1
###################################################

nc.data = dplyr::select( nc.sids@data , NAME , FIPSNO , BIR74 , SID74 , NWBIR74 , 
                         BIR79 , SID79 , NWBIR79 , x , y , lon , lat )

nc.data$SIDrate74 = nc.data$SID74/nc.data$BIR74
nc.data$SIDrate79 = nc.data$SID79/nc.data$BIR79
nc.sids@data = nc.data    ## This replaces the data with the smaller set.

# windows( 14 , 8 )
tm_shape( nc.sids ) +
  tm_polygons( col="SIDrate79" , border.col="black" ) +
  tm_text( "NAME" , "lon" , "lat" , size=0.8 , col="black" ) +
  tm_layout("SIDS Raw Rates", title.size=2)

data.frame( nc.data$NAME , nc.data$SIDrate79 )

###################################################
####  Model 2:
####  Hierarchical Model
####  non spatial
###################################################

nc.sids.Code = nimbleCode(
{
  for (i in 1:k) {
     theta[i] ~ dnorm( mu , tau )      
     y[i] ~ dpois( n[i]*exp(theta[i]) )
  }
  # prior for hyperparameters
  mu ~  dnorm( -6 , 0.001 )
  tau ~ dgamma( 1, .001)
}
)

k = 100      ## There are 100 counties
n  = nc.sids$BIR79
y  = nc.sids$SID79

nc.sids.Consts = list( k = k , n = n )
nc.sids.Data = list( y = y )
nc.sids.Inits = list( mu = -2 , tau = 12 , theta=rep(0.002,k) )

nc.sids.Model = nimbleModel( nc.sids.Code, 
                             data = nc.sids.Data, 
                             constants = nc.sids.Consts, 
                             inits = nc.sids.Inits )

compile.nc.sids.Model = compileNimble(nc.sids.Model)
str( compile.nc.sids.Model )

nc.sids.Conf = configureMCMC(nc.sids.Model, print = TRUE)
nc.sids.Conf$addMonitors(c("mu","tau","theta"))

nc.sids.MCMC = buildMCMC( nc.sids.Conf )
compile.nc.sids.MCMC = compileNimble(nc.sids.MCMC, project = nc.sids.Model)

niter = 100000
nburn =  10000
set.seed(1)

inits = list( mu=-2 , tau=12 , theta=rep(0.002,k) )
start.time = proc.time()
samples = runMCMC( compile.nc.sids.MCMC, niter = niter, nburnin = nburn,
                   inits = inits, nchains = 1, samplesAsCodaMCMC = TRUE )
stop.time = proc.time()
time.elapsed = stop.time - start.time
print( time.elapsed )

samples[1:10,]

# windows( 12 , 9 )
samples1 = samples
par(mfrow = c(2, 2), mai = c(.6, .5, .4, .1), mgp = c(1.8, 0.7, 0))
ts.plot(samples1[ , 'mu'], xlab = 'iteration', col="red" , lwd=1.5 ,
        ylab = expression(mu), main = expression(mu))
ts.plot(samples1[ , 'tau'], xlab = 'iteration', col="red" , lwd=1.5 ,
        ylab = expression(tau), main = expression(tau))
ts.plot(samples1[ , 'theta[1]'], xlab = 'iteration', col="blue" , lwd=1.5 ,
        ylab = expression(theta[1]), main = expression(theta[1]))
ts.plot(samples1[ , 'theta[2]'], xlab = 'iteration', col="blue" , lwd=1.5 ,
        ylab = expression(theta[2]), main = expression(theta[2]))

eta = exp( apply( samples[,3:102] , 2 , mean ) )

nc.sids$HB.est.rate = eta 

# windows( 14 , 8 )
tm_shape( nc.sids ) +
  tm_polygons( col="SIDrate79" , border.col="black" , 
               breaks=seq(0,0.007,0.001) ) +
  tm_text( "NAME" , "lon" , "lat" , size=0.8 , col="black" )

# windows( 14 , 8 )
tm_shape( nc.sids ) +
  tm_polygons( col="HB.est.rate" , border.col="black" , 
                breaks=seq(0,0.007,0.001)) +
  tm_text( "NAME" , "lon" , "lat" , size=0.8 , col="black" )  +
  tm_layout("SIDS Rates from \nHierarchical Model", title.size=1.8)
           
# windows( 9 , 9 )
plot( nc.sids$SIDrate79 , nc.sids$HB.est.rate , xlim=c(0,0.0062) , 
      ylim=c(0,0.005) , cex=0.2+nc.sids$BIR79/20000 , pch=19 , col="red" ,
      asp=1 , xlab="Raw Rates" , ylab="Hierarchical Model")
abline( a=0 , b=1 ,col="gray")
    
data.frame( nc.sids$NAME , nc.sids$SIDrate79 , nc.sids$HB.est.rate )

###################################################
####  Model 3a:
####  Hierarchical Model
####  SPATIAL 
####  only correlated heterogeneity
###################################################
class( nc.sids )
nc.sids.nb = poly2nb( nc.sids )
nc.sids.net = nb2lines( nc.sids.nb , coords=coordinates(nc.sids) )

# windows( 14 , 8 )
tm_shape( nc.sids ) +
  tm_borders( col="darkgray" ) +
  tm_shape( nc.sids.net ) +
  tm_lines( col="darkgreen" , lwd=2 )

nc.sids.lw = nb2listw( nc.sids.nb )

k = nrow( nc.sids@data )
num = rep(0,k)
for (i in 1:k) num[i] = length( nc.sids.lw$neighbours[[i]] )
adj = c()
for (i in 1:k) adj = c(adj,nc.sids.lw$neighbours[[i]] )
L = length(adj)

spatialCode = nimbleCode({
  mu ~ dflat()
  tau ~ dgamma( 1 , 0.001 )
  for (i in 1:L)
    weights[i]  <-  1
  s[1:k] ~ dcar_normal(adj[1:L],weights[1:L],num[1:k],tau,zero_mean=1)
  for (i in 1:k) {
    log(theta[i])  <-  mu + s[i]
    y[i] ~ dpois( BIR79[i]*theta[i] )
  }  
})

constants = list( k=k , L=L , num=num , adj=adj , BIR79=nc.sids$BIR79 )
data = list( y=nc.sids$SID79 )
inits = list( mu=0 , tau=1 , s=rep(0,k) )

nc.sids.spatial.model = nimbleModel( code=spatialCode , 
                                     constants=constants , 
                                     data=data ,
                                     inits=inits )

compile.nc.sids.spatial.model = compileNimble( nc.sids.spatial.model )

nc.sids.spatial.model.Conf = configureMCMC( nc.sids.spatial.model , print = TRUE )

nc.sids.spatial.model.Conf$addMonitors(c("mu","tau","theta"))

nc.sids.spatial.model.MCMC = buildMCMC( nc.sids.spatial.model.Conf )

compile.nc.sids.spatial.MCMC = compileNimble( nc.sids.spatial.model.MCMC )

niter = 110000
nburn =  10000
set.seed(1)

start.time = proc.time()
samples.spatial = runMCMC( compile.nc.sids.spatial.MCMC, niter = niter, 
                           nburnin = nburn, inits = inits, nchains = 1, 
                           samplesAsCodaMCMC = TRUE )
stop.time = proc.time()
time.elapsed = stop.time - start.time
print( time.elapsed )

windows( 12 , 9 )
samples1 = samples.spatial
par(mfrow = c(2, 2), mai = c(.6, .5, .4, .1), mgp = c(1.8, 0.7, 0))
ts.plot(samples1[ , 'mu'], xlab = 'iteration', col="red" , lwd=1.5 ,
        ylab = expression(mu), main = expression(mu))
ts.plot(samples1[ , 'tau'], xlab = 'iteration', col="darkgreen" , lwd=1.5 ,
        ylab = expression(tau), main = expression(tau))
ts.plot(samples1[ , 'theta[1]'], xlab = 'iteration', col="blue" , lwd=1.5 ,
        ylab = expression(theta[1]), main = expression(theta[1]))
ts.plot(samples1[ , 'theta[2]'], xlab = 'iteration', col="blue" , lwd=1.5 ,
        ylab = expression(theta[1]), main = expression(theta[1]))

theta.HB.spatial = apply( samples1[,3:102] , 2 , mean )
nc.sids$HB.spatial.est.rate = theta.HB.spatial

windows( 9 , 9 )
plot( nc.sids$HB.spatial.est.rate , nc.sids$HB.est.rate , xlim=c(0,0.0062) , 
      ylim=c(0,0.005) , cex=0.1+nc.sids$BIR79/20000 , pch=19 , col="red" ,
      asp=1 , xlab="Hierarchical Model" , ylab="Hierarchical Spatial Model")
abline( a=0 , b=1 )

windows( 9 , 9 )
plot( nc.sids$SIDrate79 ,  nc.sids$HB.spatial.est.rate , xlim=c(0,0.0062) , 
      ylim=c(0,0.005) , cex=0.1+nc.sids$BIR79/20000 , pch=19 , col="red" ,
      asp=1 , xlab="Raw Rates" , ylab="Hierarchical Spatial Model")
abline( a=0 , b=1 )

windows( 14 , 8 )
tm_shape( nc.sids ) +
  tm_polygons( col="HB.spatial.est.rate" , border.col="black" , 
               breaks=seq(0,0.007,0.001)) +
  tm_text( "NAME" , "lon" , "lat" , size=0.8 , col="black" )  +
  tm_layout("SIDS Rates from \nCAR Model", title.size=1.8)


ThreeRateEstimates = 
  data.frame( nc.sids@data$NAME , nc.sids$SIDrate79 , nc.sids$HB.est.rate , nc.sids$HB.spatial.est.rate )

###################################################
####  Model 3b:
####  Hierarchical Model
####  SPATIAL 
####  both correlated and uncorrelated heterogeneity
###################################################
class( nc.sids )
nc.sids.nb = poly2nb( nc.sids )
nc.sids.net = nb2lines( nc.sids.nb , coords=coordinates(nc.sids) )

windows( 14 , 8 )
tm_shape( nc.sids ) +
  tm_borders( col="darkgray" ) +
  tm_shape( nc.sids.net ) +
  tm_lines( col="darkgreen" , lwd=2 )

nc.sids.lw = nb2listw( nc.sids.nb )

k = nrow( nc.sids@data )
num = rep(0,k)
for (i in 1:k) num[i] = length( nc.sids.lw$neighbours[[i]] )
adj = c()
for (i in 1:k) adj = c(adj,nc.sids.lw$neighbours[[i]] )
L = length(adj)

spatialCode = nimbleCode({
  mu ~ dflat()
  tau ~ dgamma( 1 , 0.001 )
  tau1 ~ dgamma( 1 , 0.001 )
  for (i in 1:L)
    weights[i]  <-  1
  s[1:k] ~ dcar_normal(adj[1:L],weights[1:L],num[1:k],tau,zero_mean=1)
  for (i in 1:k) {
    log(theta[i])  <-  mu + s[i] + v[i]
    y[i] ~ dpois( BIR79[i]*theta[i] )
    v[i] ~ dnorm( 0 , tau1 )
  }  
})

constants = list( k=k , L=L , num=num , adj=adj , BIR79=nc.sids$BIR79 )
data = list( y=nc.sids$SID79 )
inits = list( mu=0 , tau=1 , s=rep(0,k) , tau1=1 )

nc.sids.spatial.model = nimbleModel( code=spatialCode , 
                                     constants=constants , 
                                     data=data ,
                                     inits=inits )

compile.nc.sids.spatial.model = compileNimble( nc.sids.spatial.model )

nc.sids.spatial.model.Conf = configureMCMC( nc.sids.spatial.model , 
                                            thin=100 , print = TRUE )

nc.sids.spatial.model.Conf$addMonitors(c("mu","tau","tau1","theta"))

nc.sids.spatial.model.MCMC = buildMCMC( nc.sids.spatial.model.Conf )

compile.nc.sids.spatial.MCMC = compileNimble( nc.sids.spatial.model.MCMC )

set.seed(123)
niter = 10100000
nburn =   100000


start.time = proc.time()
samples.spatial = runMCMC( compile.nc.sids.spatial.MCMC, niter = niter, 
                           nburnin = nburn, inits = inits, nchains = 1, 
                           samplesAsCodaMCMC = TRUE )
stop.time = proc.time()
time.elapsed = stop.time - start.time
print( time.elapsed )

samples1 = samples.spatial

windows( 12 , 9 )

par(mfrow = c(2, 2), mai = c(.6, .5, .4, .1), mgp = c(1.8, 0.7, 0))
ts.plot(samples1[ , 'mu'], xlab = 'iteration', col="red" , lwd=1.5 ,
        ylab = expression(mu), main = expression(mu))
ts.plot(samples1[ , 'tau'], xlab = 'iteration', col="darkgreen" , lwd=1.5 ,
        ylab = expression(tau), main = expression(tau))
ts.plot(samples1[ , 'theta[1]'], xlab = 'iteration', col="blue" , lwd=1.5 ,
        ylab = expression(theta[1]), main = expression(theta[1]))
ts.plot(samples1[ , 'tau1'], xlab = 'iteration', col="blue" , lwd=1.5 ,
        ylab = expression(theta[1]), main = expression(tau[1]))


theta.HB.spatialB = apply( samples1[,4:103] , 2 , mean )
nc.sids$HB.spatialB.est.rate = theta.HB.spatialB

windows( 9 , 9 )
plot( nc.sids$HB.spatial.est.rate , nc.sids$HB.spatialB.est.rate , xlim=c(0,0.0062) , 
      ylim=c(0,0.005) , cex=0.1+nc.sids$BIR79/20000 , pch=19 , col="red" ,
      asp=1 , xlab="Hierarchical Spatial Model (s only)" , 
      ylab="Hierarchical Spatial Model (both s and v)" )
abline( a=0 , b=1 )

windows( 9 , 9 )
plot( nc.sids$SIDrate79 ,  nc.sids$HB.spatial.est.rate , xlim=c(0,0.0062) , 
      ylim=c(0,0.005) , cex=0.1+nc.sids$BIR79/20000 , pch=19 , col="red" ,
      asp=1 , xlab="Raw Rates" , ylab="Hierarchical Spatial Model")
abline( a=0 , b=1 )


###########################################################################
####    C O N T O U R    P L O T S
###########################################################################

windows( 14 , 8 )
tm_shape( nc.sids ) +
  tm_polygons( col="HB.spatialB.est.rate" , border.col="black" , 
               breaks=seq(0,0.007,0.001)) +
  tm_text( "NAME" , "lon" , "lat" , size=0.8 , col="black" )  +
  tm_layout("SIDS Rates from \nCAR Model (s and v)", title.size=1.8)

windows( 10 , 10 )
tauSample = samples1[,"tau"]
tau1Sample = samples1[,"tau1"]
z = kde2d( tauSample , tau1Sample)
plot( as.numeric(tauSample) , as.numeric(tau1Sample) , pch="." , log="xy" )
contour( z , drawlabels=FALSE , nlevels=8, add=TRUE)

windows( 9 , 16 )
par( mfrow=c(2,1) , mar=c(4,4,1,2) )
logtauSample = log10( samples1[,"tau"] )
logtau1Sample = log10( samples1[,"tau1"] )
tauhat = 10^mean( logtauSample )
tau1hat = 10^mean( logtau1Sample )
tauhat
tau1hat
bothtaus = data.frame( samples1[,"tau"] , samples1[,"tau1"] )
#write.csv( bothtaus , "NC-SIDS-tau-tau1-posterior.csv"  )
zlogtau = kde2d( logtauSample , logtau1Sample , n=50 )
plot( as.numeric(tauSample) , as.numeric(tau1Sample) , asp=1 , pch="." ,
      xlim=c(0,5000) , ylim=c(0,5000) , xlab="tau" , ylab="tau1" , cex=0.7 )
plot( as.numeric(logtauSample) , as.numeric(logtau1Sample) , pch="." , asp=1 , cex=0.7 , 
      xaxt="n" , yaxt="n" , ylim=c(0.5,4) , xlim=c(0.5,4) , xlab="tau" , ylab="tau1" )
  ax = c(1,2,5,10,20,50,100,200,500,1000,2000,5000,10000)
  ax = c(1 , 10 , 100 , 1000 , 10000 )
  axis( 1, at=log10(ax) , labels=ax )
  axis( 2, at=log10(ax) , labels=ax )
contour( zlogtau , drawlabels=FALSE , nlevels=50 , lwd=1 , col="red" , add=TRUE)

#zvar = kde2d( 1/logtauSample , 1/logtau1Sample , n=50 )

windows( 9 , 16 )
par( mfrow=c(2,1) , mar=c(4,4,1,2) )
sigmasq = 1/as.numeric( tauSample )
sigmasq1 = 1/as.numeric( tau1Sample )
sd0 = sqrt( sigmasq )
sd1 = sqrt( sigmasq1 ) 
#zsd = kde2d( sd0 , sd1 , n=50 )
plot( sd0 , sd1 , asp=1 , pch="." ,
      xlim=c(0,0.8) , ylim=c(0,0.8) , xlab="sigma" , ylab="sigma1" , cex=0.7 )
  contour( zsd , drawlabels=FALSE , nlevels=40 , lwd=1 , col="green" , add=TRUE )

#zlsd = kde2d( log(sd0) , log(sd1) , n=50 )
plot( as.numeric( log( sd0 ) ), as.numeric( log( sd1 ) ) , xlab="sigma" , ylab="sigma2" , 
      pch="." , cex=0.7 , xlim=c(-5,0) , ylim=c(-5,0) , xaxt="n" , yaxt="n" )
  ax = c( 0.0001 , 0.001 , 0.01 , 0.1 , 1 , 10 , 100 , 1000 , 10000 )
  axis( 1, at=log10(ax) , labels=ax )
  axis( 2, at=log10(ax) , labels=ax )
  contour( zlsd , drawlabels=FALSE , nlevels=40 , lwd=1 , col="green" , add=TRUE )
  
windows( 9 , 12 )
par( mfrow=c(2,2) , byrow=TRUE )
hist( sd0 , col="green" )
hist( sd1 , col="green" )
hist( logtauSample , col="red" )
hist( logtau1Sample , col="red" )
  
  
  
  
FourRateEstimates = 
  data.frame( nc.sids@data$NAME , nc.sids$SIDrate79 , nc.sids$HB.est.rate , nc.sids$HB.spatial.est.rate )

windows( 9 , 9 )
plot( nc.sids$SIDrate79 ,  nc.sids$HB.spatial.est.rate , xlim=c(0,0.0062) , 
      ylim=c(0,0.005) , cex=0.1+nc.sids$BIR79/20000 , pch=19 , col="red" ,
      asp=1 , xlab="Raw Rates" , ylab="Hierarchical Spatial Model")
abline( a=0 , b=1 )






p <- ggscatter( data.frame(sd0,sd1) )
ggMarginal(p, type = "density")

####  TEST

library("ggExtra")
library("ggpubr")
p <- ggscatter(iris, x = "Sepal.Length", y = "Sepal.Width",
               color = "Species", palette = "jco",
               size = 3, alpha = 0.6)
ggMarginal(p, type = "density")
# Change marginal plot type
ggMarginal(p, type = "boxplot")
p

