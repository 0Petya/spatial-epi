#---------------------------------------------------------------------
#-
#-  Code from Textbook: Intro to R for Spatial Analysis & Mapping, 2e
#-
#-  Chapter 7: Spatial Attribute Analysis with R
#-
#---------------------------------------------------------------------
#-
#-  Review Lines 1 - 200
#-
#---------------------------------------------------------------------

# setwd( "C:/Users/serig/Dropbox/MyCourses/BST_5610_Fall2022" )
#library( GISTools )
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

data( pennLC_sf )   ## Pennsylvania lung cancer data - part of SpatielEpi package
class( pennLC_sf )

## Map the counties in penn.state.utm
tm_shape( pennLC_sf ) + tm_borders()

pop.df = pennLC_sf %>% 
           group_by( county ) %>% 
           summarize( pop = sum(population) , smk = mean(smoking) , cases=sum(cases) )
pop = pop.df$pop
county = pop.df$county
smk = pop.df$smk
cases = pop.df$cases
LCrate = cases/pop
geometry = pennLC_sf$geometry

PA = st_sf( county , pop , cases , smk , LCrate , geometry[seq(from=1,by=16,length=67)] )

names(PA) = c("county","pop","cases","smk","LCrate","geometry")
st_geometry(PA) <- "geometry" 
  
# windows( 9 , 6 )
tm_shape( PA ) +
  tm_polygons( col="smk" , title="% of Population" )

# windows( 9 , 6 )
tm_shape( PA ) +
  tm_polygons( col="LCrate" , title="Lung Cancer Rate" )

##  Are neighboring counties smoking percentages independent?
##  Or correlated?   
##  Correlation in space is called autocorrelation.

##  "Shuffle" the smoking percentages across counties

sample( 1:6 )
PA$smk_rand1 = sample( PA$smk )
PA$smk_rand2 = sample( PA$smk )
PA$smk_rand3 = sample( PA$smk )
PA$smk_rand4 = sample( PA$smk )
PA$smk_rand5 = sample( PA$smk )

vars = sample( c("smk" , "smk_rand1" , "smk_rand2" , "smk_rand3" ,
                 "smk_rand4" , "smk_rand5" ) )
real.data.i = which( vars == "smk" )

# windows( 9 , 9 )
tm_shape( PA ) +
  tm_polygons( col=vars , legend.show=FALSE ) +
  tm_layout(title=1:6 , title.position=c("right","top") )

PA_sp = as( PA , "Spatial" )
class( PA_sp )
PA.nb = poly2nb( PA_sp )
# Queens
PA.net = nb2lines( PA.nb , coords=coordinates(PA_sp) )
class( PA.net )

PA.nb2 = poly2nb( PA_sp , queen=FALSE )
PA.net2 = nb2lines( PA.nb2 , coords=coordinates(PA_sp) )

##  Plot both networks together in different colors
# windows( 12 , 7 )
tm_shape( PA_sp ) +
  tm_borders( col="gray" ) +
  tm_shape( PA.net ) +
    tm_lines( col="blue" , lwd=2 ) +
  tm_shape( PA.net2 ) +
    tm_lines( col="orange" , lwd=3 )

k = nrow( PA ) 
num = rep(0,k)
for (i in 1:k) num[i] = length( PA.nb[[i]] )
adj = c()
for (i in 1:k) adj = c(adj,PA.nb[[i]] )
L = length(adj)

spatialCodeR = nimbleCode({
  for ( i in 1:k )
  {
    y[i] ~ dpois( pop[i]*theta[i] )
  }
  for ( i in 1:k )
  {
    log(theta[i])  <-  beta0 + beta1*x[i] + s[i] + v[i]
  }
  beta0 ~ dnorm( 0 , 10 )
  beta1 ~ dnorm( 0 , 10 )
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


spatialCode = nimbleCode({
  for ( i in 1:k )
  {
    y[i] ~ dpois( pop[i]*theta[i] )
  }
  for ( i in 1:k )
  {
    log(theta[i])  <-  beta0 + s[i] + v[i]
  }
  beta0 ~ dnorm( 0 , 10 )
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

constants = list( k=k , L=L , num=num , adj=adj , pop=pop , x=smk )
data = list( y=cases )
inits = list( beta0=0 , beta1=0 , taus=0.1 , tauv=0.1 , s=rep(0,k) , v=rep(0,k) )

PA.spatial.model = nimbleModel( code=spatialCodeR , 
                                constants=constants , 
                                data=data ,
                                inits=inits )

compile.PA.spatial.model = compileNimble( PA.spatial.model )

PA.spatial.model.Conf = configureMCMC( PA.spatial.model , print = TRUE , thin=100 )

PA.spatial.model.Conf$addMonitors(c("beta0","beta1","taus","tauv","theta"))

PA.spatial.model.MCMC = buildMCMC( PA.spatial.model.Conf )

compile.PA.spatial.MCMC = compileNimble( PA.spatial.model.MCMC )

niter = 1100000
nburn =  100000
set.seed(1)

start.time = proc.time()
samples.spatial = runMCMC( compile.PA.spatial.MCMC, niter = niter, 
                           nburnin = nburn, inits = inits, nchains = 1, 
                           samplesAsCodaMCMC = TRUE )
stop.time = proc.time()
time.elapsed = stop.time - start.time
print( time.elapsed )

# windows( 12 , 9 )
samples1 = samples.spatial
par(mfrow = c(2, 2), mai = c(.6, .5, .4, .1), mgp = c(1.8, 0.7, 0))
ts.plot(samples1[ , 'beta0'], xlab = 'iteration', col="red" , lwd=1.5 ,
        ylab = expression(beta0), main = expression(beta0))
ts.plot(samples1[ , 'taus'], xlab = 'iteration', col="blue" , lwd=1.5 ,
        ylab = expression(taus), main = expression(taus))
ts.plot(samples1[ , 'tauv'], xlab = 'iteration', col="darkgreen" , lwd=1.5 ,
        ylab = expression(tauv), main = expression(tauv))
ts.plot(samples1[ , 'beta1'], xlab = 'iteration', col="blue" , lwd=1.5 ,
        ylab = expression(beta1), main = expression(beta1))

str( samples1 )
theta.HB.spatial = apply( samples1[,5:71] , 2 , mean )

plot( LCrate , theta.HB.spatial )
abline( a=0 , b=1 )
PA$HB.spatial.est.rate = theta.HB.spatial

# windows( 9 , 6 )
tm_shape( PA ) +
  tm_polygons( col="LCrate" , title="Lung Cancer Rate" )

# windows( 9 , 6 )
tm_shape( PA ) + 
  tm_polygons( col="HB.spatial.est.rate" , title="Lung Cancer Spatial Model")
  
source("./src/class_code/plotPost.R")
source("./src/class_code/HDIofMCMC.R")
plotPost( samples1[,"beta1"] )




##  Moran's I
##  See equation (7.4) on p. 257 & PDF Notes

##  Lagged means
##  See equation (7.3) on p. 254.  Note error in book.

PA.lw = nb2listw( PA.nb2 )
class( PA.lw )
str( PA.lw )
PA.lw$neighbours
PA.lw$weights

PA$smk.lagged.means = lag.listw( PA.lw , PA$smk )


View( PA )
windows( 12 , 12 )
tm_shape( PA ) +
  tm_polygons( col=c("smk","smk.lagged.means") , title="% of Population" )

with( data.frame( PA ) ,
      {
        plot( smk , smk.lagged.means , asp=1 , 
              xlim=range(smk) , ylim=range(smk) )
        abline( a=0 , b=1 ) 
        abline( v=mean(smk) , lty=2 )
        abline( h=mean(smk.lagged.means) , lty=2 )
      })

moran.range = function( lw )
{
  wmat = listw2mat( lw )
  return( range( eigen( (wmat + t(wmat) )/2 ) $values ) )
}

moran.range( PA.lw )

moran.test( PA$smk , PA.lw , randomisation=FALSE )
##  Note: "randomisation" not "randomization"

moran.test( PA$smk , PA.lw )
##  Note: randomisation is the default

moran.mc.output = moran.mc( PA$smk , PA.lw , 100000 )
str( moran.mc.output )

windows( 9 , 6 )
moran.test.output = moran.test( PA$smk , PA.lw )
hist( moran.mc.output$res , xlim=moran.range(PA.lw ) )
 abline( v = moran.test.output$estimate[1] )

 
 

 
 
 
#------------------------------------------------------------------
#-  
#-  moran.test  and  moran.mc  require a neighborhood list
#-
#-  If you want to apply a general weight matrix  W , you need to
#-  use the "ape" package and the Moran.I() function
#-
#------------------------------------------------------------------

##  install.packages( "ape" )
library( ape )

Moran.I( penn.state.utm$smk , listw2mat( penn.state.lw ) )

moran.test( penn.state.utm$smk , penn.state.lw )

dist.mat.penn = as.matrix( dist( coordinates( penn.state.utm ) ) )
head( round( dist.mat.penn , 3 ) )

diag( dist.mat.penn ) = 1   ## Arbitrary number, but not 0
round( inv.dist.mat.penn , 3 )

W.penn = 1/dist.mat.penn
diag( W.penn ) = 1
round( W.penn , 3 )

W1.penn = W.penn/apply( W.penn , 1 , sum )  ## normalized so that rows sum to 1

Moran.I( penn.state.utm$smk , W1.penn )
Moran.I( penn.state.utm$smk , listw2mat( penn.state.lw ) )

#--------------------------------------------------------------------
#-
#-  Compute Moran's I and give choropleth map for shuffled smoking
#-  rates.
#-
#--------------------------------------------------------------------

penn.state.utm$smk_rand = sample( penn.state.utm$smk )
moran = moran.test( penn.state.utm$smk_rand , penn.state.lw )$estimate[1]
plot1 = tm_shape( penn.state.utm ) +
  tm_polygons( col="smk_rand" , legend.show=FALSE ) +
  tm_layout( paste0( "Moran's I = " , round(moran,3) ) , 
             title.position=c("right","top") )

penn.state.utm$smk.lagged.means.rand = lag.listw( penn.state.lw , 
                                                  penn.state.utm$smk_rand )
xbar = mean( penn.state.utm$smk_rand )
ybar = mean( penn.state.utm$smk.lagged.means.rand )
plot2 = ggplot( penn.state.utm@data , aes( x=smk_rand , 
                                           y=smk.lagged.means.rand ) ) +
  xlim( 15 , 30 ) + 
  ylim( 15 , 30 ) +
  geom_point() +
  geom_hline( yintercept=ybar ) +
  geom_vline( xintercept=xbar )
  
library( grid )
windows( 12 , 6 )
grid.newpage()

pushViewport( viewport( layout=grid.layout(1,2) ) )  
print( plot1 , vp=viewport(layout.pos.col=1 , height=5 ) )
print( plot2 , vp=viewport(layout.pos.col=2 , height=4 ) )


