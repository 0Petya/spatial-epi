setwd("C:/Users/srigd/Dropbox/MyCourses/BST_5610_Fall2022")

install.packages( "nimble" )


####  The following is from Section 4.4 of the nimble manual
####  https://r-nimble.org/html_manual/cha-installing-nimble.html
path <- Sys.getenv('PATH')
newPath <- paste("C:\\Rtools\\bin;C:\\Rtools\\mingw_64\\bin;",
                 path, sep = "")
Sys.setenv(PATH = newPath) 

####  The following is an alternative to the above ...
# Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";"))
# Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")

install.packages("Rcpp")

library( nimble )
library( coda )
library( igraph )
library( R6 )


myCode = nimbleCode({
  for (i in 1:N) {
    x[i] ~ dnorm( mu , tau )    #### Note: nimble uses N(mean,precision)
  }
  mu ~ dnorm( 0 , 0.0001 )
  tau ~ dgamma( 0.1 , 0.001 )
})

####
####  Define all the data and constants as R objects
####

N = 16      ## 16 observations from N(mu,tau)
x = c( 111 , 153 , 238 , 123 , 177 , 155 , 138 , 105 ,
        98 , 155 , 149 , 122 , 170 , 142 , 149 , 133 )

myConsts = list( N = N )
myData = list( x = x )
myInits = list( mu = 100 , tau = 1 )

myModel = nimbleModel( myCode, 
                       data = myData, 
                       constants = myConsts, 
                       inits = myInits)
####
####  Compile the model code defined above.
####
compileMyModel = compileNimble(myModel)

####
####  Now, configure the simulation.  This means
####  to find a good plan to simulate the Markov
####  chain so it samples from the posterior.
####  

myConf = configureMCMC( compileMyModel , print = TRUE )
myConf$addMonitors( c("mu","tau" ) )

myMCMC = buildMCMC(myConf)
compileMyMCMC = compileNimble( myMCMC, project = myModel)

niter =  10000
nburn =   1000
set.seed(1)

start.time = proc.time()
samples = runMCMC( compileMyMCMC, niter = niter, nburnin = nburn,
                   inits = myInits, nchains = 1, samplesAsCodaMCMC = TRUE )
stop.time = proc.time()
time.elapsed = stop.time - start.time
print( time.elapsed )

head( samples )
windows( 9 , 7 )
par(mfrow = c(3, 1), mai = c(.6, .5, .4, .1), mgp = c(1.8, 0.7, 0))
ts.plot(samples[ , 'mu'], xlab = 'iteration', col="red" , lwd=1.5 ,
        ylab = expression(mu), main = expression(mu))
ts.plot(samples[ , 'tau'], xlab = 'iteration', col="red" , lwd=1.5 ,
        ylab = expression(tau), main = expression(tau))
ts.plot( sqrt(1/samples[ , 'tau']), xlab = 'iteration', col="red" , lwd=1.5 ,
        ylab = expression(sigma), main = expression(sigma))

windows( 9 , 12 )
par(mfrow = c(5, 1), mai = c(.6, .5, .4, .1), mgp = c(1.8, 0.7, 0))
xvec = seq(0,200,1)
yvec = dnorm( xvec , 0 , sqrt(1/0.0001) )
plot( xvec , yvec , type="l" )
hist( x , main="Original Data" , breaks=12 , xlim=c(0,200) )
hist( samples[ , 'mu' ] , breaks=50 , xlim=c(0,200) )
hist( samples[ , 'tau' ] , breaks=50 )
hist( 1/sqrt(samples[ , 'tau' ]) , breaks=50 )

source( "HDIofMCMC.R" )
source( "plotPost.R" )

windows( 9 , 9 )
par(mfrow = c(2, 1), mai = c(.6, .5, .4, .1), mgp = c(1.8, 0.7, 0))
plotPost( samples[,'mu'] )
plotPost( 1/sqrt(samples[,'tau']) )

