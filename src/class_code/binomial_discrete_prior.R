
theta = seq( 0 , 1 , 0.01 )              ## all possible values of theta
p = rep(1/length(theta),length(theta))   ## prior probabilities (1/101)

n = 10                                   ## of trials
x = 7                                    ## of obsereved successes in n trials
den = sum( p*dbinom(x,n,theta) )         ## denominator of Bayes theorem
num = p*dbinom(x,n,theta)                ## numerator in Bayes theorem
posterior = num/den                      ## posterior distribution
plot( theta , posterior )                ## plot the posterior

sum( posterior )                         ## should be 1
theta[ which.max( posterior ) ]          ## maximum a posteriori
sum( theta*posterior )                   ## posterior mean
cum.prob = cumsum(posterior)             ## CDF of posterior
theta[ which.max( cum.prob > 0.5 ) ]     ## posterior median

