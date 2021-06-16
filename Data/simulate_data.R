######################################
#### Generation of simulated data ####
######################################


## load packages
library(mixdist)


## set parameter values
pars <- data.frame(sim=paste('sim', seq(1:200)), dist0=NA, dist1=NA, foi=NA, 
                   mu0=NA, mu1=NA, sd0=NA, sd1=NA, seroprev=NA)

pars$foi <- runif(200, 0.001, 0.18)
pars$mu0 <- runif(200, 0.05, 1.2)
pars$mu1 <- runif(200, 2.2, 4)
pars$sd0 <- runif(200, 0.1, 0.8)
pars$sd1 <- runif(200, 0.1, 0.8)

dists <- c('norm','gamma','weibull')


## simulate function
simulate <- function(age, Nsims){
  
  Sims <- list()
  titres <- vector()
  status <- vector()
  
  for(j in 1:Nsims){
    
    pars$dist0[j] <- dists[runif(1,1,4)]
    pars$dist1[j] <- dists[runif(1,1,4)]
    
    
    for(i in 1:length(age)){
      
      status[i] <- rbinom(1, 1, 1-exp(-pars$foi[j]*age[i]))
      
      if(status[i]==0){
        if(pars$dist0[j]=='norm') titres[i] <- rnorm(1, pars$mu0[j], pars$sd0[j])
        if(pars$dist0[j]=='gamma') titres[i] <- rgamma(1, shape=(pars$mu0[j]^2/pars$sd0[j]^2), scale=(pars$sd0[j]^2/pars$mu0[j]))
        if(pars$dist0[j]=='weibull') titres[i] <- rweibull(1, weibullpar(pars$mu0[j], pars$sd0[j])$shape, weibullpar(pars$mu0[j], pars$sd0[j])$scale)
      } 
      else{
        if(pars$dist1[j]=='norm') titres[i] <- rnorm(1, pars$mu1[j], pars$sd1[j])
        if(pars$dist1[j]=='gamma') titres[i] <- rgamma(1, shape=(pars$mu1[j]^2/pars$sd1[j]^2), scale=(pars$sd1[j]^2/pars$mu1[j]))
        if(pars$dist1[j]=='weibull') titres[i] <- rweibull(1, weibullpar(pars$mu1[j], pars$sd1[j])$shape, weibullpar(pars$mu1[j], pars$sd1[j])$scale)
      } 
    }
    pars$seroprev[j] <- sum(status)/length(status)
    Sims[[j]] <- data.frame(seroval=titres, age=age, status=status)
  }
  
  return(list(Sims=Sims, pars=pars, status=status, titres=titres))
}


## set the age-distribution. 
      ## we set age-distribution equal to that from our Indonesian serology data. Ref. Prayitno et al., 2017

age_df <- data.frame(Ages =1:18, Tally = c(121,166,202,177,164,176,175,167,166,190,172,159,166,194,253,244,210,64))
a<-rep(age_df$Ages, times=age_df$Tally)

sss <- simulate(age=a, 200)
sims <- sss$Sims
pars <- sss$pars 

## save outputs
saveRDS(sims, 'MixtureSims.RDS')
write.csv(pars, 'MixtureSimParams.csv', row.names=F)

