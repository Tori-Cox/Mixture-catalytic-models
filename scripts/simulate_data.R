######################################
#### Generation of simulated data ####
######################################


## set parameter values
pars <- data.frame(sim=paste('sim', seq(1:200)), dist0=NA, dist1=NA, foi=NA, 
                   mu0=NA, mu1=NA, sd0=NA, sd1=NA, seroprev=NA)

pars$foi <- runif(200, 0.001, 0.18)
pars$mu0 <- runif(200, 0.05, 1.2)
pars$mu1 <- runif(200, 2.2, 4)
pars$sd0 <- runif(200, 0.1, 0.8)
pars$sd1 <- runif(200, 0.1, 0.8)

dists <- c('norm','gamma','weibull')


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

