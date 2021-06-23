########################################################################################
## Functions for simulating data and estimating serostatus using optimised thresholds ##
#######################################################################################


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


## threshold estimation function to be minimised
est_serostatus <-function(new_data, par){
  
  if(par[1]<par[2]){
    for(i in 1:nrow(new_data)){
      if(new_data$seroval[i] <= par[1]){
        new_data$sero[i] <- 0
      }
      if(new_data$seroval[i] >= par[2]){
        new_data$sero[i] <- 1
      }
      if(new_data$seroval[i] > par[1] & new_data$seroval[i] < par[2]){
        new_data$sero[i] <- "Equivical"
      }
      
      ifelse(new_data$sero[i]==new_data$status[i], new_data$match[i]<-0, new_data$match[i]<-1)
    }
    sum<-sum(new_data$match)
  }else{sum<-Inf}
  return(sum)
}