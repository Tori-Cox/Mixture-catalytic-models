#####################################
## MIXTURE MODEL PLOTTING - STEP 1 ##
#####################################


#mixtures fitted
mixfit$mixdata <- mixfit$mixdata[-nrow(mixfit$mixdata), ]
binwidth <-  abs(mixfit$mixdata$X[2] - mixfit$mixdata$X[1])
dist0 <- mixfit$distribution[[1]]
dist1	<- mixfit$distribution[[2]]

thresholds <- readRDS("MixtureSims_thresholds.RDS")
thresholds <- thresholds[[sim_num]]

## norm norm
if(dist0=="norm" & dist1=="norm"){    
  ggplot() + 
    geom_histogram(aes(x=mixfit$mixdata$X, y=mixfit$mixdata$count), stat="identity") +
    theme_bw()+
    
    stat_function(fun = function(x){
      dnorm(x, mean = mixfit$parameters$mu[1], sd = mixfit$parameters$sigma[1]) * binwidth * (sum(mixfit$mixdata$count)*mixfit$parameters$pi[1])
    }, colour = "blue", size = 1)+
    
    stat_function(fun = function(x){
      dnorm(x, mean = mixfit$parameters$mu[2], sd = mixfit$parameters$sigma[2]) * binwidth * (sum(mixfit$mixdata$count)*mixfit$parameters$pi[2])
    }, colour = "blue", size = 1)+
    
    
    geom_vline(aes(xintercept = mixfit$parameters$mu), colour = "red", linetype ="longdash", size = .8) +
    geom_vline(aes(xintercept = thresholds), colour = "black", linetype ="dashed", size = .8) +
    labs(y="Antibody titer count", x=paste("Log(antibody titre + 1)  \n seronegative distribution = ", dist0, "\n seropositive distribution = ", dist1, sep=""))
}

# norm gamma
if(dist0=="norm" & dist1=="gamma"){    
  ggplot() + 
    geom_histogram(aes(x=mixfit$mixdata$X, y=mixfit$mixdata$count), stat="identity") +
    theme_bw()+
    
    stat_function(fun = function(x){
      dnorm(x, mean = mixfit$parameters$mu[1], sd = mixfit$parameters$sigma[1]) * binwidth * (sum(mixfit$mixdata$count)*mixfit$parameters$pi[1])
    }, colour = "blue", size = 1)+
    
    stat_function(fun = function(x){
      dgamma(x, shape = ((mixfit$parameters$mu[2] /mixfit$parameters$sigma[2])^2), 
             rate = (mixfit$parameters$mu[2] / (mixfit$parameters$sigma[2] ^2))) * binwidth * (sum(mixfit$mixdata$count)*mixfit$parameters$pi[2])
    },  colour = "blue", size = 1)+
    
    
    geom_vline(aes(xintercept = mixfit$parameters$mu), colour = "red", linetype ="longdash", size = .8) +
    geom_vline(aes(xintercept = thresholds), colour = "black", linetype ="dashed", size = .8) +
    labs(y="Antibody titer count", x=paste("Log(antibody titre + 1)  \n seronegative distribution = ", dist0, "\n seropositive distribution = ", dist1, sep=""))
}
# norm weibull
if(dist0=="norm" & dist1=="weibull"){    
  ggplot() + 
    geom_histogram(aes(x=mixfit$mixdata$X, y=mixfit$mixdata$count), stat="identity") +
    theme_bw()+
    
    stat_function(fun = function(x){
      dnorm(x, mean = mixfit$parameters$mu[1], sd = mixfit$parameters$sigma[1]) * binwidth * (sum(mixfit$mixdata$count)*mixfit$parameters$pi[1])
    }, colour = "blue", size = 1)+
    
    stat_function(fun = function(x){
      dweibull(x, shape = weibullpar(mu = mixfit$parameters$mu[2], sigma = mixfit$parameters$sigma[2])$shape, 
               scale = weibullpar(mu = mixfit$parameters$mu[2], 
                                  sigma = mixfit$parameters$sigma[2])$scale) * binwidth * (sum(mixfit$mixdata$count)*mixfit$parameters$pi[2]) }, 
      colour = "blue", size = 1)+
    
    
    geom_vline(aes(xintercept = mixfit$parameters$mu), colour = "red", linetype ="longdash", size = .8) +
    geom_vline(aes(xintercept = thresholds), colour = "black", linetype ="dashed", size = .8) +
    labs(y="Antibody titer count", x=paste("Log(antibody titre + 1)  \n seronegative distribution = ", dist0, "\n seropositive distribution = ", dist1, sep=""))
}
#  weibull norm
if(dist0=="weibull" & dist1=="norm"){    
  ggplot() + 
    geom_histogram(aes(x=mixfit$mixdata$X, y=mixfit$mixdata$count), stat="identity") +
    theme_bw()+
    
    stat_function(fun = function(x){
      dweibull(x, shape = weibullpar(mu = mixfit$parameters$mu[1], sigma = mixfit$parameters$sigma[1])$shape, 
               scale = weibullpar(mu = mixfit$parameters$mu[1], 
                                  sigma = mixfit$parameters$sigma[1])$scale) * binwidth * (sum(mixfit$mixdata$count)*mixfit$parameters$pi[1])
      
    }, colour = "blue", size = 1)+
    
    stat_function(fun = function(x){
      dnorm(x, mean = mixfit$parameters$mu[2], sd = mixfit$parameters$sigma[2]) * binwidth * (sum(mixfit$mixdata$count)*mixfit$parameters$pi[2])
    }, colour = "blue", size = 1)+
    
    
    geom_vline(aes(xintercept = mixfit$parameters$mu), colour = "red", linetype ="longdash", size = .8) +
    geom_vline(aes(xintercept = thresholds), colour = "black", linetype ="dashed", size = .8) +
    labs(y="Antibody titer count", x=paste("Log(antibody titre + 1)  \n seronegative distribution = ", dist0, "\n seropositive distribution = ", dist1, sep=""))
}
#  weibull gamma
if(dist0=="weibull" & dist1=="gamma"){    
  ggplot() + 
    geom_histogram(aes(x=mixfit$mixdata$X, y=mixfit$mixdata$count), stat="identity") +
    theme_bw()+
    
    stat_function(fun = function(x){
      dweibull(x, shape = weibullpar(mu = mixfit$parameters$mu[1], sigma = mixfit$parameters$sigma[1])$shape, 
               scale = weibullpar(mu = mixfit$parameters$mu[1], 
                                  sigma = mixfit$parameters$sigma[1])$scale) * binwidth * (sum(mixfit$mixdata$count)*mixfit$parameters$pi[1])
      
    }, colour = "blue", size = 1)+
    
    stat_function(fun = function(x){
      dgamma(x, shape = ((mixfit$parameters$mu[2] /mixfit$parameters$sigma[2])^2), 
             rate = (mixfit$parameters$mu[2] / (mixfit$parameters$sigma[2] ^2))) * binwidth * (sum(mixfit$mixdata$count)*mixfit$parameters$pi[2])
    },  colour = "blue", size = 1)+
    
    
    geom_vline(aes(xintercept = mixfit$parameters$mu), colour = "red", linetype ="longdash", size = .8) +
    geom_vline(aes(xintercept = thresholds), colour = "black", linetype ="dashed", size = .8) +
    labs(y="Antibody titer count", x=paste("Log(antibody titre + 1)  \n seronegative distribution = ", dist0, "\n seropositive distribution = ", dist1, sep=""))
}
# weibull weibull
if(dist0=="weibull" & dist1=="weibull"){    
  ggplot() + 
    geom_histogram(aes(x=mixfit$mixdata$X, y=mixfit$mixdata$count), stat="identity") +
    theme_bw()+
    
    stat_function(fun = function(x){
      dweibull(x, shape = weibullpar(mu = mixfit$parameters$mu[1], sigma = mixfit$parameters$sigma[1])$shape, 
               scale = weibullpar(mu = mixfit$parameters$mu[1], 
                                  sigma = mixfit$parameters$sigma[1])$scale) * binwidth * (sum(mixfit$mixdata$count)*mixfit$parameters$pi[1])
      
    }, colour = "blue", size = 1)+
    
    stat_function(fun = function(x){
      dweibull(x, shape = weibullpar(mu = mixfit$parameters$mu[2], sigma = mixfit$parameters$sigma[2])$shape, 
               scale = weibullpar(mu = mixfit$parameters$mu[2], 
                                  sigma = mixfit$parameters$sigma[2])$scale) * binwidth * (sum(mixfit$mixdata$count)*mixfit$parameters$pi[2]) }, colour = "blue", size = 1)+
    
    
    geom_vline(aes(xintercept = mixfit$parameters$mu), colour = "red", linetype ="longdash", size = .8) +
    geom_vline(aes(xintercept = thresholds), colour = "black", linetype ="dashed", size = .8) +
    labs(y="Antibody titer count", x=paste("Log(antibody titre + 1)  \n seronegative distribution = ", dist0, "\n seropositive distribution = ", dist1, sep=""))
  
}
# gamma gamma
if(dist0=="gamma" & dist1=="gamma"){    
  ggplot() + 
    geom_histogram(aes(x=mixfit$mixdata$X, y=mixfit$mixdata$count), stat="identity") +
    theme_bw()+
    
    stat_function(fun = function(x){
      dgamma(x, shape = ((mixfit$parameters$mu[1] ^2) /(mixfit$parameters$sigma[1]^2)), 
             rate = (mixfit$parameters$mu[1] / (mixfit$parameters$sigma[1] ^2))) * binwidth * (sum(mixfit$mixdata$count)*mixfit$parameters$pi[1])
    },  colour = "blue", size = 1)+
    
    stat_function(fun = function(x){
      dgamma(x, shape = ((mixfit$parameters$mu[2] /mixfit$parameters$sigma[2])^2), 
             rate = (mixfit$parameters$mu[2] / (mixfit$parameters$sigma[2] ^2))) * binwidth * (sum(mixfit$mixdata$count)*mixfit$parameters$pi[2])
    },  colour = "blue", size = 1)+
    
    
    geom_vline(aes(xintercept = mixfit$parameters$mu), colour = "red", linetype ="longdash", size = .8) +
    geom_vline(aes(xintercept = thresholds), colour = "black", linetype ="dashed", size = .8) +
    labs(y="Antibody titer count", x=paste("Log(antibody titre + 1)  \n seronegative distribution = ", dist0, "\n seropositive distribution = ", dist1, sep=""))
}
# gamma weibull
if(dist0=="gamma" & dist1=="weibull"){    
  ggplot() + 
    geom_histogram(aes(x=mixfit$mixdata$X, y=mixfit$mixdata$count), stat="identity") +
    theme_bw()+
    
    stat_function(fun = function(x){
      dgamma(x, shape = ((mixfit$parameters$mu[1] ^2) /(mixfit$parameters$sigma[1]^2)), 
             rate = (mixfit$parameters$mu[1] / (mixfit$parameters$sigma[1] ^2))) * binwidth * (sum(mixfit$mixdata$count)*mixfit$parameters$pi[1])
    },  colour = "blue", size = 1)+
    
    
    stat_function(fun = function(x){
      dweibull(x, shape = weibullpar(mu = mixfit$parameters$mu[2], sigma = mixfit$parameters$sigma[2])$shape, 
               scale = weibullpar(mu = mixfit$parameters$mu[2], 
                                  sigma = mixfit$parameters$sigma[2])$scale) * binwidth * (sum(mixfit$mixdata$count)*mixfit$parameters$pi[2]) }, colour = "blue", size = 1)+
    
    
    geom_vline(aes(xintercept = mixfit$parameters$mu), colour = "red", linetype ="longdash", size = .8) +
    geom_vline(aes(xintercept = thresholds), colour = "black", linetype ="dashed", size = .8) +
    labs(y="Antibody titer count", x=paste("Log(antibody titre + 1)  \n seronegative distribution = ", dist0, "\n seropositive distribution = ", dist1, sep=""))
}
# gamma norm 
if(dist0=="gamma" & dist1=="norm"){    
  ggplot() + 
    geom_histogram(aes(x=mixfit$mixdata$X, y=mixfit$mixdata$count), stat="identity") +
    theme_bw()+
    
    stat_function(fun = function(x){
      dgamma(x, shape = ((mixfit$parameters$mu[1] ^2) /(mixfit$parameters$sigma[1]^2)), 
             rate = (mixfit$parameters$mu[1] / (mixfit$parameters$sigma[1] ^2))) * binwidth * (sum(mixfit$mixdata$count)*mixfit$parameters$pi[1])
    },  colour = "blue", size = 1)+
    
    stat_function(fun = function(x){
      dnorm(x, mean = mixfit$parameters$mu[2], sd = mixfit$parameters$sigma[2]) * binwidth * (sum(mixfit$mixdata$count)*mixfit$parameters$pi[2])
    }, colour = "blue", size = 1)+
    
    
    geom_vline(aes(xintercept = mixfit$parameters$mu), colour = "red", linetype ="longdash", size = .8) +
    geom_vline(aes(xintercept = thresholds), colour = "black", linetype ="dashed", size = .8) +
    labs(y="Antibody titer count", x=paste("Log(antibody titre + 1)  \n seronegative distribution = ", dist0, "\n seropositive distribution = ", dist1, sep=""))
}

ggsave(filename = paste("model_output/mixture/mix_fit_", sim_num, ".png", sep=""))




#save table of parameters

stored_params$dist0[sim_num] <- mixfit$distribution[[1]]
stored_params$dist1[sim_num]	<- mixfit$distribution[[2]]
stored_params$mu0[sim_num]	<- mixfit$parameters$mu[1]
stored_params$mu1[sim_num]	<- mixfit$parameters$mu[2]
stored_params$sd0[sim_num]	<- mixfit$parameters$sigma[1]
stored_params$sd1[sim_num]	<- mixfit$parameters$sigma[2]

stored_params$mu0_low[sim_num]	<- mixfit$parameters$mu[1] - 1.96*mixfit$se$mu.se[1]
stored_params$mu1_low[sim_num]	<- mixfit$parameters$mu[2] - 1.96*mixfit$se$mu.se[2]
stored_params$sd0_low[sim_num]	<- mixfit$parameters$sigma[1] - 1.96*mixfit$se$sigma.se[1]
stored_params$sd1_low[sim_num]	<- mixfit$parameters$sigma[2] - 1.96*mixfit$se$sigma.se[2]

stored_params$mu0_upp[sim_num]	<- mixfit$parameters$mu[1] + 1.96*mixfit$se$mu.se[1]
stored_params$mu1_upp[sim_num]	<- mixfit$parameters$mu[2] + 1.96*mixfit$se$mu.se[2]
stored_params$sd0_upp[sim_num]	<- mixfit$parameters$sigma[1] + 1.96*mixfit$se$sigma.se[1]
stored_params$sd1_upp[sim_num]	<- mixfit$parameters$sigma[2] + 1.96*mixfit$se$sigma.se[2]

stored_params$chisq[sim_num] <- mixfit$chisq
stored_params$p[sim_num] <- mixfit$P
stored_params$AIC[sim_num] <- mixfit$AIC
stored_params$df[sim_num] <- mixfit$df
