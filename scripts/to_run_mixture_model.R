#######################
#### Mixture Model ####
#######################


##################################
## STEP 1 - fitting the mixture ##
##################################

## prepare storage dataframe for params 
    stored_params <- data.frame(sim=1:540, dist0=NA, dist1=NA,
                                mu0=NA, mu0_low=NA, mu0_upp=NA,
                                mu1=NA, mu1_low=NA, mu1_upp=NA,
                                sd0=NA, sd0_low=NA, sd0_upp=NA,
                                sd1=NA, sd1_low=NA, sd1_upp=NA,
                                chisq=NA, p=NA, AIC=NA, df=NA,
                                seroprev=NA, seroprev_low=NA, seroprev_upp=NA,
                                foi =NA, foi_low =NA, foi_upp =NA)
    

## read in simulated data
    data <- readRDS("MixtureSims.rds")
    
## loop
    for(sim_num in 1:540){
     
    data_sub <- data[[sim_num]]
    z<-data_sub$seroval[data_sub$seroval<10] #constrain to remove outliers 
    breaks = 50
    mixdat <- mixgroup(z, breaks) # functions from mixdist package to group data
    mixpar <- mixparam(pi=c(0.5, 0.5), mu = c(0.5,3), sigma = c(0.5, 0.5))
    constr = list(conpi = "NONE", conmu = "NONE", consigma = "NONE", 
                  fixpi = NULL, fixmu = NULL, fixsigma = NULL, cov = NULL, size = NULL)
    
    weibull_norm <- list()
    weibull_norm$AIC <- NA 
    weibull_gamma <- list()
    weibull_gamma$AIC <- NA
    weibull_weibull <- list()
    weibull_weibull$AIC <- NA
    gamma_weibull <- list()
    gamma_weibull$AIC <- NA
    gamma_gamma <- list()
    gamma_gamma$AIC <- NA
    gamma_norm <- list()
    gamma_norm$AIC <- NA
    norm_weibull <- list()
    norm_weibull$AIC <- NA
    norm_norm <- list()
    norm_norm$AIC <- NA
    norm_gamma <- list()
    norm_gamma$AIC <- NA
    
    tryCatch({
      weibull_weibull <- mix(mixdat, mixpar,  dist1="weibull",  dist2="weibull", emsteps = 1,  exptol = 5e-06, constr=constr)
      }, error=function(e){
        print(paste("Unable to fit WW mixture distribution to Sim", sim_num, sep=' '))
      })
    tryCatch({
      weibull_norm <- mix(mixdat, mixpar, dist1="weibull", dist2="norm", emsteps = 1,  exptol = 5e-06, constr=constr)
    }, error=function(e){
      print(paste("Unable to fit WN mixture distribution to Sim", sim_num, sep=' '))
    })
    tryCatch({
      weibull_gamma <- mix(mixdat, mixpar,  dist1="weibull",  dist2="gamma", emsteps = 1,  exptol = 5e-06, constr=constr)
    }, error=function(e){
      print(paste("Unable to fit WG mixture distribution to Sim", sim_num, sep=' '))
    })
    tryCatch({
      norm_weibull <- mix(mixdat, mixpar,  dist1="norm",  dist2="weibull",emsteps = 1,  exptol = 5e-06, constr=constr)
    }, error=function(e){
      print(paste("Unable to fit NW mixture distribution to Sim", sim_num, sep=' '))
    })
    tryCatch({
      gamma_weibull <- mix(mixdat, mixpar,  dist1="gamma",  dist2="weibull", emsteps = 1,  exptol = 5e-06, constr=constr)
    }, error=function(e){
      print(paste("Unable to fit GW mixture distribution to Sim", sim_num, sep=' '))
    })
    tryCatch({
    norm_norm <- mix(mixdat, mixpar,  dist1="norm",  dist2="norm", emsteps = 1,  exptol = 5e-06, constr=constr)
    }, error=function(e){
      print(paste("Unable to fit NN mixture distribution to Sim", sim_num, sep=' '))
    })
    tryCatch({
      norm_gamma <- mix(mixdat, mixpar,  dist1="norm",  dist2="gamma", emsteps = 1,  exptol = 5e-06, constr=constr)
    }, error=function(e){
      print(paste("Unable to fit NG mixture distribution to Sim", sim_num, sep=' '))
    })
    tryCatch({
      gamma_norm <- mix(mixdat, mixpar,  dist1="gamma",  dist2="norm", emsteps = 1,  exptol = 5e-06, constr=constr)
    }, error=function(e){
      print(paste("Unable to fit GN mixture distribution to Sim", sim_num, sep=' '))
    })
    tryCatch({
      gamma_gamma <- mix(mixdat, mixpar,  dist1="gamma",  dist2="gamma",emsteps = 1,  exptol = 5e-06, constr=constr)
    }, error=function(e){
      print(paste("Unable to fit GG mixture distribution to Sim", sim_num, sep=' '))
    })

    #choosing optimal distribution pair based on AIC
    df <- data.frame(AIC <- c(weibull_norm$AIC, weibull_gamma$AIC, weibull_weibull$AIC,
                              norm_norm$AIC, norm_gamma$AIC, norm_weibull$AIC,
                              gamma_norm$AIC, gamma_gamma$AIC, gamma_weibull$AIC))
    colnames(df) <- "AIC"
    c_min <- min(df$AIC, na.rm = TRUE)
    
    if(is.na(weibull_weibull$AIC)==TRUE | weibull_weibull$AIC != c_min){
      if(is.na(weibull_norm$AIC)==TRUE | weibull_norm$AIC != c_min){
        if(is.na(weibull_gamma$AIC)==TRUE | weibull_gamma$AIC != c_min){
          if(is.na(norm_weibull$AIC)==TRUE | norm_weibull$AIC != c_min){
            if(is.na(norm_norm$AIC)==TRUE | norm_norm$AIC != c_min){
              if(is.na(norm_gamma$AIC)==TRUE | norm_gamma$AIC != c_min){
                if(is.na(gamma_gamma$AIC)==TRUE | gamma_gamma$AIC != c_min){
                  if(is.na(gamma_norm$AIC)==TRUE | gamma_norm$AIC != c_min){
                    choice <- gamma_weibull}
                  
                  else{choice <- gamma_norm}
                } else{choice <- gamma_gamma}
              } else{choice <- norm_gamma}
            } else{choice <- norm_norm}
          } else{choice <- norm_weibull}
        } else{choice <- weibull_gamma}
      } else{choice <- weibull_norm}
    } else{choice <- weibull_weibull}
    
    mixfit <- choice
    
     
    # source R code files for plotting and storing params
    source("scripts/plotting_mixture_model_1.R") 
   
    }
    write.csv(stored_params, file = "model_output/mixture/Mix_fit_params.csv") #save params for the simulated datasets
    
    
 
    
###############################################
## STEP 2 - calculating seroprevalence & FOI ##
###############################################
    
    stored_params <- read.csv("model_output/mixture/Mix_fit_params.csv")
    
    ## loop
    data <- readRDS("MixtureSims.rds")
    for(sim_num in 1:540){
      data_sub <- data[[sim_num]]
      
      z<-data_sub$seroval[data_sub$seroval<10]  
      a<-data_sub$age[data_sub$seroval<10] +0.5
      z_and_a <- data.frame(z=z, a=a)
      agemin <- min(data_sub$age)+0.5
      agemax <- max(data_sub$age)+0.5
      
    
    ## optimise parameters for spline fitting function
    opt_spline_param <- spline_optimisation_fun(z=z, a=a)
    
    degree_fin <- 3
    alpha_fin <- opt_spline_param[[1]]$alpha_fin
    knots_fin <- opt_spline_param[[1]]$knots_fin
    plot(opt_spline_param[[2]]$x, opt_spline_param[[2]]$yhat)

    
    ## fit spline 
    spline_mu_output <- spline.cis(B = 500, original = z_and_a) 
    
    spline_mu <- spline_mu_output[[1]]
    spline_deriv <- spline_mu_output[[2]]
    mu_a_plotting <- data.frame(age= agemin:agemax, mu_a=NA, mu_a_lower = NA, mu_a_upper = NA)
    mu_a_plotting$mu_a <- spline_mu_output[[3]]$main.curve
    mu_a_plotting$mu_a_lower <- spline_mu_output[[3]]$lower.ci
    mu_a_plotting$mu_a_upper <- spline_mu_output[[3]]$upper.ci
    
    ## estimate seroprevalence
        muS	<- stored_params$mu0[sim_num]
        muI	<- stored_params$mu1[sim_num]
        SigS	<- stored_params$sd0[sim_num]
        SigI	<- stored_params$sd1[sim_num]
        
        seroprevalence <- seroprev_cis(B=5000, mu_a_data = spline_mu, agemin=agemin, agemax=agemax, muS=muS, muI=muI, SigS=SigS, SigI=SigI)
        
        z_and_a$tally <- 1
        z_and_a_agg <- merge(data.frame(Group.1=min(z_and_a$a):max(z_and_a$a)),
                             aggregate(z_and_a, by = list(a), sum), all.x=TRUE)
        z_and_a_agg <- z_and_a_agg[, c(1,4)]
        
        for(age in min(seroprevalence$Age):max(seroprevalence$Age)){
          if(is.na(z_and_a_agg$tally[z_and_a_agg$Group.1==age]) == FALSE){
            seroprevalence$N[seroprevalence$Age==age] <- z_and_a_agg$tally[z_and_a_agg$Group.1==age]
          }}
        
        
        seroprevalence$tot <- seroprevalence$pi_a*seroprevalence$N
        seroprevalence$tot_low <- seroprevalence$pi_a_low*seroprevalence$N 
        seroprevalence$tot_upp <- seroprevalence$pi_a_upp*seroprevalence$N

        tot_seroprev <- sum(seroprevalence$tot, na.rm=TRUE)/sum(seroprevalence$N, na.rm=TRUE)
        tot_seroprev_low <- sum(seroprevalence$tot_low, na.rm=TRUE)/sum(seroprevalence$N, na.rm=TRUE)
        tot_seroprev_upp <- sum(seroprevalence$tot_upp, na.rm=TRUE)/sum(seroprevalence$N, na.rm=TRUE)
        
    
    ## estimate FOI
        FOI <- FOI_cis(B=5000, mu_a_data = spline_mu, mu_deriv_data = spline_deriv, agemin=agemin, agemax=agemax, muI=muI)
        
        FOI_values  <- FOI[[1]]   
        FOI_var <- FOI[[2]]
        
        FOI_total <- FOI_combine(data=FOI_values)
      
    
    ## plotting    
    source("scripts/plotting_mixture_model_2.R")
        
    }
    stored_params[c(21:23)][stored_params$seroprev>1,]<-NA
    
    write.csv(stored_params, file = "model_output/mixture/Mix_fit_and_est_params.csv")
    