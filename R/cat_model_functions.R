###################################
#### Catalytic Model functions ####
###################################


### data processing
        data_processing <- function(data){
          
          # new dataframe for proportions seropositive by age
          new_data <- data.frame(Age=seq(1:18), AgeMid=seq(1:18)+0.5,
                                 PosData=NA, TotData=NA, serop=NA, lower_ci=NA, upper_ci=NA)
          
          # for each age 
          for(a in unique(data$age)){
            
            age <- data[data$age==a, ] # subset by age
            new_data$PosData[a] <- sum(age$sero) # number positive aged a
            new_data$TotData[a] <- nrow(age) # total number aged a
            
            # if there are no samples from a single year of age replace Na with 0
            new_data$TotData[which(is.na(new_data$TotData)==TRUE)] <- 0
            new_data$PosData[which(is.na(new_data$PosData)==TRUE)] <- 0
            
            # seroprevalence
            new_data$serop <- new_data$PosData/new_data$TotData
            
            # confidence intervals for age-specific seroprevalence
            if(new_data$TotData[a]>0){
              ci <- binom.exact(new_data$PosData[a], new_data$TotData[a], conf.level=0.95)
              new_data$lower_ci[a] <- ci$lower
              new_data$upper_ci[a] <- ci$upper
            }
          }
          
          # remove rows where there were no samples (there are only a couple)
          new_data <- new_data[(is.na(new_data$serop)==FALSE), ]
          loc.all<- data.frame("0"=(new_data$TotData-new_data$PosData), "1"=new_data$PosData,  "Age"=new_data$Age+0.5)
          return(loc.all)
          
        }



### resampler for bootstrap method uncertainty intervals
        resampler <- function(original) {
          n <- nrow(original)
          
          repeat{
            resample.rows <- sample(1:n,size=n,replace=TRUE)
            r <- original[resample.rows,]
            if (unique(r$age)==unique(original$age)) break   
          }
          resampled_data <- data_processing(r)
          return(resampled_data)
        }


        
### model 
        ### Ref for model: Rodriguez-Barraquer et al., 2011 
        ### see https://github.com/isabelrodbar/catalytic_model_mexico/blob/master/catalytic_model_ruth.R
        

        # liklihood
        likeli.cons<- function(theta, data, npar) {
          
          ##### Get relevant data
          shift <- min(data$Age)
          age<-(1:(max(as.numeric(dimnames(data)[[1]]))))+0.5
          age.dat<-age
          
          nyk<-data[,which(dimnames(data)[[2]]=="X1")] #SEROPOSITIVE
          N<-data[,which(dimnames(data)[[2]]=="X0")]+nyk
          nxk<-N-nyk #SERONEGATIVE
          
          #### Create a vector (or matrix) for your lambdas. This npar=1 for constant foi, npar>1 for time varying foi (divided equally in npar parts)	
          if (npar==1) {
            lambda<-matrix(theta, nrow=length(age), ncol=1)
            lambda[1,]<-theta*shift #To correct for the fact that this is the floor of the ages 
          }
          
          if (npar>1) {
            cut.ages<-1:max(age)
            lambda<-matrix(NA, nrow=length(age), ncol=1)
            
            lambda[,1] <-theta[cut.ages]
            lambda[1,]<-theta[1]*shift
          }
          
          #####Compute cumulative FOI experienced up to age a (cum.lambda)
          xaprod2<-cumsum(lambda)
          
          x<-exp(-xaprod2)	 #Proportion that remains susceptible at age a
          
          # Get proportions susceptible for observed age groups
          x.lik<-x
          y.lik<-1-x.lik
          
          #Add the log likelihood terms. This is a binomial likelihood of the form x^nxk * (y)^(N-nxk)
          lognxk<-nxk*log(x.lik) 
          lognyk<-(N-nxk)*log(y.lik)
          logterms<-lognxk+lognyk
          #Sum the two log likelihoods
          minuslogl<- -sum(logterms)
          return(minuslogl)
        }


        #Fit the model	
        fit.data <- function (data, npar, Hes=FALSE,...) {
          initial.guess<-abs(rnorm(npar, .02, 0.01))
          
          parscale <- 1e-5 + abs(initial.guess)
          
          opt.vals <- optim(par=initial.guess, 
                            fn=likeli.cons,   
                            data=data,    #dataset
                            npar=npar,	#Number of parameters to fit
                            method='L-BFGS-B', 
                            lower = rep(0.000001, length(initial.guess)),  #Constrain all lambdas to be greater than 0.000001
                            control=list(trace=TRUE, parscale=parscale,  maxit=100000),  #Some tuning
                            hessian=Hes, #Whether to estimate the hessian Matrix or not. Default set to FALSE
                            ...)
          
          opt.vals$init <- initial.guess
          
          return(opt.vals)
          
        }
        
        
        
### confidence intervals onto above functions 
        fit.data.cis <- function(B, data, original, npar) {
          
          fit.main<-fit.data(data, npar)
          
          # Draw B bootstrap samples
          fit.boots <- replicate(B, fit.data(data = resampler(original = original), npar))
          fit.boots.par<-fit.boots[1,]
          
          #formatting
          if(npar==1){
            fit.boots.df <- data.frame(Age="All")
          }
          if(npar>1){
            fit.boots.df <- data.frame(Age=rownames(data))
          }
          for(i in 1:length(fit.boots.par)){
            fit.boots.df <- cbind(fit.boots.df, fit.boots.par[[i]])
          }
          fit.boots.df<-fit.boots.df[,-1]
          fit.boots.df$upp<-NA
          fit.boots.df$low<-NA
          
          for(j in 1:nrow(fit.boots.df)){
            fit.boots.df$upp[j] <- quantile(fit.boots.df[j,], probs=0.975, na.rm=TRUE)
            fit.boots.df$low[j] <- quantile(fit.boots.df[j,], probs=0.025, na.rm=TRUE)
          }
          
          cis.lower <-  fit.boots.df$low
          cis.upper <-  fit.boots.df$upp
          
          
          foi_plotting <- list(fit.main$par,cis.lower,cis.upper, fit.boots.df)
          return(foi_plotting)
        }
        
        
        seroprev.cis <- function(fit, data){
          
          bootstrap <- fit[[4]]
          shift<-min(data$Age)
          npar<-nrow(bootstrap)
          age<-1:(max(as.numeric(dimnames(data)[[1]])))
          age.dat<-as.numeric(dimnames(data)[[1]])
          nyk<-data[,which(dimnames(data)[[2]]=="X1")]
          N<-data[,which(dimnames(data)[[2]]=="X0")]+nyk
          nxk<-N-nyk
          
          y <- list()
          for(i in 1:ncol(bootstrap)){
            theta <- bootstrap[,i]
            
            if (npar==1) {
              lambda<-matrix(theta, nrow=length(age), ncol=1)
              lambda[1,]<-theta*shift #To correct for the fact that this is the floor of the ages 
            }
            if (npar>1) {
              cut.ages<-1:max(age)
              lambda<-matrix(NA, nrow=length(age), ncol=1)
              lambda[,1] <-theta[cut.ages]
              lambda[1,]<-theta[1]*shift
            }  
            #Cumulative FOI upt to age a (cum.lambda)
            xaprod2<-cumsum(lambda)
            
            x<-exp(-xaprod2)	 #Proportion susceptible
            y[[i]] <- 1-x   #Proportion immune/seropositive
            
          }
          
          #formatting
          fit.boots.df <- data.frame(Age=rownames(data)) 
          for(i in 1:length(y)){
            fit.boots.df <- cbind(fit.boots.df, y[[i]])
          }
          fit.boots.df<-fit.boots.df[,-1]
          fit.boots.df$upp<-NA
          fit.boots.df$low<-NA
          
          for(j in 1:nrow(fit.boots.df)){
            fit.boots.df$upp[j] <- quantile(fit.boots.df[j,], probs=0.975, na.rm=TRUE)
            fit.boots.df$low[j] <- quantile(fit.boots.df[j,], probs=0.025, na.rm=TRUE)
            fit.boots.df$central[j] <- quantile(fit.boots.df[j,], probs=0.5, na.rm=TRUE)
          }
          
          cis.central <- fit.boots.df$central
          cis.lower <-  fit.boots.df$low
          cis.upper <-  fit.boots.df$upp
          
          return(list(Central=cis.central,Lower=cis.lower,Upper=cis.upper))
          
        }
        
        
### calculate population level seroprevalence
          total_seroprev<- function(original, est){
          
          nages <- nrow(original)
          original$serop_model<-est$Central
          original$serop_model_low<-est$Lower
          original$serop_model_upp<-est$Upper
          original$n <- original$X0 + original$X1
          
          original$n_serop_est <- original$serop_model * original$n
          original$n_serop_est_low<- original$serop_model_low * original$n
          original$n_serop_est_upp <- original$serop_model_upp * original$n
          
          
          est_total_serop <- sum(original$n_serop_est) / sum(original$n)
          est_total_serop_low <- sum(original$n_serop_est_low) / sum(original$n)
          est_total_serop_upp <- sum(original$n_serop_est_upp) / sum(original$n)
          
          RES<- c(est_total_serop, est_total_serop_low, est_total_serop_upp)
          return(RES)
        }
        
        
