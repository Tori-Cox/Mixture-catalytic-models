######################################
## functions for calculating SP/FOI ##
######################################


## Fitting the P spline 

    ## optimising spline parameters
              spline_optimisation_fun <- function(z,a){
                
                ## set up parameters to optimise
                alpha_grid <- c(0.001, 0.01, 0.1, 0.5, 1, 5, 10, 50, 100, 500)
                knots_grid <- 5:length(agemin:agemax)
                res_knots <- matrix(ncol=2,nrow=length(knots_grid))
                res_alpha <- matrix(ncol=2,nrow=length(alpha_grid))
                
                #to work out starting alpha for the optimisation
                for (q in (1:length(alpha_grid))) {
                  
                  for (y in (1:length(knots_grid))) {
                    
                    fitC <- mpspline.fit(response=z, x.var=a, ps.intervals=knots_grid[y], degree=3,
                                         order=2, alpha = 1, kappa=1e8)
                    
                    res_knots[y, 1] <- knots_grid[y]
                    res_knots[y, 2] <- fitC$bic
                    knots_fin <- res_knots[res_knots[,2] == min(res_knots[,2]),1]
                  }
                  
                  fitC <- mpspline.fit(response=z, x.var=a, ps.intervals=knots_fin, degree=3,
                                       order=2, alpha = alpha_grid[q], kappa=1e8)
                  
                  res_alpha[q,1] <- alpha_grid[q]
                  res_alpha[q,2] <- fitC$bic
                  alpha_fin <- res_alpha[res_alpha[,2] == min(res_alpha[,2]),1]
                }
                
                opt_param <- data.frame(alpha_fin=alpha_fin, knots_fin=knots_fin)
                return(list(opt_param, fitC))
              }
              
              
              

      ## use optimised parameters in function for spline 
            spline.estimator <- function(original) {
              z <- original$z
              a <- original$a
              
              fit <- mpspline.fit(response=z, x.var=a, ps.intervals=knots_fin, degree=degree_fin,
                                  order=2, alpha = alpha_fin, kappa=1e8)
              mu_a <- fit$yhat
              return(mu_a)
            }
            
            ### for Uncertainty intervals
            resampler <- function(original) {
              n <- nrow(original)
              
              repeat{
                resample.rows <- sample(1:n,size=n,replace=TRUE)
                r <- original[resample.rows,]
                z_new <-r$z[order(r$a)]
                a_new <-r$a[order(r$a)]
                in_age_order <- data.frame(z=z_new, a=a_new)
                if (min(in_age_order$a)==min(original$a) & max(in_age_order$a)==max(original$a)) break   
                ###this is so the age lengths are the same for the spline cis function 
              }
              return(in_age_order)
            }
            
            spline.cis <- function(B, original) {
              
              spline.main <- spline.estimator(original = original)
              
              # Draw B bootstrap samples, fit the spline to each
              spline.boots <- replicate(B, spline.estimator(original = resampler(original = original)))
              
              cis.upper <- apply(spline.boots,1,quantile,probs=0.975)
              cis.lower <- apply(spline.boots,1,quantile,probs=0.025)
              spline_plotting <- list(main.curve=spline.main,lower.ci=cis.lower,upper.ci=cis.upper)
              
              df <- matrix(nrow=B, ncol=(agemax+1-agemin))
              
              for(i in 1:B){
                df[i, ] <- spline.boots[,,i]
              }
              df <- as.data.frame(df)
              colnames(df)<- as.character(agemin:agemax)
              
              if(agemin!=1.5){                             
                fix <- matrix(ncol= (agemin-1.5))
                fix <- as.data.frame(fix)
                df <- cbind(fix, df)
                colnames(df) <- 1:agemax
                agemin <- 1.5
              }
              
              df_a <- data.frame(Age = agemin:agemax, Mean = NA, SD = NA)
              for(a in (agemin-0.5):(agemax-0.5)){
                df_a[a,]$Mean <- mean(df[,a])
                df_a[a,]$SD <- sqrt(var(df[,a]))
              }
              df_a <- na.omit(df_a)
              
              
              new <- matrix(nrow=B, ncol=(agemax+1-agemin)) 
              new <- as.data.frame(new)
              colnames(new) <- agemin:agemax
              
              #for the derivative (will be one less column)
              for(REPEAT in 1:B){
                for(i in (agemin-0.5):(agemax-0.5-1)){
                  new[REPEAT,i] <- df[REPEAT,i+1] - df[REPEAT,i]
                }
              }
              
              df_a_deriv <- data.frame(Age = agemin:(agemax-1), Mean = NA, SD = NA)
              for(a in (agemin-0.5):(agemax-0.5-1)){
                df_a_deriv[a,]$Mean <- mean(new[,a])
                df_a_deriv[a,]$SD <- sqrt(var(new[,a]))
              }
              
              df_a_deriv <- na.omit(df_a_deriv)
              agemin <- as.numeric(min(original$a)) # reset agemin
              
              return(list(df_a, df_a_deriv, spline_plotting))
            }


## estimating seroprevalence using the spline
            
            seroprev_draw <- function(agemin, agemax, muS, muI, SigS, SigI, mu_a_data2){
              
              df <- matrix(nrow = (agemax+1-agemin))
              
              if(agemin!=1.5){                             
                fix <- matrix(nrow= (agemin-1.5))
                fix <- as.data.frame(fix)
                df <- rbind(fix, df)
                
                
                fix2 <- matrix(nrow= (agemin-1), ncol = 3)
                fix2 <- as.data.frame(fix2)
                colnames(fix2) <- c("Age", "Mean", "SD")
                mu_a_data2 <- rbind(fix2, mu_a_data2)
                
              }
              
              for(a in 1:(agemax-0.5)){
                mu_I <- muI
                mu_S <- muS
                
                if(is.na(mu_a_data2[a, ]$Mean)!=TRUE){ 
                
                mu_a <- rnorm(1, mean = (mu_a_data2[a, ]$Mean), sd = (mu_a_data2[a, ]$SD))
                
                if(mu_a > mu_S){
                  df[a, ] <- (mu_a - mu_S) / (mu_I-mu_S)
                }else{
                  df[a, ] <- NA}
                }
                else{
                  df[a, ] <- NA}
                }
            
              
              
              df <- as.matrix(df)
              return(df)
            }
            
            
            # confidence intervals on seroprev estimates
            seroprev_cis <- function(B, agemin, agemax, muS, muI, SigS, SigI, mu_a_data){
              
              seroprev_boots <- replicate(B, seroprev_draw(agemin=agemin, agemax=agemax, muS=muS, 
                                                           muI=muI, SigS=SigS, SigI=SigI, mu_a_data2=mu_a_data)) 
              
              df <- matrix(nrow=B, ncol=(agemax+1-agemin))
              
              if(agemin!=1.5){                             
                fix <- matrix(ncol= (agemin-1.5), nrow=B)
                fix <- as.data.frame(fix)
                df <- cbind(fix, df)
                colnames(df) <- 1:agemax
              }
              
              for(i in 1:B){
                df[i, ] <- seroprev_boots[,,i]
              }
              
              df <- as.data.frame(df)
              colnames(df)<- as.character(1.5:agemax)
              
              new_df <- list()
              for(a in (agemin-0.5):(agemax-0.5)){
                new_df[[a]] <- df[, a]
                new_df[[a]] <- na.omit(new_df[[a]])
              }
              
              serop_df <- data.frame(Age=(1.5:agemax), pi_a = NA, pi_a_low = NA, pi_a_upp = NA)
              
              for(a in (agemin-0.5):(agemax-0.5)){
                
                pi_a <- quantile(new_df[[a]], probs=0.5, type=1)
                serop_df$pi_a[serop_df$Age==a+0.5] <- pi_a
                
                pi_a_low <- quantile(new_df[[a]], probs=0.025, type=1)
                serop_df$pi_a_low[serop_df$Age==a+0.5] <- pi_a_low
                
                pi_a_upp <- quantile(new_df[[a]], probs=0.975, type=1)
                serop_df$pi_a_upp[serop_df$Age==a+0.5] <- pi_a_upp
                
              }
              serop_df <- na.omit(serop_df)
              return(serop_df)
            }



## estimating FOI using the spline
            
            FOI_draw <- function(agemin, agemax, muI, mu_a_data2, mu_deriv_data){
              
              muI_mean <- muI
              
              #distribution for mu_a term 
              mu_a_matrix <- matrix(nrow = (agemax+1-agemin))
              
              if(agemin!=1.5){ 
                
                fix2 <- matrix(nrow= (agemin-1))
                mu_a_matrix <- rbind(fix2, mu_a_matrix)
                
                fix2 <- matrix(nrow= (agemin-1), ncol = 3)
                fix2 <- as.data.frame(fix2)
                colnames(fix2) <- c("Age", "Mean", "SD")
                mu_a_data2 <- rbind(fix2, mu_a_data2)
                
              }
              
              for(a in (agemin-0.5):(agemax-0.5)){
                mu_a <- rnorm(1, mean = (mu_a_data2[a, ]$Mean), sd = (mu_a_data2[a, ]$SD))
                
                mu_a_matrix[a, ] <- mu_a
              }
              
              
              #distribution for mu_a_deriv term
              mu_deriv_matrix <- matrix(nrow = (agemax-agemin))
              
              if(agemin!=1.5){ 
                
                fix2 <- matrix(nrow= (agemin-1))
                mu_deriv_matrix <- rbind(fix2, mu_deriv_matrix)
                
                fix2 <- matrix(nrow= (agemin-1), ncol = 3)
                fix2 <- as.data.frame(fix2)
                colnames(fix2) <- c("Age", "Mean", "SD")
                mu_deriv_data<- rbind(fix2, mu_deriv_data)
                
              }
              
              for(a in (agemin-0.5):(agemax-0.5-1)){
                mu_deriv <- rnorm(1, mean = (mu_deriv_data[a, ]$Mean), sd = (mu_deriv_data[a, ]$SD))
                mu_deriv <- pmax(mu_deriv,0) #under assumption of monotonicity of spline (increasing titre mean with age)
                mu_deriv <- pmin(mu_deriv,1)
                mu_deriv_matrix[a, ] <- mu_deriv
              }
              
 
              df <- matrix(nrow = (agemax-agemin))
              if(agemin!=1.5){ 
                fix2 <- matrix(nrow= (agemin-1))
                df <- rbind(fix2, df)
                
                
              }
              
              for(a in (agemin:(agemax-1))){
                
                if((muI_mean > mu_a_matrix[a, ]) & (mu_deriv_matrix[a, ]) < (muI_mean - mu_a_matrix[a, ])) { 
                  
                  df[a, ] <- (mu_deriv_matrix[a, ]) / (muI_mean - mu_a_matrix[a, ]) #muI>mu_a always, see Bollaerts et al, 2012
                
                  }else{
                  df[a, ] <- NA
                }
                
              }
              
              return(df)
            }
            
            
            # confidence intervals on FOI estimates
            FOI_cis <- function(B, agemin, agemax, muI, mu_a_data, mu_deriv_data){
              
              
              FOI_boots <- replicate(B, FOI_draw(agemin=agemin, agemax=agemax, muI=muI, mu_a_data2=mu_a_data, mu_deriv_data=mu_deriv_data))
              
              df <- matrix(nrow=B, ncol=(agemax-agemin))
              
              if(agemin!=1.5){
                fix <- matrix(ncol= (agemin-1), nrow = B)
                fix <- as.data.frame(fix)
                df <- cbind(fix, df)
                colnames(df) <- 1:(agemax-1)
              }
              
              for(i in 1:B){
                df[i, ] <- FOI_boots[,,i]
              }
              
              
              df <- as.data.frame(df)
              colnames(df)<- as.character(1.5:(agemax-1))
              
              new_df <- list()
              for(a in (agemin-0.5):(agemax-0.5-1)){
                new_df[[a]] <- df[, a]
                new_df[[a]] <- na.omit(new_df[[a]])
              }
              
              
              FOI_df <- data.frame(Age=(1.5:(agemax-1)), FOI = NA, FOI_low = NA, FOI_upp = NA)
              FOI_sd <- data.frame(Age =(1.5:(agemax-1)), Mean = NA, SD = NA)
              
              for(a in (agemin-0.5):(agemax-0.5-1)){
                
                FOI <- quantile(new_df[[a]], probs=0.5, type=1)
                FOI_df$FOI[FOI_df$Age==a+0.5] <- FOI
                
                FOI_low <- quantile(new_df[[a]], probs=0.025, type=1)
                FOI_df$FOI_low[FOI_df$Age==a+0.5] <- FOI_low
                
                FOI_upp <- quantile(new_df[[a]], probs=0.975, type=1)
                FOI_df$FOI_upp[FOI_df$Age==a+0.5] <- FOI_upp
                
                FOI_sd[a,]$Mean <- mean(new_df[[a]])
                FOI_sd[a,]$SD <- sqrt(var(new_df[[a]]))
                
              }
              FOI_sd <- na.omit(FOI_sd)
              FOI_df <- na.omit(FOI_df)
              
              return(list(FOI_df, FOI_sd))
            }
            
            
            
            FOI_combine <- function(data){
              data <- na.omit(data)
              FOI_total <- mean(data$FOI)
              FOI_total_ci_l <- mean(data$FOI_low)
              FOI_total_ci_u <- mean(data$FOI_upp)
              FOI_total <- data.frame(FOI = FOI_total, lower_ci = FOI_total_ci_l, upper_ci = FOI_total_ci_u)
              return(FOI_total)
            }
            
            
