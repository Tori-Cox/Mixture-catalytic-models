######################################################
## Calculate titre thresholds per simulated dataset ##
######################################################


## read in simulated data
original_all <- readRDS("MixtureSims.RDS")


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


## estimate thresholds for each of the 200 simulated datasets 
stored_thresholds<-list()
for(sim in 1:200){
  original<-original_all[[sim]]
  
  #status is the true serostatus, we will define 'sero' as the estimated serostatus based on optimised thresholds
  data<-original
  data$sero <- NA
  data$match <-NA
  
  initial_threshold <- c(mean(data$seroval), mean(data$seroval)+1)
  optim_res<- optim(par=initial_threshold, fn=est_serostatus, new_data=data)
  
  print(optim_res$par)
  
  stored_thresholds[[sim]]<- optim_res$par
}

saveRDS(stored_thresholds, "MixtureSims_thresholds.RDS")



## calculate misclassifcation rates

new_all<-list()
for(sim in 1:200){
  original<-original_all[[sim]]
  original$sero<-NA
  par <- stored_thresholds[[sim]]
  
  for(i in 1:nrow(original)){
    
    if(original$seroval[i] <= par[1]){
      original$sero[i] <- 0
    }
    if(original$seroval[i] >= par[2]){
      original$sero[i] <- 1
    }
    if(original$seroval[i] > par[1] & original$seroval[i] < par[2]){
      original$sero[i] <- "Equivical"
    }
  }
  new_all[[sim]] <- original
}

correct_classification<-NA
for(sim in 1:200){
  new_single<-new_all[[sim]]
  new_single$sum<-NA
  for(i in 1:nrow(new_single)){
    ifelse(new_single$sero[i]==new_single$status[i], new_single$sum[i]<-1, new_single$sum[i]<-0)
  }
  value  <- sum(new_single$sum)
  correct_classification[sim] <- print(value/nrow(new_single)) 
}

saveRDS(correct_classification, "MixtureSims_correct_classification_vector.RDS")
saveRDS(new_all, "MixtureSims_serostatus_categorised.RDS")



## calculate number of equivical samples

eq_classification<-NA
for(sim in 1:200){
  new_single<-new_all[[sim]]
  new_single$sum<-NA
  for(i in 1:nrow(new_single)){
    ifelse(new_single$sero[i]=="Equivical", new_single$sum[i]<-1, new_single$sum[i]<-0)
  }
  value  <- sum(new_single$sum)
  eq_classification[sim] <- print(value) 
}
saveRDS(eq_classification, "MixtureSims_number_equivical_classification_vector.RDS")
