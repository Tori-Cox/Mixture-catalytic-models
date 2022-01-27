######################################################
## Calculate titre thresholds per simulated dataset ##
######################################################


## read in simulated data
original_all <- readRDS("MixtureSims.RDS")


## estimate thresholds for each of the 200 simulated datasets 
stored_thresholds<-list()
for(sim in 1:540){
  original<-original_all[[sim]]
  
  #'status' is the true serostatus, 'sero' is the estimated serostatus based on optimised thresholds
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
for(sim in 1:540){
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
      original$sero[i] <- "Equivocal"
    }
  }
  new_all[[sim]] <- original
}

correct_classification<-NA
for(sim in 1:540){
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
