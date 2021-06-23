#########################
#### Catalytic Model ####
#########################


## read in data
original_all <- readRDS("MixtureSims_serostatus_categorised.RDS")


## set up storage lists (V for time-varying, C for time-constant)
store_foiV<-list()
store_spV<-list()
store_foiC<-list()
store_spC<-list()


## loop
for(sim in 1:200){
  
  original <- original_all[[sim]]
  original <- original[original$sero!="Equivical", ]
  original$sero <- as.numeric(original$sero)
  
  loc.all <- data_processing(original) #function in cat_model_functions.R
  nages <- nrow(loc.all)
  
  
  ## model fitting 
  
        # fit initial models without CIs to check values    
        fit.c<-fit.data(data=loc.all, npar=1) ## npar= number of parameters to fit. Use 1 for a time constant model.
        fit.c$par #age invariant FOI estimate
        
        fit.v<-fit.data(data=loc.all, npar=nages)
        mean(fit.v$par) #mean of age specific FOI
        
  # fit the models with CIs
  fit.c_CI <- fit.data.cis(B=500, data=loc.all, original=original, npar=1)
  fit.v_CI <- fit.data.cis(B=500, data=loc.all, original=original, npar=nages)
  
  # age-specific seroprevalence estimates
  SP_c_FOI <- seroprev.cis(fit= fit.c_CI, data = loc.all)
  SP_v_FOI <- seroprev.cis(fit= fit.v_CI , data = loc.all)
  
  # total population seroprevalence
  SP_total_c <- total_seroprev(original=loc.all, est=SP_c_FOI)
  SP_total_v <- total_seroprev(original=loc.all, est=SP_v_FOI)
  
  
  ## plotting
  source("scripts/plotting_catalytic_models_1.R") 
  
  }

saveRDS(store_foiV, "model_output/catalytic/Cat_MixtureSims_est_FOI_variable.RDS")
saveRDS(store_spV, "model_output/catalytic/Cat_MixtureSims_est_SP_variable.RDS")
saveRDS(store_foiC, "model_output/catalytic/Cat_MixtureSims_est_FOI_constant.RDS")
saveRDS(store_spC, "model_output/catalytic/Cat_MixtureSims_est_SP_constant.RDS")

