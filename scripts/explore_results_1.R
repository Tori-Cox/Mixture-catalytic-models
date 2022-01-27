###############################################
## Simulation study - exploration of results ##
## 1. Mixture model ###########################
###############################################


## read in data
est <- read.csv("model_output/mixture/Mix_fit_and_est_params.csv") #estimated by mixture model
true<- read.csv("MixtureSimParams.csv") #simulated 'true' parameters
true_foi <- true$foi
true_serop <- true$seroprev


##################################
# calculating bias and uncertainty
##################################

        
        data <- data.frame(sim=1:540, True=true_foi, Est= est$foi, Est_low = est$foi_low,
                           Est_upp = est$foi_upp, Parameter = "FOI")
        
        data2 <- data.frame(sim=1:540, True=true_serop, Est= est$seroprev, Est_low = est$seroprev_low,
                            Est_upp = est$seroprev_upp, Parameter = "Seroprevalence")
        
        newdata <- rbind(data,data2)
        
        
        newdata$within <-0
        for(i in 1:(2*540)){
          if(is.na(newdata$Est[i])!=T & newdata$True[i] <= newdata$Est_upp[i] & newdata$True[i] >= newdata$Est_low[i]){
            newdata$within[i]<-1
          }
        }
        newdata$uncertainty <-NA
        for(i in 1:(2*540)){
          newdata$uncertainty[i]<-newdata$Est_upp[i]-newdata$Est_low[i]
        }
        newdata$bias <-NA
        for(i in 1:(2*540)){
          newdata$bias[i]<-newdata$Est[i]-newdata$True[i]
        }
        
        mix_plot <- ggplot(newdata) + facet_wrap("Parameter", scales="free")+
          geom_point(aes(x=True, y=Est, col=within), size=2, alpha=0.4) +
          geom_errorbar(aes(x=True, ymin=Est_low, ymax=Est_upp, col=within), alpha=0.4)+
          theme_bw() +
          labs(x="True values", y="Estimated values") +
          geom_abline(slope=1, linetype="dashed", size=1) +
          theme(axis.text = element_text(size=15), axis.title = element_text(size=20),
                strip.text = element_text(size=15), legend.position = "none")
        
        #ggsave(plot = mix_plot, filename = "Mix_results_sp_foi_est_vs_true.png")
        
        MIXDATA <- newdata
        saveRDS(MIXDATA, "results_MIXDATA.RDS")


        
        
#########################################################
# explore estimated versus true distribution parameters # 
#########################################################
        

      ## distribution family specification by mixture model
        
        colnames(true) <- c("sim", "dist0_true", "dist1_true",
                            "foi_true" ,"mu0_true", "mu1_true",
                            "sd0_true", "sd1_true", "seroprev_true")
        true$sim <- 1:540
        est2 <- est[,c(3:16)] #just distribution params 
        data<-cbind(true,est2)
        
        data$matchdist0 <- 0
        data$matchdist1 <- 0
        for(i in 1:540){
          if(data$dist0_true[i] == data$dist0[i]){
            data$matchdist0[i] <- 1
          }
        }
        for(i in 1:540){
          if(data$dist1_true[i] == data$dist1[i]){
            data$matchdist1[i] <- 1
          }
        }
        
        data$matchdistboth <- 0
        data$matchdistneither <- 0
        for(i in 1:540){
          if(data$matchdist0[i] == 1 & data$matchdist1[i] == 1){
            data$matchdistboth[i] <- 1
          }
          if(data$matchdist0[i] != 1 & data$matchdist1[i] != 1){
            data$matchdistneither[i] <- 1
          }
          
        }
        
        
        x<-  sum(data$matchdist0)
        y<-  sum(data$matchdist1)
        z<-  sum(data$matchdistboth)
        w<- sum(data$matchdistneither)
        dist_vect <- c(x,y,z,w) 
        
        results <- data.frame(Dist = c("Correct dist0", "Correct dist1",
                                       "Both dist correct", "Neither dist correct"),
                              Numb_correct = dist_vect, Perc = dist_vect/540)
        write.csv(results, "Mix_results_dist_perc_correct_220113.csv")
        
        
        # per distribution combination
        #dist0<-"gamma"
        #dist1<-"weibull"
        #subset <- data[data$dist0_true==dist0 & data$dist1_true==dist1, ]
        #x<-  sum(subset$matchdist0[subset$matchdist1==0])
        #y<-  sum(subset$matchdist1[subset$matchdist0==0])
        #z<-  sum(subset$matchdistboth)
        #w<- sum(subset$matchdistneither)
        #dist_vect <- c(x,y,z,w) 
        
        #results <- data.frame(Dist = c("Correct dist0", "Correct dist1",
        #         "Both dist correct", "Neither dist correct"),
        #Numb_correct = dist_vect, Perc = dist_vect/nrow(subset))
        
        
        
   ## distribution parameter specification by mixture model

        data<-merge(true,est[,-c(1:2)], by="sim")
        data$within <- 0
        
        dataFOI <- data[, c(1,4,31:34)]
        dataFOI$Parameter <- "FOI"
        colnames(dataFOI)[2:5] <- c("True", "Est", "low", "upp") 
        for(i in 1:540){
          if(is.na(dataFOI$Est[i])!=TRUE){
            if(dataFOI$True[i] > dataFOI$low[i] & dataFOI$True[i] < dataFOI$upp[i]) dataFOI$within[i] <- 1
          }
        }
        
        dataseroprev <- data[, c(1,9,28:30,34)]
        dataseroprev$Parameter <- "Seroprevalence"
        colnames(dataseroprev)[2:5] <- c("True", "Est", "low", "upp") 
        for(i in 1:540){
          if(is.na(dataseroprev$Est[i])!=TRUE){
            if(dataseroprev$True[i] > dataseroprev$low[i] & dataseroprev$True[i] < dataseroprev$upp[i]) dataseroprev$within[i] <- 1
          }
        }
        
        datamu0 <- data[, c(1,5,12:14,34)] 
        datamu0$Parameter <- "muS"
        colnames(datamu0)[2:5] <- c("True", "Est", "low", "upp")
        for(i in 1:540){
          if(is.na(datamu0$Est[i])!=T&is.na(datamu0$low[i])!=T&is.na(datamu0$upp[i])!=T){
            if(datamu0$True[i] > datamu0$low[i] & 
               datamu0$True[i] < datamu0$upp[i])datamu0$within[i] <- 1
          }
        }
        
        datamu1 <- data[, c(1,6,15:17,34)]
        datamu1$Parameter <- "muI"
        colnames(datamu1)[2:5] <- c("True", "Est", "low", "upp") 
        for(i in 1:540){
          if(is.na(datamu1$Est[i])!=T&is.na(datamu0$low[i])!=T&is.na(datamu0$upp[i])!=T){
            if(datamu1$True[i] > datamu1$low[i] & datamu1$True[i] < datamu1$upp[i]) datamu1$within[i] <- 1
          }
        }
        
        datasd0 <- data[, c(1,7,18:20,34)]
        datasd0$Parameter <- "sdS"
        colnames(datasd0)[2:5] <- c("True", "Est", "low", "upp") 
        for(i in 1:540){
          if(is.na(datasd0$Est[i])!=T&is.na(datamu0$low[i])!=T&is.na(datamu0$upp[i])!=T){
            if(datasd0$True[i] > datasd0$low[i] & datasd0$True[i] < datasd0$upp[i]) datasd0$within[i] <- 1
          }
        }
        
        datasd1 <- data[, c(1,8,21:23,34)]
        datasd1$Parameter <- "sdI"
        colnames(datasd1)[2:5] <- c("True", "Est", "low", "upp") 
        for(i in 1:540){
          if(is.na(datasd1$Est[i])!=T&is.na(datamu0$low[i])!=T&is.na(datamu0$upp[i])!=T){
            if(datasd1$True[i] > datasd1$low[i] & datasd1$True[i] < datasd1$upp[i]) datasd1$within[i] <- 1
          }
        }
        
        x<- sum(dataFOI$within)
        y<-sum(dataseroprev$within)
        z<-sum(datasd1$within)
        w<-sum(datasd0$within)
        v<-sum(datamu1$within)
        u<-sum(datamu0$within)
        
        param_vect <- c(x,y,z,w, v, u) 
        
        results <- data.frame(Param = c("FOI", "SP","sd1", "sd0", "mu1", "mu0"),
                              Numb_correct = param_vect, Perc = param_vect/540)
        write.csv(results, "Mix_results_param_perc_correct.csv")
        
      
        newdata <- rbind(datasd1, datasd0, datamu1, datamu0)#, dataseroprev, dataFOI)
        newdata$within <- as.character(newdata$within)
        newdata$Parameter <- factor(newdata$Parameter, 
                                    levels = c("muI","muS","sdI","sdS"),
                                    labels=c(expression(paste(mu,"I")),expression(paste(mu,"S")),
                                             expression(paste(sigma,"I")),expression(paste(sigma,"S"))))
        mix_plot2 <- ggplot(transform(newdata,
                                      Parameter=factor(Parameter,
                                                       levels=c(expression(paste(mu,"S")),expression(paste(mu,"I")),
                                                                expression(paste(sigma,"S")),expression(paste(sigma,"I")))))) + 
          facet_wrap("Parameter", scales="free", label = "label_parsed")+
          geom_point(aes(x=True, y=Est, col=within), size=3, alpha=0.4) +
          geom_errorbar(aes(x=True, ymin=low, ymax=upp, col=within), alpha=0.4)+
          theme_bw() +
          labs(x="True values", y="Estimated values", ) +
          geom_abline(slope=1, linetype="dashed", size=1) +
          theme(axis.text = element_text(size=15), axis.title = element_text(size=20),
                strip.text = element_text(size=15), legend.position = "none") +
          scale_y_continuous(breaks = scales::pretty_breaks(7), limits = c(0, NA)) +
          scale_x_continuous(breaks = scales::pretty_breaks(7), limits = c(NA, NA)) 
        
        ggsave(plot = mix_plot2, filename = "Mix_results_param_est_vs_true2.png")
        