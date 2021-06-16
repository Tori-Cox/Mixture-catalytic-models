###############################################
## Simulation study - exploration of results ##
## 1. Mixture model ###########################
###############################################

## load packages
library(ggplot2)
library(scales)

## read in data
est <- read.csv("Mix_fit_and_est_params.csv") #estimated by mixture model
true<- read.csv("MixtureSimParams.csv") #simulated 'true' parameters
true_foi <- true$foi
true_serop <- true$seroprev


##################################
# calculating bias and uncertainty
##################################

        
        data <- data.frame(sim=1:200, True=true_foi, Est= est$foi, Est_low = est$foi_low,
                           Est_upp = est$foi_upp, Parameter = "FOI")
        
        data2 <- data.frame(sim=1:200, True=true_serop, Est= est$seroprev, Est_low = est$seroprev_low,
                            Est_upp = est$seroprev_upp, Parameter = "Seroprevalence")
        
        newdata <- rbind(data,data2)
        
        
        newdata$within <-0
        for(i in 1:400){
          if(newdata$True[i] <= newdata$Est_upp[i] & newdata$True[i] >= newdata$Est_low[i]){
            newdata$within[i]<-1
          }
        }
        newdata$uncertainty <-0
        for(i in 1:400){
          newdata$uncertainty[i]<-newdata$Est_upp[i]-newdata$Est_low[i]
        }
        newdata$bias <-0
        for(i in 1:400){
          newdata$bias[i]<-newdata$Est[i]-newdata$True[i]
        }
        
        
        newdata2<-newdata[newdata$Parameter=="Seroprevalence",]
        sum(newdata2$within)/200
        newdata2<-newdata[newdata$Parameter=="FOI",]
        sum(newdata2$within)/200
        
        
        ggplot(newdata) + facet_wrap("Parameter", scales="free")+
          geom_point(aes(x=True, y=Est, col=within), size=2, alpha=0.4) +
          geom_errorbar(aes(x=True, ymin=Est_low, ymax=Est_upp, col=within), alpha=0.4)+
          theme_bw() +
          labs(x="True values", y="Estimated values") +
          geom_abline(slope=1, linetype="dashed", size=1) +
          theme(axis.text = element_text(size=15), axis.title = element_text(size=20),
                strip.text = element_text(size=15), legend.position = "none")
        
        MIXDATA <- newdata
        saveRDS(MIXDATA, "results_MIXDATA.RDS")


        
        
#########################################################
# explore estimated versus true distribution parameters # 
#########################################################
        

      ## distribution family specification by mixture model
        
        colnames(true) <- c("sim", "dist0_true", "dist1_true","foi_true" ,"mu0_true", "mu1_true", "sd0_true", "sd1_true", "seroprev_true")
        true$sim <- 1:200
        est2 <- est[,c(2:16)] #just distribution params 
        data<-cbind(true,est2)

        data$matchdist0 <- 0
        data$matchdist1 <- 0
        for(i in 1:200){
          if(data$dist0_true[i] == data$dist0[i]){
            data$matchdist0[i] <- 1
          }
          }
        for(i in 1:200){
          if(data$dist1_true[i] == data$dist1[i]){
            data$matchdist1[i] <- 1
          }
        }

        data$matchdistboth <- 0
        data$matchdistneither <- 0
        for(i in 1:200){
          if(data$matchdist0[i] == 1 & data$matchdist1[i] == 1){
            data$matchdistboth[i] <- 1
          }
          if(data$matchdist0[i] != 1 & data$matchdist1[i] != 1){
            data$matchdistneither[i] <- 1
          }
          
        }
        sum(data$matchdist0)
        sum(data$matchdist1)
        sum(data$matchdistboth)
        sum(data$matchdistneither)
        
        
        
   ## distribution parameter specification by mixture model

      data<-merge(true,est[,-1], by="sim")
      data$within <- 0
      data$mudiff<- data$mu1_true-data$mu0_true

      dataFOI <- data[, c(1,4,28:31)]
      dataFOI$Parameter <- "FOI"
      colnames(dataFOI)[2:5] <- c("True", "Est", "low", "upp") 
      for(i in 1:200){
        if(is.na(dataFOI$Est[i])!=TRUE){
         if(dataFOI$True[i] > dataFOI$low[i] & dataFOI$True[i] < dataFOI$upp[i]) dataFOI$within[i] <- 1
        }
        }
      
      dataseroprev <- data[, c(1,9,25:27,31)]
      dataseroprev$Parameter <- "Seroprevalence"
      colnames(dataseroprev)[2:5] <- c("True", "Est", "low", "upp") 
      for(i in 1:200){
       if(is.na(dataseroprev$Est[i])!=TRUE){
        if(dataseroprev$True[i] > dataseroprev$low[i] & dataseroprev$True[i] < dataseroprev$upp[i]) dataseroprev$within[i] <- 1
       }
        }
      
      datamu0 <- data[, c(1,5,12:14,31)] 
      datamu0$Parameter <- "muS"
      colnames(datamu0)[2:5] <- c("True", "Est", "low", "upp")
      for(i in 1:200){
        if(is.na(datamu0$Est[i])!=TRUE){
          if(datamu0$True[i] > datamu0$low[i] & datamu0$True[i] < datamu0$upp[i]) datamu0$within[i] <- 1
        }
        }
      
      datamu1 <- data[, c(1,6,15:17,31)]
      datamu1$Parameter <- "muI"
      colnames(datamu1)[2:5] <- c("True", "Est", "low", "upp") 
      for(i in 1:200){
        if(is.na(datamu1$Est[i])!=TRUE){
          if(datamu1$True[i] > datamu1$low[i] & datamu1$True[i] < datamu1$upp[i]) datamu1$within[i] <- 1
        }
        }
      
      datasd0 <- data[, c(1,7,18:20,31)]
      datasd0$Parameter <- "sdS"
      colnames(datasd0)[2:5] <- c("True", "Est", "low", "upp") 
      for(i in 1:200){
        if(is.na(datasd0$Est[i])!=TRUE){
          if(datasd0$True[i] > datasd0$low[i] & datasd0$True[i] < datasd0$upp[i]) datasd0$within[i] <- 1
        }
        }
      
      datasd1 <- data[, c(1,8,21:23,31)]
      datasd1$Parameter <- "sdI"
      colnames(datasd1)[2:5] <- c("True", "Est", "low", "upp") 
      for(i in 1:200){
        if(is.na(datasd1$Est[i])!=TRUE){
          if(datasd1$True[i] > datasd1$low[i] & datasd1$True[i] < datasd1$upp[i]) datasd1$within[i] <- 1
        }
        }

      sum(dataFOI$within)/200
      sum(dataseroprev$within)/200
      sum(datasd1$within)/200
      sum(datasd0$within)/200
      sum(datamu1$within)/200
      sum(datamu0$within)/200
      
      
      newdata <- rbind(datasd1, datasd0, datamu1, datamu0, dataseroprev, dataFOI)
      newdata$within <- as.character(newdata$within)
      ggplot(newdata) + facet_wrap("Parameter", scales="free")+
        geom_point(aes(x=True, y=Est, col=within), size=3, alpha=0.4) +
        geom_errorbar(aes(x=True, ymin=low, ymax=upp, col=within), alpha=0.4)+
        theme_bw() +
        labs(x="\nTrue values", y="Estimated values\n") +
        geom_abline(slope=1, linetype="dashed", size=1) +
        theme(axis.text = element_text(size=15), axis.title = element_text(size=20),
              strip.text = element_text(size=15), legend.position = "none") +
        scale_y_continuous(breaks = scales::pretty_breaks(7), limits = c(0, NA)) +
        scale_x_continuous(breaks = scales::pretty_breaks(7), limits = c(NA, NA)) 
