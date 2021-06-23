###############################################
## Simulation study - exploration of results ##
## 2. Catalytic models ########################
###############################################


##################################
# calculating bias and uncertainty
##################################


## calculating bias and uncertainty for time-varying FOI catalytic model:
            est_foi_variable<-list()
            est_sp_variable<-list()
            for(i in 1:200){
              summary <- read.csv(paste("model_output/catalytic/Cat_", i, "_SummaryResults.csv", sep=""))
              
              foi<- c(summary$Value[4], summary$Value[5], summary$Value[6])
              sp<- c(summary$Value[10], summary$Value[11], summary$Value[12])
              
              est_foi_variable[[i]]<- foi
              est_sp_variable[[i]]<-sp
            }
            
            est_foi_central<-NA
            est_foi_lower<-NA
            est_foi_upper<-NA
            est_serop_central<-NA
            est_serop_lower<-NA
            est_serop_upper<-NA
            for(i in 1:200){
              est_foi_central[i]<-est_foi_variable[[i]][1]
              est_foi_lower[i]<-est_foi_variable[[i]][2]
              est_foi_upper[i]<-est_foi_variable[[i]][3]
              est_serop_central[i]<-est_sp_variable[[i]][1]
              est_serop_lower[i]<-est_sp_variable[[i]][2]
              est_serop_upper[i]<-est_sp_variable[[i]][3]
            }
            
            true_foi<- read.csv("MixtureSimParams.csv")
            data <- true_foi
            data <- data[, c(1, 4)]
            colnames(data)[2]<-"True"
            data$Est <- est_foi_central
            data$Est_low <- est_foi_lower
            data$Est_upp <- est_foi_upper
            data$Parameter<-"FOI"
            
            data2 <- true_foi
            data2 <- data2[, c(1, 9)]
            colnames(data2)[2]<-"True"
            data2$Est <- est_serop_central
            data2$Est_low <- est_serop_lower
            data2$Est_upp <- est_serop_upper
            data2$Parameter<-"Seroprevalence"
            
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
            
            newdata2<-newdata[newdata$Parameter=="FOI",]
            x <- sum(newdata2$within)/200
            newdata2<-newdata[newdata$Parameter!="FOI",]
            y <- sum(newdata2$within)/200
            results <- data.frame(Param = c("FOI", "Serop"), Perc = c(x,y))
            write.csv(results, "Cat_V_results_est_perc_correct.csv")
            
            
            # store for later
            CATVDATA <- newdata
            saveRDS(CATVDATA, "Results_CATVDATA.RDS")
            
             

            
            
            
## calculating bias and uncertainty for time constant catalytic model:
            est_foi_constant<-list()
            est_sp_constant<-list()
            for(i in 1:200){
              summary <- read.csv(paste("model_output/catalytic/Cat_", i, "_SummaryResults.csv", sep=""))
              
              foi<- c(summary$Value[1], summary$Value[2], summary$Value[3])
              sp<- c(summary$Value[7], summary$Value[8], summary$Value[9])
              
              est_foi_constant[[i]]<- foi
              est_sp_constant[[i]]<-sp
            }
            
            est_foi_central<-NA
            est_foi_lower<-NA
            est_foi_upper<-NA
            est_serop_central<-NA
            est_serop_lower<-NA
            est_serop_upper<-NA
            for(i in 1:200){
              est_foi_central[i]<-est_foi_constant[[i]][1]
              est_foi_lower[i]<-est_foi_constant[[i]][2]
              est_foi_upper[i]<-est_foi_constant[[i]][3]
              est_serop_central[i]<-est_sp_constant[[i]][1]
              est_serop_lower[i]<-est_sp_constant[[i]][2]
              est_serop_upper[i]<-est_sp_constant[[i]][3]
            }
            
            # true values processed above for time-varying model
            newdata <- newdata[, 1:6]
            
            
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
            
            newdata2<-newdata[newdata$Parameter=="FOI",]
            x<- sum(newdata2$within)/200
            newdata2<-newdata[newdata$Parameter!="FOI",]
            y <- sum(newdata2$within)/200
            
            results <- data.frame(Param = c("FOI", "Serop"), Perc = c(x,y))
            write.csv(results, "Cat_C_results_est_perc_correct.csv")
            
            # store for later
            CATCDATA <- newdata
            saveRDS(CATCDATA, "Results_CATCDATA.RDS")



############################################
# plotting estimated versus true estimates #
############################################


# time-constant
con_plot <- ggplot(CATCDATA) + facet_wrap("Parameter", scales="free")+
  geom_point(aes(x=True, y=Est, col=within), size=2, alpha=0.4) +
  geom_errorbar(aes(x=True, ymin=Est_low, ymax=Est_upp, col=within), alpha=0.4)+
  theme_bw() +
  labs(x="True values", y="Estimated values") +
  geom_abline(slope=1, linetype="dashed", size=1) +
  theme(axis.text = element_text(size=15), axis.title = element_text(size=20),
        strip.text = element_text(size=15), legend.position = "none")

ggsave(plot = con_plot, filename = "Cat_C_results_est_vs_true.png")

# time-varying
var_plot <- ggplot(CATVDATA) + facet_wrap("Parameter", scales="free")+
  geom_point(aes(x=True, y=Est, col=within), size=2, alpha=0.4) +
  geom_errorbar(aes(x=True, ymin=Est_low, ymax=Est_upp, col=within), alpha=0.4)+
  theme_bw() +
  labs(x="True values", y="Estimated values") +
  geom_abline(slope=1, linetype="dashed", size=1) +
  theme(axis.text = element_text(size=15), axis.title = element_text(size=20),
        strip.text = element_text(size=15), legend.position = "none")

ggsave(plot = var_plot, filename = "Cat_V_results_est_vs_true.png")

