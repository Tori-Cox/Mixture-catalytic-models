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
            for(i in 1:540){
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
            for(i in 1:540){
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
            for(i in 1:(2*540)){
              if(newdata$True[i] <= newdata$Est_upp[i] & newdata$True[i] >= newdata$Est_low[i]){
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
            
            newdata <- newdata[! substring(newdata$sim, 5,8) %in% to_remove$sim,]
            
            newdata2<-newdata[newdata$Parameter=="FOI",]
            x <- sum(newdata2$within)/(nrow(newdata)/2)
            newdata2<-newdata[newdata$Parameter!="FOI",]
            y <- sum(newdata2$within)/(nrow(newdata)/2)
            results <- data.frame(Param = c("FOI", "Serop"), Perc = c(x,y))
            write.csv(results, "Cat_V_results_est_perc_correct.csv")
            
            
            # store for later
            CATVDATA <- newdata
            saveRDS(CATVDATA, "Results_CATVDATA.RDS")
            
            
            
## calculating bias and uncertainty for time constant catalytic model:
            est_foi_constant<-list()
            est_sp_constant<-list()
            for(i in 1:540){
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
            for(i in 1:540){
              est_foi_central[i]<-est_foi_constant[[i]][1]
              est_foi_lower[i]<-est_foi_constant[[i]][2]
              est_foi_upper[i]<-est_foi_constant[[i]][3]
              est_serop_central[i]<-est_sp_constant[[i]][1]
              est_serop_lower[i]<-est_sp_constant[[i]][2]
              est_serop_upper[i]<-est_sp_constant[[i]][3]
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
            for(i in 1:(2*540)){
              if(newdata$True[i] <= newdata$Est_upp[i] & newdata$True[i] >= newdata$Est_low[i]){
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
            
            newdata <- newdata[! substring(newdata$sim, 5,8) %in% to_remove$sim,]
            
            newdata2<-newdata[newdata$Parameter=="FOI",]
            x <- sum(newdata2$within)/(nrow(newdata)/2)
            newdata2<-newdata[newdata$Parameter!="FOI",]
            y <- sum(newdata2$within)/(nrow(newdata)/2)
            results <- data.frame(Param = c("FOI", "Serop"), Perc = c(x,y))
            write.csv(results, "Cat_C_results_est_perc_correct.csv")
            
            
            # store for later
            CATCDATA <- newdata
            saveRDS(CATCDATA, "Results_CATCDATA.RDS")
            

############################################
# plotting estimated versus true estimates #
############################################
            
            to_remove<-readRDS("MixtureSims_sim_remove.RDS")
            CATCDATA2 <- CATCDATA[! substring(CATCDATA$sim, 5,8) %in% to_remove$sim,]
            CATVDATA2 <- CATVDATA[! substring(CATVDATA$sim, 5,8) %in% to_remove$sim,]
            MIXDATA2 <- MIXDATA[! substring(MIXDATA$sim, 5,8) %in% to_remove$sim,]
            

            # time-constant
            con_plot <- ggplot(CATCDATA2) + facet_wrap("Parameter", scales="free")+
              geom_point(aes(x=True, y=Est, col=within), size=2, alpha=0.4) +
              geom_errorbar(aes(x=True, ymin=Est_low, ymax=Est_upp, col=within), alpha=0.4)+
              theme_bw() +
              labs(x="True values", y="Estimated values") +
              geom_abline(slope=1, linetype="dashed", size=1) +
              coord_cartesian(ylim=c(0,1))+
              theme(axis.text = element_text(size=15), axis.title = element_text(size=20),
                    strip.text = element_text(size=15), legend.position = "none")
            
            ggsave(plot = con_plot, filename = "Cat_C_results_est_vs_true.png")
            
            # time-varying
            var_plot <- ggplot(CATVDATA2) + facet_wrap("Parameter", scales="free")+
              geom_point(aes(x=True, y=Est, col=within), size=2, alpha=0.4) +
              geom_errorbar(aes(x=True, ymin=Est_low, ymax=Est_upp, col=within), alpha=0.4)+
              theme_bw() +
              labs(x="True values", y="Estimated values") +
              geom_abline(slope=1, linetype="dashed", size=1) +
              coord_cartesian(ylim=c(0,1))+
              theme(axis.text = element_text(size=15), axis.title = element_text(size=20),
                    strip.text = element_text(size=15), legend.position = "none")
            
            ggsave(plot = var_plot, filename = "Cat_V_results_est_vs_true.png")
            
            
############################################
### plotting misclassification of titres ###
############################################
            
            misclass <- readRDS("MixtureSims_correct_classification_vector.RDS")
            CATVDATA$misclass<- c((1-misclass), (1-misclass))
            CATCDATA$misclass<- c((1-misclass), (1-misclass))
            
            true<- read.csv("MixtureSimParams.csv")
            
            data <- data.frame(sim= 1:540,misclass=(1-misclass), Diff = (true$mu1 - true$mu0))
            quantile(data$misclass,0.9)
            data2 <- data[! data$sim %in% to_remove$sim,]
            quantile(data2$misclass,0.9)
            
            ggplot(data2) +
              geom_point(aes(x=Diff,y=misclass))+
              labs(x="Difference between true mean titres\n of seronegative and seropositive components", 
                   y="Percentage of titres\n misclassified")+
              theme_bw() +
              geom_smooth(method="loess",aes(x=Diff,y=misclass), se=T)+
              ylim(0,0.3)+
              theme(axis.text = element_text(size=15), 
                    axis.title = element_text(size=17),
                    strip.text = element_text(size=15))
            
            ggsave(filename = "Cat_misclass_versus_titre_separation.png")
            
            CATCDATA$Model2 <- rep(c("FOI (Time-constant)","Seroprevalence (Time-constant)"), each=540)
            CATVDATA$Model2 <- rep(c("FOI (Time-varying)","Seroprevalence (Time-varying)"), each=540)
            CATCDATA2 <- CATCDATA[! substring(CATCDATA$sim, 5,8) %in% to_remove$sim,]
            CATVDATA2 <- CATVDATA[! substring(CATVDATA$sim, 5,8) %in% to_remove$sim,]
            comb <- rbind(CATCDATA2, CATVDATA2)
            
            comb$absbias <- abs(comb$bias)
            #comb<-comb[comb$absbias<0.5, ]
            
            #ggplot(comb) + facet_wrap("Model2", scales="free")+
            # geom_jitter(aes(x=misclass, y=absbias))+
            # theme_bw()+
            # labs(x="Absolute Bias", y="Serostatus misclassification rate")+
            # geom_smooth(method = "lm",aes(x=misclass, y=absbias))
            
            ggplot(comb) + facet_wrap("Model2", scales="free")+
              geom_jitter(aes(x=bias, y=misclass))+
              theme_bw()+
              labs(x="Bias (Estimated value - True value)", y="Serostatus misclassification rate")+
              geom_vline(xintercept = 0, col="red", size=1.5)
            
            comb2<-comb[comb$Est<0.4|comb$Parameter!="FOI",]
            ggplot(comb2) + facet_wrap("Model2", scales="free")+
              geom_jitter(aes(x=absbias, y=misclass))+
              theme_bw()+
              labs(x="Absolute Bias", y="Serostatus misclassification rate")+
              stat_cor(aes(y=misclass, x=absbias, label = ..r.label..), label.y= c(0.55,0.95)) +
              geom_smooth(method = "lm",aes(y=misclass, x=absbias))+ylim(0,0.6)
            
            ggsave(filename = "Cat_misclass_versus_bias_abs.png")
            
            