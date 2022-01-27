###############################################
## Simulation study - exploration of results ##
## 3. Both models  ############################
###############################################


## read in processed model output results (see explore_results_1.R and explore_results_2.R)
      MIXDATA <- readRDS("results_MIXDATA.RDS") # mixture model
      CATCDATA <- readRDS("Results_CATCDATA.RDS") # time-constant FOI catalytic model
      CATVDATA <- readRDS("Results_CATVDATA.RDS") # time-variable FOI catalytic model
      
      MIXDATA$model<-"Mixture model"
      CATCDATA$model<-"Time-constant catalytic model"
      CATVDATA$model<-"Time-varying catalytic model"
      
      newdata <- rbind(CATCDATA,CATVDATA,MIXDATA)
      newdata2<-subset(newdata, newdata$model!="Time-constant catalytic model")
      newdata3<-subset(newdata, newdata$model!="Time-varying catalytic model")
      
      
#######################################
## plot comparing est vs true values ##
#######################################
      to_remove<-readRDS("MixtureSims_sim_remove.RDS")
      newdata <- rbind(CATCDATA[! substring(CATCDATA$sim, 5,8) %in% to_remove$sim,],
                       CATVDATA[! substring(CATVDATA$sim, 5,8) %in% to_remove$sim,],MIXDATA)
      newdata2<-subset(newdata, newdata$model!="Time-constant catalytic model")
      newdata3<-subset(newdata, newdata$model!="Time-varying catalytic model")
      
      newdata2<-newdata2[newdata2$Est<=1, ]
      newdata3<-newdata3[newdata3$Est<=1, ]
      newdata2<-newdata2[is.na(newdata2$Est)==F, ]
      newdata3<-newdata3[is.na(newdata3$Est)==F, ]
      
      
      ## plot comparing estimated versus true FOI and seroprevalence values
      #time varying catalytic model
      a<- 
        ggscatter(newdata2, x = "True", y = "Est",
                  add = "reg.line",                         
                  color = "model", size=1.5, alpha=0.7)+
        stat_cor(aes(color = model,label = ..r.label..), label.y= c(0.85,0.95)) +
        facet_wrap("Parameter", scales="free")+
        geom_errorbar(aes(x=True, ymin=Est_low, ymax=Est_upp, col=model), alpha=0.4)+
        theme_bw()+
        labs(x="True values", y="Estimated values", col=NULL) +
        geom_abline(slope=1, linetype="dashed", size=1) +
        coord_cartesian(ylim=c(0,1))+
        theme(axis.text = element_text(size=14),
              axis.title = element_text(size=15),
              legend.text=element_text(size=14),
              strip.text = element_text(size=15),
              legend.position = "bottom")+
        scale_color_manual(breaks = c("Mixture model",
                                      "Time-varying catalytic model"),
                           values=c("#F8766D","#00BFC4" ))+
        guides(color = guide_legend(override.aes = list(size=3,linetype=0)))
      
      
      #time constant catalytic model
      b<-
        ggscatter(newdata3, x = "True", y = "Est",
                  add = "reg.line",                         
                  color = "model", size=1.5, alpha=0.7)+
        stat_cor(aes(color = model,label = ..r.label..), label.y= c(0.85,0.95)) +
        facet_wrap("Parameter", scales="free")+
        geom_errorbar(aes(x=True, ymin=Est_low, ymax=Est_upp, col=model), alpha=0.4)+
        theme_bw()+
        labs(x="True values", y="Estimated values", col=NULL) +
        geom_abline(slope=1, linetype="dashed", size=1) +
        coord_cartesian(ylim=c(0,1))+
        theme(axis.text = element_text(size=14),
              axis.title = element_text(size=15),
              legend.text=element_text(size=14),
              strip.text = element_text(size=15),
              legend.position = "bottom")+
        scale_color_manual(breaks = c("Mixture model",
                                      "Time-constant catalytic model"),
                           values=c("#F8766D","#00BA38"))+
        guides(color = guide_legend(override.aes = list(size=3,linetype=0)))
      
      png(filename = "Comp_results_est_vs_true.png")
      z <- grid.arrange(b,a)
      dev.off()
      
###################################
## plots of bias and uncertainty ##
###################################
      
      summary_stats<- data.frame(Parameter = rep(c("FOI", "Seroprevalence"), each = 9),
                                 Model = rep(c("Mixture model", 
                                               "Time-constant Catalytic model",
                                               "Time-varying Catalytic model"), times=3),
                                 Stat = rep(c("Bias", "Uncertainty", "Coverage"), each=3),
                                 Mean=NA, 
                                 Low=NA,Upp=NA, 
                                 Count=NA)
      CATCDATA2 <- CATCDATA[! substring(CATCDATA$sim, 5,8) %in% to_remove$sim,]
      CATVDATA2 <- CATVDATA[! substring(CATVDATA$sim, 5,8) %in% to_remove$sim,]
      MIXDATA2 <- MIXDATA[! MIXDATA$sim %in% to_remove$sim,]
      
      
      summary_stats$Count[summary_stats$Parameter=="FOI" & 
                            summary_stats$Model=="Mixture model" &
                            summary_stats$Stat=="Coverage"] <-  
        sum(MIXDATA2$within[MIXDATA2$Parameter=="FOI"])
      summary_stats$Count[summary_stats$Parameter=="FOI" & 
                            summary_stats$Model=="Time-constant Catalytic model" &
                            summary_stats$Stat=="Coverage"] <-  
        sum(CATCDATA2$within[CATCDATA2$Parameter=="FOI"])
      summary_stats$Count[summary_stats$Parameter=="FOI" & 
                            summary_stats$Model=="Time-varying Catalytic model" &
                            summary_stats$Stat=="Coverage"] <-  
        sum(CATVDATA2$within[CATVDATA2$Parameter=="FOI"])
      summary_stats$Count[summary_stats$Parameter!="FOI" & 
                            summary_stats$Model=="Mixture model" &
                            summary_stats$Stat=="Coverage"] <-  
        sum(MIXDATA2$within[MIXDATA2$Parameter!="FOI"])
      summary_stats$Count[summary_stats$Parameter!="FOI" & 
                            summary_stats$Model=="Time-constant Catalytic model" &
                            summary_stats$Stat=="Coverage"] <-  
        sum(CATCDATA2$within[CATCDATA2$Parameter!="FOI"])
      summary_stats$Count[summary_stats$Parameter!="FOI" & 
                            summary_stats$Model=="Time-varying Catalytic model" &
                            summary_stats$Stat=="Coverage"] <-  
        sum(CATVDATA2$within[CATVDATA2$Parameter!="FOI"])
      
      
      summary_stats$Mean[summary_stats$Parameter=="FOI" & 
                           summary_stats$Model=="Mixture model" &
                           summary_stats$Stat=="Uncertainty"] <-  
        mean(MIXDATA2$uncertainty[MIXDATA2$Parameter=="FOI"],na.rm=T)
      summary_stats$Mean[summary_stats$Parameter=="FOI" & 
                           summary_stats$Model=="Time-constant Catalytic model" &
                           summary_stats$Stat=="Uncertainty"] <-  
        mean(CATCDATA2$uncertainty[CATCDATA2$Parameter=="FOI"],na.rm=T) 
      summary_stats$Mean[summary_stats$Parameter=="FOI" & 
                           summary_stats$Model=="Time-varying Catalytic model" &
                           summary_stats$Stat=="Uncertainty"] <-  
        mean(CATVDATA2$uncertainty[CATVDATA2$Parameter=="FOI"],na.rm=T) 
      summary_stats$Mean[summary_stats$Parameter!="FOI" & 
                           summary_stats$Model=="Mixture model" &
                           summary_stats$Stat=="Uncertainty"] <-  
        mean(MIXDATA2$uncertainty[MIXDATA2$Parameter!="FOI"],na.rm=T) 
      summary_stats$Mean[summary_stats$Parameter!="FOI" & 
                           summary_stats$Model=="Time-constant Catalytic model" &
                           summary_stats$Stat=="Uncertainty"] <-  
        mean(CATCDATA2$uncertainty[CATCDATA2$Parameter!="FOI"],na.rm=T) 
      summary_stats$Mean[summary_stats$Parameter!="FOI" & 
                           summary_stats$Model=="Time-varying Catalytic model" &
                           summary_stats$Stat=="Uncertainty"] <-  
        mean(CATVDATA2$uncertainty[CATVDATA2$Parameter!="FOI"],na.rm=T) 
      
      
      summary_stats$Mean[summary_stats$Parameter=="FOI" & 
                           summary_stats$Model=="Mixture model" &
                           summary_stats$Stat=="Bias"] <-  
        mean(MIXDATA2$bias[MIXDATA2$Parameter=="FOI"],na.rm=T)
      summary_stats$Mean[summary_stats$Parameter=="FOI" & 
                           summary_stats$Model=="Time-constant Catalytic model" &
                           summary_stats$Stat=="Bias"] <-  
        mean(CATCDATA2$bias[CATCDATA2$Parameter=="FOI"],na.rm=T) 
      summary_stats$Mean[summary_stats$Parameter=="FOI" & 
                           summary_stats$Model=="Time-varying Catalytic model" &
                           summary_stats$Stat=="Bias"] <-  
        mean(CATVDATA2$bias[CATVDATA2$Parameter=="FOI"],na.rm=T) 
      summary_stats$Mean[summary_stats$Parameter!="FOI" & 
                           summary_stats$Model=="Mixture model" &
                           summary_stats$Stat=="Bias"] <-  
        mean(MIXDATA2$bias[MIXDATA2$Parameter!="FOI"],na.rm=T) 
      summary_stats$Mean[summary_stats$Parameter!="FOI" & 
                           summary_stats$Model=="Time-constant Catalytic model" &
                           summary_stats$Stat=="Bias"] <-  
        mean(CATCDATA2$bias[CATCDATA2$Parameter!="FOI"],na.rm=T) 
      summary_stats$Mean[summary_stats$Parameter!="FOI" & 
                           summary_stats$Model=="Time-varying Catalytic model" &
                           summary_stats$Stat=="Bias"] <-  
        mean(CATVDATA2$bias[CATVDATA2$Parameter!="FOI"],na.rm=T) 
      
      
      summary_stats[summary_stats$Parameter=="FOI" & 
                      summary_stats$Model=="Mixture model" &
                      summary_stats$Stat=="Bias",][,c(5:6)] <-  
        quantile((MIXDATA2[MIXDATA2$Parameter=="FOI",]$bias), probs=c(0.025,0.975), na.rm=T)
      summary_stats[summary_stats$Parameter=="FOI" & 
                      summary_stats$Model=="Time-constant Catalytic model" &
                      summary_stats$Stat=="Bias",][,c(5:6)] <-
        quantile((CATCDATA2[CATCDATA2$Parameter=="FOI",]$bias), probs=c(0.025,0.975), na.rm=T)
      summary_stats[summary_stats$Parameter=="FOI" & 
                      summary_stats$Model=="Time-varying Catalytic model" &
                      summary_stats$Stat=="Bias",][,c(5:6)] <-
        quantile((CATVDATA2[CATVDATA2$Parameter=="FOI",]$bias), probs=c(0.025,0.975), na.rm=T)
      summary_stats[summary_stats$Parameter!="FOI" & 
                      summary_stats$Model=="Mixture model" &
                      summary_stats$Stat=="Bias",][,c(5:6)] <-  
        quantile((MIXDATA2[MIXDATA2$Parameter!="FOI",]$bias), probs=c(0.025,0.975), na.rm=T)
      summary_stats[summary_stats$Parameter!="FOI" & 
                      summary_stats$Model=="Time-constant Catalytic model" &
                      summary_stats$Stat=="Bias",][,c(5:6)] <-
        quantile((CATCDATA2[CATCDATA2$Parameter!="FOI",]$bias), probs=c(0.025,0.975), na.rm=T)
      summary_stats[summary_stats$Parameter!="FOI" & 
                      summary_stats$Model=="Time-varying Catalytic model" &
                      summary_stats$Stat=="Bias",][,c(5:6)] <-
        quantile((CATVDATA2[CATVDATA2$Parameter!="FOI",]$bias), probs=c(0.025,0.975), na.rm=T)
      
      summary_stats[summary_stats$Parameter=="FOI" & 
                      summary_stats$Model=="Mixture model" &
                      summary_stats$Stat=="Uncertainty",][,c(5:6)] <-  
        quantile((MIXDATA2[MIXDATA2$Parameter=="FOI",]$uncertainty), probs=c(0.025,0.975), na.rm=T)
      summary_stats[summary_stats$Parameter=="FOI" & 
                      summary_stats$Model=="Time-constant Catalytic model" &
                      summary_stats$Stat=="Uncertainty",][,c(5:6)] <-
        quantile((CATCDATA2[CATCDATA2$Parameter=="FOI",]$uncertainty), probs=c(0.025,0.975), na.rm=T)
      summary_stats[summary_stats$Parameter=="FOI" & 
                      summary_stats$Model=="Time-varying Catalytic model" &
                      summary_stats$Stat=="Uncertainty",][,c(5:6)] <-
        quantile((CATVDATA2[CATVDATA2$Parameter=="FOI",]$uncertainty), probs=c(0.025,0.975), na.rm=T)
      summary_stats[summary_stats$Parameter!="FOI" & 
                      summary_stats$Model=="Mixture model" &
                      summary_stats$Stat=="Uncertainty",][,c(5:6)] <-  
        quantile((MIXDATA2[MIXDATA2$Parameter!="FOI",]$uncertainty), probs=c(0.025,0.975), na.rm=T)
      summary_stats[summary_stats$Parameter!="FOI" & 
                      summary_stats$Model=="Time-constant Catalytic model" &
                      summary_stats$Stat=="Uncertainty",][,c(5:6)] <-
        quantile((CATCDATA2[CATCDATA2$Parameter!="FOI",]$uncertainty), probs=c(0.025,0.975), na.rm=T)
      summary_stats[summary_stats$Parameter!="FOI" & 
                      summary_stats$Model=="Time-varying Catalytic model" &
                      summary_stats$Stat=="Uncertainty",][,c(5:6)] <-
        quantile((CATVDATA2[CATVDATA2$Parameter!="FOI",]$uncertainty), probs=c(0.025,0.975), na.rm=T)
      
      
      summary_stats$Model <- factor(summary_stats$Model, 
                                    levels=c("Mixture model", 
                                             "Time-constant Catalytic model", 
                                             "Time-varying Catalytic model"))
      summary_stats_cov<-summary_stats[summary_stats$Stat=="Coverage",]
      for(i in 1:6){
        summary_stats_cov$Upp[i] <- (binom.exact(summary_stats_cov$Count[i], 509, conf.level = .95)$upper) *100
        summary_stats_cov$Low[i] <- (binom.exact(summary_stats_cov$Count[i], 509, conf.level = .95)$lower) *100
      }
      
      write.csv(rbind(summary_stats[summary_stats$Stat!="Coverage",], 
                      summary_stats_cov), "Comp_results_bias_CIs.csv")
      
      
      not_cov <- ggplot(summary_stats[summary_stats$Parameter=="FOI" & summary_stats$Stat!="Coverage",]) + 
        facet_wrap("Parameter", scales="free")+
        geom_point(aes(x=Stat,y=Mean,col=Model),size=2, pch = 15,position = position_dodge(width = 0.6)) +
        geom_errorbar(aes(x=Stat,ymin=Low, ymax=Upp,col=Model),size=1,position = position_dodge(width = 0.6), width=0.4) +
        theme_bw() +
        labs(x=NULL, y="Value", col=NULL) +
        theme(axis.text = element_text(size=15),
              axis.title = element_text(size=15),
              legend.text=element_text(size=13),
              strip.text = element_text(size=15),
              legend.position =  c(0.3,0.8))+
        ylim(-0.5,1)
      not_cov2 <- ggplot(summary_stats[summary_stats$Parameter!="FOI" & summary_stats$Stat!="Coverage",]) + 
        facet_wrap("Parameter", scales="free")+
        geom_point(aes(x=Stat,y=Mean,col=Model),size=2, pch = 15,position = position_dodge(width = 0.6)) +
        geom_errorbar(aes(x=Stat,ymin=Low, ymax=Upp,col=Model),size=1,position = position_dodge(width = 0.6), width=0.4) +
        theme_bw() +
        labs(x=NULL, y="Value", col=NULL) +
        theme(axis.text = element_text(size=15),
              axis.title = element_text(size=15),
              legend.text=element_text(size=13),
              strip.text = element_text(size=15),
              legend.position = "none")+
        ylim(-0.5,1)
      cov <- ggplot(summary_stats_cov[summary_starts_cov$Parameter=="FOI",]) + facet_wrap("Parameter")+
        geom_point(aes(x=Stat,y=(Count/509)*100,col=Model),size=2,pch=15, position = position_dodge(width = 0.6))+
        geom_errorbar(aes(x=Stat,ymin=Low, ymax=Upp,col=Model),size=1,position = position_dodge(width = 0.6), width=0.4) +
        theme_bw() +
        labs(x=NULL, y="Percentage of \nsimulations (%)", col=NULL) +
        theme(axis.text = element_text(size=15),
              axis.title = element_text(size=15),
              legend.text=element_text(size=15),
              legend.title=element_text(size=15),
              strip.text = element_text(size=15),
              legend.position = "none")+
        geom_hline(yintercept=95, linetype="dashed")+
        ylim(0,100)
      cov2 <- ggplot(summary_stats_cov[summary_starts_cov$Parameter!="FOI",]) + facet_wrap("Parameter")+
        geom_point(aes(x=Stat,y=(Count/509)*100,col=Model),
                   size=2,pch=15, 
                   position = position_dodge(width = 0.6))+
        geom_errorbar(aes(x=Stat,ymin=Low, ymax=Upp,col=Model),size=1,position = position_dodge(width = 0.6), width=0.4) +
        theme_bw() +
        labs(x=NULL, y="Percentage of \nsimulations (%)", col=NULL) +
        theme(axis.text = element_text(size=15),
              axis.title = element_text(size=15),
              legend.text=element_text(size=15),
              legend.title=element_text(size=15),
              strip.text = element_text(size=15),
              legend.position = "none")+
        geom_hline(yintercept=95, linetype="dashed")+
        ylim(0,100)
      
      plot <- grid.arrange(not_cov, cov, not_cov2, cov2, ncol=2,nrow=2, widths=c(4,2))
      ggsave(plot, filename = "Sim_bias_CI_cov.png")
      