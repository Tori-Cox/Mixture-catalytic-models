###############################################
## Simulation study - exploration of results ##
## 3. Both models  ############################
###############################################


## load packages
library(ggplot2)
library(gridExtra)
library(grid)


## read in processed model output results (see explore_results_1.R and explore_results_2.R)
      MIXDATA <- readRDS("results_MIXDATA.RDS") # mixture model
      CATCDATA <- readRDS("results_CATCDATA.RDS") # time-constant FOI catalytic model
      CATVDATA <- readRDS("results_CATVDATA.RDS") # time-variable FOI catalytic model
      
      MIXDATA$model<-"Mixture model"
      CATCDATA$model<-"Time-constant catalytic model"
      CATVDATA$model<-"Time-varying catalytic model"
      
      newdata <- rbind(CATCDATA,CATVDATA,MIXDATA)
      newdata2<-subset(newdata, newdata$model!="Time-constant catalytic model")
      newdata3<-subset(newdata, newdata$model!="Time-varying catalytic model")
      
      

###################################
## plots of bias and uncertainty ##
###################################
      
      ## plot comparing estimated versus true FOI and seroprevalence values
      a<- ggplot(newdata2) + facet_wrap("Parameter", scales="free")+
        geom_errorbar(aes(x=True, ymin=Est_low, ymax=Est_upp, col=model), alpha=0.6)+
        geom_point(aes(x=True, y=Est, col=model), size=1.5, alpha=0.4) +
        theme_bw() +
        labs(x="True values", y="Estimated values", col=NULL) +
        geom_abline(slope=1, linetype="dashed", size=1) +
        theme(axis.text = element_text(size=15),
              axis.title = element_text(size=16),
              legend.text=element_text(size=14),
              strip.text = element_text(size=15),
              legend.position = c(0.2,0.85))
      b<-ggplot(newdata3) + facet_wrap("Parameter", scales="free")+
        geom_errorbar(aes(x=True, ymin=Est_low, ymax=Est_upp, col=model), alpha=0.6)+
        geom_point(aes(x=True, y=Est, col=model), size=1.5, alpha=0.4) +
        theme_bw() +
        labs(x=NULL, y="Estimated values", col=NULL) +
        geom_abline(slope=1, linetype="dashed", size=1) +
        theme(axis.text = element_text(size=15),
              axis.title = element_text(size=16),
              legend.text=element_text(size=14),
              strip.text = element_text(size=15),
              legend.position = c(0.2,0.85))
      
      z <- grid.arrange(b,a)
      
      
      ## bias and uncertainty plot
      p1<-ggplot(newdata) + facet_wrap("Parameter", scales="free")+
        geom_freqpoly(aes(x=uncertainty,col=model),size=1) +
        theme_bw() +
        labs(x="\nWidth of confidence interval\n\n", y=NULL, col="Model") +
        theme(axis.text = element_text(size=15),
              axis.title = element_text(size=15),
              legend.text=element_text(size=15),
              legend.title=element_text(size=15),
              strip.text = element_text(size=15))
      
      p2<-ggplot(newdata) + facet_wrap("Parameter", scales="free")+
        geom_freqpoly(aes(x=bias,col=model), size=1) +
        theme_bw() +
        labs(x="\nBias (Estimate-True value)\n", y=NULL, col="Model") +
        geom_vline(linetype="dashed", size=1, xintercept = 0, col="darkgrey") +
        theme(axis.text = element_text(size=15),
              axis.title = element_text(size=15),
              #legend.position = "none",
              legend.text=element_text(size=15),
              legend.title=element_text(size=15),
              strip.text = element_text(size=15))
      
      z2 <- grid.arrange(p1,p2)
      
  
          
###########################
## Calculate CIs on bias ##
###########################

      # FOI
      MIXDATA_f<- MIXDATA[MIXDATA$Parameter=="FOI",]
      CATCDATA_f<- CATCDATA[CATCDATA$Parameter=="FOI",]
      CATVDATA_f<- CATVDATA[CATVDATA$Parameter=="FOI",]
      
      quantile((abs(MIXDATA_f$bias)*100), probs=c(0.5,0.025,0.975))
      quantile((abs(CATCDATA_f$bias)*100), probs=c(0.5,0.025,0.975))
      quantile((abs(CATVDATA_f$bias)*100), probs=c(0.5,0.025,0.975))
      
      # Seroprevalence
      MIXDATA_S<- MIXDATA[MIXDATA$Parameter!="FOI",]
      CATCDATA_S<- CATCDATA[CATCDATA$Parameter!="FOI",]
      CATVDATA_S<- CATVDATA[CATVDATA$Parameter!="FOI",]
      
      quantile((abs(MIXDATA_S$bias)*100), probs=c(0.5,0.025,0.975))
      quantile((abs(CATCDATA_S$bias)*100), probs=c(0.5,0.025,0.975))
      quantile((abs(CATVDATA_S$bias)*100), probs=c(0.5,0.025,0.975))





