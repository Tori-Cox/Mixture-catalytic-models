##############################
## CATALYTIC MODEL PLOTTING ##
##############################


## summary results
    FOI_constant<-fit.c_CI[[1]]
    FOI_constant_LOW<-fit.c_CI[[2]]
    FOI_constant_UPP<-fit.c_CI[[3]]
    
    FOI_variable <-mean(fit.v_CI[[1]])
    FOI_variable_LOW <- mean(fit.v_CI[[2]])
    FOI_variable_UPP <- mean(fit.v_CI[[3]])
    
    SP_constant<-SP_total_c[1]
    SP_constant_LOW<-SP_total_c[2]
    SP_constant_UPP<-SP_total_c[3]
    
    SP_variable <-SP_total_v[1]
    SP_variable_LOW <- SP_total_v[2]
    SP_variable_UPP <- SP_total_v[3]
    
    
    Results <- c(FOI_constant,FOI_constant_LOW, FOI_constant_UPP, 
                 FOI_variable, FOI_variable_LOW, FOI_variable_UPP,
                 SP_constant,SP_constant_LOW, SP_constant_UPP, 
                 SP_variable, SP_variable_LOW, SP_variable_UPP)
    Results <- data.frame(Variable = c("FOI_c", "FOI_c_LOW", "FOI_c_UPP", 
                                       "FOI_v","FOI_v_LOW", "FOI_v_UPP",
                                       "SP_c", "SP_c_LOW", "SP_c_UPP", 
                                       "SP_v","SP_v_LOW", "SP_v_UPP"), Value = Results)
    write.csv(Results, file=paste("model_output/catalytic/Cat_", sim, "_SummaryResults.csv", sep=""))
    
    
    ## enter into storage list for all 200 datasets
    store_foiC[[sim]] <- c(fit.c_CI[[1]], fit.c_CI[[2]],fit.c_CI[[3]])
    store_foiV[[sim]] <- c(mean(fit.v_CI[[1]]), mean(fit.v_CI[[2]]), mean(fit.v_CI[[3]]))
    store_spC[[sim]] <- c(SP_total_c[1], SP_total_c[2], SP_total_c[3])
    store_spV[[sim]] <- c(SP_total_v[1], SP_total_v[2], SP_total_v[3])
    


## save FOI by age for time-varying model 
    FOI_by_age <- data.frame(Age = loc.all$Age, FOI = fit.v_CI[[1]], FOI_low = fit.v_CI[[2]], FOI_upp = fit.v_CI[[3]])
    write.csv(FOI_by_age, paste("model_output/catalytic/Cat_", sim, "_FOI_by_age.csv", sep=""))

    
   
## save seroprevalence by age csv file and plots
    shift <- min(loc.all$Age)
    age<-loc.all$Age
    nyk<-loc.all[,which(dimnames(loc.all)[[2]]=="X1")]
    N<-loc.all[,which(dimnames(loc.all)[[2]]=="X0")]+nyk
    prop.seropos<-nyk/N
    
    ## constant FOI
    y<- SP_c_FOI$Central
    y_low<-SP_c_FOI$Lower
    y_upp<-SP_c_FOI$Upper
    
    Seropositive_by_age <- data.frame(Age = age, Serop_count = prop.seropos, Serop_model = y, Serop_model_low = y_low, Serop_model_upp = y_upp)
    write.csv(Seropositive_by_age, paste("model_output/catalytic/Cat_", sim, "_Seropositive_by_age_CONSTANT_FOI.csv", sep=""))
    
    seroplot <- ggplot(Seropositive_by_age, aes(x=Age, y=Serop_model)) +
        geom_line(col = "black") +
        geom_point(aes(x=Age, y=Serop_count), col="black")+
        labs(x="Age (Years)", y="Seroprevalence", 
             title = paste("Catalytic model estimates of seroprevalence", "\n(constant force of infection)", sep="")) +
        geom_ribbon(aes(x=Age, ymin = Serop_model_low, ymax = Serop_model_upp), fill = "blue", alpha = 0.1) +
        theme(plot.title = element_text(hjust=0.5), 
              panel.background = element_rect(fill = NA), 
              panel.border = element_rect(colour = "black", fill = NA)) +
       coord_cartesian(ylim=0:1) +
        scale_x_continuous(breaks=min(age):max(age))
    
    seroplot
    ggsave(filename = paste("model_output/catalytic/Cat_", sim, "_Seroprevalence_constant_FOI.png", sep=""), plot = seroplot)
    
     
    ## variable FOI
    y<- SP_v_FOI$Central
    y_low<-SP_v_FOI$Lower
    y_upp<-SP_v_FOI$Upper
    
    Seropositive_by_age <- data.frame(Age = age, Serop_count = prop.seropos, Serop_model = y, Serop_model_low = y_low, Serop_model_upp = y_upp)
    write.csv(Seropositive_by_age, paste("model_output/catalytic/Cat_", sim, "_Seropositive_by_age_VARIABLE_FOI.csv", sep=""))
    
    seroplot2 <- ggplot(Seropositive_by_age, aes(x=Age, y=Serop_model)) +
      geom_line(col = "black") +
      geom_point(aes(x=Age, y=Serop_count), col="black")+
      labs(x="Age (Years)", y="Seroprevalence", 
           title = paste("Catalytic model estimates of seroprevalence", "\n(variable force of infection)", sep="")) +
      geom_ribbon(aes(x=Age, ymin = Serop_model_low, ymax = Serop_model_upp), fill = "blue", alpha = 0.1) +
      theme(plot.title = element_text(hjust=0.5), 
            panel.background = element_rect(fill = NA), 
            panel.border = element_rect(colour = "black", fill = NA)) +
      coord_cartesian(ylim=0:1) +
      scale_x_continuous(breaks=min(age):max(age))
    seroplot2
    ggsave(filename = paste("model_output/catalytic/Cat_", sim, "_Seroprevalence_variable_FOI.png", sep=""), plot = seroplot2)
    
    
    
    