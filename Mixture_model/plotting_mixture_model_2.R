#####################################
## MIXTURE MODEL PLOTTING - STEP 2 ##
#####################################


#mu_a spline
mu_plot <- ggplot(z_and_a, aes(x=a, y=z)) + geom_point(col="grey") + 
  labs(x="\n Age (Years)",y="log(antibody titre + 1) \n", 
       title = paste("P spline for the mean log(titre+1) score per age \n", sim_num, "\n", sep="")) +
  geom_line(data=mu_a_plotting, aes(x=age, y=mu_a)) + 
  geom_line(data=mu_a_plotting, aes(x=age, y=mu_a_upper), col = "red") + 
  geom_line(data=mu_a_plotting, aes(x=age, y=mu_a_lower), col = "red") +
  theme(axis.title = element_text(size=15), axis.text = element_text(size=12), 
        plot.title = element_text(hjust = 0.5, size=15), panel.border = element_rect(colour = "black", fill = NA), panel.background = element_rect(fill = NA)) +
  scale_x_continuous(breaks=agemin:agemax)
mu_plot
ggsave(filename = paste("Mix_",sim_num, "_spline.png", sep=""), plot = mu_plot)


#seroprevalence
sero_plot <- ggplot(seroprevalence, aes(x = Age, y = pi_a)) + 
  geom_point(col = "black") + geom_line(col = "black") +
  labs(x="Age (Years)", y="Seroprevalence", title = paste("Mixture model estimates of seroprevalence \n", sim_num, sep="")) +
  geom_ribbon(aes(x=Age, ymin = pi_a_low, ymax = pi_a_upp), fill = "blue", alpha = 0.1) +
  theme(plot.title = element_text(hjust=0.5), panel.background = element_rect(fill = NA), panel.border = element_rect(colour = "black", fill = NA)) +
  coord_cartesian(ylim=0:1) +
  scale_x_continuous(breaks=agemin:agemax)
sero_plot
ggsave(filename = paste("Mix_",sim_num, "_Seroprevalence.png", sep=""), plot = sero_plot)


#save table of seroprevalence
write.csv(seroprevalence, file=paste("Mix_",sim_num, "_Seroprevalence_by_age", ".csv", sep=""))


#FOI
FOI_plot <- ggplot(FOI_values, aes(x = Age, y = FOI)) +
  geom_point(col="black") + geom_line(col="black")+
  labs(x="Age", y="FOI", title = paste("Mixture model estimates of FOI \n",sim_num, sep="")) +
  geom_ribbon(aes(x=Age, ymin=FOI_low, ymax = FOI_upp), fill = "blue", alpha = 0.1) +
  theme(plot.title = element_text(hjust=0.5), panel.background = element_rect(fill = NA), panel.border = element_rect(colour = "black", fill = NA)) +
  coord_cartesian(ylim=0:1) +
  scale_x_continuous(breaks=agemin:(agemax-1))

FOI_plot
ggsave(filename = paste("Mix_",sim_num, "_FOI.png", sep=""), plot = FOI_plot)

write.csv(FOI_values, file=paste("Mix_",sim_num, "_FOI_by_age", ".csv", sep=""))


#save table of parameters
stored_params$foi[sim_num] <- FOI_total$FOI
stored_params$foi_low[sim_num] <- FOI_total$lower_ci
stored_params$foi_upp[sim_num] <- FOI_total$upper_ci
stored_params$seroprev[sim_num] <- tot_seroprev
stored_params$seroprev_low[sim_num] <- tot_seroprev_low
stored_params$seroprev_upp[sim_num] <- tot_seroprev_upp




