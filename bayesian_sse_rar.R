set.seed(919)
sample_size_calculation <-
  function(eta = 0.95,
           zeta = 0.90,
           known_sd=1,
           r = c(1 /3, 2 / 3),
           q_prior = c(1, 1),
           delta_star = 0.3) {
    
    # eta point on normal distribution
    z_eta <- qnorm(eta)
    # zeta point on normal distribution
    z_zeta <- qnorm(zeta)
    #precision
    nu<-1/known_sd^2
    V=((z_eta+z_zeta)/delta_star)^2
    sample<-round(sum((4)*V/nu-q_prior))
    n_treatment=sample*r
    results <- round(n_treatment)
    names(results) <- c(paste0("treatment", seq_along(n_treatment)))
    return(as.data.frame(t(results)))
      } 
    
treatment_effect_statistics<-function(known_sd=1,mu_0_1=0,mu_0_2=0,sigma_0_1=1,
                                      sigma_0_2=1,treatment_1,
                                      treatment_2,n1,n2,delta_star){
  #get the posterior values for standard deviation
  sigma_1<-sqrt(1/sigma_0_1^2+n1/known_sd^2)^(-1)
  sigma_2<-sqrt(1/sigma_0_2^2+n2/known_sd^2)^(-1)
  #posterior values for means
  mu_star1<-(mu_0_1/sigma_0_1^2+sum(treatment_1)/known_sd^2)*sigma_1^2
  mu_star2<-(mu_0_2/sigma_0_2^2+sum(treatment_2)/known_sd^2)*sigma_1^2
  #treatment difference
  delta_mean<-mu_star2-mu_star1
  #treatment difference sd
  delta_sd<-sqrt(sigma_1^2+sigma_2^2)
  #standardized difference
  statistics<-(delta_mean-delta_star)/delta_sd
  #function output
  values<-c(mu_star1,mu_star2,sigma_1,sigma_2,delta_mean,delta_sd,statistics)
  names(values)<-c("posterior_mean_treatment1",
                   "posterior_mean_treatment2",
                   "posterior_sd_treatment1",
                   "posterior_sd_treatment2",
                   "delta_mean",
                   "delta_sd",
                   "statistics")
  return(values)
}


allocation_calculation <- function(delta_scaled, n = 1, N = 2) {
  #in normal distribution, probability of max is F(y)=P(y>X1,y>X2,y>X3,y>X4)=
  #=P(y>X1)*P(y>X2)*P(y>X2)*P(y>X3)*P(y>X4),and since they all come from the 
  #same distribution, this turns into P(y>X_i)^number of treatments
  #P(y>x_i) in R is pnorm(x_i)
  prob_2_max <- pnorm(delta_scaled) ^ length(delta_scaled)
  #allocation ratio is scaled based on number of patients recruited at
  #interim analysis, and total number to be recruited
  r <- prob_2_max ^ (n / (2 * N)) / sum(prob_2_max ^ (n / (2 * N)))
  return(r)
}



conditional_power<-function(theta=0.5,
                            mean_treatment_effects=c(0,1),
                            interim_sds=c(1,1),
                            interim_sample_sizes=c(20,20),
                            final_sample_sizes=c(40,40),
                            alpha=0.05){
  #final and interim informations, based on the formula
  final_information<-sum(interim_sds^2/final_sample_sizes)^(-1)
  interim_information<-sum(interim_sds^2/interim_sample_sizes)^(-1)
  #Z statistics
  Z_k<-(mean_treatment_effects[-1]-mean_treatment_effects[1])*sqrt(interim_information)
  #upper and lower ends of the equations
  upper_sided<- (Z_k*sqrt(interim_information)-
    qnorm(1-alpha/2)*sqrt(final_information)+
      theta*(final_information-interim_information))/sqrt(final_information-interim_information)
  lower_sided<-(-Z_k*sqrt(interim_information)-
                  qnorm(1-alpha/2)*sqrt(final_information)-
                  theta*(final_information-interim_information))/sqrt(final_information-interim_information)
  #conditional power
  conditional_power<-pnorm(upper_sided)+pnorm(lower_sided)
  return(conditional_power)
}



predictive_power<-function(mean_treatment_effects=c(0,1),
                           interim_sds=c(1,1),
                           interim_sample_sizes=c(20,20),
                           final_sample_sizes=c(40,40),
                           alpha=0.05){
  #same as in conditional power, based on the formula
  final_information<-sum(interim_sds^2/final_sample_sizes)^(-1)
  interim_information<-sum(interim_sds^2/interim_sample_sizes)^(-1)
  Z_k<-abs((mean_treatment_effects[-1]-mean_treatment_effects[1])*sqrt(interim_information))
  upper_sided<- (Z_k*sqrt(final_information)-
                   qnorm(1-alpha/2)*sqrt(interim_information))/sqrt(final_information-interim_information)
  lower_sided<-(-Z_k*sqrt(final_information)-
                  qnorm(1-alpha/2)*sqrt(interim_information))/sqrt(final_information-interim_information)
  predictive_power<-pnorm(upper_sided)+pnorm(lower_sided)
  return(predictive_power)
}

obrien_fleming_boundary<-function(interim_sample_sizes=c(20,20),
                                  final_sample_sizes=c(40,40),
                                  alpha=0.025){
  #based on OBF spending function
  boundary<-2-2*pnorm((qnorm(1-alpha/2)/sqrt(sum(interim_sample_sizes)/sum(final_sample_sizes))))
  return(abs(qnorm(boundary)))
}
                 
pocock_boundary<-function(interim_sample_sizes=c(20,20),
                          final_sample_sizes=c(40,40),
                          alpha=0.05){
  #based on pocock spending function
  boundary<-alpha*log(1+exp(1)*sum(interim_sample_sizes)/sum(final_sample_sizes))
  return(abs(qnorm(boundary)))
}


run_simulation<-function(efficacy="O'Brien-Fleming",futility="CP",
                         active_treatment_effect=0,futility_boundary=0.1,eta = 0.95,
                         zeta = 0.90, alpha=0.05){
  #initialize vectors for stopping reasons and for the final sample size
  stopped_for_efficacy_1<-c()
  stopped_for_futility_1<-c()
  stopped_for_efficacy_2<-c()
  stopped_for_futility_2<-c()
  final_sample_size<-tibble()
  #perform 5000 simulations, change if needed
  for(i in 1:5000){
    q_prior = c(1, 1)
    r=c(1/2,1/2)
    #calculate initial sample size
    number_to_be_recruited<-sample_size_calculation(eta=eta,zeta=zeta,r=r,q_prior = q_prior)
    #interim_analysis_at_25%
    interim_sample_sizes<-round(number_to_be_recruited*0.25)
    treatment1<-rnorm(n=interim_sample_sizes$treatment1,0,1)
    treatment2<-rnorm(n=interim_sample_sizes$treatment2,active_treatment_effect,1)
    posteriors<-treatment_effect_statistics(known_sd = 1,treatment_1 = treatment1,
                                            treatment_2 = treatment2,
                                            sigma_0_1 = sd(treatment1),
                                            sigma_0_2 = sd(treatment2),
                                            n1=interim_sample_sizes$treatment1,
                                            n2=interim_sample_sizes$treatment2,
                                            delta_star = 0.3)
    delta_scaled<-c(0,posteriors["delta_mean"]/posteriors["delta_sd"])
    #get new allocation
    new_r<-allocation_calculation(delta_scaled = delta_scaled,
                                  n = sum(interim_sample_sizes),
                                  N=sum(number_to_be_recruited))
    q_prior<-q_prior+interim_sample_sizes
    number_to_be_recruited<-sample_size_calculation(eta=eta,zeta=zeta,r=new_r,q_prior = q_prior)
    total_sample_size<-number_to_be_recruited+interim_sample_sizes
    #get efficacy boundaries
    if(efficacy=="O'Brien-Fleming"){
      efficacy_bound<-obrien_fleming_boundary(interim_sample_sizes = interim_sample_sizes,
                                        final_sample_sizes = total_sample_size,alpha = alpha)  
    }else{
      efficacy_bound<-pocock_boundary(interim_sample_sizes = interim_sample_sizes,
                                      final_sample_sizes = total_sample_size,alpha =alpha)
    }
    #check if stopped for efficacy
    eff_criteria<-posteriors["statistics"]>=efficacy_bound
    #get futility boundaries
    if(futility=="CP"){
      fut_power<-conditional_power(theta=0.3,
                                    mean_treatment_effects=c(posteriors["posterior_mean_treatment1"],
                                                             posteriors["posterior_mean_treatment2"]),
                                    interim_sds=c(sd(treatment1),sd(treatment2)),
                                    interim_sample_sizes=interim_sample_sizes,
                                    final_sample_sizes=total_sample_size,
                                    alpha=alpha)  
    }else{
      fut_power<-predictive_power(mean_treatment_effects=c(posteriors["posterior_mean_treatment1"],
                                                           posteriors["posterior_mean_treatment2"]),
                                  interim_sds=c(sd(treatment1),sd(treatment2)),
                                  interim_sample_sizes=interim_sample_sizes,
                                  final_sample_sizes=total_sample_size,
                                  alpha=alpha)
    }
    #check if stopped for futility
    fut_criteria<-fut_power<=futility_boundary
    stopped_for_efficacy_1<-c(stopped_for_efficacy_1,eff_criteria)
    stopped_for_futility_1<-c(stopped_for_futility_1,fut_criteria)
    #if stopped, output sample size and move on
    if(any(eff_criteria,fut_criteria)){
      final_sample_size<-final_sample_size%>%bind_rows(interim_sample_sizes)
      stopped_for_efficacy_2<-c(stopped_for_efficacy_2,NA)
      stopped_for_futility_2<-c(stopped_for_futility_2,NA)
    }else{
      #if not stopped at first interim analysis, repeated step 1
      interim_sample_sizes<-round(total_sample_size*0.5)
      treatment1<-c(treatment1, rnorm(n=interim_sample_sizes$treatment1,0,1))
      treatment2<-c(treatment2, rnorm(n=interim_sample_sizes$treatment2,active_treatment_effect,1))
      posteriors<-treatment_effect_statistics(known_sd = 1,treatment_1 = treatment1,
                                              treatment_2 = treatment2,
                                              mu_0_1=posteriors["posterior_mean_treatment1"],
                                              mu_0_2=posteriors["posterior_mean_treatment2"],
                                              sigma_0_1=posteriors["posterior_sd_treatment1"],
                                              sigma_0_2=posteriors["posterior_sd_treatment2"],
                                              n1=interim_sample_sizes$treatment1,
                                              n2=interim_sample_sizes$treatment2,
                                              delta_star = 0.3)
      
      delta_scaled<-c(0,posteriors["delta_mean"]/posteriors["delta_sd"])
      
      new_r<-allocation_calculation(delta_scaled = delta_scaled,n = sum(interim_sample_sizes),
                                    N=sum(total_sample_size))
      q_prior<-q_prior+interim_sample_sizes
      
      number_to_be_recruited<-sample_size_calculation(eta=eta,zeta=zeta,r=new_r,q_prior = q_prior)
      total_sample_size<-number_to_be_recruited+interim_sample_sizes
      
      if(efficacy=="O'Brien-Fleming"){
        efficacy_bound<-obrien_fleming_boundary(interim_sample_sizes = interim_sample_sizes,
                                                final_sample_sizes = total_sample_size,alpha = alpha)  
      }else{
        efficacy_bound<-pocock_boundary(interim_sample_sizes = interim_sample_sizes,
                                        final_sample_sizes = total_sample_size,alpha = alpha)
      }
      
      eff_criteria<-posteriors["statistics"]>=efficacy_bound
      if(futility=="CP"){
        fut_power<-conditional_power(theta=0.3,
                                     mean_treatment_effects=c(posteriors["posterior_mean_treatment1"],
                                                              posteriors["posterior_mean_treatment2"]),
                                     interim_sds=c(sd(treatment1),sd(treatment2)),
                                     interim_sample_sizes=interim_sample_sizes,
                                     final_sample_sizes=total_sample_size,
                                     alpha=alpha)  
      }else{
        fut_power<-predictive_power(mean_treatment_effects=c(posteriors["posterior_mean_treatment1"],
                                                             posteriors["posterior_mean_treatment2"]),
                                    interim_sds=c(sd(treatment1),sd(treatment2)),
                                    interim_sample_sizes=interim_sample_sizes,
                                    final_sample_sizes=total_sample_size,
                                    alpha=alpha)
      }
      fut_criteria<-fut_power<=futility_boundary
      stopped_for_efficacy_2<-c(stopped_for_efficacy_2,eff_criteria)
      stopped_for_futility_2<-c(stopped_for_futility_2,fut_criteria)
      #if stopped for futility or efficacy at the second interim analysis
      #output sample size at that point, if not, output the final sample size
      if(any(eff_criteria,fut_criteria)){
        final_sample_size<-final_sample_size%>%bind_rows(interim_sample_sizes)
      }else{
        final_sample_size<-final_sample_size%>%bind_rows(total_sample_size)
      }
    }
    final_sample_size$stopped_for_efficacy_1<-stopped_for_efficacy_1
    final_sample_size$stopped_for_futility_1<-stopped_for_futility_1
    final_sample_size$stopped_for_efficacy_2<-stopped_for_efficacy_2
    final_sample_size$stopped_for_futility_2<-stopped_for_futility_2
  }
  #put everything in a tibble for output 
  simulation_summary<-tibble()
  simulation_summary<-simulation_summary%>%bind_rows(colSums(final_sample_size,na.rm=T)/5000)
  simulation_summary$sd_treatment1<-sd(final_sample_size$treatment1,na.rm = T)
  simulation_summary$sd_treatment2<-sd(final_sample_size$treatment2,na.rm = T)
  simulation_summary$efficacy<-efficacy
  simulation_summary$futility<-futility
  simulation_summary$active_treatment<-active_treatment_effect
  return(simulation_summary)
  
}

#code for tables 1 and 2 
df1<-run_simulation(efficacy="O'Brien-Fleming",futility="CP",
                    active_treatment_effect=0)

df1<-df1%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="PP",
                                    active_treatment_effect=0))


df1<-df1%>%bind_rows(run_simulation(efficacy="Pocock",futility="CP",
                                    active_treatment_effect=0))


df1<-df1%>%bind_rows(run_simulation(efficacy="Pocock",futility="PP",
                                    active_treatment_effect=0))

df1<-df1%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="CP",
                    active_treatment_effect=0.1))

df1<-df1%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="PP",
                                    active_treatment_effect=0.1))


df1<-df1%>%bind_rows(run_simulation(efficacy="Pocock",futility="CP",
                                    active_treatment_effect=0.1))


df1<-df1%>%bind_rows(run_simulation(efficacy="Pocock",futility="PP",
                                    active_treatment_effect=0.1))

df1<-df1%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="CP",
                                    active_treatment_effect=0.2))

df1<-df1%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="PP",
                                    active_treatment_effect=0.2))


df1<-df1%>%bind_rows(run_simulation(efficacy="Pocock",futility="CP",
                                    active_treatment_effect=0.2))


df1<-df1%>%bind_rows(run_simulation(efficacy="Pocock",futility="PP",
                                    active_treatment_effect=0.2))



df1<-df1%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="CP",
                                    active_treatment_effect=0.3))

df1<-df1%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="PP",
                                    active_treatment_effect=0.3))


df1<-df1%>%bind_rows(run_simulation(efficacy="Pocock",futility="CP",
                                    active_treatment_effect=0.3))


df1<-df1%>%bind_rows(run_simulation(efficacy="Pocock",futility="PP",
                                    active_treatment_effect=0.3))


df1<-df1%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="CP",
                                    active_treatment_effect=0.4))

df1<-df1%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="PP",
                                    active_treatment_effect=0.4))


df1<-df1%>%bind_rows(run_simulation(efficacy="Pocock",futility="CP",
                                    active_treatment_effect=0.4))


df1<-df1%>%bind_rows(run_simulation(efficacy="Pocock",futility="PP",
                                    active_treatment_effect=0.4))

df1<-df1%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="CP",
                                    active_treatment_effect=0.5))

df1<-df1%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="PP",
                                    active_treatment_effect=0.5))


df1<-df1%>%bind_rows(run_simulation(efficacy="Pocock",futility="CP",
                                    active_treatment_effect=0.5))


df1<-df1%>%bind_rows(run_simulation(efficacy="Pocock",futility="PP",
                                    active_treatment_effect=0.5))


df1<-df1%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="CP",
                                    active_treatment_effect=0.6))

df1<-df1%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="PP",
                                    active_treatment_effect=0.6))


df1<-df1%>%bind_rows(run_simulation(efficacy="Pocock",futility="CP",
                                    active_treatment_effect=0.6))


df1<-df1%>%bind_rows(run_simulation(efficacy="Pocock",futility="PP",
                                    active_treatment_effect=0.6))


#write_csv(df1,"Table1_2.csv")


#code for tables 3 and 4.


df2<-run_simulation(efficacy="O'Brien-Fleming",futility="CP",
                    active_treatment_effect=0,futility_boundary = 0.2)

df2<-df2%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="PP",
                                    active_treatment_effect=0,futility_boundary = 0.2))


df2<-df2%>%bind_rows(run_simulation(efficacy="Pocock",futility="CP",
                                    active_treatment_effect=0,futility_boundary = 0.2))


df2<-df2%>%bind_rows(run_simulation(efficacy="Pocock",futility="PP",
                                    active_treatment_effect=0,futility_boundary = 0.2))

df2<-df2%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="CP",
                                    active_treatment_effect=0.1,futility_boundary = 0.2))

df2<-df2%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="PP",
                                    active_treatment_effect=0.1,futility_boundary = 0.2))


df2<-df2%>%bind_rows(run_simulation(efficacy="Pocock",futility="CP",
                                    active_treatment_effect=0.1,futility_boundary = 0.2))


df2<-df2%>%bind_rows(run_simulation(efficacy="Pocock",futility="PP",
                                    active_treatment_effect=0.1,futility_boundary = 0.2))

df2<-df2%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="CP",
                                    active_treatment_effect=0.2,futility_boundary = 0.2))

df2<-df2%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="PP",
                                    active_treatment_effect=0.2,futility_boundary = 0.2))


df2<-df2%>%bind_rows(run_simulation(efficacy="Pocock",futility="CP",
                                    active_treatment_effect=0.2,futility_boundary = 0.2))


df2<-df2%>%bind_rows(run_simulation(efficacy="Pocock",futility="PP",
                                    active_treatment_effect=0.2,futility_boundary = 0.2))



df2<-df2%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="CP",
                                    active_treatment_effect=0.3,futility_boundary = 0.2))

df2<-df2%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="PP",
                                    active_treatment_effect=0.3,futility_boundary = 0.2))


df2<-df2%>%bind_rows(run_simulation(efficacy="Pocock",futility="CP",
                                    active_treatment_effect=0.3,futility_boundary = 0.2))


df2<-df2%>%bind_rows(run_simulation(efficacy="Pocock",futility="PP",
                                    active_treatment_effect=0.3,futility_boundary = 0.2))


df2<-df2%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="CP",
                                    active_treatment_effect=0.4,futility_boundary = 0.2))

df2<-df2%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="PP",
                                    active_treatment_effect=0.4,futility_boundary = 0.2))


df2<-df2%>%bind_rows(run_simulation(efficacy="Pocock",futility="CP",
                                    active_treatment_effect=0.4,futility_boundary = 0.2))


df2<-df2%>%bind_rows(run_simulation(efficacy="Pocock",futility="PP",
                                    active_treatment_effect=0.4,futility_boundary = 0.2))

df2<-df2%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="CP",
                                    active_treatment_effect=0.5,futility_boundary = 0.2))

df2<-df2%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="PP",
                                    active_treatment_effect=0.5,futility_boundary = 0.2))


df2<-df2%>%bind_rows(run_simulation(efficacy="Pocock",futility="CP",
                                    active_treatment_effect=0.5,futility_boundary = 0.2))


df2<-df2%>%bind_rows(run_simulation(efficacy="Pocock",futility="PP",
                                    active_treatment_effect=0.5,futility_boundary = 0.2))


df2<-df2%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="CP",
                                    active_treatment_effect=0.6,futility_boundary = 0.2))

df2<-df2%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="PP",
                                    active_treatment_effect=0.6,futility_boundary = 0.2))


df2<-df2%>%bind_rows(run_simulation(efficacy="Pocock",futility="CP",
                                    active_treatment_effect=0.6,futility_boundary = 0.2))


df2<-df2%>%bind_rows(run_simulation(efficacy="Pocock",futility="PP",
                                    active_treatment_effect=0.6,futility_boundary = 0.2))





#code for tables 7 and 8


df3<-run_simulation(efficacy="O'Brien-Fleming",futility="CP",
                    active_treatment_effect=0,eta = 0.9,
                    zeta = 0.80,alpha=0.1)

df3<-df3%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="PP",
                                    active_treatment_effect=0,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))


df3<-df3%>%bind_rows(run_simulation(efficacy="Pocock",futility="CP",
                                    active_treatment_effect=0,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))


df3<-df3%>%bind_rows(run_simulation(efficacy="Pocock",futility="PP",
                                    active_treatment_effect=0,eta = 0.9,
                     zeta = 0.80,alpha=0.1))

df3<-df3%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="CP",
                                    active_treatment_effect=0.1,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))

df3<-df3%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="PP",
                                    active_treatment_effect=0.1,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))


df3<-df3%>%bind_rows(run_simulation(efficacy="Pocock",futility="CP",
                                    active_treatment_effect=0.1,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))


df3<-df3%>%bind_rows(run_simulation(efficacy="Pocock",futility="PP",
                                    active_treatment_effect=0.1,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))

df3<-df3%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="CP",
                                    active_treatment_effect=0.2,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))

df3<-df3%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="PP",
                                    active_treatment_effect=0.2,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))


df3<-df3%>%bind_rows(run_simulation(efficacy="Pocock",futility="CP",
                                    active_treatment_effect=0.2,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))


df3<-df3%>%bind_rows(run_simulation(efficacy="Pocock",futility="PP",
                                    active_treatment_effect=0.2,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))



df3<-df3%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="CP",
                                    active_treatment_effect=0.3,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))

df3<-df3%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="PP",
                                    active_treatment_effect=0.3,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))


df3<-df3%>%bind_rows(run_simulation(efficacy="Pocock",futility="CP",
                                    active_treatment_effect=0.3,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))


df3<-df3%>%bind_rows(run_simulation(efficacy="Pocock",futility="PP",
                                    active_treatment_effect=0.3,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))


df3<-df3%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="CP",
                                    active_treatment_effect=0.4,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))

df3<-df3%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="PP",
                                    active_treatment_effect=0.4,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))


df3<-df3%>%bind_rows(run_simulation(efficacy="Pocock",futility="CP",
                                    active_treatment_effect=0.4,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))


df3<-df3%>%bind_rows(run_simulation(efficacy="Pocock",futility="PP",
                                    active_treatment_effect=0.4,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))

df3<-df3%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="CP",
                                    active_treatment_effect=0.5,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))

df3<-df3%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="PP",
                                    active_treatment_effect=0.5,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))


df3<-df3%>%bind_rows(run_simulation(efficacy="Pocock",futility="CP",
                                    active_treatment_effect=0.5,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))


df3<-df3%>%bind_rows(run_simulation(efficacy="Pocock",futility="PP",
                                    active_treatment_effect=0.5,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))


df3<-df3%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="CP",
                                    active_treatment_effect=0.6,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))

df3<-df3%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="PP",
                                    active_treatment_effect=0.6,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))


df3<-df3%>%bind_rows(run_simulation(efficacy="Pocock",futility="CP",
                                    active_treatment_effect=0.6,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))


df3<-df3%>%bind_rows(run_simulation(efficacy="Pocock",futility="PP",
                                    active_treatment_effect=0.6,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))


#code for tables 5-6
df4<-run_simulation(efficacy="O'Brien-Fleming",futility="CP",
                    active_treatment_effect=0,futility_boundary = 0.2,eta = 0.9,
                    zeta = 0.80,alpha=0.1)

df4<-df4%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="PP",
                                    active_treatment_effect=0,futility_boundary = 0.2,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))


df4<-df4%>%bind_rows(run_simulation(efficacy="Pocock",futility="CP",
                                    active_treatment_effect=0,futility_boundary = 0.2,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))


df4<-df4%>%bind_rows(run_simulation(efficacy="Pocock",futility="PP",
                                    active_treatment_effect=0,futility_boundary = 0.2,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))

df4<-df4%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="CP",
                                    active_treatment_effect=0.1,futility_boundary = 0.2,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))

df4<-df4%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="PP",
                                    active_treatment_effect=0.1,futility_boundary = 0.2,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))


df4<-df4%>%bind_rows(run_simulation(efficacy="Pocock",futility="CP",
                                    active_treatment_effect=0.1,futility_boundary = 0.2,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))


df4<-df4%>%bind_rows(run_simulation(efficacy="Pocock",futility="PP",
                                    active_treatment_effect=0.1,futility_boundary = 0.2,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))

df4<-df4%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="CP",
                                    active_treatment_effect=0.2,futility_boundary = 0.2,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))

df4<-df4%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="PP",
                                    active_treatment_effect=0.2,futility_boundary = 0.2,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))


df4<-df4%>%bind_rows(run_simulation(efficacy="Pocock",futility="CP",
                                    active_treatment_effect=0.2,futility_boundary = 0.2,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))


df4<-df4%>%bind_rows(run_simulation(efficacy="Pocock",futility="PP",
                                    active_treatment_effect=0.2,futility_boundary = 0.2,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))



df4<-df4%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="CP",
                                    active_treatment_effect=0.3,futility_boundary = 0.2,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))

df4<-df4%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="PP",
                                    active_treatment_effect=0.3,futility_boundary = 0.2,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))


df4<-df4%>%bind_rows(run_simulation(efficacy="Pocock",futility="CP",
                                    active_treatment_effect=0.3,futility_boundary = 0.2,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))


df4<-df4%>%bind_rows(run_simulation(efficacy="Pocock",futility="PP",
                                    active_treatment_effect=0.3,futility_boundary = 0.2,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))


df4<-df4%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="CP",
                                    active_treatment_effect=0.4,futility_boundary = 0.2,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))

df4<-df4%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="PP",
                                    active_treatment_effect=0.4,futility_boundary = 0.2,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))


df4<-df4%>%bind_rows(run_simulation(efficacy="Pocock",futility="CP",
                                    active_treatment_effect=0.4,futility_boundary = 0.2,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))


df4<-df4%>%bind_rows(run_simulation(efficacy="Pocock",futility="PP",
                                    active_treatment_effect=0.4,futility_boundary = 0.2,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))

df4<-df4%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="CP",
                                    active_treatment_effect=0.5,futility_boundary = 0.2,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))

df4<-df4%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="PP",
                                    active_treatment_effect=0.5,futility_boundary = 0.2,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))


df4<-df4%>%bind_rows(run_simulation(efficacy="Pocock",futility="CP",
                                    active_treatment_effect=0.5,futility_boundary = 0.2,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))


df4<-df4%>%bind_rows(run_simulation(efficacy="Pocock",futility="PP",
                                    active_treatment_effect=0.5,futility_boundary = 0.2,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))


df4<-df4%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="CP",
                                    active_treatment_effect=0.6,futility_boundary = 0.2,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))

df4<-df4%>%bind_rows(run_simulation(efficacy="O'Brien-Fleming",futility="PP",
                                    active_treatment_effect=0.6,futility_boundary = 0.2,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))


df4<-df4%>%bind_rows(run_simulation(efficacy="Pocock",futility="CP",
                                    active_treatment_effect=0.6,futility_boundary = 0.2,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))


df4<-df4%>%bind_rows(run_simulation(efficacy="Pocock",futility="PP",
                                    active_treatment_effect=0.6,futility_boundary = 0.2,eta = 0.9,
                                    zeta = 0.80,alpha=0.1))

#write_csv(df1,"confirmatory_futility_01.csv")

#write_csv(df2,"confirmatory_futility_02.csv")

#write_csv(df3,"exploratory_futility_01.csv")

#write_csv(df4,"exploratory_futility_02.csv")
