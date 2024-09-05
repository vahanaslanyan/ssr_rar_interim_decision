sample_size_calculation <-
  function(eta = 0.95,
           zeta = 0.90,
           known_sd=1,
           r = c(1 /3, 2 / 3),
           q_prior = c(1, 1),
           delta_star = 0.2) {
    
    # eta point on normal distribution
    z_eta <- qnorm(eta)
    # zeta point on normal distribution
    z_zeta <- qnorm(zeta)
    #precision
    nu<-1/known_sd^2
    V=((z_eta+z_zeta)/delta_star)^2
    sample<-round(sum((1+length(r)^(-0.5))*V/nu-q_prior))
    n_treatment=sample*r
    results <- round(n_treatment)
    names(results) <- c(paste0("treatment", seq_along(n_treatment)))
    return(as.data.frame(t(results)))
      } 
    
treatment_effect_statistics<-function(known_sd=1,mu_0=0,sigma_0=1,treatment_1,
                                      treatment_2,n1,n2,delta_star){
  sigma_1<-sqrt(1/sigma_0^2+n1/known_sd^2)^(-1)
  sigma_2<-sqrt(1/sigma_0^2+n2/known_sd^2)^(-1)
  mu_star1<-(mu_0/sigma_0^2+sum(treatment_1)/known_sd^2)*sigma_1^2
  mu_star2<-(mu_0/sigma_0^2+sum(treatment_2)/known_sd^2)*sigma_1^2
  delta_mean<-mu_star2-mu_star1
  delta_sd<-sqrt(sigma_1^2+sigma_2^2)
  statistics<-(delta_mean-delta_star)/delta_sd
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
  
  final_information<-sum(interim_sds^2/final_sample_sizes)^(-1)
  interim_information<-sum(interim_sds^2/interim_sample_sizes)^(-1)
  Z_k<-(mean_treatment_effects[-1]-mean_treatment_effects[1])*sqrt(interim_information)
  upper_sided<- (Z_k*sqrt(interim_information)-
    qnorm(1-alpha/2)*sqrt(final_information)+
      theta*(final_information-interim_information))/sqrt(final_information-interim_information)
  lower_sided<-(-Z_k*sqrt(interim_information)-
                  qnorm(1-alpha/2)*sqrt(final_information)-
                  theta*(final_information-interim_information))/sqrt(final_information-interim_information)
  conditional_power<-pnorm(upper_sided)+pnorm(lower_sided)
  return(conditional_power)
}



predictive_power<-function(mean_treatment_effects=c(0,1),
                           interim_sds=c(1,1),
                           interim_sample_sizes=c(20,20),
                           final_sample_sizes=c(40,40),
                           alpha=0.05){
  final_information<-sum(interim_sds^2/final_sample_sizes)^(-1)
  interim_information<-sum(interim_sds^2/interim_sample_sizes)^(-1)
  Z_k<-abs((mean_treatment_effects[-1]-mean_treatment_effects[1])*sqrt(interim_information))
  upper_sided<- (Z_k*sqrt(interim_information)-
                   qnorm(1-alpha/2)*sqrt(final_information))/sqrt(final_information-interim_information)
  lower_sided<-(-Z_k*sqrt(interim_information)-
                  qnorm(1-alpha/2)*sqrt(final_information))/sqrt(final_information-interim_information)
  predictive_power<-pnorm(upper_sided)+pnorm(lower_sided)
  return(predictive_power)
}

obrien_fleming_boundary<-function(interim_sample_sizes=c(20,20),
                                  final_sample_sizes=c(40,40),
                                  alpha=0.025){
  boundary<-2-2*pnorm((qnorm(1-alpha/2)/sqrt(sum(interim_sample_sizes)/sum(final_sample_sizes))))
  return(abs(qnorm(boundary)))
}
                 
pocock_boundary<-function(interim_sample_sizes=c(20,20),
                          final_sample_sizes=c(40,40),
                          alpha=0.05){
  boundary<-alpha*log(1+exp(1)*sum(interim_sample_sizes)/sum(final_sample_sizes))
  return(abs(qnorm(boundary)))
}




        