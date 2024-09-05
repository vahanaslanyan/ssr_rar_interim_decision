q_prior = c(1, 1)
number_to_be_recruited<-sample_size_calculation(r=c(1/2,1/2))
#interim_analysis_at_25%
interim_sample_sizes<-number_to_be_recruited*0.25
ob_boundary<-obrien_fleming_boundary(interim_sample_sizes = interim_sample_sizes,
                                  final_sample_sizes = number_to_be_recruited,alpha = 0.05)
pck_bound<-pocock_boundary(interim_sample_sizes = interim_sample_sizes,
                final_sample_sizes = number_to_be_recruited,alpha = 0.05)


treatment_1<-rnorm(n=interim_sample_sizes$treatment1,0,1)
treatment_2<-rnorm(n=interim_sample_sizes$treatment2,0.6,1)

posteriors<-treatment_effect_statistics(known_sd = 1,treatment_1 = treatment_1,
                                        treatment_2 = treatment_2,
                                        n1=interim_sample_sizes$treatment1,
                                        n2=interim_sample_sizes$treatment2,
                                        delta_star = 0.2)

delta_scaled<-c(0,posteriors["delta_mean"]/posteriors["delta_sd"])

new_r<-allocation_calculation(delta_scaled = delta_scaled,n = sum(interim_sample_sizes),
                       N=sum(number_to_be_recruited))
q_prior<-q_prior+interim_sample_sizes

number_to_be_recruited<-sample_size_calculation(r=new_r,q_prior = q_prior)
total_sample_size<-number_to_be_recruited+interim_sample_sizes
cond_power<-conditional_power(theta=0.2,mean_treatment_effects=c(posteriors["posterior_mean_treatment1"],
                                                     posteriors["posterior_mean_treatment2"]),
                            interim_sds=c(1,1),
                            interim_sample_sizes=interim_sample_sizes,
                            final_sample_sizes=total_sample_size,
                            alpha=0.05)

ppos<-predictive_power(mean_treatment_effects=c(posteriors["posterior_mean_treatment1"],
                                                    posteriors["posterior_mean_treatment2"]),
                 interim_sds=c(1,1),
                 interim_sample_sizes=interim_sample_sizes,
                 final_sample_sizes=total_sample_size,
                 alpha=0.05)

