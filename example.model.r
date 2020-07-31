
model

{

	for( k in 1 : ncounts )

	{

		log(lambda[k]) <- (beta_space_time[species[k],space_time[k]] * (p_forest[k]-0.5)) + int[region[k],species[k]] + noise[k] #noise is the overdispersion parameter

	 	noise[k] ~ dnorm(0, taunoise)

    #noise[k] ~ dt(0, taunoise, nu) #alternative t-distributed noise = heavy-tailed overdispersion



		count[k] ~ dpois(lambda[k])

	}
	
	#	nu ~ dgamma(2, 0.1) #degrees of freedom (i.e., the heavy-tail component of the t-distribution), if nu is large (infinite) this is equivalent to the normal distribution noise

	#taunoise ~ dscaled.gamma(0.5,10) #weak prior on the overdispersion
taunoise ~ dgamma(0.001,0.001)
	sdnoise <- 1 / pow(taunoise, 0.5)

for(s in 1:nspecies){
### intercepts
  alpha[s] ~ dnorm(0,0.1) #weakly informative prior on the species intercept (mean across regions, hyperparameter)
  #tau_regs_s[s] ~ dscaled.gamma(0.1,10) # slightly informative prior on the regional variation in species-level abundance - defines a prior on the sd = (1/sqrt(tau_regs_s[s])) = 0.1*half-t(df = 10) (half-t is just the positive values of a t-distribution) - this puts 95% of the prior distribution < 0.22, which is still higher than likely
  tau_regs_s[s] ~ dgamma(0.001,0.001) #
  for(g in 1:nregions){
  int[g,s] ~ dnorm(alpha[s],tau_regs_s[s]) #int[g,s] = species-level intercept in region-g (random effect centered on species-level mean)
}
## betas (effects of p_forest)


#modifications for time (==1) and space (==2)
beta_space_time[s,1] ~ dnorm(0,0.01) #treats the species-level time slopes as fixed effects
### difference between time and space
beta_space_time[s,2] <- beta_space_time[s,1]*beta_mod[s]

#exponential model to estimate a multiplicative change in slope
# allows for consistent differences that ignore the direction of the effect for a given species
log(beta_mod[s]) <- e_beta_mod[s]
e_beta_mod[s] ~ dnorm(B_mod,tau_beta_mod)


beta_dif[s] <- beta_space_time[s,1]-beta_space_time[s,2]


  }#s
 


 # hyperparameter for beta_mod treated as random effect
 B_mod ~ dnorm(0,0.1)
 tau_beta_mod ~ dgamma(0.001,0.001)
 #tau_beta_space_time ~ dscaled.gamma(0.1,10)
 
 sd_beta_mod <- 1/pow(tau_beta_mod,0.5)


}


