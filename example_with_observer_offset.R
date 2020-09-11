
# JAGS model description --------------------------------------------------
# this generates a character vector of the entire model description
# which is then written to a text file using the cat argument below.
# here I've just included the model description
# the simulated data in "example2.R" could work if it had observer-structures added in

modl <- "
model

{

	for( k in 1 : ncounts )

	{
# the p_forest[k]*0.5 centers the percent forest values so that the intercepts represent the mean abundance (when p_forest == 0.5, 50% forest), instead of the abundance when p_forest == 0
# this is relatively simple to change, depending on the range of your p_forest values
# obs_offset - entered as data, vector of length number of observers - log-scale difference in mean count of all species for this observer
# obs[k] - entered as data, vector of length ncounts - numerical indicators for each observer
log(lambda[k]) <- (beta_time_space[species[k],space_time[k]] * (p_forest[k]-0.5)) + int[region[k],species[k]] + obs_offset[obs[k]] + noise[k] #noise is the overdispersion parameter

	 	noise[k] ~ dnorm(0, taunoise)

    #noise[k] ~ dt(0, taunoise, nu) #alternative t-distributed noise = heavy-tailed overdispersion
    # for now, I suggest you stick with the normal overdispersion, the t-distribution requires more time to converge



		count[k] ~ dpois(lambda[k])

	}
	
	#	nu ~ dgamma(2, 0.1) #degrees of freedom (i.e., the heavy-tail component of the t-distribution), if nu is large (infinite) this is equivalent to the normal distribution noise


# taunoise ~ dgamma(0.001,0.001) #prior on the precision (inverse of the variance) of the overdispersion
# 	sdnoise <- 1 / pow(taunoise, 0.5) #transformation of the precision into a standard deviation scale that is easier to interpret

 sdnoise ~ dt(0, 1, 4) T(0,)     #half-t prior on the standard deviation
 taunoise <- 1/pow(sdnoise,2)
 


for(s in 1:nspecies){
### species-level intercepts
  alpha[s] ~ dnorm(0,1) #weakly informative prior on the species intercept (mean across regions, hyperparameter)
  #tau_regs_s[s] ~ dgamma(0.001,0.001) #prior on the precision for the regional variation in species mean abundance
  ### alternative half-t prior

 sd_regs_s[s] ~ dt(0, 1, 4) T(0,)     #half-t prior on the standard deviation
 
 tau_regs_s[s] <- 1/pow(sd_regs_s[s],2)
  
  for(g in 1:nregions){
  int[g,s] ~ dnorm(alpha[s],tau_regs_s[s]) #int[g,s] = species-level intercept in region-g (random effect centered on species-level mean)
}


########### slopes (betas = effects of p_forest)
### beta_time_space is a matrix with nspecies rows and 2 columns. 
### the first column is the time-slope, the second colum is the space-slope
#modifications for time (==1) and space (==2)

beta_time_space[s,1] ~ dnorm(0,0.01) #treats the species-level time slopes as fixed effects (so makes no assumption about the direction or magnitude of the species response to forest)

## alternatively, one could estimate the slopes as random effects, which would reduce the influence of data-sparse species
# beta_time_space[s,1] ~ dnorm(0,tau_beta_time_space1) #treats the species-level time slopes as fixed effects (so makes no assumption about the direction or magnitude of the species response to forest)
#beta_time_space_alt[s] ~ dnorm(0,0.01) #alternate calcluation of a slope for the space effect, that only gets informed when psi[s] is low, and I[s] == 0 - see below

################################### Alternative 1
### multiplicative modification of the time-slope to generate the space-slope
### using a strictly positive, multiplicative modification (instead of a smple additive effect) so that
### the sharing of information among species only assumes some change in the magnitude of the slope, instead of a change in direction
### if beta_mod[s] == 1, then there is no difference in the space and time slopes
### if beta_mod[s] > 1 then there is a stronger effect of space than time
### if beta_mod[s] < 1 then there is a weaker effect of space than time
#beta_time_space[s,2] <- beta_time_space[s,1]*beta_mod[s]

#exponential model to estimate a strictly positive, multiplicative change in slope
## allows for consistent differences in the strength of the slope, while ignoring the direction of the effect for a given species
# log(beta_mod[s]) <- e_beta_mod[s]
# e_beta_mod[s] ~ dnorm(B_mod,tau_beta_mod)
## this multiplicative approach may fall apart if there is no time-slope
## that is, if beta_time_space[s,1] == 0, then beta_mod[s] would be imposible to estimate
## hopefully, the sharing of infor among species, and the prior would then shrink that species estimates for both slopes to 1
## it could work well in the opposite direction too, if space was the base-slope and the modification was to estimate time-slopes
#######################################################


################################### Alternative 1.1 - commented out due to email discussion
### multiplicative modification of the time-slope to generate the space-slope, plus an indicator variable that removes the relationship between the time and space slope when the time slope is lower than a particular threshold
### using a strictly positive, multiplicative modification (instead of a simple additive effect) so that
### the sharing of information among species only assumes some change in the magnitude of the slope, instead of a change in direction
### if beta_mod[s] == 1, then there is no difference in the space and time slopes
### if beta_mod[s] > 1 then there is a stronger effect of space than time
### if beta_mod[s] < 1 then there is a weaker effect of space than time
#beta_time_space[s,2] <- (beta_time_space[s,1]*I[s])*beta_mod[s] + (beta_time_space_alt[s]*((I[s]-1)*-1))

#exponential model to estimate a strictly positive, multiplicative change in slope
## allows for consistent differences in the strength of the slope, while ignoring the direction of the effect for a given species
#log(beta_mod[s]) <- e_beta_mod[s]
#e_beta_mod[s] ~ dnorm(B_mod,tau_beta_mod)
## this multiplicative approach may fall apart if there is no time-slope
## that is, if beta_time_space[s,1] == 0, then beta_mod[s] would be imposible to estimate
## hopefully, the sharing of infor among species, and the prior would then shrink that species estimates for both slopes to 1
## it could work well in the opposite direction too, if space was the base-slope and the modification was to estimate time-slopes
## could also estimate the probability that the beta_time_space[s,1]
# I[s] ~ dbern(psi[s])
#    alpha_psi[s] ~ dunif(2,3)
#    beta_psi[s] ~ dunif(1,2)
#    psi[s] ~ dbeta(alpha_psi[s],beta_psi[s])
    
#######################################################

####################################################### Alternative 2 - commented out
## alternative additive model (currently commented out) that shares information among species but does not require the time-slope to be non-zero
## however, this approach then assumes that there is a consistent difference in direction of the slopes
## that is, time-slopes are always more negative (if B_mod < 0)
## this seems like an arbitrary assumption, with no clear basis in reality, unless we're talking about only forest birds?
## this is what I was imagining while we were talking, but then as I coded it out, I started to doubt this assumption
#beta_time_space[s,2] <- beta_time_space[s,1]+beta_mod[s]
#beta_mod[s] ~ dnorm(B_mod,tau_beta_mod)
#######################################################


####################################################### Alternative 3
## alternative random-effect model that does not share information among species on the direction or magnitude of the difference
## but does share information on the varaition among species in the difference
## but has a regularizing prior to reduce differences for data-sparse species to 0
## this is tidy, but it assumes that the difference between the two kinds of slopes == 0
## if there are solid data to the contrary, this will show it
## this will also reduce the chances of seeing wild variation among species that are largely the result of sampling error
## i.e., it reduces the influence of poor data-species, but allows them to stay in the model.

beta_time_space[s,2] <- beta_time_space[s,1]+beta_mod[s]
beta_mod[s] ~ dnorm(0,tau_beta_mod)

## this same mean-zero random effect approach could also be used for the time-slope estimates to provide the same benefits for data-sparse species there
## e.g., see note above 
## 
#######################################################


beta_dif[s] <- beta_time_space[s,1]-beta_time_space[s,2]


  }#s
 


 # hyperparameter for beta_mod treated as random effect
 #B_mod ~ dnorm(0,7) #informative prior placing 95% of the prior between 0.5 and 2.0 for the multiplicative difference averaged across species - possibly too strong a prior
 #tau_beta_mod ~ dgamma(0.001,0.001)
 #tau_beta_time_space1 ~ dgamma(0.001,0.001) # alternative random effects precision for the time-slopes
 #sd_beta_mod <- 1/pow(tau_beta_mod,0.5)

### alternative half-t prior

 sd_beta_modt ~ dt(0, 1, 4) T(0,)     # prior that is effectively uninformative
 tau_beta_mod <- 1/pow(sd_beta_mod,2)
}


"
cat(modl,file = "example_model_obs_offset.r")


