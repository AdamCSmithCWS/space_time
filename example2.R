## ideas for space-time comparison model
library(jagsUI)
library(tidyverse)
library(ggmcmc)

#number of space_time replicates (i.e., number of regions)
nreps = 21


# setting up fake data structure ------------------------------------------

dat <- expand.grid(species = c("Red-eyed Vireo","American Robin","Ovenbird","Wood Thrush","Hermit Thrush","Black-capped Chickadee"),
                   region = c("BCR12","BCR13","BCR23","BCR34"),
                   space_time = rep(c("time","space"),nreps),
                   p_forest = seq(0.1,0.8,length.out = 13))


#changing the factor-level information in dat to unique integers
#this is required for JAGS, it can't handle factors or character variables, only numeric and integer
# so for example, this changes "Red-eyed Vireo" to 1, "American Robin" to 2, etc.
for(i in c("species","region","space_time")){
dat[,paste0(i,"_f")] = as.integer(dat[,i])
}

nspecies = length(unique(dat$species)) #number of species
nregions = length(unique(dat$region)) #number of regions


## True trends for time 
B = c(1,0.5,2,-2,-1.5,1) #true species slopes for p_forest (change in log-scale counts for a unit change in p_forest) note: a unit change in proportional variable = change from no forest to complete forest, so these values are relatively large

## True proportional change in slope for space i.e., here space-slopes are 90% of the value for time - on average across all species
True.mean.modifier.space <- 0.9
True.mod.B <- rnorm(nspecies,True.mean.modifier.space,0.05) #random variation among species in the modifier

#table of the True slopes for each species in time (column-1) and space (column-2)
B.space.time = matrix(c(B,True.mod.B*B),nrow = nspecies,ncol = 2,byrow = FALSE) #random 1%/year variance around the specie-leve mean response based on time or space


## species abundances = True intercepts
a.sp = runif(nspecies,0.5,2) #coarse mean counts overall

a.reg = rnorm(nregions,0,0.1)#random variation in abundance among regions

#relatively small amoung of overdispersion
noise <- 0.001 #itnroduced overdispersion

#NOTE: this all ignore observers, but that shouldn't be too hard to add in afterwards

for(i in 1:nrow(dat)){
  forst.eff <- (B.space.time[dat[i,"species_f"],dat[i,"space_time_f"]]*(dat[i,"p_forest"]-0.5)) #this -0.5 just centers the values of percent forest without changing the scale of the variable, it also means the intercepts represent the mean values of abundance and not some hypoethical mean count when forest == 0
  int <- a.sp[dat[i,"species_f"]]+a.reg[dat[i,"region_f"]]
  lmbd <- exp(int+forst.eff + rnorm(1,0,noise))
  dat[i,"count"] <- rpois(1,lambda = lmbd)
}

ncounts = nrow(dat)



# JAGS model description --------------------------------------------------
# this generates a character vector of the entire model description
# which is then written to a text file using the cat argument below.

modl <- "
model

{

	for( k in 1 : ncounts )

	{
# the p_forest[k]*0.5 centers the percent forest values so that the intercepts represent the mean abundance (when p_forest == 0.5, 50% forest), instead of the abundance when p_forest == 0
# this is relatively simple to change, depending on the range of your p_forest values

log(lambda[k]) <- (beta_time_space[species[k],space_time[k]] * (p_forest[k]-0.5)) + int[region[k],species[k]] + noise[k] #noise is the overdispersion parameter

	 	noise[k] ~ dnorm(0, taunoise)

    #noise[k] ~ dt(0, taunoise, nu) #alternative t-distributed noise = heavy-tailed overdispersion
    # for now, I suggest you stick with the normal overdispersion, the t-distribution requires more time to converge



		count[k] ~ dpois(lambda[k])

	}
	
	#	nu ~ dgamma(2, 0.1) #degrees of freedom (i.e., the heavy-tail component of the t-distribution), if nu is large (infinite) this is equivalent to the normal distribution noise

taunoise ~ dgamma(0.001,0.001) #prior on the precision (inverse of the variance) of the overdispersion
	sdnoise <- 1 / pow(taunoise, 0.5) #transformation of the precision into a standard deviation scale that is easier to interpret

for(s in 1:nspecies){
### species-level intercepts
  alpha[s] ~ dnorm(0,0.1) #weakly informative prior on the species intercept (mean across regions, hyperparameter)
  tau_regs_s[s] ~ dgamma(0.001,0.001) #prior on the precision for the regional variation in species mean abundance
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


################################### Alternative 1
### multiplicative modification of the time-slope to generate the space-slope
### using a strictly positive, multiplicative modification (instead of a smple additive effect) so that
### the sharing of information among species only assumes some change in the magnitude of the slope, instead of a change in direction
### if beta_mod[s] == 1, then there is no difference in the space and time slopes
### if beta_mod[s] > 1 then there is a stronger effect of space than time
### if beta_mod[s] < 1 then there is a weaker effect of space than time
beta_time_space[s,2] <- beta_time_space[s,1]*beta_mod[s]

#exponential model to estimate a strictly positive, multiplicative change in slope
## allows for consistent differences in the strength of the slope, while ignoring the direction of the effect for a given species
log(beta_mod[s]) <- e_beta_mod[s]
e_beta_mod[s] ~ dnorm(B_mod,tau_beta_mod)
## this multiplicative approach may fall apart if there is no time-slope
## that is, if beta_time_space[s,1] == 0, then beta_mod[s] would be imposible to estimate
## hopefully, the sharing of infor among species, and the prior would then shrink that species estimates for both slopes to 1
## it could work well in the opposite direction too, if space was the base-slope and the modification was to estimate time-slopes
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


####################################################### Alternative 3 - commented out
## alternative random-effect model that does not share information among species
## but has a regularizing prior to reduce differences for data-sparse species to 0
## this is tidy, but it assumes that the difference between the two kinds of slopes == 0
## if there are solid data to the contrary, this will show it
## this will also reduce the chances of seeing wild variation among species that are largely the result of sampling error
## i.e., it reduces the influence of poor data-species, but allows them to stay in the model.
#beta_time_space[s,2] <- beta_time_space[s,1]+beta_mod[s]
#beta_mod[s] ~ dnorm(0,tau_beta_mod)
## this same mean-zero random effect approach could also be used for the time-slope estimates to provide the same benefits for data-sparse species there
## e.g., 
## 
#######################################################


beta_dif[s] <- beta_time_space[s,1]-beta_time_space[s,2]


  }#s
 


 # hyperparameter for beta_mod treated as random effect
 B_mod ~ dnorm(0,0.1)
 tau_beta_mod ~ dgamma(0.001,0.001)
 #tau_beta_time_space1 ~ dgamma(0.001,0.001) # alternative random effects precision for the time-slopes
 
 sd_beta_mod <- 1/pow(tau_beta_mod,0.5)


}


"
cat(modl,file = "example.model.r")

jags_dat <- list(count = dat$count,
                 species = dat$species_f,
                 space_time = dat$space_time_f,
                 p_forest = dat$p_forest,
                 region = dat$region_f,
                 nspecies = nspecies,
                 nregions = nregions,
                 ncounts = ncounts)



parms <- c("sd_beta_mod",
           "B_mod",
           "beta_mod",
           "beta_time_space",
           "beta_dif",
           "alpha") 
  burnInSteps = 500            # Number of steps to "burn-in" the samplers. this is sufficient for testing, but you'll want to increase this
nChains = 3                   # Number of chains to run.
numSavedSteps=300         # Total number of steps in each chain to save. this is sufficient for testing, but you'll want to increase this
thinSteps=10                   # Number of steps to "thin" (1=keep every step).
nIter = ceiling( ( (numSavedSteps * thinSteps )+burnInSteps)) # Steps per chain.



out = jagsUI(data = jags_dat,
                  parameters.to.save = parms,
                  n.chains = 3,
                  n.burnin = burnInSteps,
                  n.thin = thinSteps,
                  n.iter = nIter,
                  parallel = T,
                  modules = NULL,
                  model.file = "example.model.r")



# out = jagsUI object, see the package help to explore
out$mean$beta_time_space #posterior means of the slope parameters

100*((B.space.time - out$mean$beta_time_space)/B.space.time) #% difference in posterior means and the true simulated values


out$summary #table showing a summary of the posterior distribution and some basic convergence stats for all the monitored parameters

out_ggs = ggs(out$samples)
ggmcmc(out_ggs,file = "convergence summaries.pdf", paparam_page = 8)

#### 
# library(bbsBayes)
# model_to_file(model = "slope",filename = "bbs_slope_model.r")
# model_to_file(model = "slope",filename = "bbs_slope_model_heavy_tails.r",heavy_tailed = TRUE)




