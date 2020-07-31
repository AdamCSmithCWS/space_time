## ideas for space-time comparison model
library(jagsUI)
library(tidyverse)
library(ggmcmc)


nreps = 21

dat <- expand.grid(species = c("Red-eyed Vireo","American Robin","Ovenbird","Wood Thrush","Hermit Thrush","Black-capped Chickadee"),
                   region = c("BCR12","BCR13","BCR23","BCR34"),
                   space_time = rep(c("time","space"),nreps),
                   p_forest = seq(0.1,0.8,length.out = 13))


for(i in c("species","region","space_time")){
dat[,paste0(i,"_f")] = as.integer(dat[,i])
}

nspecies = length(unique(dat$species))
nregions = length(unique(dat$region))
B = c(1,0.5,2,-2,-1.5,1) #mean species slopes for p_forest
True.mean.dif <- 0.9
True.dif.B <- rnorm(nspecies,True.mean.dif,0.05)

B.space.time = matrix(c(B,True.dif.B*B),nrow = nspecies,ncol = 2,byrow = FALSE) #random 1%/year variance around the specie-leve mean response based on time or space

a.sp = runif(nspecies,0.5,2) #coarse mean counts overall

a.reg = rnorm(nregions,0,0.1)#random variation in abundance among regions

noise <- 0.001 #itnroduced overdispersion

for(i in 1:nrow(dat)){
  forst.eff <- (B.space.time[dat[i,"species_f"],dat[i,"space_time_f"]]*(dat[i,"p_forest"]-0.5)) #this -0.5 just centers the values of percent forest without changing the scale of the variable, it also means the intercepts represent the mean values of abundance and not some hypoethical mean count when forest == 0
  int <- a.sp[dat[i,"species_f"]]+a.reg[dat[i,"region_f"]]
  lmbd <- exp(int+forst.eff + rnorm(1,0,noise))
  dat[i,"count"] <- rpois(1,lambda = lmbd)
}

ncounts = nrow(dat)

modl <- "
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
           "beta_space_time",
           "beta_dif",
           "alpha") 
  burnInSteps = 500            # Number of steps to "burn-in" the samplers.
nChains = 3                   # Number of chains to run.
numSavedSteps=300         # Total number of steps in each chain to save.
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



out$mean$beta_space_time


# library(bbsBayes)
# model_to_file(model = "slope",filename = "bbs_slope_model.r")
# model_to_file(model = "slope",filename = "bbs_slope_model_heavy_tails.r",heavy_tailed = TRUE)




