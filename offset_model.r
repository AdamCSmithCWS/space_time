
model
{
 
  for( k in 1:ncounts_obs ){
    log(lambda[k]) <- obs_offset[obs[k]] + species_effect[species[k]] + region_effect[region[k]] + route_effect[route[k]] + noise[k] #noise is the overdispersion parameter
    
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
  
  
  ############### observer-offsets
  for( o in 1:nobs){
    obs_offset[o] ~ dnorm(0,tau_obs)
  }
  sd_obs ~ dt(0, 1, 4) T(0,)     #half-t prior on the standard deviation
  tau_obs <- 1/pow(sd_obs,2)

  ############### route-offsets
  for( t in 1:nroutes){
    route_effect[t] ~ dnorm(0,tau_route) #random effect sharing information among routes
  }
  sd_route ~ dt(0, 1, 4) T(0,)     #half-t prior on the standard deviation
  tau_route <- 1/pow(sd_route,2)
  
  
  ############### species_effects
  for( s in 1:nspecies){
    species_effect[s] ~ dnorm(0,1) #fixed-effects so  no information shared among species
  }
  
  ############### region_effects
  for( r in 1:nregions){
    region_effect[r] ~ dnorm(0,1) #fixed-effects so  no information shared among species
  }
  
  
  
  
  } 
