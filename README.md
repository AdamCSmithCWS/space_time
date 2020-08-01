# space_time
toy example of a space for time substitution model

Example of some alternative model structures for an over-dispersed Poisson regression that would compare the slopes for a time-series effect of changing forest cover to the slopes for a spatial-series effect of changing forest cover
A few alternatives for how to structure the comparison including a multiplicative difference between temporal slopes and spatial slopes that shares information among species without requiring any assumption about the direction of the species response to forest cover.

The main script sets up some toy data. Defines the model in JAGS. and runs the model in parallel on 3 cores and provides some initial summaries and convergence checks using ggmcmc package




