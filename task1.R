##############Simulation conditions###################################

t_start = 0
t_end = 2000
t_step = 0.05
t_store = 500
norm = FALSE

simul.cond = list(t_start = t_start, t_end = t_end, t_step = t_step, t_store = t_store, norm = norm)

###############Common Parameters#####################################
n_species = 10
migration_p = 0.05 
metacommunity.eveness = 1 #vary from 1 (uneven) to 1000 (highly even)
metacommunity.p = rdirichlet(1, rep(.5, n_species)*metacommunity.eveness)
growth_rates <- rbeta(n_species, shape1 = 0.5, shape2 = 0.5)
#growth_rates <- rep(1, n_species)

###############For the per capita based models#####################################
n.individuals = 1000 
community.initial = rmultinom(1, n.individuals, rdirichlet(1, rep(1, n_species))) [,]
k_events = 1

######################Params Hubbel Model ########################



hub_params = c(simul.cond, list(community.initial = community.initial, migration_p = migration_p,
                                metacommunity.p = metacommunity.p, k_events = k_events))

hub <- parse(text = "do.call(simulateHubbell, hub_params)")


hub_model <- eval(hub)

makePlot(hub_model$matrix)


hub_moments <- generateMoments(hub, n.instances = 50, t_store = 500, is.perCapita = TRUE)


########
evenness_gradient <- seq(1,20, 1)
shannon_diversity <- rep(0, length(evenness_gradient))




basisM <- matrix(0, ncol=n_species, nrow = length(evenness_gradient))

for (i in seq_along(evenness_gradient)){
  metacommunity.eveness = evenness_gradient[i]
  metacommunity.p = rdirichlet(1, rep(.5, n_species)*metacommunity.eveness)
  hub_params = c(simul.cond, list(community.initial = community.initial, migration_p = migration_p,
                                  metacommunity.p = metacommunity.p, k_events = k_events))
  
  hub <- parse(text = "do.call(simulateHubbell, hub_params)")
  
  moment <- generateMoments(hub, n.instances = 3, t_store = 500, is.perCapita = TRUE)
  shannon <- diversity(moment$basis)
  shannon_diversity[i]<-mean(shannon)
  
  basisM[i,] <- colMeans(moment$basis)
  
  
}
makeUMAP(basisM, gradient = evenness_gradient)
