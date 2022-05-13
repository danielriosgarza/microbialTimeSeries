##############Simulation conditions###################################

t_start = 0
t_end = 2000
t_step = 0.05
t_store = 500
norm = TRUE

simul.cond = list(t_start = t_start, t_end = t_end, t_step = t_step, t_store = t_store, norm = norm)


n.instances = 10

###############Common parameters#####################################
n_species = 25
migration_p = 0
metacommunity.eveness = 1 #vary from 1 (uneven) to 1000 (highly even)
metacommunity.p = rdirichlet(1, rep(.5, n_species)*metacommunity.eveness)
growth_rates <- rbeta(n_species, shape1 = 0.5, shape2 = 0.5)
#growth_rates <- rep(1, n_species) 
x0 <- rep(0.001, n_species)
self.interactions <- runif(n_species, max(growth_rates)/3, max(growth_rates))

##############Stochastic ODE Models###################################
stochastic = FALSE
sigma_drift = 0.01
sigma_epoch = 0.011
sigma_external = 0
sigma.migration = 0
p.epoch = 0
n_external.perturbation = 0
t_external_events = getPerturbT(t_end, n_external.perturbation)
t_external_durations = rep(1, n_external.perturbation)

stoch.cond <- list(stochastic=stochastic, sigma_drift = sigma_drift,
                   sigma_epoch = sigma_epoch, sigma_external = sigma_external,
                   sigma.migration = sigma.migration, p.epoch = p.epoch, 
                   t_external_events = t_external_events,
                   t_external_durations = t_external_durations)


######################Params Logistic Model ########################
carrying.k <- growth_rates/self.interactions
death_rates = rep(0, n_species)


slm_params <- c(stoch.cond, simul.cond, list(n_species = n_species, growth_rates = growth_rates,
                                             carrying.k = carrying.k, death_rates = death_rates,
                                             x0 = x0, migration_p=migration_p ))

slm <- parse(text = "do.call(simulateStochasticLogistic, slm_params)")

slm_det <- generateMoments(slm, n.instances = n.instances, t_store = 500, is.perCapita = FALSE)

slm_params$stochastic = TRUE

slm_stoch <- generateMoments(slm, n.instances = n.instances, t_store = 500, is.perCapita = FALSE)

###############

alpha_det <- diversity(slm_det$basis)
alpha_stoch <- diversity(slm_stoch$basis)

ra.det <- rankabundance(as.data.frame(slm_det$basis))
ra.stoch <- rankabundance(as.data.frame(slm_stoch$basis))

rankabunplot(ra.det)
rankabunplot(ra.stoch, addit = TRUE)

##################
