##############Simulation conditions###################################

t_start = 0
t_end = 2000
t_step = 0.05
t_store = 500
norm = TRUE

simul.cond = list(t_start = t_start, t_end = t_end, t_step = t_step, t_store = t_store, norm = norm)



###############Common parameters#####################################
n_species = 20
migration_p = 0.5 
metacommunity.eveness = 8 #vary from 1 (uneven) to 1000 (highly even)
metacommunity.p = rdirichlet(1, rep(.5, n_species)*metacommunity.eveness)
growth_rates <- rbeta(n_species, shape1 = 0.5, shape2 = 0.5)
#growth_rates <- rep(1, n_species) 
x0 <- rep(0.001, n_species)
self.interactions <- runif(n_species, max(growth_rates)/3, max(growth_rates))

##############Stochastic ODE Models###################################
stochastic = FALSE
sigma_drift = 0.01
sigma_epoch = 0.1
sigma_external = 0.3
sigma.migration = 0.049
p.epoch = .1
n_external.perturbation = 1
t_external_events = getPerturbT(t_end, n_external.perturbation)
t_external_durations = rep(1, n_external.perturbation)

stoch.cond <- list(stochastic=stochastic, sigma_drift = sigma_drift,
                   sigma_epoch = sigma_epoch, sigma_external = sigma_external,
                   sigma.migration = sigma.migration, p.epoch = p.epoch, 
                   t_external_events = t_external_events,
                   t_external_durations = t_external_durations)


######################Params Logistic Model ########################
carrying.k <- growth_rates/self.interactions
death_rates = rep(0.023, n_species)


slm_params <- c(stoch.cond, simul.cond, list(n_species = n_species, growth_rates = growth_rates,
                                             carrying.k = carrying.k, death_rates = death_rates,
                                             x0 = x0, migration_p=migration_p ))

slm <- parse(text = "do.call(simulateStochasticLogistic, slm_params)")


#################
slm_m <- eval(slm)

########
slm_moments <- generateMoments(slm, n.instances = 25, t_store = 500, is.perCapita = FALSE)
-

###############For the per capita based models#####################################
n.individuals = 1000 
community.initial = rmultinom(1, n.individuals, rdirichlet(1, rep(1, n_species))) [,]
k_events = 1

######################Params Hubbel Model ########################



hub_params = c(simul.cond, list(community.initial = community.initial, migration_p = migration_p,
                                metacommunity.p = metacommunity.p, k_events = k_events))

hub <- parse(text = "do.call(simulateHubbell, hub_params)")

###################
hub_moments <- generateMoments(hub, n.instances = 25, t_store = 500, is.perCapita = TRUE)

##################

gr <- c(rep(1, 25), rep(2,25))

m <- rbind(slm_moments$basis, hub_moments$basis)

makeUMAP(m, gradient = gr)


#############################

growth_rates[1] = 5*max(growth_rates)
slm_params$growth_rates = growth_rates
slm_moments2 <- generateMoments(slm, n.instances = 25, t_store = 500, is.perCapita = FALSE)

##############################

gr <- c(rep(1, 25), rep(2,25), rep(3,25))

m <- rbind(slm_moments$basis, slm_moments2$basis, hub_moments$basis)

makeUMAP(m, gradient = gr)
