##############Simulation conditions###################################

t.start = 0
t.end = 2000
t.step = 0.05
t.store = 500
norm = TRUE

simul.cond = list(t.start = t.start, t.end = t.end, t.step = t.step, t.store = t.store, norm = norm)


n.instances = 10

###############Common parameters#####################################
n.species = 25
migration.p = 0
metacommunity.eveness = 1 #vary from 1 (uneven) to 1000 (highly even)
metacommunity.p = rdirichlet(1, rep(.5, n.species)*metacommunity.eveness)
growth.rates <- rbeta(n.species, shape1 = 0.5, shape2 = 0.5)
#growth.rates <- rep(1, n.species) 
x0 <- rep(0.001, n.species)
self.interactions <- runif(n.species, max(growth.rates)/3, max(growth.rates))

##############Stochastic ODE Models###################################
stochastic = FALSE
sigma.drift = 0.01
sigma.epoch = 0.011
sigma.external = 0
sigma.migration = 0
p.epoch = 0
n.external.perturbation = 0
t.external_events = getPerturbT(t.end, n.external.perturbation)
t.external_durations = rep(1, n.external.perturbation)

stoch.cond <- list(stochastic=stochastic, sigma.drift = sigma.drift,
                   sigma.epoch = sigma.epoch, sigma.external = sigma.external,
                   sigma.migration = sigma.migration, p.epoch = p.epoch, 
                   t.external_events = t.external_events,
                   t.external_durations = t.external_durations)


######################Params Logistic Model ########################
carrying.k <- growth.rates/self.interactions
death.rates = rep(0, n.species)


slm_params <- c(stoch.cond, simul.cond, list(n.species = n.species, growth.rates = growth.rates,
                                             carrying.k = carrying.k, death.rates = death.rates,
                                             x0 = x0, migration.p=migration.p ))

slm <- parse(text = "do.call(simulateStochasticLogistic, slm_params)")

slm_det <- generateMoments(slm, n.instances = n.instances, t.store = 500, is.perCapita = FALSE)

slm_params$stochastic = TRUE

slm_stoch <- generateMoments(slm, n.instances = n.instances, t.store = 500, is.perCapita = FALSE)

###############

alpha_det <- diversity(slm_det$basis)
alpha_stoch <- diversity(slm_stoch$basis)

ra.det <- rankabundance(as.data.frame(slm_det$basis))
ra.stoch <- rankabundance(as.data.frame(slm_stoch$basis))

rankabunplot(ra.det)
rankabunplot(ra.stoch, addit = TRUE)

##################
