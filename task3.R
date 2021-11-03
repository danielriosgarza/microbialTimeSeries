##############Simulation conditions###################################

t.start = 0
t.end = 2000
t.step = 0.05
t.store = 500
norm = TRUE

simul.cond = list(t.start = t.start, t.end = t.end, t.step = t.step, t.store = t.store, norm = norm)



###############Common parameters#####################################
n.species = 20
migration.p = 0.5 
metacommunity.eveness = 8 #vary from 1 (uneven) to 1000 (highly even)
metacommunity.p = rdirichlet(1, rep(.5, n.species)*metacommunity.eveness)
growth.rates <- rbeta(n.species, shape1 = 0.5, shape2 = 0.5)
#growth.rates <- rep(1, n.species) 
x0 <- rep(0.001, n.species)
self.interactions <- runif(n.species, max(growth.rates)/3, max(growth.rates))

##############Stochastic ODE Models###################################
stochastic = FALSE
sigma.drift = 0.01
sigma.epoch = 0.1
sigma.external = 0.3
sigma.migration = 0.049
p.epoch = .1
n.external.perturbation = 1
t.external_events = getPerturbT(t.end, n.external.perturbation)
t.external_durations = rep(1, n.external.perturbation)

stoch.cond <- list(stochastic=stochastic, sigma.drift = sigma.drift,
                   sigma.epoch = sigma.epoch, sigma.external = sigma.external,
                   sigma.migration = sigma.migration, p.epoch = p.epoch, 
                   t.external_events = t.external_events,
                   t.external_durations = t.external_durations)


######################Params Logistic Model ########################
carrying.k <- growth.rates/self.interactions
death.rates = rep(0.023, n.species)


slm_params <- c(stoch.cond, simul.cond, list(n.species = n.species, growth.rates = growth.rates,
                                             carrying.k = carrying.k, death.rates = death.rates,
                                             x0 = x0, migration.p=migration.p ))

slm <- parse(text = "do.call(simulateStochasticLogistic, slm_params)")


#################
slm_m <- eval(slm)

########
slm_moments <- generateMoments(slm, n.instances = 25, t.store = 500, is.perCapita = FALSE)
-

###############For the per capita based models#####################################
n.individuals = 1000 
community.initial = rmultinom(1, n.individuals, rdirichlet(1, rep(1, n.species))) [,]
k.events = 1

######################Params Hubbel Model ########################



hub_params = c(simul.cond, list(community.initial = community.initial, migration.p = migration.p,
                                metacommunity.p = metacommunity.p, k.events = k.events))

hub <- parse(text = "do.call(simulateHubbell, hub_params)")

###################
hub_moments <- generateMoments(hub, n.instances = 25, t.store = 500, is.perCapita = TRUE)

##################

gr <- c(rep(1, 25), rep(2,25))

m <- rbind(slm_moments$basis, hub_moments$basis)

makeUMAP(m, gradient = gr)


#############################

growth.rates[1] = 5*max(growth.rates)
slm_params$growth.rates = growth.rates
slm_moments2 <- generateMoments(slm, n.instances = 25, t.store = 500, is.perCapita = FALSE)

##############################

gr <- c(rep(1, 25), rep(2,25), rep(3,25))

m <- rbind(slm_moments$basis, slm_moments2$basis, hub_moments$basis)

makeUMAP(m, gradient = gr)
