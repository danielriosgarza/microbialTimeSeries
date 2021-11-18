
##############Simulation conditions###################################

t.start = 0
t.end = 2000
t.step = 0.05
t.store = 500
norm = FALSE

simul.cond = list(t.start = t.start, 
                  t.end = t.end, 
                  t.step = t.step, 
                  t.store = t.store, 
                  norm = norm)



###############Common parameters#####################################
n.species = 25
migration.p = 0.1 
metacommunity.eveness = 1 #vary from 1 (uneven) to 1000 (highly even)
metacommunity.probability = rdirichlet(1, rep(.5, n.species)*metacommunity.eveness)
growth.rates <- rbeta(n.species, shape1 = 0.5, shape2 = 0.5)
#growth.rates <- rep(1, n.species) 
x0 <- rep(0.001, n.species)
self.interactions <- runif(n.species, max(growth.rates)/3, max(growth.rates))

##############Stochastic ODE Models###################################
stochastic = FALSE
sigma.drift = 0.01
sigma.epoch = 0.1
sigma.external = 0.3
sigma.migration = 0.3
epoch.p = .1
n.external.perturbation = 1
t.external_events = getPerturbT(t.end, n.external.perturbation)
t.external_durations = rep(1, n.external.perturbation)

stoch.cond <- list(stochastic=stochastic, sigma.drift = sigma.drift,
    sigma.epoch = sigma.epoch, 
    sigma.external = sigma.external,
    sigma.migration = sigma.migration, 
    epoch.p = epoch.p, 
    t.external_events = t.external_events,
    t.external_durations = t.external_durations)


###############For the per capita based models#####################################
n.individuals = 1000 
community.initial = rmultinom(1, n.individuals, rdirichlet(1, rep(1, n.species)))
k.events = 1

######################Params Hubbel Model ########################



hub_params = c(simul.cond, 
    list(community.initial = community.initial, 
        migration.p = migration.p,
        metacommunity.probability = metacommunity.probability, 
        k.events = k.events))

hub <- parse(text = "do.call(simulateHubbell, hub_params)")



######################Params Hubbel Rates Model ########################


hub_rates_params = c(simul.cond, 
    list(community.initial = community.initial,
        migration.p = migration.p,
        metacommunity.probability = metacommunity.probability, 
        k.events = k.events,
        growth.rates = growth.rates))


hub_rates <- parse(text = "do.call(simulateHubbellRates, hub_rates_params)")


######################Params Logistic Model ########################
carrying.capacities <- growth.rates/self.interactions
death.rates = rep(0.001, n.species)


slm_params <- c(stoch.cond, simul.cond, 
    list(n.species = n.species, 
        growth.rates = growth.rates,
        carrying.capacities = carrying.capacities, 
        death.rates = death.rates,
        x0 = x0, 
        migration.p = migration.p))

slm <- parse(text = "do.call(simulateStochasticLogistic, slm_params)")

######################Params GLV Model ########################

####interaction matrix####
symmetric = FALSE
connectance = 0.5
diagonal = -1*self.interactions
mutualism.w = 1
commensalism.w = 1
parasitism.w = 1
amensalism.w = 1
competition.w = 1
interaction.w <- c(mutualism.w, commensalism.w, parasitism.w, amensalism.w,  competition.w)
scale = 0.5
sigma.A = 1
interactions = rnorm(n=n.species^2, sd = sigma.A)


rA_params <- list(n.species = n.species, connectance = connectance, 
                  diagonal = diagonal, 
                  mutualism = mutualism.w,
                  commensalism = commensalism.w,
                  parasitism = parasitism.w,
                  amensalism = amensalism.w,
                  competition = competition.w,
                  scale = scale, interactions = interactions, symmetric = symmetric)

A = eval(parse(text = "do.call(randomA, rA_params)"))
##########################

glv_params <-  c(stoch.cond, simul.cond, list(n.species = n.species, A=A, 
                                              growth.rates = growth.rates,
                                              x0 = x0, migration.p=migration.p ))


glv <- parse(text = "do.call(simulateGLV, glv_params)")

######################Params Consumer Resource Model ########################

######Consumption/production matrix######
n.resources = as.integer(1.3*n.species)
consumption.w = 0 
production.w = 1
maintenance = 0

mean.production = production.w * n.resources
mean.consumption = consumption.w * n.resources


rE_params <- list(n.species = n.species, 
                  n.resources = n.resources, 
                  mean.production = mean.production, 
                  mean.consumption = mean.consumption, 
                  maintenance = maintenance)

E = eval(parse(text= 'do.call(randomE, rE_params)'))

###########################################
resource.eveness = 1
resource_dist = rdirichlet(1, rep(.5, n.resources)*resource.eveness)
resource_concentration = 1000
resources = resource_dist*resource_concentration
k.mul = 50
monod.k <- matrix(rgamma(n=n.species*n.resources, shape = k.mul*max(resources),
              rate = 1), nrow = n.species)


crm_params <- c(simul.cond, list(n.species = n.species, n.resources = n.resources, 
                                eff=E, x0=x0, resources = resources, 
                                growth.rates = growth.rates, monod.k = monod.k))

crm <- parse(text = "do.call(simulateConsumerResource, crm_params)")

