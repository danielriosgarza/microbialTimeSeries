
##############Simulation conditions###################################

t_start = 0
t_end = 2000
t_step = 0.05
t_store = 500
norm = FALSE

simul.cond = list(t_start = t_start, 
                  t_end = t_end, 
                  t_step = t_step, 
                  t_store = t_store, 
                  norm = norm)



###############Common parameters#####################################
n_species = 25
migration_p = 0.1 
metacommunity.eveness = 1 #vary from 1 (uneven) to 1000 (highly even)
metacommunity_probability = rdirichlet(1, rep(.5, n_species)*metacommunity.eveness)
growth_rates <- rbeta(n_species, shape1 = 0.5, shape2 = 0.5)
#growth_rates <- rep(1, n_species) 
x0 <- rep(0.001, n_species)
self.interactions <- runif(n_species, max(growth_rates)/3, max(growth_rates))

##############Stochastic ODE Models###################################
stochastic = FALSE
sigma_drift = 0.01
sigma_epoch = 0.1
sigma_external = 0.3
sigma_migration = 0.3
epoch_p = .1
n_external.perturbation = 1
t_external_events = getPerturbT(t_end, n_external.perturbation)
t_external_durations = rep(1, n_external.perturbation)

stoch.cond <- list(stochastic=stochastic, sigma_drift = sigma_drift,
    sigma_epoch = sigma_epoch, 
    sigma_external = sigma_external,
    sigma_migration = sigma_migration, 
    epoch_p = epoch_p, 
    t_external_events = t_external_events,
    t_external_durations = t_external_durations)


###############For the per capita based models#####################################
n.individuals = 1000 
community.initial = rmultinom(1, n.individuals, rdirichlet(1, rep(1, n_species)))
k_events = 1

######################Params Hubbel Model ########################



hub_params = c(simul.cond, 
    list(community.initial = community.initial, 
        migration_p = migration_p,
        metacommunity_probability = metacommunity_probability, 
        k_events = k_events))

hub <- parse(text = "do.call(simulateHubbell, hub_params)")



######################Params Hubbel Rates Model ########################


hub_rates_params = c(simul.cond, 
    list(community.initial = community.initial,
        migration_p = migration_p,
        metacommunity_probability = metacommunity_probability, 
        k_events = k_events,
        growth_rates = growth_rates))


hub_rates <- parse(text = "do.call(simulateHubbellRates, hub_rates_params)")


######################Params Logistic Model ########################
carrying_capacities <- growth_rates/self.interactions
death_rates = rep(0.001, n_species)


slm_params <- c(stoch.cond, simul.cond, 
    list(n_species = n_species, 
        growth_rates = growth_rates,
        carrying_capacities = carrying_capacities, 
        death_rates = death_rates,
        x0 = x0, 
        migration_p = migration_p))

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
interactions = rnorm(n=n_species^2, sd = sigma.A)


rA_params <- list(n_species = n_species, connectance = connectance, 
                  diagonal = diagonal, 
                  mutualism = mutualism.w,
                  commensalism = commensalism.w,
                  parasitism = parasitism.w,
                  amensalism = amensalism.w,
                  competition = competition.w,
                  scale = scale, interactions = interactions, symmetric = symmetric)

A = eval(parse(text = "do.call(randomA, rA_params)"))
##########################

glv_params <-  c(stoch.cond, simul.cond, list(n_species = n_species, A=A, 
                                              growth_rates = growth_rates,
                                              x0 = x0, migration_p=migration_p ))


glv <- parse(text = "do.call(simulateGLV, glv_params)")

######################Params Consumer Resource Model ########################

######Consumption/production matrix######
n_resources = as.integer(1.3*n_species)
consumption.w = 0 
production.w = 1
maintenance = 0

mean_production = production.w * n_resources
mean_consumption = consumption.w * n_resources


rE_params <- list(n_species = n_species, 
                  n_resources = n_resources, 
                  mean_production = mean_production, 
                  mean_consumption = mean_consumption, 
                  maintenance = maintenance)

E = eval(parse(text= 'do.call(randomE, rE_params)'))

###########################################
resource.eveness = 1
resource_dist = rdirichlet(1, rep(.5, n_resources)*resource.eveness)
resource_concentration = 1000
resources = resource_dist*resource_concentration
k.mul = 50
monod.k <- matrix(rgamma(n=n_species*n_resources, shape = k.mul*max(resources),
              rate = 1), nrow = n_species)


crm_params <- c(simul.cond, list(n_species = n_species, n_resources = n_resources, 
                                eff=E, x0=x0, resources = resources, 
                                growth_rates = growth_rates, monod.k = monod.k))

crm <- parse(text = "do.call(simulateConsumerResource, crm_params)")

