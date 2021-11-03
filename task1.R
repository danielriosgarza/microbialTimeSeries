##############Simulation conditions###################################

t.start = 0
t.end = 2000
t.step = 0.05
t.store = 500
norm = FALSE

simul.cond = list(t.start = t.start, t.end = t.end, t.step = t.step, t.store = t.store, norm = norm)

###############Common Parameters#####################################
n.species = 10
migration.p = 0.05 
metacommunity.eveness = 1 #vary from 1 (uneven) to 1000 (highly even)
metacommunity.p = rdirichlet(1, rep(.5, n.species)*metacommunity.eveness)
growth.rates <- rbeta(n.species, shape1 = 0.5, shape2 = 0.5)
#growth.rates <- rep(1, n.species)

###############For the per capita based models#####################################
n.individuals = 1000 
community.initial = rmultinom(1, n.individuals, rdirichlet(1, rep(1, n.species))) [,]
k.events = 1

######################Params Hubbel Model ########################



hub_params = c(simul.cond, list(community.initial = community.initial, migration.p = migration.p,
                                metacommunity.p = metacommunity.p, k.events = k.events))

hub <- parse(text = "do.call(simulateHubbell, hub_params)")


hub_model <- eval(hub)

makePlot(hub_model$matrix)