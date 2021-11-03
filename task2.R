##############Simulation conditions###################################

t.start = 0
t.end = 2000
t.step = 0.01
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
n.individuals = 500 
community.initial = rmultinom(1, n.individuals, rdirichlet(1, rep(1, n.species))) [,]
k.events = 10

######################Params Hubbel Rates Model ########################

n.instances = 10

hub_rates_params = c(simul.cond, list(community.initial = community.initial, migration.p = migration.p,
                                      metacommunity.p = metacommunity.p, k.events = k.events,
                                      growth.rates = growth.rates))


hub_rates <- parse(text = "do.call(simulateHubbellRates, hub_rates_params)")


hubR_moments <- generateMoments(hub_rates, n.instances = n.instances, t.store = 500, is.perCapita = TRUE)


growth.rates <- rep(1, n.species)
hub_rates_params$growth.rates = growth.rates = growth.rates
hub_rates_neutral <- parse(text = "do.call(simulateHubbellRates, hub_rates_params)")

hubR_neutral_moments <- generateMoments(hub_rates, n.instances = n.instances, t.store = 500, is.perCapita = TRUE)

hub_params = c(simul.cond, list(community.initial = community.initial, migration.p = migration.p,
                                metacommunity.p = metacommunity.p, k.events = k.events))

hub <- parse(text = "do.call(simulateHubbell, hub_params)")

hub_moments <- generateMoments(hub, n.instances = n.instances, t.store = 500, is.perCapita = TRUE)

#########UMAP##########

gr <- c(rep(1, n.instances),rep(2,n.instances), rep(3,n.instances))

m<-rbind(hubR_moments$basis, hubR_neutral_moments$basis, hub_moments$basis)

makeUMAP(m, gradient = gr)
