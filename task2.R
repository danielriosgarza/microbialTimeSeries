##############Simulation conditions###################################

t_start = 0
t_end = 2000
t_step = 0.01
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
n.individuals = 500 
community.initial = rmultinom(1, n.individuals, rdirichlet(1, rep(1, n_species))) [,]
k_events = 10

######################Params Hubbel Rates Model ########################

n.instances = 10

hub_rates_params = c(simul.cond, list(community.initial = community.initial, migration_p = migration_p,
                                      metacommunity.p = metacommunity.p, k_events = k_events,
                                      growth_rates = growth_rates))


hub_rates <- parse(text = "do.call(simulateHubbellRates, hub_rates_params)")


hubR_moments <- generateMoments(hub_rates, n.instances = n.instances, t_store = 500, is.perCapita = TRUE)


growth_rates <- rep(1, n_species)
hub_rates_params$growth_rates = growth_rates = growth_rates
hub_rates_neutral <- parse(text = "do.call(simulateHubbellRates, hub_rates_params)")

hubR_neutral_moments <- generateMoments(hub_rates, n.instances = n.instances, t_store = 500, is.perCapita = TRUE)

hub_params = c(simul.cond, list(community.initial = community.initial, migration_p = migration_p,
                                metacommunity.p = metacommunity.p, k_events = k_events))

hub <- parse(text = "do.call(simulateHubbell, hub_params)")

hub_moments <- generateMoments(hub, n.instances = n.instances, t_store = 500, is.perCapita = TRUE)

#########UMAP##########

gr <- c(rep(1, n.instances),rep(2,n.instances), rep(3,n.instances))

m<-rbind(hubR_moments$basis, hubR_neutral_moments$basis, hub_moments$basis)

makeUMAP(m, gradient = gr)
