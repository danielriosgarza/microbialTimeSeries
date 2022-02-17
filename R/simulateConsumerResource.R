#' Consumer-resource model simulation
#'
#' Simulates time series with the consumer-resource model
#'
#' @param n.species Integer: number of species
#' @param n.resources Interger: number of resources
#' @param names.species Character: names of species. If NULL,
#' `paste0("sp", seq_len(n.species))` is used.
#' (default: \code{names.species = NULL})
#' @param names.resources Character: names of resources. If NULL,
#' `paste0("res", seq_len(n.resources))` is used.
#' @param E matrix: matrix of efficiency. A matrix defining the efficiency of
#' resource-to-biomass conversion (positive values) and the relative conversion
#' of metabolic by-products (negative values). If NULL, 
#' `randomE(n.species, n.resources)` is used.
#' (default: \code{E = NULL})
#' @param x0 Numeric: initial abundances of simulated species. If NULL,
#' `runif(n = n.species, min = 0.1, max = 10)` is used.
#' (default: \code{x0 = NULL})
#' @param resources Numeric: initial concentrations of resources. If NULL,
#' `runif(n = n.resources, min = 1, max = 100)` is used.
#' (default: \code{resources = NULL})
#' @param growth.rates Numeric: vector of maximum growth rates(mu) of species.
#' If NULL, `rep(1, n.species)` is used.
#' (default: \code{growth.rates = NULL})
#' @param monod.constant matrix: the constant of additive monod growth of
#' n.species consuming n.resources. If NULL,
#' `matrix(rgamma(n = n.species*n.resources, shape = 50*max(resources), rate = 1), nrow = n.species)`
#' is used.
#' (default: \code{monod.constant = NULL})
#' @param sigma.drift Numeric: standard deviation of a normally distributed 
#' noise applied in each time step (t.step)
#' (default: \code{sigma.drift = 0.001})
#' @param sigma.epoch Numeric: standard deviation of a normally distributed 
#' noise applied to random periods of the community composition with frequency 
#' defined by the epoch.p parameter
#' (default: \code{sigma.epoch = 0.1})
#' @param sigma.external Numeric: standard deviation of a normally distributed 
#' noise applied to user-defined external events/disturbances
#' (default: \code{sigma.external = 0.3})
#' @param sigma.migration Numeric: standard deviation of a normally distributed 
#' variable that defines the intensity of migration at each time step (t.step)
#' (default: \code{sigma.migration = 0.01})
#' @param epoch.p Numeric: the probability/frequency of random periodic 
#' changes introduced to the community composition
#' (default: \code{epoch.p = 0.001})
#' @param t.external_events Numeric: the starting time points of defined 
#' external events that introduce random changes to the community composition
#' (default: \code{t.external_events = NULL})
#' @param t.external_durations Numeric: respective duration of the external 
#' events that are defined in the 't.external_events' (times) and 
#' sigma.external (std).
#' (default: \code{t.external_durations = NULL})
#' @param stochastic Logical: whether to introduce noise in the simulation.
#' If False, sigma.drift, sigma.epoch, and sigma.external are ignored.
#' (default: \code{stochastic = FALSE})
#' @param migration.p Numeric: the probability/frequency of migration from a 
#' metacommunity.
#' (default: \code{migration.p = 0.01})
#' @param metacommunity.probability Numeric: Normalized probability distribution
#' of the likelihood that species from the metacommunity can enter the community
#' during the simulation. If NULL, `rdirichlet(1, alpha = rep(1,n.species))` is 
#' used.
#' (default: \code{metacommunity.probability = NULL})
#' @param error.variance Numeric: the variance of measurement error.
#' By default it equals to 0, indicating that the result won't contain any 
#' measurement error. This value should be non-negative.
#' (default: \code{error.variance = 0})
#' @param norm Logical scalar: returns normalized abundances (proportions
#' in each generation) 
#' (default: \code{norm = FALSE})
#' @param t.end Numeric scalar indicating the final time of the simulation
#' (default: \code{t.end = 1000})
#' @param trophic.priority Matrix: a matrix defining the orders of resources to 
#' be consumed by each species. If NULL, by default, this feature won't be 
#' turned on, and species will consume all resources simultaneously to grow. 
#' The dimension should be identical to matrix E.
#' (default: \code{trophic.priority = NULL})
#' @param ... additional parameters, see \code{\link{utils}} to know more.
#'
#' @examples
#' ExampleConsumerResource <- simulateConsumerResource(n.species = 2, 
#'     n.resources = 4)
#' makePlot(ExampleConsumerResource$matrix)
#' makePlot(ExampleConsumerResource$resources)
#' 
#' # example to get relative abundance and relative proportion of resources
#' ExampleConsumerResource <- simulateConsumerResource(n.species = 2, 
#'     n.resources = 4, norm = TRUE)
#' makePlot(ExampleConsumerResource$matrix)
#' makePlot(ExampleConsumerResource$resources)
#'
#' # example with user-defined values (names.species, names.resources, E, x0, 
#' # resources, growth.rates, error.variance, t.end, t.step)
#' n.resources <- 6
#' ExampleE <- randomE(n.species = 4, n.resources = n.resources, mean.con = 3, 
#'     mean.prod = 1, maintenance = 0.4)
#' ExampleResources <- rep(100, n.resources)
#' ExampleConsumerResource <- simulateConsumerResource(n.species = 4, 
#'     n.resources = 6, names.species = letters[1:4], 
#'     names.resources = paste0("res",LETTERS[1:6]), E = ExampleE, 
#'     x0 = rep(0.001, 4), resources = ExampleResources, 
#'     growth.rates <- runif(4),
#'     error.variance = 1,
#'     t.end = 5000, t.step = 1)
#' makePlot(ExampleConsumerResource$matrix)
#' makePlot(ExampleConsumerResource$resources)
#' 
#' # example with trophic levels
#' n.species <- 10
#' n.resources <- 15
#' 
#' ExampleEfficiencyMatrix <- randomE(n.species = 10, n.resources = 15,
#'                                    trophic.levels = c(6,3,1),
#'                                    trophic.preferences = list(c(rep(1,5), rep(-1, 5), rep(0, 5)), 
#'                                                               c(rep(0,5), rep(1, 5), rep(-1, 5)), 
#'                                                               c(rep(0,10), rep(1, 5))))
#' makeHeatmap(ExampleEfficiencyMatrix)
#' 
#' # ExampleResources <- rep(100, n.resources)
#' ExampleResources <- c(rep(500, 5), rep(200, 5), rep(50, 5))
#' ExampleConsumerResource <- simulateConsumerResource(n.species = n.species, 
#'                                                     n.resources = n.resources, names.species = letters[1:n.species], 
#'                                                     names.resources = paste0("res",LETTERS[1:n.resources]), E = ExampleEfficiencyMatrix, 
#'                                                     x0 = rep(0.001, n.species), resources = ExampleResources, 
#'                                                     growth.rates = rep(1, n.species),
#'                                                     error.variance = 0.001,
#'                                                     t.end = 5000, t.step = 1)
#' makePlot(ExampleConsumerResource$matrix)
#' makePlotRes(ExampleConsumerResource$resources)
#' 
#' 
#' # example with trophic priority
#' n.species <- 4
#' n.resources <- 6
#' ExampleE <- randomE(n.species = n.species, n.resources = n.resources, mean.consumption = n.resources, mean.production = 0)
#' ExampleTrophicPriority <- t(apply(matrix(sample(n.species * n.resources), nrow = n.species), 1, order))
#' # make sure that for non-consumables resources for each species, 
#' # the priority is 0 (smaller than any given priority) 
#' ExampleTrophicPriority <- (ExampleE > 0) * ExampleTrophicPriority
#' ExampleConsumerResource <- simulateConsumerResource(n.species = n.species, 
#'                                                     n.resources = n.resources, 
#'                                                     E = ExampleE, 
#'                                                     trophic.priority = ExampleTrophicPriority,
#'                                                     t.end = 2000
#'                                                     )
#' makePlot(ExampleConsumerResource$matrix)
#' makePlotRes(ExampleConsumerResource$resources)
#' @docType methods
#' @aliases simulateConsumerResource-numeric
#' 
#' @importFrom deSolve ode
#'
#' @export
simulateConsumerResource <- function(n.species, n.resources,
    names.species = NULL,
    names.resources = NULL,
    E = NULL,
    x0 = NULL,
    resources = NULL,
    resources.dilution = NULL,
    growth.rates = NULL,
    monod.constant = NULL,
    dilution.rate = 0,
    sigma.drift = 0.001,
    sigma.epoch = 0.1,
    sigma.external = 0.3,
    sigma.migration = 0.01,
    epoch.p = 0.001,
    t.external_events = NULL,
    t.external_durations = NULL,
    stochastic = FALSE,
    migration.p = 0.01, 
    metacommunity.probability = NULL,
    error.variance = 0,
    norm = FALSE,
    t.end = 1000,
    trophic.priority = NULL,
    ...){
    
    t.dyn <- SimulationTimes(t.end = t.end,...)
    
    # calculate the time points influenced by the disturbances
    tEvent = eventTimes(
        t.events = t.external_events,
        t.duration = t.external_durations,
        t.end = t.end,
        ... = ...)
    
    
    # define the consumer-resource model
    
    consumerResourceModel <- function(t, state, params){
        with(as.list(c(state, params)),{
            
            x0 <- pmax(0, state[startsWith(names(state), "consumer")])
            resources <- pmax(0, state[startsWith(names(state), "resource")])
            growth.rates <- params[['growth.rates']]
            E <- params[['E']]
            monod.constant <- params[['monod.constant']]
            
            resources.dilution <- params[['resources.dilution']]
            dilution.rate <- params[['dilution.rate']]
            trophic.priority <- params[['trophic.priority']]
            
            if(!is.null(trophic.priority)){
                if (!identical(dim(ExampleE), dim(ExampleTrophicPriority))){
                    stop("The dimension of trophic priority is not correct.")
                }
                # modify E in each step
                Emod <- trophic.priority
                ## resources <= a relatively small number(instead of 0) ######
                Emod[, resources <= 0.1] <- 0 
                Emod <- E*(t(apply(Emod, 1, getRowMax))>0) + E*(E<0)
                E <- Emod
            }

            B <- matrix(rep(resources, length(x0)),
                        ncol = length(resources), byrow = TRUE) + monod.constant
            growth <- ((E*(E>0)/B) %*% resources)*x0
            
            
            consumption <- -resources*(t(1/B*(E>0))%*%x0)
            
            production <- -t(E*(E<0)) %*% growth
            
            
            #consumption <- (t(growth) %*% ((E>0)/B))*resources
            #production <- -(t(growth) %*% (E*(E<0)/B))*resources
            dResources <- consumption + production - dilution.rate*(resources - resources.dilution)
            dConsumers <- growth.rates*growth - dilution.rate*x0
            dxdt <- list(c(dConsumers, dResources))
            return(dxdt)
        })
    }
    
    # set the default values
    if (is.null(names.species)) {
        names.species <- paste0("sp", seq_len(n.species))
    }
    if (is.null(names.resources)) {
        names.resources <- paste0("res", seq_len(n.resources))
    }
    if (is.null(E)) {
        E <- randomE(n.species, n.resources)
    }
    if (is.null(x0)) {
        x0 <- runif(n = n.species, min = 0.1, max = 10)
    }
    if (is.null(resources)) {
        resources <- runif(n = n.resources, min = 1, max = 100)
    }
    if (is.null(resources.dilution)) {
        resources.dilution <- resources
    }
    if (is.null(growth.rates)) {
        growth.rates <- rep(1, n.species)
    }
    if (is.null(monod.constant)) {
        monod.constant <- matrix(rgamma(n = n.species*n.resources,
                                        shape = 50*max(resources), rate = 1), nrow = n.species)
    }
    if (is.null(metacommunity.probability)) {
        metacommunity.probability <- rdirichlet(1, alpha = rep(1,n.species))
    }
    
    #normalize metacommunity.probability
    metacommunity.probability <- metacommunity.probability/
        sum(metacommunity.probability)
    
    # define the perturbation event
    perturb <- function(t, y, parameters){
        with(as.list(y),{
            #continuous or episodic perturbation
            epoch.rN <- 0
            external.rN <- 0
            migration.rN <- 0
            if (rbinom(1,1, parameters$epoch.p)){
                epoch.rN <- rnorm(parameters$n.species, sd=parameters$sigma.epoch)
                epoch.rN <- parameters$stochastic*epoch.rN
            }
            
            if (rbinom(1,1, parameters$migration.p)){
                migration.rN <- rmultinom(n = 1, size = 1, 
                                          prob = parameters$metacommunity.probability)[,]*abs(rnorm(n=1, 
                                                                                                    mean=0, sd = parameters$sigma.migration))
                
            }
            
            if (t %in% parameters$tEvent){
                external.rN <- rnorm(parameters$n.species,
                                     sd=parameters$sigma.external)
                external.rN <- parameters$stochastic*external.rN
            }
            drift.rN <- rnorm(parameters$n.species, sd=parameters$sigma.drift)
            drift.rN <- parameters$stochastic*drift.rN
            
            
            #perturbation is applied to the current population
            consumer <- pmax(y[grepl('consumer', names(y))], 0)
            
            consumer <- consumer*(1+drift.rN)*(1+epoch.rN)*(1+external.rN)+ migration.rN
            resource <- y[grepl('resource', names(y))]

            return(c(consumer, resource))})
    }
    
    
    state.init <- c(x0, resources)
    names(state.init) <- c(paste0("consumer", seq(length(x0))),
                           paste0("resource", seq(length(resources))))
    
    parameters <- list(growth.rates = growth.rates, E = E, monod.constant = monod.constant,  
                       n.species = n.species, sigma.drift = sigma.drift, stochastic = stochastic,
                       sigma.epoch = sigma.epoch, epoch.p = epoch.p,
                       sigma.external = sigma.external, tEvent = tEvent,
                       migration.p = migration.p, metacommunity.probability = metacommunity.probability,
                       sigma.migration = sigma.migration, 
                       resources.dilution = resources.dilution,
                       dilution.rate = dilution.rate,
                       trophic.priority = trophic.priority)
    
    out <- as.data.frame(ode(y = state.init, times = t.dyn$t.sys,
                             func = consumerResourceModel, parms = parameters, 
                             events = list(func = perturb, time = t.dyn$t.sys)))
    
    
    species.index <- grepl('consumer', names(out))
    resource.index <- grepl('resource', names(out))
    
    out.species.matrix <- as.matrix(out[,species.index])
    out.species.matrix <- as.matrix(out.species.matrix[t.dyn$t.index,])
    
    out.resource.matrix <- as.matrix(out[,resource.index])
    out.resource.matrix <- as.matrix(out.resource.matrix[t.dyn$t.index,])
    
    if(error.variance > 0){
        measurement.error <- rgamma(length(out.species.matrix), 1/error.variance, 1/error.variance)
        out.species.matrix <- out.species.matrix * matrix(measurement.error, ncol = n.species)
    }
    
    if(norm){
        out.species.matrix <- out.species.matrix/rowSums(out.species.matrix)
        out.resource.matrix <- out.resource.matrix/rowSums(out.resource.matrix)
    }
    
    colnames(out.species.matrix) <- names.species
    colnames(out.resource.matrix) <- names.resources
    
    out.species.matrix <- cbind(out.species.matrix, 
        time = t.dyn$t.sys[t.dyn$t.index])
    out.resource.matrix <- cbind(out.resource.matrix, 
        time = t.dyn$t.sys[t.dyn$t.index])
    
    out.list <- list(matrix = out.species.matrix, 
        resources = out.resource.matrix,
        E.Mat = E,
        monod.constant = monod.constant,
        error.variance = error.variance)
    return(out.list)
}
