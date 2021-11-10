#' Stochastic Logistic simulation
#'
#' Simulates time series with the (stochastic) logistic model and forms a list 
#' as a result.
#'
#' Simulates a community time series using the logistic model.
#' The change rate of the species was defined as
#' `dx/dt = growth.rates*x0*(1-(x0/carrying.capacities))*rN - death.rates*x0`, 
#' where `growth.rates` is the growth rates of each species,
#' `x0` is the initial abundances of each species,
#' `carrying.capacities` is the vector of maximum carrying capacities,
#' `rN` is a random number ranged from 0 to 1 which changes in each time step,
#' `death.rates` is the vector of constant death rates.
#' Also, the vectors of initial dead species abundances can be set.
#' The number of species will be set to 0 if the dead species abundances
#' surpass the alive species abundances.
#'
#' @param n.species Integer: number of species
#' @param growth.rates Numeric: growth rates of simulated species. If NULL, the 
#' default, `runif(n = n.species, min = 0.1, max = 0.2)` is used.
#' (default: \code{growth.rates = NULL})
#' @param carrying.capacities Numeric: The max population of species supported 
#' in the community. If NULL, by default, carrying.capacities is
#' `runif(n = n.species, min = 1000, max = 2000)`.
#' (default: \code{carrying.capacities = NULL})
#' @param death.rates Numeric: death rates of each species. If NULL, the 
#' default death.rates is `runif(n = n.species, min = 0.0005, max = 0.0025)`.
#' (default: \code{death.rates = NULL})
#' @param x0 Numeric: initial abundances of simulated species. If NULL, the 
#' default, then x0 is set to `runif(n = n.species, min = 0.1, max = 10)`.
#' (default: \code{x0 = NULL})
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
#' @param epoch.p Numeric: the probability/frequency of random periodic (epoch) 
#' changes introduced to the community composition
#' (default: \code{epoch.p = 0.001})
#' @param t.external_events Numeric: the starting time points of defined 
#' external events that introduce random changes to the community composition
#' (default: \code{t.external_events = c(0, 240, 480)})
#' @param t.external_durations Numeric: the durations of the defined external 
#' events.
#' (default: \code{t.external_durations = c(0, 1, 1)})
#' @param migration.p Numeric: the probability/frequency of migration from a 
#' metacommunity.
#' (default: \code{migration.p = 0.01})
#' @param metacommunity.probability Numeric: Normalized probability distribution
#' of the likelihood that species from the metacommunity can enter the community
#' during the simulation. If NULL as default, then 
#' `rdirichlet(1, alpha = rep(1,n.species))` is used.
#' (default: \code{metacommunity.probability = NULL})
#' @param stochastic Logical: whether to introduce noise in the simulation
#' (default: \code{stochastic = TRUE})
#' @param t.end Numeric: the end time of the simulationTimes, defining the 
#' modeled time length of the community. 
#' (default: \code{t.end = 1000})
#' @param ... additional parameters, see \code{\link{utils}} to know more.
#'
#' @docType methods
#' @examples
#' ## ATTENTION: Don't set a large value to t.step, otherwise the computer won't
#' #give a correct solution to the logistic ODE(ordinary differential equation).
#' #Keeping t.step under 0.05 or 0.01 is a good practice.
#'
#' # Example of logistic model without stochasticity or external disturbances
#' ExampleLogistic <- simulateStochasticLogistic(n.species = 5, stochastic = FALSE)
#' makePlot(ExampleLogistic$matrix)
#' 
#' # Example of stochastic logistic model
#' ExampleLogistic <- simulateStochasticLogistic(n.species = 5)
#' makePlot(ExampleLogistic$matrix)
#' 
#' # example with given initial parameters
#' ExampleLogistic <- simulateStochasticLogistic(n.species = 2,
#'     growth.rates = c(0.2, 0.1), 
#'     carrying.capacities = c(1000, 2000), 
#'     death.rates = c(0.001, 0.0015), 
#'     x0 = c(3, 0.1),
#'     sigma.drift = 0.001, 
#'     sigma.epoch = 0.3, 
#'     sigma.external = 0.5,
#'     sigma.migration = 0.002,
#'     epoch.p = 0.001,
#'     t.external_events = c(100, 200, 300), 
#'     t.external_durations = c(0.1, 0.2, 0.3),
#'     migration.p = 0.01,
#'     metacommunity.probability = rdirichlet(1, alpha = rep(1, 2)),
#'     stochastic = TRUE,
#'     norm = TRUE, 
#'     t.end = 400, 
#'     t.start = 0, t.step = 0.01,
#'     t.store = 1500)
#' makePlot(ExampleLogistic$matrix)
#' 
#' @return \code{simulateStochasticLogistic} returns a list of community dynamic
#' matrix(species abundance as rows and time points as columns) and its 
#' inputs(including metacommunity.probability and migration.p)
#'
#' @docType methods
#' @aliases simulateStochasticLogistic-numeric
#' @aliases simulateStochasticLogistic,numeric-method
#'
#' @importFrom deSolve ode
#'
#' @export
simulateStochasticLogistic <- function(n.species,
    growth.rates = NULL,
    carrying.capacities = NULL,
    death.rates = NULL,
    x0 = NULL,
    sigma.drift = 0.001,
    sigma.epoch = 0.1,
    sigma.external = 0.3,
    sigma.migration = 0.01,
    epoch.p = 0.001,
    t.external_events = c(0, 240, 480), 
    t.external_durations = c(0, 1, 1),  
    migration.p = 0.01,
    metacommunity.probability = NULL,
    stochastic = TRUE,
    norm = FALSE,
    t.end = 1000,...){

    # define the stochastic logistic model
    stochasticLogisticModel <- function (t, state, parameters){
        current <- pmax(0,state[names(state) == 'current'])
        live <- state[names(state) == 'live']
        dead <- state[names(state) == 'dead']
        growth.rates <- parameters$growth.rates
        carrying.capacities <- parameters$carrying.capacities
        death.rates <- parameters$death.rates

        dlive <- growth.rates*live*(1-(live/(carrying.capacities)))
        ddead <- death.rates*current
        dcurrent <- (dlive-ddead)
        dxdt <- list(c(dcurrent, dlive, ddead))
        return(dxdt)
    }
    
    # set the default values
    if (is.null(growth.rates)) {
        growth.rates <- runif(n = n.species, min = 0.1, max = 0.2)
    }
    if (is.null(carrying.capacities)) {
        carrying.capacities <- runif(n = n.species, min = 1000, max = 2000)
    }
    if (is.null(death.rates)) {
        death.rates <- runif(n = n.species, min = 0.0005, max = 0.0025)
    }
    if (is.null(x0)) {
        x0 <- runif(n = n.species, min = 0.1, max = 10)
    }
    if (is.null(metacommunity.probability)){
        metacommunity.probability <- rdirichlet(1, alpha = rep(1,n.species))
    }

    # select the time points to simulate and to store
    t.dyn <- SimulationTimes(t.end = t.end, ...)

    #continuous or episodic perturbation
    perturb <- function(t, y, parms){with(as.list(y), {
        epoch.rN <- 0
        external.rN <- 0
        migration.rN <- 0
        if (rbinom(1,1, parms$epoch.p)){
            epoch.rN <- rnorm(parms$n.species, sd=parms$sigma.epoch)
            epoch.rN <- parameters$stochastic*epoch.rN
        }
        if (rbinom(1,1, parameters$migration.p)){
            migration.rN <- rmultinom(n = 1, size = 1, prob = parameters$metacommunity.probability)[,]*abs(rnorm(n=1, mean=0, sd = parameters$sigma.migration))
            
            
        }
        if (t %in% parms$tEvent){
            external.rN <- rnorm(parms$n.species, sd=parms$sigma.external)
            external.rN <- parameters$stochastic*external.rN
        }
        drift.rN <- rnorm(parms$n.species, sd=parms$sigma.drift)
        drift.rN <- parameters$stochastic*drift.rN

        #perturbation is applied to the current population
        current <- pmax(y[names(y)=='current'], 0)
        current <- current*(1+drift.rN)*(1+epoch.rN)*(1+external.rN)+ migration.rN
        live <- y[names(y)=='live']
        dead <- y[names(y)=='dead']
        return(c(current, live, dead))})
    }

    tEvent = eventTimes(t.events = t.external_events,
        t.duration = t.external_durations,
        t.end = t.end, ...)

    parameters <- list(growth.rates=growth.rates, carrying.capacities=carrying.capacities, death.rates=death.rates, n.species = n.species,
        sigma.drift = sigma.drift, stochastic = stochastic,
        sigma.epoch = sigma.epoch, epoch.p = epoch.p,
        sigma.external = sigma.external, tEvent = tEvent,
        migration.p = migration.p, metacommunity.probability = metacommunity.probability,
        sigma.migration = sigma.migration)
    yinit <- c(x0, x0, numeric(n.species))
    names(yinit) <- rep(c("current", "live", "dead"), each = n.species)

    out <- as.data.frame(ode(func=stochasticLogisticModel,
        y=yinit, times=t.dyn$t.sys, parms=parameters,
        events = list(func = perturb, time = t.dyn$t.sys)))

    out.matrix <- as.matrix(out[,names(out)=='current'])
    out.matrix <- out.matrix[t.dyn$t.index,]
    
    if(norm){
        out.matrix <- out.matrix/rowSums(out.matrix)
    }
    
    # TODO: if is.null(names), do this
    colnames(out.matrix) <- seq_len(n.species)
    
    out.matrix <- cbind(out.matrix, time = t.dyn$t.sys[t.dyn$t.index])
    
    out.list <- list(matrix = out.matrix, 
                     metacommunity.probability = metacommunity.probability,
                     migration.p = migration.p)
    
    return(out.list)
}
