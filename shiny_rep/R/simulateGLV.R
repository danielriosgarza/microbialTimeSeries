#' Generalized Lotka-Volterra (gLV) simulation
#'
#' Simulates time series with the generalized Lotka-Volterra model
#'
#' @param n.species Integer: number of species
#' @param names.species Character: names of species. If NULL,
#' `paste0("sp", seq_len(n.species))` is used.
#' (default: \code{names.species = NULL})
#' @param A matrix: interaction matrix defining the positive and negative 
#' interactions between n.species. If NULL, `randomA(n.species)` is used.
#' (default: \code{A = NULL})
#' @param x0 Numeric: initial abundances of simulated species. If NULL, 
#' `runif(n = n.species, min = 0, max = 1)` is used.
#' (default: \code{x0 = NULL})
#' @param growth.rates Numeric: growth rates of simulated species. If NULL,
#' `runif(n = n.species, min = 0, max = 1)` is used.
#' (default: \code{growth.rates = NULL})
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
#' (default: \code{stochastic = TRUE})
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
#' @param ... additional parameters, see \code{\link{utils}} to know more.
#' 
#' @examples
#' 
#' # generate a random interaction matrix
#' ExampleA <- randomA(n.species = 4, diagonal = -1)
#' 
#' # run the model with default values (only stochastic migration considered)
#' ExampleGLV <- simulateGLV(n.species = 4, A = ExampleA)
#' # visualize the result
#' makePlot(ExampleGLV$matrix)
#' 
#' # run the model with two external disturbances at time points 240 and 480 
#' # with durations equal to 1 (10 time steps when t.step by default is 0.1).
#' ExampleGLV <- simulateGLV(n.species = 4, A = ExampleA,
#'     t.external_events = c(0, 240, 480), t.external_durations = c(0, 1, 1))
#' # visualize the result
#' makePlot(ExampleGLV$matrix)
#' 
#' # run the model with no pertubation nor migration
#' set.seed(42)
#' ExampleGLV <- simulateGLV(n.species = 4, A = ExampleA, stochastic = FALSE, 
#'     sigma.migration = 0)
#' # visualize the result
#' makePlot(ExampleGLV$matrix)
#' 
#' # run the model with no pertubation nor migration but with measurement error
#' set.seed(42)
#' ExampleGLV <- simulateGLV(n.species = 4, A = ExampleA, stochastic = FALSE, 
#'     error.variance = 0.001, sigma.migration = 0)
#' # visualize the result
#' makePlot(ExampleGLV$matrix)
#'
#' @docType methods
#' @aliases simulateGLV-numeric
#' @aliases simulateGLV,numeric-method
#'
#' @importFrom deSolve ode
#'
#' @export

# define the GLV Model
glvModel <- function(t, x0, parameters){
    x0[x0 < 10^-8] <- 0 
    growth.rates <- parameters$growth.rates
    A <- parameters$A
    # rate of change
    dx <- x0*(growth.rates+A %*% x0)
    # return rate of change
    list(dx)
}

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
            # TODO: is migration also stochastic? if so, add the following:####
            # migration.rN <- parameters$stochastic*migration.rN
        }
        
        if (t %in% parameters$tEvent){
            external.rN <- rnorm(parameters$n.species,
                sd=parameters$sigma.external)
            external.rN <- parameters$stochastic*external.rN
        }
        drift.rN <- rnorm(parameters$n.species, sd=parameters$sigma.drift)
        drift.rN <- parameters$stochastic*drift.rN

        #perturbation is applied to the current population
        y <- y * (1 + drift.rN)*(1 + epoch.rN)*(1 + external.rN) + migration.rN 
        return(y*(y>0))
    })
}

simulateGLV <- function(n.species, 
    names.species = NULL,
    A = NULL, 
    x0 = NULL,
    growth.rates = NULL,
    sigma.drift = 0.001,
    sigma.epoch = 0.1,
    sigma.external = 0.3,
    sigma.migration = 0.01,
    epoch.p = 0.001,
    t.external_events = NULL,
    t.external_durations = NULL,
    stochastic = TRUE,
    migration.p = 0.01, 
    metacommunity.probability = NULL,
    error.variance = 0,
    norm = FALSE, 
    t.end = 1000, ...){
    
    # set the default values
    if (is.null(names.species)) {
        names.species <- paste0("sp", seq_len(n.species))
    }
    if (is.null(A)) {
        A <- randomA(n.species)
    }
    if (is.null(x0)) {
        x0 <- runif(n = n.species, min = 0, max = 1)
    }
    if (is.null(growth.rates)) {
        growth.rates <- runif(n = n.species, min = 0, max = 1)
    }
    if (is.null(metacommunity.probability)) {
        metacommunity.probability <- rdirichlet(1, alpha = rep(1,n.species))
    }
    
    #normalize metacommunity.probability
    metacommunity.probability <- metacommunity.probability/
        sum(metacommunity.probability)
    
    # select the time points to simulate and to store
    t.dyn <- SimulationTimes(t.end = t.end, ...)
    
    # calculate the time points influenced by the disturbances
    tEvent = eventTimes(
        t.events = t.external_events,
        t.duration = t.external_durations,
        t.end = t.end,
        ... = ...)
    
    parameters <- list(growth.rates=growth.rates, A = A, n.species = n.species,
        sigma.drift = sigma.drift, stochastic = stochastic,
        sigma.epoch = sigma.epoch, epoch.p = epoch.p,
        sigma.external = sigma.external, tEvent = tEvent,
        migration.p = migration.p, 
        metacommunity.probability = metacommunity.probability, 
        sigma.migration= sigma.migration)
    
    out <- ode(
        y = x0,
        times = t.dyn$t.sys,
        func = glvModel,
        parms = parameters,
        events = list(func = perturb, time = t.dyn$t.sys),
        maxsteps=10^9, method = "ode45")
    
    
    out.matrix <- out[,2:ncol(out)]
    
    out.matrix <- out.matrix[t.dyn$t.index,]
    
    if(error.variance > 0){
        measurement.error <- rnorm(n = length(t.dyn$t.index)*n.species, 
                                   mean = 0, sd = sqrt(error.variance))
        measurement.error <- matrix(measurement.error, 
                                    nrow = length(t.dyn$t.index))
        out.matrix <- out.matrix + measurement.error
    }
    
    if(norm){
        out.matrix <- out.matrix/rowSums(out.matrix)
    }
    
    colnames(out.matrix) <- names.species
    
    out.matrix <- cbind(out.matrix, time = t.dyn$t.sys[t.dyn$t.index])
    
    out.list <- list(matrix = out.matrix, 
        metacommunity.probability = metacommunity.probability,
        migration.p = migration.p,
        error.variance = error.variance)
    
    return(out.list)
}
