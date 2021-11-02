#' Stochastic Logistic simulation
#'
#' Simulates time series with the (stochastic) logistic model and
#' forms a \linkS4class{SummarizedExperiment} object or a matrix.
#'
#' Simulates a community time series using the logistic model.
#' The change rate of the species was defined as
#' `dx/dt = growth.rates*x0*(1-(x0/carrying.k))*rN - death.rates*x0`, where
#' growth.rates is the vector of growth rates,
#' x0 is the vector of initial species abundances,
#' carrying.k is the vector of maximum carrying capacities,
#' rN is a random number ranged from 0 to 1 which changes in each time step,
#' death.rates is the vector of constant death rates.
#' Also, the vectors of initial dead species abundances can be set.
#' The number of species will be set to 0 if the dead species abundances
#' surpass the alive species abundances.
#'
#' The resulting abundance matrix model is used to construct
#' \linkS4class{SummarizedExperiment} object.
#'
#' @param n.species Integer: number of species
#' @param growth.rates Numeric: growth rates
#' (default: \code{growth.rates = runif(n = n.species, min = 0.1, max = 0.2)})
#' @param carrying.k Numeric: carrying capacities
#' (default: \code{carrying.k = runif(n = n.species, min = 1000, max = 2000)})
#' @param death.rates Numeric: death rates
#' (default: \code{death.rates = runif(n = n.species, min = 0.0005, max = 0.0025)})
#' @param x0 Numeric: initial abundances
#' (default: \code{x0 = runif(n = n.species, min = 0.1, max = 10)})
#' @param sigma.drift Numeric: degree of drift (turnover of species) in each
#' time step.
#' (default: \code{sigma.drift = 0.001})
#' @param sigma.epoch Numeric: degree of epoch change of community
#' (default: \code{sigma.epoch = 0.1})
#' @param sigma.external Numeric: degree of external events/disturbances
#' (default: \code{sigma.external = 0.3})
#' @param p.epoch Numeric: probability/frequency of inherit periodic changes of
#' community
#' (default: \code{p.epoch = 0.001})
#' @param t.external_events Numeric: starting times of external events
#' (default: \code{t.external_events = c(0, 240, 480)})
#' @param t.external_durations Numeric: durations of external events
#' (default: \code{t.external_durations = c(0, 1, 1)})
#' @param t.end Numeric scalar indicating the final time of the dimulation
#' (default: \code{t.end = 2000})
#' @param return.matrix Logical: whether to export only the stored time points
#' in a matrix
#' (default: \code{return.matrix = FALSE})
#' @param stochastic Logical: whether the logistic model should be stochastic
#' (controlled by multiplying the growth rate by a random number)
#' (default: \code{stochastic = TRUE})
#' @param ... additional parameters including 't.start', 't.step', and 't.store'
#' see \code{\link{SimulationTimes}} for more information
#'
#' @docType methods
#' @examples
#' ## ATTENTION: Don't set a large value to t.step, otherwise the computer won't
#' #give a correct solution to the logistic ODE(ordinary differential equation).
#' #Keeping t.step under 0.05 or 0.01 is a good practice.
#'
#' #while (!exists("ExampleLogistic"))
#' ExampleLogistic <- simulateStochasticLogistic(n.species = 5)
#' #plot the calculated points
#' matplot(t(assays(ExampleLogistic)[["counts"]]), type = "l")
#'
#' #calculation by setting initial parameters explicitly
#' ExampleLogistic <- simulateStochasticLogistic(
#' n.species = 2,
#' growth.rates = c(0.2, 0.1), carrying.k = c(1000, 2000), death.rates = c(0.001, 0.0015), x0 = c(3, 0.1),
#' sigma.drift = 0.001, sigma.epoch = 0.3, sigma.external = 0.5,p.epoch = 0.001,
#' t.external_events = c(100, 200, 300), t.external_durations = c(1, 2, 3),
#' t.start = 0, t.end = 1500, t.step = 0.01,
#' t.store = 1500, return.matrix = FALSE, stochastic = TRUE)
#'
#' @return \code{simulateStochasticLogistic} returns a
#' \linkS4class{SummarizedExperiment} class object containing matrix with
#' species abundance as rows and time points as columns
#'
#' @docType methods
#' @aliases simulateStochasticLogistic-numeric
#' @aliases simulateStochasticLogistic,numeric-method
#'
#' @importFrom deSolve ode
#' @importFrom stats runif
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom methods setGeneric
#'
#' @export

setGeneric("simulateStochasticLogistic",signature = "n.species",
    function(n.species,
        growth.rates = runif(n = n.species, min = 0.1, max = 0.2),
        carrying.k = runif(n = n.species, min = 1000, max = 2000),
        death.rates = runif(n = n.species, min = 0.0005, max = 0.0025),
        x0 = runif(n = n.species, min = 0.1, max = 10),
        sigma.drift = 0.001,
        sigma.epoch = 0.1,
        sigma.external = 0.3,
        sigma.migration = 0.01,
        p.epoch = 0.001,
        t.external_events = c(0, 240, 480),
        t.external_durations = c(0, 1, 1),
        migration.p = 0.01, 
        metacommunity.p = NULL,
        stochastic = TRUE,
        norm = FALSE,
        t.end = 1000,...)
        standardGeneric("simulateStochasticLogistic"))

setMethod("simulateStochasticLogistic", signature = c(n.species="numeric"),
    function(n.species,
             growth.rates = runif(n = n.species, min = 0.1, max = 0.2),
             carrying.k = runif(n = n.species, min = 1000, max = 2000),
             death.rates = runif(n = n.species, min = 0.0005, max = 0.0025),
             x0 = runif(n = n.species, min = 0.1, max = 10),
             sigma.drift = 0.001,
             sigma.epoch = 0.1,
             sigma.external = 0.3,
             sigma.migration = 0.01,
             p.epoch = 0.001,
             t.external_events = c(0, 240, 480),
             t.external_durations = c(0, 1, 1),
             migration.p = 0.01, 
             metacommunity.p = NULL,
             stochastic = TRUE,
             norm = FALSE,
             t.end = 1000,...){

        # define the stochastic logistic model
        stochasticLogisticModel <- function (t, state, parameters){
            current <- pmax(0,state[names(state) == 'current'])
            live <- state[names(state) == 'live']
            dead <- state[names(state) == 'dead']
            growth.rates <- parameters$growth.rates
            carrying.k <- parameters$carrying.k
            death.rates <- parameters$death.rates

            dlive <- growth.rates*live*(1-(live/(carrying.k)))
            ddead <- death.rates*current
            dcurrent <- (dlive-ddead)
            dxdt <- list(c(dcurrent, dlive, ddead))
            return(dxdt)
        }

        if (is.null(metacommunity.p)){
            metacommunity.p <- rdirichlet(1, alpha = rep(1,n.species))
        }

        # select the time points to simulate and to store
        t.dyn <- SimulationTimes(t.end = t.end,...)

        #continuous or episodic perturbation
        perturb <- function(t, y, parms){with(as.list(y), {
            epoch.rN <- 0
            external.rN <- 0
            migration.rN <- 0
            if (rbinom(1,1, parms$p.epoch)){
                epoch.rN <- rnorm(parms$n.species, sd=parms$sigma.epoch)
                epoch.rN <- parameters$stochastic*epoch.rN
            }
            if (rbinom(1,1, parameters$migration.p)){
                migration.rN <- rmultinom(n = 1, size = 1, prob = parameters$metacommunity.p)[,]*abs(rnorm(n=1, mean=0, sd = parameters$sigma.migration))
                
                
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

        parameters <- list(growth.rates=growth.rates, carrying.k=carrying.k, death.rates=death.rates, n.species = n.species,
            sigma.drift = sigma.drift, stochastic = stochastic,
            sigma.epoch = sigma.epoch, p.epoch = p.epoch,
            sigma.external = sigma.external, tEvent = tEvent,
            migration.p = migration.p, metacommunity.p = metacommunity.p,
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
        
        
        colnames(out.matrix) <- seq_len(n.species)
        
        
        
        out.matrix <- cbind(out.matrix, time = t.dyn$t.sys[t.dyn$t.index])
        
        
       
        out.list <- list(matrix = out.matrix, 
                         metacommunity.p = metacommunity.p,
                         migration.p = migration.p)
        
        return(out.list)
    }
)
