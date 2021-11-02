#' Generalized Lotka-Volterra (gLV) simulation
#'
#' Simulates time series with the generalized Lotka-Volterra model and forms a
#' \linkS4class{SummarizedExperiment} object.
#'
#' Simulates a community time series using the generalized Lotka-Volterra model,
#' defined as dx/dt = diag(x0)(growth.rates+Ax), where x0 is the vector of species abundances
#' ,diag(x0) is a diagonal matrix with the diagonal values set to x0.
#' A is the interaction matrix and growth.rates is the vector of growth rates.
#'
#' The resulting abundance matrix model is used to construct
#'\linkS4class{SummarizedExperiment} object.
#'
#' @param n.species Integer: number of species
#' @param A interaction matrix
#' @param x0 Numeric: initial abundances
#' @param growth.rates Numeric: growth rates
#' @param t.start Numeric scalar indicating the initial time of the simulation.
#' (default: \code{t.start = 0})
#' @param t.end Numeric scalar indicating the final time of the dimulation
#' (default: \code{t.end = 1000})
#' @param t.store Integer scalar indicating the number of evenly distributed
#' time points to keep
#' @param t.step Numeric scalar indicating the interval between simulation steps
#' (default: \code{t.step = 0.01})
#' @param sigma.drift Numeric: degree of drift (turnover of species) in each
#' time step.
#' (default: \code{sigma.drift = 0.01})
#' @param sigma.epoch Numeric: degree of epoch change of community
#' (default: \code{sigma.epoch = 0.3})
#' @param sigma.external Numeric: degree of external events/disturbances
#' (default: \code{sigma.external = 0.3})
#' @param p.epoch Numeric: probability/frequency of inherit periodic changes of
#' community
#' (default: \code{p.epoch = 0.01})
#' @param t.external_events Numeric: starting times of external events
#' (default: \code{t.external_events = c(12, 36, 48)})
#' @param t.external_durations Numeric: durations of external events
#' (default: \code{t.external_durations = c(3, 10, 99)})
#' @param stochastic Logical: whether the gLV model should be stochastic
#' (default: \code{stochastic = FALSE})
#' @param norm Logical scalar: returns normalised abundances (proportions
#' in each generation) (default: \code{norm = FALSE})
#' @param t.start Numeric scalar indicating the initial time of the simulation.
#' (default: \code{t.start = 0})
#' @param t.store Integer scalar indicating the number of evenly distributed
#' time points to keep (default: \code{t.store = 100})
#' @param arguments Logical: decides whether the generated matrix parameters
#' are included in \linkS4class{SummarizedExperiment} object
#' @param ... additional arguments that can be called from miaSim::tDyn
#'
#' @return
#' \code{simulateGLV} returns a \linkS4class{SummarizedExperiment} object
#' containing abundance matrix
#'
#' @examples
#' row_data <- data.frame(Kingdom = "Animalia",
#'                 Phylum = rep(c("Chordata", "Echinodermata"), c(500, 500)),
#'                 Class = rep(c("Mammalia", "Asteroidea"), each = 500),
#'                 ASV = paste0("X", seq_len(1000)),
#'                 row.names = rownames(paste0("species", seq_len(1000))),
#'                 stringsAsFactors = FALSE)
#'
#' row_data <- t(row_data)
#'
#' col_data <- DataFrame(sampleID = seq_len(1000),
#'                     time = as.Date(1000, origin = "2000-01-01"),
#'                     row.names = colnames(paste0("sample", seq_len(1000))))
#'
#' A <- miaSim::powerlawA(4, alpha = 1.01)
#'
#' SEobject <- simulateGLV(n.species = 4, A, t.start = 0, t.store = 1000)
#' rowData(SEobject) <- row_data
#' colData(SEobject) <- col_data
#'
#' @docType methods
#' @aliases simulateGLV-numeric
#' @aliases simulateGLV,numeric-method
#'
#' @importFrom utils str
#' @importFrom deSolve ode
#' @importFrom stats runif
#' @importFrom SummarizedExperiment metadata<-
#' @importFrom methods setGeneric
#'
#' @export
setGeneric("simulateGLV", signature = "n.species",
    function(n.species, A,
             x0 = runif(n.species),
             growth.rates = runif(n.species),
             t.step = 0.1,
             sigma.drift = 0.01,
             sigma.epoch = 0.3,
             sigma.external = 0.3,
             sigma.migration = 0.01,
             p.epoch = 0.01,
             t.external_events = c(12, 36, 48),
             t.external_durations = c(3, 10, 99),
             stochastic = FALSE,
             migration.p = 0.01, 
             metacommunity.p = NULL,
             norm = FALSE,
             arguments = TRUE,
             t.end=1000, ...)
        standardGeneric("simulateGLV"))

dxdt <- function(t, x0, parameters){
    x0[x0 < 10^-8] <- 0 
    growth.rates <- parameters$growth.rates
    A <- parameters$A
    # rate of change
    dx <- x0*(growth.rates+A %*% x0)
    # return rate of change
    list(dx)
}

perturb <- function(t, y, parameters){
    with(as.list(y),{
        #continuous or episodic perturbation
        epoch.rN <- 0
        external.rN <- 0
        migration.rN <- 0
        if (rbinom(1,1, parameters$p.epoch)){
            epoch.rN <- rnorm(parameters$n.species, sd=parameters$sigma.epoch)
            epoch.rN <- parameters$stochastic*epoch.rN
        }
        
        if (rbinom(1,1, parameters$migration.p)){
            migration.rN <- rmultinom(n = 1, size = 1, prob = parameters$metacommunity.p)[,]*abs(rnorm(n=1, mean=0, sd = parameters$sigma.migration))
            
            
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

setMethod("simulateGLV", signature = c(n.species="numeric"),
    function(n.species, A,
        x0 = runif(n.species),
        growth.rates = runif(n.species),
        sigma.drift = 0.01,
        sigma.epoch = 0.3,
        sigma.external = 0.3,
        sigma.migration = 0.01,
        p.epoch = 0.01,
        t.external_events = c(12, 36, 48),
        t.external_durations = c(3, 10, 99),
        stochastic = FALSE,
        migration.p = 0.01, 
        metacommunity.p = NULL,
        norm = FALSE,
        arguments = TRUE, ...){
        
        t.dyn <- SimulationTimes(t.end = t.end, ...)
        
        tEvent = eventTimes(
            t.events = t.external_events,
            t.duration = t.external_durations,
            t.start = t.start, t.end = t.end,
            t.step = t.step, t.store = t.store)
        
        if (is.null(metacommunity.p)){
            metacommunity.p <- rdirichlet(1, alpha = rep(1,n.species))[,]
        }
        #normalize metacommunity.p
        metacommunity.p <- metacommunity.p/sum(metacommunity.p)
        
        parameters <- list(growth.rates=growth.rates, A = A, n.species = n.species,
            sigma.drift = sigma.drift, stochastic = stochastic,
            sigma.epoch = sigma.epoch, p.epoch = p.epoch,
            sigma.external = sigma.external, tEvent = tEvent,
            migration.p = migration.p, metacommunity.p = metacommunity.p, 
            sigma.migration= sigma.migration)
        
        out <- ode(
            y = x0,
            times = t.dyn$t.sys,
            func = dxdt,
            parms = parameters,
            events = list(func = perturb, time = t.dyn$t.sys),
            maxsteps=10^9, method = "ode45")
        
        
        out.matrix <- out[,2:ncol(out)]
        
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
