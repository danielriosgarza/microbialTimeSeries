#' Hubbell's neutral model simulation
#'
#' Neutral species abundances simulation according to the Hubbell model.
#'
#' @param community.initial Numeric: initial species composition
#' @param names.species Character: names of species. If NULL,
#' `paste0("sp", seq_len(n.species))` is used.
#' (default: \code{names.species = NULL})
#' @param migration.p Numeric: the probability/frequency of migration from a 
#' metacommunity
#' (default: \code{migration.p = 0.01})
#' @param metacommunity.probability Numeric: normalized probability distribution
#' of the likelihood that species from the metacommunity can enter the community
#' during the simulation. If NULL, `rdirichlet(1, alpha = rep(1,n.species))` is 
#' used.
#' (default: \code{metacommunity.probability = NULL})
#' @param k.events Integer: number of events to simulate before updating the 
#' sampling distributions.
#' (default: \code{k.events = 1})
#' @param growth.rates Numeric: maximum growth rates(mu) of species.
#' If NULL, `rep(1, n.species)` is used.
#' (default: \code{growth.rates = NULL})
#' @param error.variance Numeric: the variance of measurement error.
#' By default it equals to 0, indicating that the result won't contain any 
#' measurement error. This value should be non-negative.
#' (default: \code{error.variance = 0})
#' @param norm Logical: whether the time series should be returned with
#' the abundances as proportions (\code{norm = TRUE}) or
#' the raw counts (default: \code{norm = FALSE})
#' @param t.end Numeric: simulation end time (default: \code{t.end = 1000})
#' @param ... additional parameters including 't.start', 't.step', and 't.store'
#' see \code{\link{utils}} for more information.
#'
#' @docType methods
#' @aliases simulateHubbellRates-numeric
#' @aliases simulateHubbellRates,numeric-method
#' @aliases simulateNeutralRates
#'
#' @examples
#' set.seed(42)
#' ExampleHubbellRates <- simulateHubbellRates(rep(100, 5))
#' makePlot(ExampleHubbellRates$matrix)
#' 
#' # no migration, all stochastic birth and death
#' set.seed(42)
#' ExampleHubbellRates <- simulateHubbellRates(rep(100, 5), migration.p = 0)
#' makePlot(ExampleHubbellRates$matrix)
#' 
#' # all migration, no stochastic birth and death
#' set.seed(42)
#' ExampleHubbellRates <- simulateHubbellRates(rep(100, 5), migration.p = 1, 
#'      t.end = 20, t.store = 200)
#' makePlot(ExampleHubbellRates$matrix)
#' 
#' # all migration, no stochastic birth and death, but with measurement errors
#' set.seed(42)
#' ExampleHubbellRates <- simulateHubbellRates(rep(100, 5), migration.p = 1, 
#'      t.end = 20, t.store =200, error.variance = 100)
#' makePlot(ExampleHubbellRates$matrix) 
#' 
#' # model with specified inputs
#' set.seed(42)
#' ExampleHubbellRates <- simulateHubbellRates(rep(100, 5), migration.p = 0.1,
#'     metacommunity.probability <- c(1,2,3,4,5), k.events = 5,
#'     growth.rates <- c(1,2,3,4,5))
#' makePlot(ExampleHubbellRates$matrix)
#' 
#' @return \code{simulateHubbellRates} returns a list of initial states, 
#' parameters of the model, including a matrix with species abundance as rows 
#' and time points as columns.
#' 
#' @docType methods
#' @aliases simulateHubbellRates-numeric
#' @aliases simulateHubbellRates,numeric-method
#'
#' @importFrom gtools rdirichlet
#' @importFrom stats rbinom
#' @importFrom stats rmultinom
#'
#' @references Rosindell, James et al. "The unified neutral theory of
#' biodiversity and biogeography at age ten." Trends in ecology & evolution
#' vol. 26,7 (2011).
#
#' @export
simulateHubbellRates <- function(community.initial, 
    names.species = NULL,
    migration.p = 0.01, 
    metacommunity.probability = NULL,
    k.events = 1, 
    growth.rates = NULL,
    error.variance = 0,
    norm = FALSE,
    t.end=1000,...){
    
    # set the default values
    n.species <- length(community.initial)
    if (is.null(names.species)) {
        names.species <- paste0("sp", seq_len(n.species))
    }
    
    if (is.null(metacommunity.probability)){
        metacommunity.probability <- rdirichlet(1, alpha = rep(1,n.species))
    }
    # normalize metacommunity.probability
    metacommunity.probability <- metacommunity.probability/
        sum(metacommunity.probability)
    
    if (is.null(growth.rates)){
        growth.rates <- rep(1,n.species)
    }
    
    t.dyn <- SimulationTimes(t.end = t.end,...)
    t.store <- length(t.dyn$t.index)
    
    birth.p <- 1 - migration.p
    community <- community.initial
    
    propensities <- sum(community)*(c(migration.p, 1-migration.p))
    event.probabilities <- propensities/(sum(propensities))
    
    out.matrix <- matrix(0, nrow=length(t.dyn$t.index), ncol = n.species)
    out.matrix[1,] = community.initial
    
    stored_time = t.dyn$t.sys[t.dyn$t.index]
    current_t <- stored_time[1]
    last_stored_t <- stored_time[1]
    
    while(current_t <= t.end){
        tau_events <- min(min(community[community>0]),k.events)
        tau <- rgamma(n = 1, shape = tau_events, scale = 1/(sum(propensities)))
        current_t <- current_t + tau
        
        composition.propensities <- community*growth.rates

        composition.probabilities <- composition.propensities/sum(composition.propensities)
        
        #k deaths
        community <- community -
            (rmultinom(n=1, size=tau_events, prob=composition.probabilities))
        
        n_births <- rbinom(n=1, size=tau_events, prob=event.probabilities[2])
        n_migration <- tau_events-n_births
        
        community <- community +
          (rmultinom(n=1, size=n_births, prob=composition.probabilities)) +
          (rmultinom(n=1, size=n_migration, prob=metacommunity.probability))
        
        index <- ((current_t >= stored_time )  & (last_stored_t < stored_time))
        
        if (sum(index)>0) {
            out.matrix[index,] <- t(matrix(rep(t(community),sum(index)), 
                ncol = sum(index)))
            last_stored_t <- stored_time[max(seq(t.store)[index])]
        }
    }
    
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
    
    #out.matrix$t <- t.dyn$t.sys[t.dyn$t.index]
    #SE <- SummarizedExperiment(assays = list(counts = out.matrix))
    out.list <- list(matrix = out.matrix, 
        community = community, 
        community.initial = community.initial,
        metacommunity.probability = metacommunity.probability,
        error.variance = error.variance)
    return(out.list)
}
