#' Hubbell's neutral model simulation
#'
#' Neutral species abundances simulation according to the Hubbell model.
#'
#' @param x0 Numeric: initial species composition
#' @param names.species Character: names of species. If NULL,
#' `paste0("sp", seq_len(n.species))` is used.
#' (default: \code{names.species = NULL})
#' @param migration.p Numeric: the probability/frequency of migration from a 
#' metacommunity
#' (default: \code{m = 0.02})
#' @param metacommunity.probability Numeric: normalized probability distribution
#' of the likelihood that species from the metacommunity can enter the community
#' during the simulation. If NULL, `rdirichlet(1, alpha = rep(1,n.species))` is 
#' used.
#' (default: \code{metacommunity.probability = NULL})
#' @param k.events Integer: number of events to simulate before updating the 
#' sampling distributions.
#' (default: \code{k.events = 1})
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
#' @examples
#' set.seed(42)
#' ExampleHubbell <- simulateHubbell(rep(1000, 10))
#' makePlot(ExampleHubbell$matrix)
#' 
#' set.seed(42)
#' ExampleHubbell <- simulateHubbell(rep(1000, 10), error.variance = 10)
#' makePlot(ExampleHubbell$matrix)
#' 
#' set.seed(42)
#' ExampleHubbell <- simulateHubbell(rep(1000, 10), migration.p = 0.1)
#' makePlot(ExampleHubbell$matrix)
#' 
#' set.seed(42)
#' ExampleHubbell <- simulateHubbell(rep(1000, 10), migration.p = 0.1,
#'     error.variance = 10)
#' makePlot(ExampleHubbell$matrix)
#' 
#' set.seed(42)
#' ExampleHubbell <- simulateHubbell(rep(1000, 5), 
#'     metacommunity.probability = c(0.01, 0.1, 0.19, 0.3, 0.4))
#' makePlot(ExampleHubbell$matrix)
#' 
#' set.seed(42)
#' ExampleHubbell <- simulateHubbell(rep(1000, 5), 
#'     metacommunity.probability = c(0.01, 0.1, 0.19, 0.3, 0.4),
#'     k.events = 5)
#' makePlot(ExampleHubbell$matrix)
#' 
#' @return \code{simulateHubbell} returns a list of initial states, parameters
#' of the model, including a matrix with species abundance as rows and time 
#' points as columns.
#' @docType methods
#' @aliases simulateHubbell-numeric
#' @aliases simulateHubbell,numeric-method
#' @aliases simulateNeutral
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
simulateHubbell <- function(x0, 
    names.species = NULL,
    migration.p = 0.01, 
    metacommunity.probability = NULL,
    k.events = 1,
    error.variance = 0,
    norm = FALSE, 
    t.end=1000, ...){
    
    # set the default values
    n.species <- length(x0)
    if (is.null(names.species)) {
        names.species <- paste0("sp", seq_len(n.species))
    }
    if (is.null(metacommunity.probability)){
        metacommunity.probability <- rdirichlet(1, alpha = rep(1,n.species))
    }
    # normalize metacommunity.probability
    metacommunity.probability <- metacommunity.probability/
        sum(metacommunity.probability)
    
    t.dyn <- SimulationTimes(t.end = t.end,...)
    birth.p <- 1 - migration.p
    community <- x0
    
    out.matrix <- matrix(0, nrow=length(t.dyn$t.index), ncol = n.species)
    
    counter = 1
    
    out.matrix[counter,] <- community
    
    for (i in t.dyn$t.sys[2:length(t.dyn$t.sys)]){
        probabilities <- community/sum(community)
        n.deaths <- min(min(community[community>0]),k.events)
        
        # deaths
        community = community - t(rmultinom(n = 1, size = n.deaths, 
            prob = probabilities))
        n.births = sum(rbinom(n=n.deaths, size=1, p = birth.p))
        n.migration = n.deaths-n.births
        
        community = community + t(rmultinom(n = 1, size = n.births, 
            prob = probabilities)) +
            t(rmultinom(n = 1, size = n.migration, 
                prob = metacommunity.probability))
        
        if (i %in% t.dyn$t.sys[t.dyn$t.index]){
            counter = counter + 1
            out.matrix[counter,] <- community
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
        x0 = x0,
        metacommunity.probability = metacommunity.probability,
        migration.p = migration.p, 
        birth.p = birth.p,
        error.variance = error.variance)
    return(out.list)
}
