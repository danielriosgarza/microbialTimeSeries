#' Hubbell's neutral model simulation
#'
#' Neutral species abundances simulation according to the Hubbell model.
#'
#' @param x0 Numeric: initial species composition
#' @param names_species Character: names of species. If NULL,
#' `paste0("sp", seq_len(n_species))` is used.
#' (default: \code{names_species = NULL})
#' @param migration_p Numeric: the probability/frequency of migration from a 
#' metacommunity
#' (default: \code{m = 0.02})
#' @param metacommunity_probability Numeric: normalized probability distribution
#' of the likelihood that species from the metacommunity can enter the community
#' during the simulation. If NULL, `rdirichlet(1, alpha = rep(1,n_species))` is 
#' used.
#' (default: \code{metacommunity_probability = NULL})
#' @param k_events Integer: number of events to simulate before updating the 
#' sampling distributions.
#' (default: \code{k_events = 1})
#' @param error_variance Numeric: the variance of measurement error.
#' By default it equals to 0, indicating that the result won't contain any 
#' measurement error. This value should be non-negative.
#' (default: \code{error_variance = 0})
#' @param norm Logical: whether the time series should be returned with
#' the abundances as proportions (\code{norm = TRUE}) or
#' the raw counts (default: \code{norm = FALSE})
#' @param t_end Numeric: simulation end time (default: \code{t_end = 1000})
#' @param ... additional parameters including 't_start', 't_step', and 't_store'
#' see \code{\link{utils}} for more information.
#'
#' @examples
#' set.seed(42)
#' ExampleHubbell <- simulateHubbell(rep(1000, 10))
#' makePlot(ExampleHubbell$matrix)
#' 
#' set.seed(42)
#' ExampleHubbell <- simulateHubbell(rep(1000, 10), error_variance = 10)
#' makePlot(ExampleHubbell$matrix)
#' 
#' set.seed(42)
#' ExampleHubbell <- simulateHubbell(rep(1000, 10), migration_p = 0.1)
#' makePlot(ExampleHubbell$matrix)
#' 
#' set.seed(42)
#' ExampleHubbell <- simulateHubbell(rep(1000, 10), migration_p = 0.1,
#'     error_variance = 10)
#' makePlot(ExampleHubbell$matrix)
#' 
#' set.seed(42)
#' ExampleHubbell <- simulateHubbell(rep(1000, 5), 
#'     metacommunity_probability = c(0.01, 0.1, 0.19, 0.3, 0.4))
#' makePlot(ExampleHubbell$matrix)
#' 
#' set.seed(42)
#' ExampleHubbell <- simulateHubbell(rep(1000, 5), 
#'     metacommunity_probability = c(0.01, 0.1, 0.19, 0.3, 0.4),
#'     k_events = 5)
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
    names_species = NULL,
    migration_p = 0.01, 
    metacommunity_probability = NULL,
    k_events = 1,
    error_variance = 0,
    norm = FALSE, 
    t_end=1000, ...){
    
    # set the default values
    n_species <- length(x0)
    if (is.null(names_species)) {
        names_species <- paste0("sp", seq_len(n_species))
    }
    if (is.null(metacommunity_probability)){
        metacommunity_probability <- rdirichlet(1, alpha = rep(1,n_species))
    }
    # normalize metacommunity_probability
    metacommunity_probability <- metacommunity_probability/
        sum(metacommunity_probability)
    
    t_dyn <- simulationTimes(t_end = t_end,...)
    birth_p <- 1 - migration_p
    community <- x0
    
    out_matrix <- matrix(0, nrow=length(t_dyn$t_index), ncol = n_species)
    
    counter = 1
    
    out_matrix[counter,] <- community
    
    for (i in t_dyn$t_sys[2:length(t_dyn$t_sys)]){
        probabilities <- community/sum(community)
        n_deaths <- min(min(community[community>0]),k_events)
        
        # deaths
        community = community - t(rmultinom(n = 1, size = n_deaths, 
            prob = probabilities))
        n_births = sum(rbinom(n=n_deaths, size=1, p = birth_p))
        n_migration = n_deaths-n_births
        
        community = community + t(rmultinom(n = 1, size = n_births, 
            prob = probabilities)) +
            t(rmultinom(n = 1, size = n_migration, 
                prob = metacommunity_probability))
        
        if (i %in% t_dyn$t_sys[t_dyn$t_index]){
            counter = counter + 1
            out_matrix[counter,] <- community
        }
    
    }
    
    if(error_variance > 0){
        measurement_error <- rnorm(n = length(t_dyn$t_index)*n_species, 
                                   mean = 0, sd = sqrt(error_variance))
        measurement_error <- matrix(measurement_error, 
                                    nrow = length(t_dyn$t_index))
        out_matrix <- out_matrix + measurement_error
    }
    
    if(norm){
        out_matrix <- out_matrix/rowSums(out_matrix)
    }
    colnames(out_matrix) <- names_species
    
    out_matrix <- cbind(out_matrix, time = t_dyn$t_sys[t_dyn$t_index])
    
    #out_matrix$t <- t_dyn$t_sys[t_dyn$t_index]
    #SE <- SummarizedExperiment(assays = list(counts = out_matrix))
    out_list <- list(matrix = out_matrix, 
        community = community, 
        x0 = x0,
        metacommunity_probability = metacommunity_probability,
        migration_p = migration_p, 
        birth_p = birth_p,
        error_variance = error_variance)
    return(out_list)
}
