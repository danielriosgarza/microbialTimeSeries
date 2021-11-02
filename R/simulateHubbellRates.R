#' Hubbell's neutral model simulation
#'
#' Neutral species abundances simulation according to the Hubbell model.
#'
#' @param community.initial Numeric: initial species composition
#' @param migration.p Numeric:
#' @param metacommunity.probability Numeric: 
#' @param k.events Integer: 
#' @param growth.rates Numeric:
#' @param norm Logical: whether the time series should be returned with
#' the abundances as proportions (\code{norm = TRUE}) or
#' the raw counts (default: \code{norm = FALSE})
#' @param t.end Numeric: simulation end time (default: \code{t.end = 1000})
#'
#' @docType methods
#' @aliases simulateHubbell-numeric
#' @aliases simulateHubbell,numeric-method
#' @aliases simulateNeutral
#'
#' @examples
#' colData <- DataFrame(sampleID = c(seq_len(100)),
#'                         time = as.Date(100, origin = "2000-01-01"))
#'
#' rowData <- data.frame(Kingdom = "Animalia",
#'                 Phylum = rep(c("Platyhelminthes", "Mollusca"), c(50, 50)),
#'                 Class = rep(c("Turbellaria", "Polyplacophora"), each = 50),
#'                 ASV1 = paste0("D", seq_len(100)),
#'                 ASV2 = paste0("E", seq_len(100)),
#'                 ASV3 = paste0("F", seq_len(100)),
#'                 ASV4 = paste0("G", seq_len(100)),
#'                 ASV5 = paste0("H", seq_len(100)),
#'                 ASV6 = paste0("J", seq_len(100)),
#'                 ASV7 = paste0("K", seq_len(100)),
#'                 row.names = rownames(paste0("species", seq_len(10))),
#'                 stringsAsFactors = FALSE)
#'
#' rowData <- t(rowData)
#'
#' ExampleHubbell <- simulateHubbell(n.species = 8, M = 10, I = 1000, d = 50,
#'                                                         m = 0.02, tend = 100)
#' rowData(ExampleHubbell) <- rowData
#' colData(ExampleHubbell) <- colData
#'
#' @return \code{simulateHubbell} returns a \linkS4class{SummarizedExperiment}
#' class object containing matrix with species abundance as rows and
#' time points as columns
#'
#' @importFrom stats rbinom
#' @importFrom stats rmultinom
#'
#' @references Rosindell, James et al. "The unified neutral theory of
#' biodiversity and biogeography at age ten." Trends in ecology & evolution
#' vol. 26,7 (2011).
#
#' @export
simulateHubbellRates <- function(community.initial, 
                                 migration.p = 0.01, 
                                 metacommunity.probability = NULL,
                                 k.events = 1, 
                                 growth.rates = NULL,
                                 norm = FALSE,
                                 t.end=1000,...)
{
  t.dyn <- SimulationTimes(t.end = t.end,...)
  
  t.store <- length(t.dyn$t.index)
  
  n.species <- length(community.initial)
  
  birth.p <- 1 - migration.p
  
  community <- community.initial
  
  
  
  if (is.null(metacommunity.probability)){
    metacommunity.probability <- rdirichlet(1, alpha = rep(1,n.species))
  }
  metacommunity.probability <- metacommunity.probability/sum(metacommunity.probability)
  
  if (is.null(growth.rates)){
    growth.rates <- rep(1,n.species)
  }
  
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
    
    
    
    composition.probabilities <- community/sum(community)
    
    #k deaths
    
    community <- community -
      (rmultinom(n = 1, size = tau_events, prob = composition.probabilities))
    
    n_births <- rbinom(n=1, size=tau_events, prob = event.probabilities[2])
    
    n_migration <- tau_events-n_births
    
    community <- community +
      (rmultinom(n = 1, size = n_births, prob = composition.probabilities)) +
      (rmultinom(n = 1, size = n_migration, prob = metacommunity.probability))
    
    index <- ((current_t >= stored_time )  & (last_stored_t < stored_time))
    
    
    
    if (sum(index)>0) {
      
      
      out.matrix[index,] <- t(matrix(rep(t(community),sum(index)), ncol = sum(index)))
      last_stored_t <- stored_time[max(seq(t.store)[index])]
      
      
    }
    
    
  }
  

  
  if(norm){
    out.matrix <- out.matrix/rowSums(out.matrix)
  }
  
  
  colnames(out.matrix) <- seq_len(n.species)
  
  out.matrix <- cbind(out.matrix, time = t.dyn$t.sys[t.dyn$t.index])
  
  
  #out.matrix$t <- t.dyn$t.sys[t.dyn$t.index]
  #SE <- SummarizedExperiment(assays = list(counts = out.matrix))
  out.list <- list(matrix = out.matrix, 
                   community = community, 
                   community.initial = community.initial,
                   metacommunity.probability = metacommunity.probability
                   
  )
  return(out.list)
}
