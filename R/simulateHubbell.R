#' Hubbell's neutral model simulation
#'
#' Neutral species abundances simulation according to the Hubbell model.
#'
#' @param community.initial Numeric: initial species composition
#' @param migration.p Numeric: 
#' (default: \code{m = 0.02})
#' @param metacommunity.p Numeric:
#' @param k.events:
#' @param norm Logical: whether the time series should be returned with
#' the abundances as proportions (\code{norm = TRUE}) or
#' the raw counts (default: \code{norm = FALSE})
#' @param t.end Numeric: simulation end time (default: \code{t.end = 1000})
#' @param ... additional parameters including 't.start', 't.step', and 't.store'
#' see \code{\link{SimulationTimes}} for more information
#'
#' @docType methods
#' @aliases simulateHubbell-numeric
#' @aliases simulateHubbell,numeric-method
#' @aliases simulateNeutral
#' @importFrom gtools rdirichlet
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
simulateHubbell <- function(community.initial, 
                            migration.p = 0.01, 
                            metacommunity.p = NULL,
                            k.events = 1,
                            norm = FALSE, 
                            t.end=1000, ...){
  
  
  t.dyn <- SimulationTimes(t.end = t.end,...)
  
  n.species <- length(community.initial)
  
  birth.p <- 1 - migration.p
  
  community <- community.initial
  
  
  
  if (is.null(metacommunity.p)){
    metacommunity.p <- rdirichlet(1, alpha = rep(1,n.species))
  }
  #normalize metacommunity.p
  metacommunity.p <- metacommunity.p/sum(metacommunity.p)
  
  out.matrix <- matrix(0, nrow=length(t.dyn$t.index), ncol = n.species)
  
  
  counter = 1

  out.matrix[counter,] <- community
  
  
  
  for (i in t.dyn$t.sys[2:length(t.dyn$t.sys)]){
    
    
    probabilities <- community/sum(community)
    n.deaths <- min(min(community[community>0]),k.events)
    
    #deaths
    community = community - t(rmultinom(n = 1, size = n.deaths, prob = probabilities))
    n.births = sum(rbinom(n=n.deaths, size=1, p = birth.p))
    n.migration = n.deaths-n.births
    
    community = community + t(rmultinom(n = 1, size = n.births, prob = probabilities)) +
      t(rmultinom(n = 1, size = n.migration, prob = metacommunity.p))
    
    if (i %in% t.dyn$t.sys[t.dyn$t.index]){
      counter = counter + 1
      out.matrix[counter,] <- community
      
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
                   metacommunity.p = metacommunity.p,
                   migration.p = migration.p, 
                   birth.p = birth.p
                   )
  return(out.list)
}
