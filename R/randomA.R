#' Generate random interaction matrix for the Generalized Lotka-Volterra model(GLV)
#' @param n.species integer number of species
#' @param diagonal values defining the strength of self-interactions. Input can be a number 
#' (will be applied to all species) or a vector of length n.species. 
#' For stability, self-interaction should be negative.
#' (default: \code{diagonal = -0.5})
#' @param connectance numeric frequency of inter-species interactions. i.e. proportion of non-zero 
#' off-diagonal terms.
#' Should be in the [0,1] interval.
#' (default: \code{connectance = 0.2})
#' @param scale numeric: scale of the off-diagonal elements compared to the diagonal. 
#'  (default: \code{scale = 0.1})
#' @param mutualism numeric: relative proportion of interactions terms consistent with mutualism
#' positive <-> positive
#' (default: \code{mutualism = 1})
#' @param commensalism numeric: relative proportion of interactions terms consistent with commensalism
#' positive <-> neutral
#' (default: \code{commensalism = 1})
#' @param parasitism numeric: relative proportion of interactions terms consistent with parasitism
#' positive <-> negative
#' (default: \code{parasitism = 1})
#' @param amensalism numeric: relative proportion of interactions terms consistent with amensalism
#' neutral <-> negative
#' (default: \code{amensalism = 1})
#' @param competition numeric: relative proportion of interactions terms consistent with competition
#' negative <-> negative
#' (default: \code{competition = 1})
#' @param interactions numeric: values of the n.species^2 pairwise interaction strengths.  
#' Diagonal terms will be replaced by the 'diagonal' parameter
#' If NULL, interactions are drawn from runif(n.species^2, min=0, max=abs(diagonal)).
#' Negative values are first converted to positive then the signs are defined by the
#' the relative weights of the biological interactions (i.e. mutualism, commensalism, 
#' parasitism, amensalism, competion) 
#' (default: \code{interactions = NULL})
#' @param symmetric logical whether the strength of mutualistic and competitive interactions are 
#' symmetric
#' (default: \code{symmetric=FALSE})
#' @examples
#' 
#' n.species = 10
#' dense_A <- randomA(n.species = n.species, scale=1, diagonal = -1.0, connectance = 0.9)
#' makeHeatmap(dense_A, lowColor = 'blue', highColor='red')
#'
#' sparse_A <- randomA(n.species = n.species, diagonal = -1.0, connectance = 0.09)
#' makeHeatmap(sparse_A, lowColor = 'blue', highColor='red')
#'
#' user_interactions <- rbeta(n = n.species^2, .5,.5)
#' user_A <- randomA(n.species, interactions = user_interactions)
#' makeHeatmap(user_A, lowColor = 'blue', highColor='red')
#' 
#' competitive_A <- randomA(n.species=n.species, mutualism=0, 
#' commensalism = 0, parasitism =0, amensalism =0, competition=1, connectance =1, scale=1 )
#' makeHeatmap(competitive_A, lowColor = 'blue', highColor='red')
#' 
#' parasitism_A <- randomA(n.species=n.species, mutualism=0, 
#' commensalism = 0, parasitism =1, amensalism =0, competition=0, connectance =1, scale=1, symmetric=TRUE )
#' makeHeatmap(parasitism_A, lowColor = 'blue', highColor='red')


#' @return
#' \code{randomA} returns a matrix A with dimensions (n.species x n.species)
#'
#' @docType methods
#' @aliases randomA-numeric
#' @aliases randomA,numeric-method
#' @export


randomA <- function(n.species, 
                    diagonal = -0.5, 
                    connectance = 0.2, 
                    scale = 0.1,
                    mutualism = 1,
                    commensalism = 1,
                    parasitism = 1,
                    amensalism = 1,
                    competition = 1,
                    interactions=NULL, 
                    symmetric = FALSE){
  if(connectance > 1 || connectance < 0) {
    stop("'connectance' should be in range [0,1]")
  }
  interaction.weights = c(mutualism, commensalism, parasitism, amensalism, competition)
        A <- interactions
      
        if (is.null(interactions)){
            A <- runif(n.species^2, min=0, max=abs(diagonal))
        }
      
        A <- matrix(A,
                        nrow = n.species,
                        ncol = n.species
            )
        
        if(symmetric){
          A[lower.tri(A)] <- t(A)[lower.tri(A)]
        }
        
        I <- getInteractions(n.species, interaction.weights, connectance)
        
        A <- I*abs(A)*(scale*min(abs(diagonal)))
        
        diag(A) <- diagonal
        
        
            
        return(A)
}
