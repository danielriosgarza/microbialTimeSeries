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
#' @param interaction.weights numeric length 5 vector with the relative proportion of random 
#' interaction pairs consistent with mutualism(1,1), commensalism (1,0), parasitism (1,-1), 
#' amensalism (0,-1), or competition (-1,-1), where 1 stands for positive, -1 for negative, 
#' and 0 for neutral relationships).
#' Neutralism (0,0) is defined by the connectance parameter. 
#' Any numerical values are accepted but will be converted to positive and divided by their sum.
#' (default: \code{interaction.weights = c(1,1,1,1,1)})
#' @param distribution numeric a n.species*n.species dimensional vector that contains a draw of
#' interaction strengths. If NULL, interactions are drawn from runif(n.species^2, min=0, max=abs(diagonal)).
#' Negative values are converted to positive. The biological interactions that are drawn from the
#' "interaction.weights" parameter define the signs. 
#' (default: \code{distribution = NULL})
#' @param symmetric logical return a symmetric interaction matrix
#' (default: \code{symmetric=FALSE})
#' @examples
#' high_inter_A <- randomA(n.species = 10, diagonal = -0.4, min.strength = -0.8,
#'                                     max.strength = 0.8, connectance = 0.5)
#'
#' low_inter_A <- randomA(n.species = 10, connectance = 0.01)
#'
#' @return
#' \code{randomA} returns a matrix A with dimensions (n.species x n.species)
#'
#' @docType methods
#' @aliases randomA-numeric
#' @aliases randomA,numeric-method
#' @export


randomA <- function(n.species, diagonal = -0.5, connectance = 0.2, interaction.weights = c(1,1,1,1,1), 
             scale=0.1, distribution=NULL, symmetric = FALSE){
  if(connectance > 1 || connectance < 0) {
    stop("'connectance' should be in range [0,1]")
  }
        A <- distribution
      
        if (is.null(distribution)){
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
