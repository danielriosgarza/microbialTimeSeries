#' Generate random efficiency matrix
#' Generate random efficiency matrix for consumer resource model from a normal
#' distribution. Positive efficiencies indicate the consumption of resources,
#' whilst negatives indicate that the species would produce the resource.
#' Efficiencies larger than 1 would be replaced by 1.
#' @param n.species Integer: number of species
#' @param n.resources Integer: number of resources
#' @param min.con Integer: minimum number of resources consumed by each species
#' @param max.con Integer: maximum number of resources consumed by each species
#' @param min.prod Integer: minimum number of resources produced by each species
#' @param max.prod Integer: maximum number of resources produced by each species
#' @param sd Numeric: standard deviation of the normal distribution
#'
#' @examples
#' # example with specific parameters
#' ExampleEfficiencyMatrix <- randomE(n.species = 3, n.resources = 6,
#' min.con = 3, max.con = 4, min.prod = 1, max.prod = 1, sd = 0.4)
#' # example with minimum parameters
#' ExampleEfficiencyMatrix2 <- randomE(n.species = 5, n.resources = 12)
#'
#' @return
#' \code{randomE} returns a matrix E with dimensions (n.species x n.resources),
#' and each row represents a species.
#'
#' @export
randomE <- function(n.species,
                    n.resources,
                    mean.con = n.resources/4,
                    mean.prod = n.resources/6,
                    maintenance = .5){
    
    efficiency.matrix <- matrix(0, nrow = n.species, ncol = n.resources)
    
    for (i in seq(n.species)){
        
        irow <- efficiency.matrix[i,]
        consumption <- irow
        production <- irow
        
        index.consumption <- sample(seq(n.resources),
                              size = min(n.resources, max(1, rpois(1, mean.con))))
        
        consumption[index.consumption] <- 1
        
        irow <- rdirichlet(1, consumption)[,]
        
        max.prod = sum(irow==0)
        index.production <- sample(seq(n.resources)[irow==0], size=min(max.prod, rpois(1, mean.prod)))
        
        production[index.production] <- 1
        
        prod <- (-1)*(1-maintenance)* rdirichlet(1, production)[,]
        
        irow[index.production] <- prod[index.production]
        
        efficiency.matrix[i,] <- irow
        
        
    }
    
    
    return(efficiency.matrix)
}
