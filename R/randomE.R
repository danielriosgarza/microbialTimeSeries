#' Generate random efficiency matrix
#' 
#' Generate random efficiency matrix for consumer resource model from a normal
#' distribution. Positive efficiencies indicate the consumption of resources,
#' whilst negatives indicate that the species would produce the resource.
#' Efficiencies larger than 1 would be replaced by 1.
#' 
#' @param n.species Integer: number of species
#' @param n.resources Integer: number of resources
#' @param mean.consumption Numeric: mean number of resources consumed by each 
#' species drawn from a poisson distribution
#' (default: \code{mean.consumption = n.resources/4})
#' @param mean.production Numeric: mean number of resources produced by each 
#' species drawn from a poisson distribution
#' (default: \code{mean.production = n.resources/6})
#' @param maintenance Numeric: proportion of resources that cannot be converted 
#' into products
#' (default: \code{maintenance = 0.5})
#' 
#' @examples
#' # example with specific parameters
#' ExampleEfficiencyMatrix <- randomE(n.species = 3, n.resources = 6,
#' mean.consumption = 3, mean.production = 1, maintenance = 0.4)
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
    mean.consumption = n.resources/4,
    mean.production = n.resources/6,
    maintenance = 0.5){
    
    efficiency.matrix <- matrix(0, nrow = n.species, ncol = n.resources)
    
    for (i in seq(n.species)){
        
        irow <- efficiency.matrix[i,]
        consumption <- irow
        production <- irow
        
        index.consumption <- sample(seq(n.resources),
            size = min(n.resources, max(1, rpois(1, mean.consumption))))
        
        consumption[index.consumption] <- 1
        
        irow <- rdirichlet(1, consumption)[,]
        
        max.prod = sum(irow==0)
        index.production <- sample(seq(n.resources)[irow==0], 
            size=min(max.prod, rpois(1, mean.production)))
        
        production[index.production] <- 1
        
        prod <- (-1)*(1-maintenance)* rdirichlet(1, production)[,]
        
        irow[index.production] <- prod[index.production]
        
        efficiency.matrix[i,] <- irow
    }
    
    return(efficiency.matrix)
}
