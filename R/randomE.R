#' Generate random efficiency matrix
#' 
#' Generate random efficiency matrix for consumer resource model from a normal
#' distribution. Positive efficiencies indicate the consumption of resources,
#' whilst negatives indicate that the species would produce the resource.
#' 
#' @param n.species Integer: number of species
#' @param n.resources Integer: number of resources
#' @param names.species Character: names of species. If NULL,
#' `paste0("sp", seq_len(n.species))` is used.
#' (default: \code{names.species = NULL})
#' @param names.resources Character: names of resources. If NULL,
#' `paste0("res", seq_len(n.resources))` is used.
#' @param mean.consumption Numeric: mean number of resources consumed by each 
#' species drawn from a poisson distribution
#' (default: \code{mean.consumption = n.resources/4})
#' @param mean.production Numeric: mean number of resources produced by each 
#' species drawn from a poisson distribution
#' (default: \code{mean.production = n.resources/6})
#' @param maintenance Numeric: proportion of resources that cannot be converted 
#' into products
#' (default: \code{maintenance = 0.5})
#' @param trophic.levels Integer: number of species in microbial trophic levels.
#' If NULL, by default, microbial trophic levels would not be considered.
#' (default: \code{trophic.levels = NULL})
#' @param trophic.preferences List: preferred resources and productions of each 
#' trophic level. Positive values indicate the consumption of resources,
#' whilst negatives indicate that the species would produce the resource.
#' If `length(trophic.preferences)` is smaller than `length(trophic.levels)`,
#' then NULL values would be appended to lower trophic levels.
#' If NULL, by default, the consumption preference will be defined randomly.
#' (default: \code{trophic.preferences = NULL})
#' 
#' @examples
#' # example with minimum parameters
#' ExampleEfficiencyMatrix <- randomE(n.species = 5, n.resources = 12)
#' 
#' # examples with specific parameters
#' ExampleEfficiencyMatrix <- randomE(n.species = 3, n.resources = 6,
#'     names.species = letters[1:3], 
#'     names.resources = paste0("res",LETTERS[1:6]),
#'     mean.consumption = 3, mean.production = 1 )
#' ExampleEfficiencyMatrix <- randomE(n.species = 3, n.resources = 6,
#'     maintenance = 0.4)
#' ExampleEfficiencyMatrix <- randomE(n.species = 3, n.resources = 6,
#'     mean.consumption = 3, mean.production = 1, maintenance = 0.4)
#' 
#' # examples with microbial trophic levels
#' ExampleEfficiencyMatrix <- randomE(n.species = 10, n.resources = 15,
#'     trophic.levels = c(6,3,1), 
#'     trophic.preferences = list(c(rep(1,5), rep(-1, 5), rep(0, 5)), 
#'         c(rep(0,5), rep(1, 5), rep(-1, 5)),
#'         c(rep(0,10), rep(1, 5))))
#' ExampleEfficiencyMatrix <- randomE(n.species = 10, n.resources = 15,
#'     trophic.levels = c(6,3,1),
#'     trophic.preferences = list(c(rep(1,5), rep(-1, 5), rep(0, 5)), NULL, NULL))
#' ExampleEfficiencyMatrix <- randomE(n.species = 10, n.resources = 15,
#'     trophic.levels = c(6,3,1))
#' makeHeatmap(ExampleEfficiencyMatrix, lowColor = "Red4", highColor = "Blue")
#' 
#' @return
#' \code{randomE} returns a matrix E with dimensions (n.species x n.resources),
#' and each row represents a species.
#'
#' @export
randomE <- function(n.species,
    n.resources,
    names.species = NULL,
    names.resources = NULL,
    mean.consumption = n.resources/4,
    mean.production = n.resources/6,
    maintenance = 0.5,
    trophic.levels = NULL,
    trophic.preferences = NULL){
    
    # set the default values
    if (is.null(names.species)) {
        names.species <- paste0("sp", seq_len(n.species))
    }
    if (is.null(names.resources)) {
        names.resources <- paste0("res", seq_len(n.resources))
    }
    if (is.null(trophic.levels)) {
        trophic.levels <- n.species
    }
    if (sum(trophic.levels) != n.species) {
        stop("Sum of 'trophic.levels' should equal to 'n.species'.")
    }
    if (!is.null(trophic.preferences)) {
        while(length(trophic.preferences) < length(trophic.levels)){
            warning("Autofilling 'trophic.preferences' with NULL")
            trophic.preferences <- c(trophic.preferences, list(NULL))
        }
    }
    efficiency.matrix <- matrix(0, nrow = n.species, ncol = n.resources,
        dimnames = list(names.species, names.resources))
    
    list.auto.trophic.preference <- list(NULL)
    for (j in seq_len(length(trophic.levels))) {
        n.species.this.level <- trophic.levels[j]
        
        for (i in seq(n.species.this.level)){
            irow <- efficiency.matrix[i+sum(trophic.levels[0:(j-1)]),]
            consumption <- irow
            production <- irow
            # calculate consumption
            consumption.pref <- trophic.preferences[[j]]*(trophic.preferences[[j]]>0)
            if (length(consumption.pref) == 0 && is.null(list.auto.trophic.preference[[j]])) {
                # no consumption preference nor auto.trophic.preference
                # consumption.pref <- NULL
                consumption.pref <- rep(1, n.resources)
                index.consumption <- sample(seq(n.resources), 
                    size = min(max(1, rpois(1, mean.consumption)), n.resources))
            } else { # with consumption preference
                if (length(consumption.pref) == 0) {
                    consumption.pref <- list.auto.trophic.preference[[j]]
                }
                index.consumption <- sample(seq(n.resources),
                    size = min(sum(consumption.pref > 0),
                        max(1, rpois(1, mean.consumption))),
                    replace = FALSE,
                    prob = consumption.pref)
            }
            consumption[index.consumption] <- 1
            irow <- rdirichlet(1, consumption * consumption.pref * 100)

            # calculate production
            production.pref <- trophic.preferences[[j]]*(trophic.preferences[[j]]<0)
            if (sum(production.pref) == 0) { # no production preference
                production.pref <- NULL
                setprod <- setdiff(seq(n.resources), index.consumption)
                if(length(setprod)>0){
                    index.production <- unique(
                        sample(setprod,
                            size = rpois(1, mean.production),
                            replace = TRUE)) 
                    index.production <- setdiff(index.production, index.consumption)
                } else{
                    index.production <- c()
                }
            } else { # with production preference
                index.production <- sample(seq(n.resources),
                    size = min(sum(production.pref < 0),
                        rpois(1, mean.production)),
                    replace = FALSE,
                    prob = abs(production.pref))
            }

            production[index.production] <- 1
            prod <- (-1)*(1-maintenance)*rdirichlet(1, production)
            irow[index.production] <- prod[index.production]
            
            
            efficiency.matrix[i+sum(trophic.levels[0:(j-1)]),] <- irow
        }
        
        # automatically generate consumption of next level according to 
        # the production of this level
        if (j < length(trophic.levels)){
            if (j+1 > length(list.auto.trophic.preference) || is.null(trophic.preferences[[j+1]])) {
                eff.mat <- efficiency.matrix[1:n.species.this.level + sum(trophic.levels[0:(j-1)]),]
                eff.mat[eff.mat > 0] <- 0
                eff.mat <- - eff.mat
                list.auto.trophic.preference[[j+1]] <- colSums(eff.mat)
            }
        }
        
    }
    return(efficiency.matrix)
}
