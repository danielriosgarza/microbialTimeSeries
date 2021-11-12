#' Consumer-resource model simulation
#'
#' Simulates time series with the consumer-resource model
#' 
#' The resulting abundance matrix model is used to construct
#' \linkS4class{SummarizedExperiment} object.
#'
#' @param n.species Integer: number of species
#' @param n.resources Interger: number of resources
#' @param names.species Character: names of species. If NULL,
#' `paste0("sp", seq_len(n.species))` is used.
#' (default: \code{names.species = NULL})
#' @param names.resources Character: names of resources. If NULL,
#' `paste0("res", seq_len(n.resources))` is used.
#' @param E matrix: matrix of efficiency. A matrix defining the efficiency of
#' resource-to-biomass conversion (positive values) and the relative conversion
#' of metabolic by-products (negative values)
#' (default: \code{E = NULL})
#' @param x0 Numeric: initial abundances of simulated species. If NULL,
#' `runif(n = n.species, min = 0.1, max = 10)` is used.
#' (default: \code{x0 = NULL})
#' @param resources Numeric: initlal concentrations of resources. If NULL,
#' `runif(n = n.resources, min = 1, max = 100)` is used.
#' (default: \code{resources = NULL})
#' @param growth.rates Numeric: vector of maximum growth rates(mu) of species.
#' If NULL, `rep(1, n.species)` is used.
#' (default: \code{growth.rates = NULL})
#' @param monod.constant matrix: the constant of additive monod growth of
#' n.species consuming n.resources. If NULL,
#' `matrix(rgamma(n = n.species*n.resources, shape = 50*max(resources), rate = 1), nrow = n.species)`
#' is used.
#' (default: \code{monod.constant = NULL})
#' @param return.matrix Logical: whether to export only the stored time points
#' in a matrix
#' (default: \code{return.matrix = FALSE})
#' @param t.end Numeric scalar indicating the final time of the dimulation
#' (default: \code{t.end = 1000})
#' @param ... additional parameters including 't.start', 't.step', and 't.store'
#'
#' @examples
#' ExampleConsumerResource <- simulateConsumerResource(n.species = 2, 
#'     n.resources = 4)
#' # visualize
#' makePlot(ExampleConsumerResource$matrix)
#' makePlot(ExampleConsumerResource$resources)
#'
#' @return \code{simulateConsumerResource} returns a
#' \linkS4class{SummarizedExperiment} class object containing matrix with
#' species and resources abundance as rows and time points as columns
#'
#' @export
simulateConsumerResource <- function(n.species, n.resources,
    names.species = NULL,
    names.resources = NULL,
    E = NULL,
    x0 = NULL,
    resources = NULL,
    growth.rates = NULL,
    monod.constant = NULL,
    norm = FALSE,
    t.end = 1000, ...){
    
    t.dyn <- SimulationTimes(t.end = t.end,...)

    # define the consumer-resource model
    consumerResourceModel <- function(t, state, params){
        with(as.list(c(state, params)),{
            x0 <- pmax(0, state[startsWith(names(state), "consumer")])
            resources <- pmax(0, state[startsWith(names(state), "resource")])
            growth.rates <- params[['growth.rates']]
            E <- params[['E']]
            monod.constant <- params[['monod.constant']]
            B <- matrix(rep(resources, length(x0)),
                ncol = length(resources), byrow = TRUE) + monod.constant
            growth <- ((E*(E>0)/B) %*% resources)*x0
            
            consumption <- (t(growth) %*% ((E>0)/B))*resources
            production <- -(t(growth) %*% (E*(E<0)/B))*resources
            dResources <- -consumption + production
            dConsumers <- growth.rates*growth
            dxdt <- list(c(dConsumers, dResources))
            return(dxdt)
        })
    }
    
    # set the default values
    if (is.null(names.species)) {
        names.species <- paste0("sp", seq_len(n.species))
    }
    if (is.null(names.resources)) {
        names.resources <- paste0("res", seq_len(n.resources))
    }
    if (is.null(E)) {
        E <- randomE(n.species, n.resources)
    }
    if (is.null(x0)) {
        x0 <- runif(n = n.species, min = 0.1, max = 10)
    }
    if (is.null(resources)) {
        resources <- runif(n = n.resources, min = 1, max = 100)
    }
    if (is.null(growth.rates)) {
        growth.rates <- rep(1, n.species)
    }
    if (is.null(monod.constant)) {
        monod.constant <- matrix(rgamma(n = n.species*n.resources,
            shape = 50*max(resources), rate = 1), nrow = n.species)
    }
    
    state.init <- c(x0, resources)
    names(state.init) <- c(paste0("consumer", seq(length(x0))),
        paste0("resource", seq(length(resources))))
    parameters <- list(growth.rates = growth.rates, E = E, monod.constant = monod.constant)

    out <- as.data.frame(ode(y = state.init, times = t.dyn$t.sys,
        func = consumerResourceModel, parms = parameters))
    
    species.index <- grepl('consumer', names(out))
    resource.index <- grepl('resource', names(out))
    
    out.species.matrix <- as.matrix(out[,species.index])
    out.species.matrix <- out.species.matrix[t.dyn$t.index,]
    
    out.resource.matrix <- as.matrix(out[,resource.index])
    out.resource.matrix <- out.resource.matrix[t.dyn$t.index,]
    
    if(norm){
        out.species.matrix <- out.species.matrix/rowSums(out.species.matrix)
        out.resource.matrix <- out.resource.matrix/rowSums(out.resource.matrix)
    }
    
    colnames(out.species.matrix) <- names.species
    colnames(out.resource.matrix) <- names.resources
    
    out.species.matrix <- cbind(out.species.matrix, time = t.dyn$t.sys[t.dyn$t.index])
    out.resource.matrix <- cbind(out.resource.matrix, time = t.dyn$t.sys[t.dyn$t.index])
    
    out.list <- list(matrix = out.species.matrix, 
                     resources = out.resource.matrix,
                     E.Mat = E,
                     monod.constant = monod.constant)
    
    return(out.list)
    
    
}
