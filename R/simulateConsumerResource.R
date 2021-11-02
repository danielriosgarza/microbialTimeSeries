#' Consumer-resource model simulation
#'
#' Simulates time series with the consumer-resource model and
#' forms a \linkS4class{SummarizedExperiment} object.
#'
#' Simulates a community time series using the consumer-resource model.
#' The change rate of each species was defined as
#' "dx/dt = growth.rates*sum(monod)*X", where
#' growth.rates is the vector of maximum growth rates for the species,
#' monod is the monod growth rate, S/(Ks+S), where S is the concentration of the
#' limiting resource, and Ks is the half-velocity constant for species X and S.
#' X is the vector of abundances of species.
#' The concentrations of resource will be set to 0 if they were calculated
#' less than 0 during the iteration.
#'
#' The resulting abundance matrix model is used to construct
#' \linkS4class{SummarizedExperiment} object.
#'
#' @param n.species Integer: number of species
#' @param n.resources Interger: number of resources
#' @param eff matrix: matrix of efficiency. How efficient are resources
#' converted into biomass, negative values represent excreted resources
#' (default: \code{eff = randomE(n.species, n.resources)})
#' @param x0 Numeric: vector of species
#' (default: \code{x0 = runif(n = n.species, min = 0.1, max = 10)})
#' @param resources Numeric: vector of resources
#' (default: \code{resources = runif(n = n.resources, min = 1, max = 100)})
#' @param growth.rates Numeric: vector of maximum mu of species
#' (default: \code{growth.rates = c(rep(1, n.species))})
#' @param monod.k matrix: matrix of K values in monod model
#' (default: \code{monod.k = matrix(rgamma(n=n.species*n.resources, shape = 50,
#' rate = 0.25), nrow = n.species)})
#'
#' @examples
#' # example1 users provide least parameters.
#' ExampleConsumerResource <- simulateConsumerResource(n.species = 2,
#' n.resources = 4)
#' # visualize the dynamics of the model
#' matplot(t(assays(ExampleConsumerResource)[["counts"]]), type = "l")
#'
#' @return \code{simulateConsumerResource} returns a
#' \linkS4class{SummarizedExperiment} class object containing matrix with
#' species and resources abundance as rows and time points as columns
#'
#' @export
simulateConsumerResource <- function(n.species, n.resources,
    eff = randomE(n.species, n.resources),
    x0 = rep(0.001, n.species),
    resources = sample(seq(50:500), n.resources),
    growth.rates = runif(n.species),
    monod.k = matrix(rgamma(n=n.species*n.resources, shape = 50*max(resources),
        rate = 1), nrow = n.species),
    norm = FALSE,
    t.end = 1000, ...){
    
    t.dyn <- SimulationTimes(t.end = t.end,...)

    # define the consumer-resource model
    consumerResourceModel <- function(t, state, params){
        with(as.list(c(state, params)),{
            x0 <- pmax(0, state[startsWith(names(state), "consumer")])
            resources <- pmax(0, state[startsWith(names(state), "resource")])
            growth.rates <- params[['growth.rates']]
            eff <- params[['eff']]
            monod.k <- params[['monod.k']]
            B <- matrix(rep(resources, length(x0)),
                ncol = length(resources), byrow = TRUE) + monod.k
            growth <- ((eff*(eff>0)/B) %*% resources)*x0
            
            consumption <- (t(growth) %*% ((eff>0)/B))*resources
            production <- -(t(growth) %*% (eff*(eff<0)/B))*resources
            dResources <- -consumption + production
            dConsumers <- growth.rates*growth
            dxdt <- list(c(dConsumers, dResources))
            return(dxdt)
        })
    }
    state.init <- c(x0, resources)
    names(state.init) <- c(paste0("consumer", seq(length(x0))),
        paste0("resource", seq(length(resources))))
    parameters <- list(growth.rates = growth.rates, eff = eff, monod.k = monod.k)

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
    
    colnames(out.species.matrix) <- seq_len(n.species)
    colnames(out.resource.matrix) <- seq_len(n.resources)
    
    out.species.matrix <- cbind(out.species.matrix, time = t.dyn$t.sys[t.dyn$t.index])
    out.resource.matrix <- cbind(out.resource.matrix, time = t.dyn$t.sys[t.dyn$t.index])
    
    out.list <- list(matrix = out.species.matrix, 
                     resources = out.resource.matrix,
                     eff.Mat = eff,
                     monod.k = monod.k)
    
    return(out.list)
    
    
}
