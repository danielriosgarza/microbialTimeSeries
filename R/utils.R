#' Generate simulation times and the indices of time points to return
#' in simulation functions.
#'
#' @param t.start Numeric scalar indicating the initial time of the simulation.
#' (default: \code{t.start = 0})
#' @param t.end Numeric scalar indicating the final time of the dimulation
#' (default: \code{t.end = 1000})
#' @param t.step Numeric scalar indicating the interval between simulation steps
#' (default: \code{t.step = 0.1})
#' @param t.store Integer scalar indicating the number of evenly distributed
#' time points to keep (default: \code{t.store = 100})
#'
#' @return lists containing simulation times (t.sys) and the indices to keep.
#' @examples
#' Time <- SimulationTimes(t.start = 0, t.end = 100, t.step = 0.5,
#'     t.store = 100)
#' DefaultTime <- SimulationTimes(t.end = 1000)
#'
#' @docType methods
#' @aliases SimulationTimes-numeric
#' @aliases SimulationTimes,numeric-method
#'
#' @keywords internal
#' @export

SimulationTimes <- function(t.start = 0, t.end = 1000, 
    t.step = 0.1, t.store = 1000){
    t.total <- t.end-t.start
    t.sys <- seq(t.start, t.end, by = t.step)
    t.index <- seq(1, length(t.sys)-1, by=floor(length(t.sys)/t.store))
    return(list("t.sys" = t.sys, "t.index" = t.index[1:t.store]))
}

isPositiveInteger <- function(x, tol = .Machine$double.eps^0.5) {
    return(abs(x - round(x)) < tol && x > 0)
}

# ExampleEventTimes <- eventTimes(t.events = c(10,20,30), t.duration = rep(3,3))
eventTimes <- function(t.events = NULL, t.duration = NULL,
                       t.end=1000, ...){
    tdyn <- SimulationTimes(t.end = t.end,...)
    t.result = c()
    for (i in seq(length(t.events))){
        p1 <- tdyn$t.sys[(tdyn$t.sys >= t.events[i]) &
            (tdyn$t.sys < (t.events[i]+t.duration[i]))]
        t.result <- c(t.result, p1)
    }
    return(t.result)
}

applyInterctionType <- function(I, pair, interType){
    if (rbinom(1,1,0.5)){
        pair <- rev(pair)
    }
    if (interType=='mutualism'){
        I[pair[1],pair[2]]=1
        I[pair[2],pair[1]]=1
        return(I)
    }else if (interType=='commensalism'){
        I[pair[1],pair[2]]=1
        I[pair[2],pair[1]]=0
        return(I)
    }else if (interType=='parasitism'){
        I[pair[1],pair[2]]=1
        I[pair[2],pair[1]]=-1
        return(I)
    }else if (interType=='amensalism'){
        I[pair[1],pair[2]]=0
        I[pair[2],pair[1]]=-1
        return(I)
    }else if (interType=='competition'){
        I[pair[1],pair[2]]=-1
        I[pair[2],pair[1]]=-1
        return(I)
    }
}

getInteractions <- function(n.species, weights, connectance){
    I <- matrix(0, n.species, n.species)
    interactions <- c('mutualism', 'commensalism', 'parasitism', 
                      'amensalism', 'competition')
    probs <- abs(weights)/sum(abs(weights))
    combinations <- combn(n.species, 2)
    for (i in sample(seq_along(combinations[1,]),as.integer(ncol(combinations)*connectance), replace=FALSE)){
        I <- applyInterctionType(I, combinations[,i], sample(interactions, 1, prob = probs))
    }
    return (I)
}
