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

getRowMax <- function(aRow){
    maxPos <- which.max(aRow)
    newRow <- rep(FALSE, length(aRow))
    newRow[maxPos] <- TRUE
    return(aRow * newRow)
}

getMoments <- function(simulaionMatrix, simulationResources = NULL, is.perCapita = FALSE, is.resource = FALSE){
    
    simul <- simulaionMatrix
    resources <- simulationResources
    
    if (is.perCapita){
        
        S = sum(simul[1,])
        
        meanH <- colMeans(simul)
        covH <- cov(simul)
        m2 = 2*(1/S)*(1/(S-1))*t(apply(covH, 1, cumsum))
        
        diag(m2) <- (1/S)*diag(covH)
        
        return(list(m1 = (1/S)*meanH, m2 = (1/S)*covH, basis = simulaionMatrix))
    }
    if (is.resource) {
        return(list(m1Matrix = colMeans(simul), 
                    m2Matrix = cov(simul), 
                    basisMatrix = simulaionMatrix,
                    m1Resources = colMeans(resources),
                    m2Resources = cov(resources),
                    basisResources = resources))
    }
    return(list(m1 = colMeans(simul), m2 = cov(simul), basis = simulaionMatrix))
    
    
}


generateMoments <- function(modelGenerateExp, n.instances, t.store, is.perCapita=FALSE, is.resource = FALSE){
    modelSimul <- eval(modelGenerateExp)
    modelMatrix <- modelSimul$matrix
    simulMatrix <- modelMatrix[,colnames(modelMatrix)!="time"]
    summaryMatrix <- matrix(0, nrow = n.instances, ncol = ncol(simulMatrix))
    
    if (is.resource){
        modelResources <- modelSimul$resources
        simulResources <- modelResources[,colnames(modelResources)!="time"]
        summaryResources <- matrix(0, nrow = n.instances, ncol = ncol(simulResources))
    }
    for (i in 1:n.instances){
        print(paste(i, "of", n.instances, "instances"))
        simul = eval(modelGenerateExp)
        modelMatrix <- simul$matrix
        summ <- modelMatrix[,colnames(modelMatrix)!="time"]
        summaryMatrix[i,] <- summ[t.store,]
        if (is.resource){
            modelResources <- simul$resources
            summr <- modelResources[,colnames(modelResources)!="time"]
            summaryResources[i,] <-summr[t.store,]
        }
        
    }
    return(getMoments(summaryMatrix, summaryResources, is.perCapita = is.perCapita, is.resource))
}


getKL <- function(mA, mB){
    mBinv = pinv(mB$m2)
    detA = as.numeric(unlist(pdeterminant(mB$m2))[1])
    detB = as.numeric(unlist(pdeterminant(mA$m2))[1])
    return (max(0, (1/2)*( (detA - detB) - 
                               length(mA$m1) + 
                               sum(diag(mBinv%*%mA$m2)) + 
                               (mB$m1 - mA$m1)%*%mBinv%*%matrix(mB$m1 - mA$m1))))
}

getPerturbT <- function(endTime, n.perturbs){
    st <- SimulationTimes(t.end = endTime, t.store = n.perturbs + 1)
    return (st$t.sys[st$t.index[2:length(st$t.index)]])
}


# plotting functions ####
makePlot <- function(out.matrix, title = "abundance of species by time", obj = "species", y.label = "x.t"){
    df <- as.data.frame(out.matrix)
    dft <-  melt(df, id="time")
    names(dft)[2] = obj
    names(dft)[3] = y.label
    lgd = ncol(df)<= 20
    ggplot(dft, aes_string(names(dft)[1], names(dft)[3], col = names(dft)[2])) +
        geom_line(show.legend = lgd, lwd=0.5) +
        ggtitle(title) + 
        theme_linedraw() +
        theme(plot.title = element_text(hjust = 0.5, size = 14))
}

makePlotRes <- function(out.matrix, title = "quantity of compounds by time"){
    df <- as.data.frame(out.matrix)
    dft <-  melt(df, id="time")
    names(dft)[2] = "resources"
    names(dft)[3] = "S.t"
    lgd = ncol(df)<= 20
    ggplot(dft, aes(time, S.t, col = resources)) + 
        geom_line(show.legend = lgd, lwd=0.5) + 
        ggtitle(title) + 
        theme_linedraw() + 
        theme(plot.title = element_text(hjust = 0.5, size = 14))
}

makePiePlot <- function(multinomdist, label = 'Meta\ncommunity', title = "Metacommunity \nspecies abundance\n"){
    df <- data.frame(group = seq(length(multinomdist)), probability = multinomdist)
    fig <- ggplot(df, aes(x=group,y=1,fill=probability, )) + 
        geom_tile(colour="#edfaf9",size=0.005) +
        theme(axis.title = element_blank()) + 
        scale_fill_gradient2(label, low = "white", high = "magenta3", midpoint = max(multinomdist)/8) +  
        theme_void() +
        coord_fixed(ratio = length(multinomdist)/4) +
        ggtitle(title) # + theme_linedraw()
    fig
}

makeHeatmap <-function(matrix.A, 
                       title = "Consumption/production matrix",
                       y.label = 'resources',
                       x.label = 'species',
                       midpoint_color = NULL, 
                       lowColor = "red", 
                       midColor = "white", 
                       highColor = "blue"){
    df <- melt(t(matrix.A))
    if (is.null(midpoint_color)) {
        midpoint_color <- 0
    }
    names(df)<- c("x", "y", "strength")
    df$y <- factor(df$y, levels=rev(unique(sort(df$y))))
    fig <- ggplot(df, aes(x,y,fill=strength)) + geom_tile() + coord_equal() +
        theme(axis.title = element_blank()) + 
        scale_fill_gradient2('strength', low = lowColor, mid = midColor, high = highColor, midpoint = midpoint_color)+
        theme_void() + ggtitle(title)
    
    if (ncol(matrix.A)<=10 & nrow(matrix.A)<=10){
        fig <- fig + geom_text(aes(label = round(strength, 2))) 
    } else if (ncol(matrix.A)<=15 & nrow(matrix.A)<=15){
        fig <- fig + geom_text(aes(label = round(strength, 1)))
    } else {
        fig <- fig
    }
    
    fig <- fig + labs(x = x.label, y = y.label)+
        theme_linedraw() + 
        theme(plot.title = element_text(hjust = 0.5, size = 14), axis.text.x = element_text(
            angle = 90))
    
    if (nrow(matrix.A) >= 20){
        # too many species 
        fig <- fig + theme(
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
        )
    }
    if (ncol(matrix.A) >= 20){
        # too many resources
        fig <- fig + theme(
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()
        )
    }
    fig
}

makeRegression <- function (fit) {
    
    ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
        geom_point() +
        stat_smooth(method = "lm", col = "red") +
        labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 2),
                           "Intercept =",signif(fit$coef[[1]],2 ),
                           " Slope =",signif(fit$coef[[2]], 2),
                           " P =",signif(summary(fit)$coef[2,4], 2)),)
}


makeUMAP <- function(matrix, n_neighbors=10, min_dist=0.1, gradient=NULL, gradient_title = 'gradient', group=NULL, group2=NULL){
    custom.config = umap.defaults
    custom.config$n_neighbors = n_neighbors
    custom.config$min_dist = min_dist
    
    df <- as.data.frame(umap(matrix,config = custom.config)$layout)
    df$gradient <- gradient
    
    if (is.null(gradient)){
        df$gradient <- 1
        
    }
    colnames(df) = c('UMAP_2', 'UMAP_1', gradient_title)
    if (is.null(group)){
        ggplot(df, aes_string('UMAP_2', 'UMAP_1', color=gradient_title)) + 
            geom_point() + 
            scale_color_gradient(low="blue", high="red")
    } else {
        if (is.null(group2)){
            ggplot(df, aes_string('UMAP_2', 'UMAP_1', color=gradient_title)) + 
                geom_point(aes(color = group)) + theme_bw()
        } else {
            ggplot(df, aes_string('UMAP_2', 'UMAP_1', color=gradient_title)) + 
                geom_point(aes(color = group, shape = group2)) + theme_bw()
        }
        
    }
    
    
}

