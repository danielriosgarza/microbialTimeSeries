installIfNeeded<-function(...) {
  libs<-unlist(list(...))
  req<-unlist(lapply(libs,require,character.only=TRUE))
  need<-libs[req==FALSE]
  if(length(need)>0){ 
    install.packages(need)
    lapply(need,require,character.only=TRUE)
  }
}

installIfNeeded("ggplot2", "deSolve", "reshape2", "Matrix", 
                "mvtnorm", "pracma", "gtools", "maotai")


getMoments <- function(simulaionMatrix, is.perCapita = FALSE){
  
  simul <- simulaionMatrix
  
  if (is.perCapita){
    
    S = sum(simul[1,])
    
    meanH <- colMeans(simul)
    covH <- cov(simul)
    m2 = 2*(1/S)*(1/(S-1))*t(apply(covH, 1, cumsum))
    
    diag(m2) <- (1/S)*diag(covH)
    
    return(list(m1 = (1/S)*meanH, m2 = (1/S)*covH, basis = simulaionMatrix))
  }
  
  return(list(m1 = colMeans(simul), m2 = cov(simul), basis = simulaionMatrix))
  
  
}


generateMoments <- function(modelGenerateExp, n.instances, t.store, is.perCapita=FALSE){

  modelMatrix <- eval(modelGenerateExp)$matrix
  
  simul <- modelMatrix[,colnames(modelMatrix)!="time"]
  
  summaryMatrix <- matrix(0, nrow = n.instances, ncol = ncol(simul))
  
  for (i in 1:n.instances){
    print(i)
    modelMatrix <- eval(modelGenerateExp)$matrix
    
    simul <- modelMatrix[,colnames(modelMatrix)!="time"]
    
    summaryMatrix[i,] <- simul[t.store,]
    
  }
  
  return(getMoments(summaryMatrix, is.perCapita = is.perCapita))
  
  
  
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
