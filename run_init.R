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
                "mvtnorm", "pracma", "gtools", "maotai", "umap", "BiodiversityR")


setwd("C:/Users/u0139894/Documents/GitHub/microbialTimeSeries")

Rfiles = gsub(" ", "", paste("./R/", list.files("./R")))
sapply(Rfiles, source)