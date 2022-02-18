# example of consumer-resource model ####
## repeat times ####
n.rep <- 50

## output dataframes ####
result.df <- data.frame(
    n.species = integer(),
    theta = numeric(),
    i = integer(),
    n.resources = integer(),
    value = numeric()
)

sorensen.df <- data.frame(
    n.species = integer(),
    theta = numeric(),
    rho.mean = numeric(),
    rho.sd = numeric()
)
## different numbers of organisms ####
for (n.species in c(13, 3, 4)){
    for (theta in c(1, 0.75, 0.5)) {
        sorensen <- c()
        for (i in seq_len(n.rep)){
            
            ### generate E ####
            Etest <- randomE(n.species = n.species, n.resources = 32, mean.consumption = theta*32, maintenance = 0)
            ### calculate rho ####
            Etest.pos <- Etest
            Etest.pos[Etest.pos<0] <- 0
            for (j in seq_len(n.species - 1)){
                for (k in 2:n.species){
                    sorensen <- c(sorensen, 
                                  sum(apply(Etest.pos[c(j,k),], 2, min)))
                }
            }
            for (n.resources in c(1,2,4,8,16,32)) {
                if (n.resources > 1){
                    Priority <- t(apply(matrix(sample(n.species * n.resources), nrow = n.species), 1, order))
                } else {
                    Priority <- NULL
                }
                sample.resources <- sample(seq_len(32), n.resources)
                print(paste(n.species, theta, i, n.resources))
                CRMtest <- simulateConsumerResource(n.species = n.species,
                                                    n.resources = n.resources,
                                                    x0 = rep(10, n.species),
                                                    resources = rep(1000/n.resources, n.resources),
                                                    E = Etest[,sample.resources],
                                                    trophic.priority = NULL,
                                                    stochastic = TRUE)
                CRMspecies <- CRMtest$matrix[1000, seq_len(n.species)]
                CRMspeciesTotal <- sum(CRMspecies)
                result.df[nrow(result.df)+1,] <- c(n.species, theta, i, n.resources, CRMspeciesTotal)
                # makePlotRes(CRMtest$resources)
                # makePlot(CRMtest$matrix)
            }
        }
        rho.mean <- mean(sorensen)
        rho.sd <- var(sorensen)
        sorensen.df[nrow(sorensen.df)+1, ] <- c(n.species, theta, rho.mean, rho.sd)
    }
}

p.result.df <- ggplot(result.df, aes(x = n.resources, y = value, group = n.resources)) + 
    geom_boxplot() + 
    facet_grid(n.species ~ theta) +
    theme_bw()+
    scale_x_continuous(trans = "log2")
p.result.df
