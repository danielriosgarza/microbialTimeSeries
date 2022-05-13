# example of consumer-resource model ####
## repeat times ####
n.rep <- 50

## output dataframes ####
result.df <- data.frame(
    n_species = integer(),
    theta = numeric(),
    i = integer(),
    n_resources = integer(),
    value = numeric()
)

sorensen.df <- data.frame(
    n_species = integer(),
    theta = numeric(),
    rho.mean = numeric(),
    rho.sd = numeric()
)
## different numbers of organisms ####
for (n_species in c(13, 3, 4)){
    for (theta in c(1, 0.75, 0.5)) {
        sorensen <- c()
        for (i in seq_len(n.rep)){
            for (n_resources in c(1,2,4,8,16,32)) {
                ### generate E ####
                Etest <- randomE(n_species = n_species, n_resources = n_resources, mean_consumption = theta*n_resources, exact = TRUE)

                ### calculate rho ####
                if (n_resources == 32){
                    Etest.pos <- Etest
                    Etest.pos[Etest.pos<0] <- 0
                    for (j in seq_len(n_species - 1)){
                        for (k in 2:n_species){
                            sorensen <- c(sorensen, 
                                          sum(apply(Etest.pos[c(j,k),], 2, min)))
                        }
                    }
                }
                
                if (n_resources > 1){
                    Priority <- t(apply(matrix(sample(n_species * n_resources), nrow = n_species), 1, order))
                    Priority <- (Etest > 0) * Priority
                } else {
                    Priority <- NULL
                }
                print(paste(n_species, theta, i, n_resources))
                CRMtest <- simulateConsumerResource(n_species = n_species,
                                                    n_resources = n_resources,
                                                    x0 = rep(10, n_species),
                                                    resources = rep(1000/n_resources, n_resources),
                                                    E = Etest,
                                                    trophic_priority = Priority,
                                                    stochastic = TRUE)
                CRMspecies <- CRMtest$matrix[1000, seq_len(n_species)]
                CRMspeciesTotal <- sum(CRMspecies)
                result.df[nrow(result.df)+1,] <- c(n_species, theta, i, n_resources, CRMspeciesTotal)
                # makePlotRes(CRMtest$resources)
                # makePlot(CRMtest$matrix)
            }
        }
        rho.mean <- mean(sorensen)
        rho.sd <- var(sorensen)
        sorensen.df[nrow(sorensen.df)+1, ] <- c(n_species, theta, rho.mean, rho.sd)
    }
}

p.fig2.result.df <- ggplot(result.df, aes(x = n_resources, y = value, group = n_resources)) +
    geom_boxplot() + 
    facet_grid(. ~ factor(n_species, levels = c(13, 3, 4))) +
    theme_bw() +
    scale_x_continuous(trans = "log2")
p.fig2.result.df

p.result.df <- ggplot(result.df, aes(x = n_resources, y = value, group = n_resources)) + 
    geom_boxplot() + 
    facet_grid(factor(n_species, levels = c(13, 3,4)) ~ factor(theta, levels = c(1, 0.75,0.5))) +
    theme_bw()+
    scale_x_continuous(trans = "log2")
p.result.df

p.result.df <- p.result.df + geom_text(data = sorensen.df,
                                       mapping = aes(x = 2^2.5, y = 250, label = paste0("ρ = ", round(rho.mean, 2), "±", round(rho.sd, 2))))
p.result.df

# paired one-sided t test
ttest.df <- data.frame(n_species = integer(),
                       theta = numeric(),
                       p = numeric())

for (n.row in seq_len(nrow(sorensen.df))) {
    n_species <- sorensen.df[n.row, "n_species"]
    theta <- sorensen.df[n.row, "theta"]
    ttestres <- t.test(x = result.df[result.df$n_species == n_species & result.df$theta == theta & result.df$n_resources == 1, "value"], 
                       y = result.df[result.df$n_species == n_species & result.df$theta == theta & result.df$n_resources == 32, "value"],
                       alternative = "less")
    ttest.df[nrow(ttest.df)+1, ] <- c(n_species, theta, ttestres$p.value)
}

p.result.df <- p.result.df + geom_text(data = ttest.df,
                                       mapping = aes(x = 2^2.5, y = 300, label = paste0("P = ", signif(p, 2))))
p.result.df
