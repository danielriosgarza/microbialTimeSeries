# example of gLV model ####
# Gibson

library(ggplot2)
library(igraph)
library(colourvalues)
library(GGally)
library(network)
library(sna)
library(dplyr)
library(philentropy)
library(cluster)
library(ape)

rpower <- function(n, alpha, norm = FALSE) {
    u <- runif(n, min = 0, max = 1 - .Machine$double.eps^0.5)
    power <- (1 - u)^(1/(1-alpha))
    if (norm) return( power / mean(power) ) else return(power)
}

params <- data.frame(
    alpha = c(7, 3, 2, 1.6, 1.01),
    t.end = c(20, 10, 5, 2, 1),
    t.step = c(0.02, 0.01, 0.005, 0.002, 0.001)
)

set.seed(42)
n.species <- 100
H <- list()
df <- list()
plot_hist <- list()
N <- list()
interactions_custom <- list()
A <- list()
g <- list()
g_plot <- list()
otu_table <- list()
jsd <- list()
PCoA <- list()
PCoA_coord <- list()
bestk <- list()
PAM <- list()
SI <- list()
pcoa_plot <- list()
groups <- list()
for (row in seq_len(nrow(params))) {
    # for each power-law distribution 
    print(paste("row =", row))
    H[[row]] <- rpower(n = n.species, alpha = params$alpha[row], norm = TRUE)
    df[[row]] <- data.frame(H[[row]])
    plot_hist[[row]] <- ggplot(df[[row]], aes(H[[row]])) + 
        geom_histogram(fill = 'steelblue', color = 'grey20', bins = 10) + 
        scale_y_log10(breaks = c(1, 10, 100), limits = c(1, 100)) + 
        theme_bw()
    plot_hist[[row]]
    
    
    N[[row]] <- matrix(rnorm(n = n.species^2, mean = 0, sd = 1), nrow = n.species)
    diag(N[[row]]) <- 0
    
    interactions_custom[[row]] <- N[[row]] %*% diag(H[[row]])
    A[[row]] <- randomA(
        n.species = 100,
        diagonal = -1,
        scale.offDiagonal = 0.07,
        connectance = 1,
        interactions = interactions_custom[[row]])
    g[[row]] <- graph_from_adjacency_matrix(A[[row]], weighted = TRUE, diag = FALSE)
    print(paste("max edge weight =", max(E(g[[row]])$weight)))
    E(g[[row]])$color <- colour_values(E(g[[row]])$weight, palette = "viridis", alpha = 63, include_alpha = TRUE)
    g_plot[[row]] <- ggnet2(
        net = g[[row]], 
        mode = "circle",
        node.color = "black",
        node.size = 1,
        edge.color = "color", edge.alpha = 0.25)
    g_plot[[row]]
    
    local_species_pool <- list()
    local_A <- list()
    simulation_GLV <- list()
    otu_table[[row]] <- data.frame(matrix(nrow = 0, ncol = 100))
    colnames(otu_table[[row]]) <- paste0("sp", 1:100)
    dfs_to_bind <- list()
    
    for (local_community in seq_len(100)){
        # for each local community
        if (local_community %% 10 == 0) print(paste("local_community =", local_community))
        local_species_pool[[local_community]] <- sample(x = 100, size = 80)
        local_A[[local_community]] <- A[[row]][local_species_pool[[local_community]], local_species_pool[[local_community]]] 
        simulation_GLV[[local_community]] <- simulateGLV(
            n.species = 80, 
            names.species = paste0("sp", local_species_pool[[local_community]]),
            A = local_A[[local_community]], 
            x0 = rep(1, 80),
            growth.rates = rep(1, 80), 
            t.end = params$t.end[row], 
            t.step = params$t.step[row], 
            stochastic = FALSE,
            norm = TRUE)
        # makePlot(simulation_GLV[[local_community]]$matrix)
        dfs_to_bind <- append(dfs_to_bind, list(simulation_GLV[[local_community]]$matrix[time = 1000,]))
    }
    otu_table[[row]] <- bind_rows(dfs_to_bind)
    otu_table[[row]] <- subset(otu_table[[row]], select = -time)
    otu_table[[row]][is.na(otu_table[[row]])] <- 0
    
    jsd[[row]] <- JSD(as.matrix(otu_table[[row]]))
    PCoA[[row]] <- ape::pcoa(as.dist(jsd[[row]]))
    
    PCoA_coord[[row]] <- PCoA[[row]]$vectors[, 1:2]
    colnames(PCoA_coord[[row]]) <- c("PCo1", "PCo2")
    
    bestk[[row]] <- 2
    PAM[[row]] <- pam(PCoA[[row]]$vectors, k = 2)
    SI[[row]] <- PAM[[row]]$silinfo$avg.width
    
    for (k in 3:10){
        PAMtemp <- pam(PCoA[[row]]$vectors, k = k)
        if (PAMtemp$silinfo$avg.width > SI[[row]]) {
            SI[[row]] <- PAMtemp$silinfo$avg.width
            bestk[[row]] <- k
            PAM[[row]] <- PAMtemp
        }
    }
    
    groups[[row]] <- factor(PAM[[row]]$clustering)
    pcoa_plot[[row]] <- ggplot(as.data.frame(PCoA_coord[[row]]), aes(x = PCo1, y = PCo2)) + 
        geom_point(aes(color = groups[[row]])) + 
        # coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1)) + 
        theme_bw() 
    pcoa_plot[[row]]
}

# export plots
for (row in seq_len(nrow(params))) {
    ggsave(paste0("./files/figures/hist", params$alpha[row],".pdf"), plot = plot_hist[[row]], dpi = 300, width = 6, height = 6, units = "cm", scale = 2)
    ggsave(paste0("./files/figures/net", params$alpha[row],".pdf"), plot = g_plot[[row]], dpi = 300, width = 6, height = 6, units = "cm", scale = 2)
    ggsave(paste0("./files/figures/pcoa_", params$alpha[row],".pdf"), plot = pcoa_plot[[row]], dpi = 300, width = 6, height = 6, units = "cm", scale = 2)
}
