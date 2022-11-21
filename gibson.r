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
library(intergraph)

rpower <- function(n, alpha, norm = FALSE) {
    u <- runif(n, min = 0, max = 1 - .Machine$double.eps^0.5)
    power <- (1 - u)^(1/(1-alpha))
    if (norm) return(power/mean(power)) else return(power)
}

params <- data.frame(
    alpha = c(7, 3, 2, 1.6, 1.01),
    t_end = c(20, 10, 5, 2, 1),
    t_step = c(0.02, 0.01, 0.005, 0.002, 0.001)
)
set.seed(42)
n_species <- 100

H <- list()
df <- list()
plot_hist <- list()
N <- list()
interactions_custom <- list()
A <- list()
g <- list()
g_plot <- list()

local_species_pool <- list()
local_A <- list()


# generate constants ####
for (row in seq_len(nrow(params))) {
    # for each power-law distribution 
    print(paste("row =", row))
    
    H[[row]] <- rpower(n = n_species, alpha = params$alpha[row], norm = TRUE)
    df[[row]] <- data.frame(H[[row]])
    plot_hist[[row]] <- ggplot(df[[row]], aes(H[[row]])) + 
        geom_histogram(fill = 'steelblue', color = 'grey20', bins = 10) + 
        scale_y_log10(breaks = c(1, 10, 100), limits = c(1, 100)) + 
        theme_bw()
    # plot_hist[[row]]
    
    N[[row]] <- matrix(rnorm(n = n_species^2, mean = 0, sd = 1), nrow = n_species)
    diag(N[[row]]) <- 0
    interactions_custom[[row]] <- N[[row]] %*% diag(H[[row]])
    A[[row]] <- randomA(
        n_species = 100,
        diagonal = -1,
        scale_off_diagonal = 0.07,
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
    
    local_species_pool[[row]] <- list()
    local_A[[row]] <- list()
    
    for (local_community in seq_len(100)){
        local_species_pool[[row]][[local_community]] <- sample(x = 100, size = 80)
        local_A[[row]][[local_community]] <- A[[row]][local_species_pool[[row]][[local_community]], local_species_pool[[row]][[local_community]]] 
    }
    # save histogram and network plots ####
    ggsave(paste0("./files/figures/hist", params$alpha[row],".pdf"), plot = plot_hist[[row]], dpi = 300, width = 6, height = 6, units = "cm", scale = 2)
    ggsave(paste0("./files/figures/net", params$alpha[row],".pdf"), plot = g_plot[[row]], dpi = 300, width = 6, height = 6, units = "cm", scale = 2)
}

# identify SIS using H[[row]] ####
SIS <- list()
group_sis <- list()
for (row in seq_len(nrow(params))) {
    SIS[[row]] <- head(order(H[[row]], decreasing = TRUE), 1)
    group_sis[[row]] <- list()
    for (local_community in seq_len(100)){
        if (paste0("sp", SIS[[row]]) %in% rownames(local_A[[row]][[local_community]])){
            group_sis[[row]][[local_community]] <- "with SIS"
        } else {
            group_sis[[row]][[local_community]] <- "without SIS"
        }
    }
    group_sis[[row]] <- unlist(group_sis[[row]])
}

# simulate GLV models using different growth rates for 4 scenarios ####

## prepare result lists ####
simulation_GLV <- list()
otu_table <- list()
jsd <- list()
PCoA <- list()
PCoA_coord <- list()
bestk <- list()
PAM <- list()
SI <- list()
groups <- list()
pcoa_plot <- list()
dfs_to_bind <- list()

scenarios <- c("original", "various growth rates", "SIS grows faster", "non-SIS grows faster")
growthRates <- vector(mode = "list", length = 4)
growthRates[[1]] <- rep(1, 80)
growthRates[[2]] <- NULL

## for each scenario ####
for (scenario in seq_along(scenarios)){
    print(paste0("scenario: ", scenarios[scenario]))
    
    ### make new lists ####
    simulation_GLV[[scenario]] <- list()
    otu_table[[scenario]] <- list()
    jsd[[scenario]] <- list()
    PCoA[[scenario]] <- list()
    PCoA_coord[[scenario]] <- list()
    bestk[[scenario]] <- list()
    PAM[[scenario]] <- list()
    SI[[scenario]] <- list()
    groups[[scenario]] <- list()
    pcoa_plot[[scenario]] <- list()
    
    ### run simulations ####
    for (row in seq_len(nrow(params))) {
        simulation_GLV[[scenario]][[row]] <- list()
        otu_table[[scenario]][[row]] <- data.frame(matrix(nrow = 0, ncol = 100))
        colnames(otu_table[[scenario]][[row]]) <- paste0("sp", 1:100)
        dfs_to_bind[[scenario]] <- list()
        
        #### overwrite growth rates ####
        for (local_community in seq_len(100)){
            if (scenario == 3){
                growthRates[[scenario]] <- 9.9 * local_species_pool[[row]][[local_community]] == SIS[[row]] + 0.1
            }
            if (scenario == 4){
                growthRates[[scenario]] <- -9.9 * local_species_pool[[row]][[local_community]] == SIS[[row]] + 10
            }
            
            simulation_GLV[[scenario]][[row]][[local_community]] <- simulateGLV(
                n_species = 80, 
                names_species = paste0("sp", local_species_pool[[row]][[local_community]]),
                A = local_A[[row]][[local_community]], 
                x0 = rep(1, 80),
                growth_rates = growthRates[[scenario]], 
                t_end = params$t_end[row], 
                t_step = params$t_step[row], 
                stochastic = FALSE,
                norm = TRUE)
            
            # makePlot(simulation_GLV[[scenario]][[row]][[local_community]]$matrix)
            dfs_to_bind[[scenario]] <- append(dfs_to_bind[[scenario]], list(simulation_GLV[[scenario]][[row]][[local_community]]$matrix[time = 1000,]))
        }
        
        ### form otu tables ####
        otu_table[[scenario]][[row]] <- bind_rows(dfs_to_bind[[scenario]])
        otu_table[[scenario]][[row]] <- subset(otu_table[[scenario]][[row]], select = -time)
        otu_table[[scenario]][[row]][is.na(otu_table[[scenario]][[row]])] <- 0
        
        ### calculate jensen-shannon distance and plot PCoA ####
        jsd[[scenario]][[row]] <- JSD(as.matrix(otu_table[[scenario]][[row]]))
        PCoA[[scenario]][[row]] <- ape::pcoa(as.dist(jsd[[scenario]][[row]]))
        PCoA_coord[[scenario]][[row]] <- PCoA[[scenario]][[row]]$vectors[, 1:2]
        colnames(PCoA_coord[[scenario]][[row]]) <- c("PCo1", "PCo2")
        
        bestk[[scenario]][[row]] <- 2
        PAM[[scenario]][[row]] <- pam(PCoA[[scenario]][[row]]$vectors, k = 2)
        SI[[scenario]][[row]] <- PAM[[scenario]][[row]]$silinfo$avg.width
        
        for (k in 3:10){
            PAMtemp <- pam(PCoA[[scenario]][[row]]$vectors, k = k)
            if (PAMtemp$silinfo$avg.width > SI[[scenario]][[row]]) {
                SI[[scenario]][[row]] <- PAMtemp$silinfo$avg.width
                bestk[[scenario]][[row]] <- k
                PAM[[scenario]][[row]] <- PAMtemp
            }
        }
        
        groups[[scenario]][[row]] <- factor(PAM[[scenario]][[row]]$clustering)
        
        dataframe_pcoa <- as.data.frame(PCoA_coord[[scenario]][[row]])
        dataframe_pcoa$group <- groups[[scenario]][[row]]
        pcoa_plot[[scenario]][[row]] <- ggplot(dataframe_pcoa,
                                               aes(x = PCo1, y = PCo2, color = group)) + 
            geom_point() +
            # coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1)) + 
            theme_bw()
        # pcoa_plot[[scenario]][[row]]
    }
}

## export plots of PCoA ####
for (row in seq_len(nrow(params))) {
    for (scenario in seq_along(scenarios)){
        ggsave(paste0("./files/figures/pcoa_", scenario, "_", params$alpha[row],".pdf"), plot = pcoa_plot[[scenario]][[row]], dpi = 300, width = 6, height = 6, units = "cm", scale = 2)
    }
}

## validate relationship between cluster and SIS ####
# with/without SIS forms distinct clusters in PCoA 
pcoa_plot_sis <- list()
for(row in seq_len(nrow(params))){
    pcoa_plot_sis[[row]] <- ggplot(as.data.frame(PCoA_coord[[1]][[row]]), aes(x = PCo1, y = PCo2)) + 
        geom_point(aes(color = group_sis[[row]])) + 
        theme_bw() 
    pcoa_plot_sis[[row]]
    
    ggsave(paste0("./files/figures/pcoa_validation_", params$alpha[row],".pdf"), plot = pcoa_plot_sis[[row]], dpi = 300, width = 6, height = 6, units = "cm", scale = 2)
}

## compare communities in different scenarios of growth rates with original ####
### cosine_similarity ####
cosine_similarity <- vector(mode = "list", length = 4)
for (scenario in seq_along(scenarios)[-1]){
    cosine_similarity[[scenario]] <- list()
    for (row in seq_len(nrow(params))){
        # in each powerlaw alpha generated interspecies interactions scenario 
        cosine_similarity[[scenario]][[row]] <- list()
        for (local_community in seq_len(100)){
            # 2 values in both calculations: without/with growth rates
            cosine_similarity[[scenario]][[row]][[local_community]] <- 
                cosine_dist(
                    P = otu_table[[1]][[row]][[local_community]],
                    Q = otu_table[[scenario]][[row]][[local_community]], 
                    testNA = FALSE)
            # cosine_dist varies from 0 to 1 for non-negative vectors, comparable to a similarity metric
        }
        avg_cos_sim <- mean(unlist(cosine_similarity[[scenario]][[row]]))
        print(paste(scenarios[[scenario]], "row =", row, ",", "average cosine similarity =", avg_cos_sim))
    }
}
## results where growth rates are randomly draw from 0 to 1 ####
## scenario 2: pcoa_2_...pdf
# Result: the less alpha is (the more SIS influence the community), the more similar are the results of final communities

## results where SIS grows faster 10 vs 0.1  ####
## scenario 3: pcoa_3_...pdf
# Result: the faster SIS grows, the more similar are the results of final communities

## results where non-SIS grows faster 10 vs 0.1 ####
## scenario 4: pcoa_4_...pdf
# Result: the faster non-SIS grows, the more similar are the results of final communities

# ## model CRM to validate results ####
# simulationCRM <- list()
# 
# for (row in seq_len(nrow(params))) {
#     # for each power-law distribution 
#     print(paste("row =", row))
#     simulationCRM[[row]] <- list()
#     for (local_community in seq_len(100)){
#         localE <- - t(local_A[[row]][[local_community]])
#         localE[rowSums(localE) < 0,] <- localE[rowSums(localE) < 0,] + sum(-localE[rowSums(localE) < 0,] + 1)/80
#         
#         simulationCRM[[row]][[local_community]] <- simulateConsumerResource(
#             n_species = 80, 
#             n_resources = 80, 
#             resources = rep(100, 80),
#             names_species = paste0("sp", local_species_pool[[row]][[local_community]]), 
#             E = localE, 
#             x0 = rep(1, 80),
#             stochastic = FALSE,
#             norm = TRUE,
#             migration_p = 0, 
#             error_variance = 0, t_end = 500
#             )
#         # makePlot(simulationCRM[[row]][[local_community]]$matrix)
#         # View(simulationCRM[[row]][[local_community]]$matrix)
#     }
# }


# SIS as fast grower : will this get clusters closer? ####
# PCoA: NO. Comparing pcoa_1_... and pcoa_3_... there's no denser or sparser trend.

# use makeUmap to see trends of clusters ####
umap <- list()
for (scenario in seq_along(scenarios)){
    umap[[scenario]] <- list()
    for (row in seq_len(nrow(params))){
        umap[[scenario]][[row]] <- makeUMAP(matrix = otu_table[[scenario]][[row]], gradient_title = "clusters", group = groups[[scenario]][[row]])
        ggsave(paste0("./files/figures/umap_", scenario, "_", params$alpha[row], ".pdf"),  plot = umap[[scenario]][[row]], dpi = 300, width = 6, height = 6, units = "cm", scale = 2)
    }
}
