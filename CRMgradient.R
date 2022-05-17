# example of consumer-resource model ####

# make a gradient of environments and a gradient of communities dissimilarity
# trying to show the gradient of dissimilarity 

library(ggplot2)
library(vegan)
# philentropy::cosine_dist cosine distance

set.seed(42)

# shared parameters ####
n_species <- 10
n_resources <- 5
E <- randomE(n_species, n_resources, mean_consumption = 1, mean_production = 3)
growth_rates <- runif(n_species)
monod_constant <- matrix(rbeta(10*5, 10,10),nrow=10, ncol=5)
t_store <- 100
n.instances <- 1 # no stochastic process: no need to repeat

# generating function #### 
gradient.df.generator <- function(n_row, n_col, density_row, max_gradient, error_interval){
    list_initial <- list()
    dissimilarity.gradient <- seq(from = 0, to = max_gradient, length.out = n_row)
    for (i in seq_len(n_row)){
        print(i)
        if (i == 1){
            row_temp <- rbeta(n_col, 1, 10)
            col_to_remove <- sample(x = seq_len(n_col), size = n_col-n_col*density_row)
            row_temp[col_to_remove] <- 0
            list_initial[[i]] <- row_temp
        } else {
            while (length(list_initial) < i) {
                row_temp <- rbeta(n_col, 1, 10)
                col_to_remove <- sample(x = seq_len(n_col), size = n_col-n_col*density_row)
                row_temp[col_to_remove] <- 0
                diff_temp <- abs(vegdist(rbind(list_initial[[1]], row_temp), method = "bray") - dissimilarity.gradient[i])
                if (diff_temp < error_interval) {
                    list_initial[[i]] <- row_temp
                }
            }
        }
    }
    dataframe_to_return <- as.data.frame(t(matrix(unlist(list_initial), ncol = n_row)))
    return(dataframe_to_return)
}

# generate communities #### 
n.community <- 10
density.community <- 0.8
community.initial.df <- gradient.df.generator(n_row = n.community, n_col = n_species, density_row = density.community, max_gradient = 0.8, error_interval = 0.05)

dist.community.initial.df <- vegdist(community.initial.df, method = "bray")
makeHeatmap(as.matrix(dist.community.initial.df), title = "dissimilarity matrix")
makeUMAP(matrix = community.initial.df, group = factor(seq_len(n.community)))


crm_params <- list(
    n_species = n_species,
    n_resources = n_resources,
    E = E,
    resources = rep(1,n_resources),
    monod_constant = monod_constant,
    migration_p = 0,
    stochastic = FALSE,
    t_start = 0,
    t_end = 20,
    t_store = t_store,
    growth_rates = growth_rates,
    norm=FALSE)

# generate resource gradients ####

resourceConcentration <- 10^seq(0,5,1) # 1 to 100000
n.medium <- 10
density.medium <- 0.8
resource.initial.df <- gradient.df.generator(n_row = n.medium, n_col = n_resources, density_row = density.medium, max_gradient = 0.8, error_interval = 0.05)

# test one run ####
crmExample <- simulateConsumerResource(n_species = n_species, n_resources = n_resources, E = E, x0 = as.numeric(community.initial.df[1,]), resources = as.numeric(resourceConcentration[3]*resource.initial.df[1,]), growth_rates = growth_rates, monod_constant = monod_constant, stochastic = FALSE, t_end = 20, t_store = 100, norm = FALSE)
makePlot(crmExample$matrix)
makePlotRes(crmExample$resources)

# generateMoments ####
set.seed(42)
# basisComposition <- matrix(0, ncol=n_species, nrow = 0)
# basisResources <- matrix(0, ncol=n_resources, nrow = 0)

community.simulation <- list()
resource.simulation <- list()

crmgradient <- parse(text = "do.call(simulateConsumerResource, crm_params)")
simulation_counter <- 0
for (resConc in resourceConcentration) {
    for (medium in seq_len(n.medium)){
        crm_params$resources <- as.numeric(resource.initial.df[medium,]*resConc)
        for (community in seq_len(n.community)){
            simulation_counter <- simulation_counter + 1
            print(paste("resConc", resConc, "medium", medium, "community", community, "simulation_counter", simulation_counter))
            crm_params$x0 <- as.numeric(community.initial.df[community,])
            # iterations
            crmMoments <- generateMoments(
                modelGenerateExp = crmgradient, 
                n.instances = n.instances,
                t_store = t_store, 
                is.resource = TRUE)
            community.simulation[[simulation_counter]] <- crmMoments$basisMatrix
            resource.simulation[[simulation_counter]] <- crmMoments$basisResources
        }
    }
}

basisComposition <- do.call(rbind, community.simulation)
basisResources <- do.call(rbind, resource.simulation)
rm(simulation_counter, community.simulation, resource.simulation)

basisComposition_prop <- basisComposition / rowSums(basisComposition)

resource_concentration_type <- as.factor(rep(resourceConcentration, each = n.medium*n.community))
medium_type <- as.factor(rep(seq_len(n.medium), each = n.community ,times = length(resourceConcentration) ))
community_type <- as.factor(rep(seq_len(n.community), times = length(resourceConcentration)*n.medium))

#plot the result in a UMAP space
makeUMAP(basisComposition, group = medium_type, group2 = resource_concentration_type, gradient_title = 'Init.Comm.')
makeUMAP(basisResources, group = community_type, group2 = resource_concentration_type, gradient_title = 'Init.Comm.')
umap_CRM_gradient <- umap(basisComposition_prop)
# umap_CRM_gradient <- umap(basisComposition)
umap_CRM_coor <- as.data.frame(umap_CRM_gradient$layout)
colnames(umap_CRM_coor) <- c("UMAP_1", "UMAP_2")
umap_CRM_coor <- cbind(umap_CRM_coor, resource_concentration_type, medium_type, community_type)
umap_CRM_gradient_plot <- ggplot(umap_CRM_coor, 
                                 aes(UMAP_1, UMAP_2, 
                                     alpha = resource_concentration_type, 
                                     color = medium_type, 
                                     shape = community_type)) +
    geom_point() + 
    scale_shape_manual(values = c(0, 1, 2, 5, 6, 8, 15, 16, 17, 18)) +
    scale_alpha_manual(values = seq(0.25, 15, 0.1)) + 
    theme_bw()
ggsave(paste0("./files/figures/CRMgradient1.pdf"), plot = umap_CRM_gradient_plot , dpi = 300, width = 12, height = 10, units = "cm", scale = 2)
ggsave(paste0("./files/figures/CRMgradient2.pdf"), plot = umap_CRM_gradient_plot + facet_grid(resource_concentration_type ~ .), dpi = 300, width = 8, height = 16, units = "cm", scale = 2)
ggsave(paste0("./files/figures/CRMgradient3.pdf"), plot = umap_CRM_gradient_plot + facet_grid(medium_type ~ resource_concentration_type), dpi = 300, width = 16, height = 16, units = "cm", scale = 2)
ggsave(paste0("./files/figures/CRMgradient4.pdf"), plot = umap_CRM_gradient_plot + facet_grid(community_type ~ resource_concentration_type), dpi = 300, width = 16, height = 16, units = "cm", scale = 2)
ggsave(paste0("./files/figures/CRMgradient5.pdf"), plot = umap_CRM_gradient_plot + facet_grid(community_type ~ medium_type), dpi = 300, width = 20, height = 16, units = "cm", scale = 2)

## TODO: for figure 3: reduce medium types and add community types, to see the 'enterotypes' always form stable cluster with abundant resources
## TODO: for figure 4: measure the average dissimilarity between communities (with community 1) to plot a 'saturation curve', to show at which concentration the communities differ totally because their intrinsic consumer-resource matrix (an implicit interspecies interactions)
## TODO: alter growth rate distributions using rbeta

# dissimilarity of initial communities vs final communities ####
dist.community.initial.jaccard <- vegdist(community.initial.df, method = "jaccard")
makeHeatmap(as.matrix(dist.community.initial.jaccard), title = "jaccard dissimilarity initial communities", x.label = "community", y.label = "community")

# basicComposition_proportional
basisComposition_prop <- cbind(basisComposition_prop, 
                               as.numeric(as.character(resource_concentration_type)), 
                               medium_type, 
                               community_type)
colnames(basisComposition_prop)[11] <- "resource_concentration_type"

## control resource concentration and different mediums #### 
## (diversity~difference correlation)
dist.community.simulation.list <- list()
simulation_counter2 <- 0 #1
# dist.community.simulation.list[[1]] <- c(unname(as.matrix(dist.community.initial.jaccard)[1,]), 0, 0)
for (resConc in resourceConcentration) {
    for (medium in seq_len(n.medium)) {
        simulation_counter2 <- simulation_counter2 + 1
        print(paste("resConc", resConc, "medium", medium, "simulation_counter2", simulation_counter2))
        relative_concentration_df <- subset(basisComposition_prop, resource_concentration_type == resConc & medium_type == medium)[,1:10]
        dist.community.simulation <- vegdist(relative_concentration_df, method = "jaccard")
        dist.community.simulation <- as.matrix(dist.community.simulation)[1,]
        dist.community.simulation <- c(unname(dist.community.simulation), resConc, medium)
        dist.community.simulation.list[[simulation_counter2]] <- dist.community.simulation
    }
}
dist.community.simulation.df <- as.data.frame(do.call(rbind, dist.community.simulation.list))
colnames(dist.community.simulation.df) <- c(paste0("com", seq_len(10)), "resource_concentration_type", "medium_type")
rm(dist.community.simulation.list, simulation_counter2)

# plot dissimilarity distribution in simulation
dist.community.simulation.df.melt <- melt(dist.community.simulation.df, id.vars = c("resource_concentration_type", "medium_type"))
dist.community.initial.jaccard.vector <- unname(as.matrix(dist.community.initial.jaccard)[,1])
levels(dist.community.simulation.df.melt$variable) <- dist.community.initial.jaccard.vector
dist.community.simulation.df.melt$variable <- as.numeric(as.character(dist.community.simulation.df.melt$variable))

dist.community.simulation.df.melt$resource_concentration_type <- as.factor(dist.community.simulation.df.melt$resource_concentration_type)
dist.community.simulation.df.melt$medium_type <- as.factor(dist.community.simulation.df.melt$medium_type)
dist.community.simulation.df.plot <- ggplot(
    dist.community.simulation.df.melt,
    aes(x = variable, y = value, color = resource_concentration_type)) +
    geom_point(alpha = 0.4) + 
    # geom_boxplot() + geom_jitter(alpha = 0.3) + 
    # facet_wrap(medium_type ~ .) + 
    geom_smooth(method = "lm", se = FALSE) + 
    xlab("initial Jaccard dissimilarity") + 
    ylab("final Jaccard dissimilarity")+
    theme_bw()
dist.community.simulation.df.plot
ggsave("./files/figures/CRMJaccardDissimilarityChange.pdf", plot = dist.community.simulation.df.plot , dpi = 300, width = 8, height = 6, units = "cm", scale = 2)


# the higher resource concentration is, the more diversity changes
dist.community.simulation.df.plot2 <- ggplot(
    dist.community.simulation.df.melt,
    aes(x = resource_concentration_type, y = value, fill = resource_concentration_type)) +
    geom_boxplot() + 
    geom_jitter(alpha = 0.3) + 
    # facet_wrap(medium_type ~ .) + 
    xlab("resource concentration gradients") + 
    ylab("final Jaccard dissimilarity")+
    theme_bw()
dist.community.simulation.df.plot2
ggsave("./files/figures/CRMJaccardDissimilarityResourceConcentration.pdf",
       plot = dist.community.simulation.df.plot2 , dpi = 300, width = 8, height = 6, units = "cm", scale = 2)

# plot changes in dissimilarity before and after simulation 
# abs(dist.community.simulation.df - dist.community.initial) with a. dissimilarity between medium and b. medium concentration
change.dist <- dist.community.simulation.df
change.dist[,1:10] <- sweep(change.dist[,1:10], 2, dist.community.initial.jaccard.vector)
change.dist.melt <- melt(change.dist, id.vars =  c("resource_concentration_type", "medium_type"))
# levels(change.dist.melt$variable) <- dist.community.initial.jaccard.vector
# change.dist.melt$variable <- as.numeric(as.character(change.dist.melt$variable))
change.dist.melt$resource_concentration_type <- as.factor(change.dist.melt$resource_concentration_type)
change.dist.melt$medium_type <- as.factor(change.dist.melt$medium_type)
# levels(change.dist.melt$variable) <- dist.community.initial.jaccard.vector
# change.dist.melt$variable <- as.numeric(as.character(change.dist.melt$variable))

change.dist.plot <- ggplot(
    change.dist.melt,
    aes(x = variable, y = value, color = medium_type)) +
    geom_point() +
    # geom_boxplot() + 
    # geom_jitter(alpha = 0.3, color = medium_type) + 
    # facet_wrap(medium_type ~ .) + 
    theme_bw()
change.dist.plot
change.dist.plot2 <- ggplot(
    change.dist.melt,
    aes(x = variable, y = medium_type, fill = value)) + 
    geom_tile() +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue")+
    facet_wrap(resource_concentration_type ~ .)+
    theme_bw()
change.dist.plot2
change.dist.plot3 <- ggplot(
    change.dist.melt,
    aes(x = variable, y = value)) +
    geom_boxplot() + geom_jitter(alpha = 0.3) + 
    #facet_wrap(resource_concentration_type ~ .) + 
    theme_bw()
change.dist.plot3

save.image("CRMgradient.RData")
