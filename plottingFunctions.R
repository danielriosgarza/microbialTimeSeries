library(ggplot2)
library(reshape2)

makePiePlot <- function(multinomdist, label = 'Meta\ncommunity', title = "Metacommunity \nspecies abundance\n"){
  df <- data.frame(group = seq(length(multinomdist)), probability = multinomdist)
  fig <- ggplot(df, aes(x=group,y=1,fill=probability, )) + 
    geom_tile(colour="#edfaf9",size=0.005) +
    theme(axis.title = element_blank()) + 
    scale_fill_gradient2(label, low = "white", high = "magenta3", midpoint = max(multinomdist)/8) +  
    theme_void() +
    coord_fixed(ratio = length(multinomdist)/4) +
    ggtitle(title)
  
  fig
}

makePlot <- function(out.matrix){
  df <- as.data.frame(out.matrix)
  dft <-  melt(df, id="time")
  names(dft)[2] = "species"
  names(dft)[3] = "x.t"
  lgd = TRUE
  if (ncol(df)>10){
    lgd = FALSE
  }
  ggplot(dft, aes(time, x.t, col = species)) + geom_line(show.legend = lgd, lwd=0.5)
  
}

makeHeatmap <-function(out.matrix, midpoint_color = NULL){
  
  
  #out.matrix = t(out.matrix[,colnames(out.matrix)!="time"])
  out.matrix = t(out.matrix)
  
  if (is.null(midpoint_color)){
    midpoint_color = max(out.matrix[1,])/8
  }
  df = melt(out.matrix)
  names(df)<- c("x", "y", "abundance")
  df$y <- factor(df$y, levels=rev(unique(sort(df$y))))
  fig <- ggplot(df, aes(x,y,fill=abundance)) + 
    geom_tile(colour="#edfaf9",size=0.005) +
    theme(axis.title = element_blank()) + 
    scale_fill_gradient2('abundance', low = "white", high = "magenta3", midpoint = midpoint_color) +  
    theme_void() +
    scale_y_discrete(expand=c(0,0))
  
  # if (ncol(out.matrix)<21){
  #   fig <- fig + geom_text(aes(label = round(abundance, 1)))
  # }
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


makeUMAP <- function(matrix, n_neighbors=10, min_dist=0.1, gradient=NULL){
  custom.config = umap.defaults
  custom.config$n_neighbors = n_neighbors
  custom.config$min_dist = min_dist
  
  df <- as.data.frame(umap(matrix,config = custom.config)$layout)
  df$gradient <- gradient
  
  if (is.null(gradient)){
    df$gradient <- 1
    
  }
  colnames(df) = c('UMAP_2', 'UMAP_1', 'gradient')
  ggplot(df, aes(UMAP_2, UMAP_1, color=gradient)) + 
    geom_point() + 
    scale_color_gradient(low="blue", high="red")
  
}