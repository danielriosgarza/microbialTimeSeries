library(shiny)
library(ggplot2)
library(deSolve)
library(reshape2)
library(gtools)

Rfiles = gsub(" ", "", paste("./R/", list.files("./R")))
sapply(Rfiles, source)

makePlot <- function(out_matrix){
  df <- as.data.frame(out_matrix)
  dft <-  melt(df, id="time")
  names(dft)[2] = "species"
  names(dft)[3] = "x.t"
  lgd = TRUE
  if (ncol(df)>10){
    lgd = FALSE
  }
  ggplot(dft, aes(time, x.t, col = species)) + geom_line(show.legend = lgd, lwd=0.5)
  
}

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

makeHeatmap <-function(matrix.A, title){
  df = melt(matrix.A)
  names(df)<- c("x", "y", "strength")
  df$y <- factor(df$y, levels=rev(unique(sort(df$y))))
  fig <- ggplot(df, aes(x,y,fill=strength)) + geom_tile() + coord_equal() +
  theme(axis.title = element_blank()) + 
  scale_fill_gradient2('strength', low = "red", mid = "white", high = "blue", midpoint = 0)+  theme_void() + ggtitle(title)
  
  if (ncol(matrix.A)<21){
    fig <- fig + geom_text(aes(label = round(strength, 1)))
  }
  
  fig
}


getPerturbT <- function(endTime, n.perturbs){
  st <- SimulationTimes(t_end = endTime, t_store = n.perturbs + 1)
  return (st$t_sys[st$t_index[2:length(st$t_index)]])
}
#constants
sp.x0 <- 0.001

ui <- fluidPage(
  titlePanel("Generalized Lotka-Volterra (GLV)"),
  sidebarLayout(
    sidebarPanel(width = 2,
                 sliderInput('n.sp', 'number of species', value=5, min=1, max = 100, step = 1),
                 br(),
                 textOutput("x0"),
                 br(),
                 sliderInput('t_end', 'simulation time', value=100, min=10, max = 2000, step = 10)),
    mainPanel(tabsetPanel(tabPanel("Metacommunity",
                                   sliderInput("multiplier", "Eveness", min = 1, max=1000, value=1.0, step = 0.1),
                                   
                                   tags$div(sliderInput("migration.r", "Migration rate", min = 0, max=1, value=0.1, step = 0.01),  style="display:inline-block"),
                                   
                                   tags$div(sliderInput("sigma_migration", "sd of migration", min = 0, max=1, value=0.01, step = 0.001),  style="display:inline-block"),
                                   
                                   
                                   plotOutput("metacommunityPlot", width = 600, height = 400)
    ),
    tabPanel("Growth rates", 
             htmlOutput("betaDist"),
             br(),
             tags$div(sliderInput("alpha", "alpha", min = 0, max=10, value=0.5, step = 0.1),  style="display:inline-block"),
             tags$div(sliderInput("beta", "beta", min = 0, max=10, value=0.5, step = 0.1),  style="display:inline-block"),
             textOutput("gr.rates1"),
             verbatimTextOutput("gr.rates2"),
             textOutput("gr.rates3"),
             plotOutput("gr.expect", width = 150, height = 150),
             br(),
             br(),
             textOutput("carryingK1"),
             verbatimTextOutput("carryingK2"),
             br()),
    tabPanel("Interactions",
             htmlOutput("interactionMat"),
             br(),
             sliderInput('s.distA', 'sd of interaction matrix', value=1, min=0, max = 10.0, step = .01),
             sliderInput('scale.A', 'scaling of interaction matrix', value=.5, min=0, max = 1.0, step = .0001),
             sliderInput('connectance', 'Connectance of interaction matrix', value=0.5, min=0, max = 1.0, step = .01),
             tags$div(sliderInput("mutualism", "mutualism (1,1)", min = 0, max=1, value=0.5, step = 0.1),  style="display:inline-block"),
             tags$div(sliderInput("commensalism", "commensalism (1,0)", min = 0, max=1, value=0.5, step = 0.1),  style="display:inline-block"),
             tags$div(sliderInput("parasitism", "parasitism (1,-1)", min = 0, max=1, value=0.5, step = 0.1),  style="display:inline-block"),
             tags$div(sliderInput("amensalism", "amensalism (0,-1)", min = 0, max=1, value=0.5, step = 0.1),  style="display:inline-block"),
             tags$div(sliderInput("competition", "competition (-1,-1)", min = 0, max=1, value=0.5, step = 0.1),  style="display:inline-block"),
             radioButtons("symmetric", "Use symmetric interactions?", c('yes', 'no'), selected ='no' ),
             textOutput("inter.w1"),
             verbatimTextOutput("inter.w2"),
             plotOutput("matrixHmap", width = 600, height = 400)),
  
  tabPanel("Stochasticity",
           radioButtons("stoch", "Use stochasticity?", c('yes', 'no'), selected ='no' ),
           sliderInput('sigma_drift', 'sd of random drift', value=0.001, min=0.0, max = 0.2, step = 0.001),
           tags$div(sliderInput('p.epoch', 'frequency of strong episodic drift', value=0.01, min=0.0, max = 0.5, step = 0.01), style="display:inline-block"),
           tags$div(sliderInput('sigma_epoch', 'sd of strong episodic drift', value=0.01, min=0.0, max = 0.5, step = 0.01), style="display:inline-block"),
           br(),
           tags$div(sliderInput('perturbations.n', 'number of perturbation events', value=1, min=0, max = 100, step = 1), style="display:inline-block"),
           tags$div(sliderInput('perturbations.t', 'duration of perturbation events (t)', value=0.1, min=0, max = 10, step = .1), style="display:inline-block"),
           tags$div(sliderInput('sigma.perturbation', 'sd of perturbation', value=0.01, min=0.0, max = 0.5, step = 0.01), style="display:inline-block"),
           textOutput("ptT")),
  
  tabPanel("Simulation",
           plotOutput("communityPlot", width = 600, height = 400))))),
  
)

server <- function(input, output, session) {
  nSpecies <- reactive(input$n.sp)
  alpha <- reactive(input$alpha)
  beta <- reactive(input$beta)
  sigma.A <- reactive(input$s.distA)
  scale.A <- reactive(input$scale.A)
  distribution.A <- reactive(rnorm(n=nSpecies()^2, sd = sigma.A()))
  connectance.p <- reactive(input$connectance)
  mutualism.w <-reactive(input$mutualism)
  commensalism.w <-reactive(input$commensalism)
  parasitism.w <-reactive(input$parasitism)
  amensalism.w <-reactive(input$amensalism)
  competition.w <-reactive(input$competition)
  
  eveness <- reactive(input$multiplier)
  migration.rate <- reactive(input$migration.r)
  migration.s <- reactive(input$sigma_migration)
  metacommunity.p <- reactive(rdirichlet(1, rep(.5, nSpecies())*eveness()))
  output$metacommunityPlot <- renderPlot(makePiePlot(metacommunity.p()[,]))
  interaction.w <- reactive(c(mutualism.w(), commensalism.w(), parasitism.w(), amensalism.w(),  competition.w() ))
  output$inter.w1 <- renderText("Interaction frequencies are: \n")
  output$inter.w2 <- renderPrint(sprintf('%.2f', interaction.w()/sum(interaction.w()))) 
  matrix.d <- reactive(runif(nSpecies(), max(g.rates())/3, max(g.rates())))
  
  g.rates <- reactive({rbeta(nSpecies(), shape1 = alpha(), shape2 = beta())})
  s.drift <- reactive(input$sigma_drift)
  prob.epoch <- reactive(input$p.epoch)
  s.epoch <- reactive(input$sigma_epoch)
  n.perturbs <- reactive(input$perturbations.n)
  t.perturbs <- reactive(input$perturbations.t)
  s.perturbs <- reactive(input$sigma.perturbation) 
  endTime <- reactive(input$t_end)
  
  output$betaDist <- renderText('<b>Distribution of growth rates</b> (
<a href="https://en.wikipedia.org/wiki/Beta_distribution">beta distribution</a>)')
  output$gr.rates1 <- renderText("Growth rates are: \n")
  output$gr.rates2 <- renderPrint(sprintf('%.2f', g.rates())) 
  output$gr.rates3 <- renderText("The expected distribution is:\n")
  output$gr.expect <- renderPlot({ggplot(data.frame(x = c(-4, 4)), aes(x = x)) + stat_function(fun = dbeta, args = list(shape1=alpha(), shape2=beta())) + xlim(0,1)}, res=60)
  k.capacity <- reactive({g.rates()/matrix.d()})
  output$carryingK1 <- renderText("Carrying capacities (K) are: ")
  output$carryingK2 <- renderPrint(sprintf('%.2f', k.capacity()))
  output$interactionMat <- renderText('<b>Interaction matrix</b>')
  symmetricM <- reactiveValues(data = NULL)
  observeEvent(input$symmetric, {
    if (input$symmetric=='yes'){
      symmetricM$data = TRUE
    }else{symmetricM$data = FALSE
    }})
  
  matrix.A <- reactive(randomA(n_species = nSpecies(), 
                                 diagonal = -1*matrix.d(), 
                                 connectance = connectance.p(), 
                                 interaction.w = interaction.w(), 
                                 scale = scale.A(), 
                                 distribution = distribution.A(),
                                 symmetric =  symmetricM$data))
  
  output$matrixHmap <- renderPlot(makeHeatmap(matrix.A(), "Interaction\nmatrix"))
  
  output$x0 <- renderText(paste('Initial populations are set to constant and equal to ', sp.x0))
  
  stochastic <- reactiveValues(data = NULL)
  
  observeEvent(input$stoch, {
    if (input$stoch=='yes'){
      stochastic$data = TRUE
    }else{stochastic$data = FALSE
    }
  })
  timePerturb <- reactive(getPerturbT(endTime(), n.perturbs()))
  
  output$ptT <- renderText({paste("Times of perturbations are: ", paste(signif(timePerturb(),2), collapse=", "))})
  simul <- reactive(simulateGLV(n_species = nSpecies(),
                                 A = matrix.A(),
                                 growth_rates = g.rates(),
                                 x0 = rep(sp.x0, nSpecies()),
                                 t_end = endTime(),
                                 stochastic = stochastic$data,
                                 sigma_drift = s.drift(),
                                 sigma_epoch = s.epoch(),
                                 t_external_events= timePerturb(),
                                 t_external_durations= rep(t.perturbs(), n.perturbs()),
                                 sigma_external = s.perturbs(),
                                 p.epoch = prob.epoch(),
                                sigma_migration = migration.s(),
                                migration_p = migration.rate(),
                                metacommunity.p = metacommunity.p())) 

  #
  output$communityPlot <- renderPlot(makePlot(simul()$matrix))

  
}

shinyApp(ui, server) 
