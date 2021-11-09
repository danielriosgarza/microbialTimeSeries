library(shiny)
library(ggplot2)
library(deSolve)
library(reshape2)
library(gtools)

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


getPerturbT <- function(endTime, n.perturbs){
  st <- SimulationTimes(t.end = endTime, t.store = n.perturbs + 1)
  return (st$t.sys[st$t.index[2:length(st$t.index)]])
}
#constants
sp.x0 <- 0.001

Rfiles = gsub(" ", "", paste("./R/", list.files("./R")))
sapply(Rfiles, source)


ui <- fluidPage(
  
  
  
    titlePanel("Logistic Model"),
    
    sidebarLayout(
      sidebarPanel(width = 2,
    
    sliderInput('n.sp', 'number of species', value=5, min=1, max = 100, step = 1),
    
    br(),
    textOutput("x0"),
    
    br(),
    
    sliderInput('n.dr', 'constant death rates', value=0.001, min=0.00, max = 0.1, step = 0.001),
    textOutput("deathRates"),
    
    
    br(),
    sliderInput('t.end', 'simulation time', value=100, min=10, max = 2000, step = 10)),
    
    mainPanel(tabsetPanel(tabPanel("Metacommunity",
                                   sliderInput("multiplier", "Eveness", min = 1, max=1000, value=1.0, step = 0.1),
                                   
                                   tags$div(sliderInput("migration.r", "Migration rate", min = 0, max=1, value=0.1, step = 0.01),  style="display:inline-block"),
                                   
                                   tags$div(sliderInput("sigma.migration", "sd of migration", min = 0, max=1, value=0.01, step = 0.001),  style="display:inline-block"),
                                   
                                   
                                   plotOutput("metacommunityPlot", width = 600, height = 400)),
                          
                          tabPanel("Growth rates", 
                                   htmlOutput("betaDist"),
                                   br(),
                                   tags$div(sliderInput("alpha", "alpha", min = 0, max=10, value=0.5, step = 0.1),  style="display:inline-block"),
                                   tags$div(sliderInput("beta", "beta", min = 0, max=10, value=0.5, step = 0.1),  style="display:inline-block"),
                                   textOutput("gr.rates"),
                                   plotOutput("gr.expect", width = 150, height = 150),
                                   br(),
                                   textOutput("carryingK")),
                          tabPanel("Stochasticity",
                                   radioButtons("stoch", "Use stochasticity?", c('yes', 'no'), selected ='no' ),
                                   sliderInput('sigma.drift', 'sd of random drift', value=0.001, min=0.0, max = 0.2, step = 0.001),
                                   tags$div(sliderInput('p.epoch', 'frequency of strong episodic drift', value=0.01, min=0.0, max = 0.5, step = 0.01), style="display:inline-block"),
                                   tags$div(sliderInput('sigma.epoch', 'sd of strong episodic drift', value=0.01, min=0.0, max = 0.5, step = 0.01), style="display:inline-block"),
                                   br(),
                                   tags$div(sliderInput('perturbations.n', 'number of perturbation events', value=1, min=0, max = 100, step = 1), style="display:inline-block"),
                                   tags$div(sliderInput('perturbations.t', 'duration of perturbation events (t)', value=0.1, min=0, max = 10, step = .1), style="display:inline-block"),
                                   tags$div(sliderInput('sigma.perturbation', 'sd of perturbation', value=0.01, min=0.0, max = 0.5, step = 0.01), style="display:inline-block"),
                                   textOutput("ptT")),
                          tabPanel("Simulation",
                                   plotOutput("communityPlot", width = 600, height = 400)
                                   )))))



server <- function(input, output, session) {
  nSpecies <- reactive(input$n.sp)
  eveness <- reactive(input$multiplier)
  migration.rate <- reactive(input$migration.r)
  
  alpha <- reactive(input$alpha)
  beta <- reactive(input$beta)
  
  d.rates <- reactive(as.numeric(rep(input$n.dr, nSpecies())))
  g.rates <- reactive({rbeta(nSpecies(), shape1 = alpha(), shape2 = beta())})
  metacommunity.p <- reactive(rdirichlet(1, rep(.5, nSpecies())*eveness()))
  output$metacommunityPlot <- renderPlot(makePiePlot(metacommunity.p()[,]))
  s.migration <- reactive(input$sigma.migration)
  
  s.drift <- reactive(input$sigma.drift)
  prob.epoch <- reactive(input$p.epoch)
  s.epoch <- reactive(input$sigma.epoch)
  
  n.perturbs <- reactive(input$perturbations.n)
  
  t.perturbs <- reactive(input$perturbations.t)
  
  s.perturbs <- reactive(input$sigma.perturbation) 
  
  endTime <- reactive(input$t.end)
  output$betaDist <- renderText('<b>Distribution of growth rates</b> (
<a href="https://en.wikipedia.org/wiki/Beta_distribution">beta distribution</a>)')
  
  output$gr.rates <- renderText({paste("Growth rates are: ", paste(signif(g.rates(),2), collapse=", "), ". The expected distribution is:")})
     
  output$gr.expect <- renderPlot({ggplot(data.frame(x = c(-4, 4)), aes(x = x)) + stat_function(fun = dbeta, args = list(shape1=alpha(), shape2=beta())) + xlim(0,1)}, res=60)
  # 
  k.capacity <- reactive({g.rates()/runif(nSpecies(), max(g.rates())/3, max(g.rates()))})
  # 
  output$carryingK <- renderText(paste("Carrying capacities (K) are: ", paste(signif(k.capacity(),2), collapse=", ")))
  #   
  

  output$deathRates <- renderText(paste("Death rates are: ", paste(d.rates(), collapse=", ")))

  
  # 
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
  simul <- reactive({simulateStochasticLogistic(n.species = nSpecies(),
                                                growth.rates = g.rates(),
                                                carrying.k = k.capacity(),
                                                x0 = rep(sp.x0, nSpecies()),
                                                death.rates = d.rates(),
                                                t.end = endTime(),
                                                t.store = endTime(),
                                                stochastic = stochastic$data,
                                                sigma.drift = s.drift(),
                                                sigma.epoch = s.epoch(),
                                                sigma.migration = s.migration(),
                                                migration.p = migration.rate(),
                                                metacommunity.p = metacommunity.p(),
                                                t.external_events= timePerturb(),
                                                t.external_durations= rep(t.perturbs(), n.perturbs()),
                                                sigma.external = s.perturbs(),
                                                p.epoch = prob.epoch())}) 
   
  output$communityPlot <- renderPlot(makePlot(simul()$matrix))
}


shinyApp(ui, server)
