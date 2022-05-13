#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

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

makePlotRes <- function(out_matrix){
    df <- as.data.frame(out_matrix)
    dft <-  melt(df, id="time")
    names(dft)[2] = "resources"
    names(dft)[3] = "S.t"
    lgd = TRUE
    if (ncol(df)>10){
        lgd = FALSE
    }
    ggplot(dft, aes(time, S.t, col = resources)) + geom_line(show.legend = lgd, lwd=0.5)
    
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
        scale_fill_gradient2('strength', low = "red", mid = "white", high = "blue", midpoint = 0)+
        theme_void() + ggtitle(title) 
    
    if (ncol(matrix.A)<21 & nrow(matrix.A)<21){
        fig <- fig + geom_text(aes(label = round(strength, 1)))
    }
    
    fig <- fig +  labs(x = "species", y = "resources")+
        theme(
            axis.title.x = element_text(size = 14, face = "bold.italic"),
            axis.title.y = element_text(size = 14, face = "bold.italic")
        )
    fig
}

#constants
sp.x0 <- 0.001

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Consumer resource model"),
    sidebarLayout(
      sidebarPanel(width = 2,
                   sliderInput('n.sp', 'number of species', value=5, min=1, max = 100, step = 1),
                   br(),
                   sliderInput('n.res', 'number of resources', value=5, min=1, max = 100, step = 1),
                   
                   textOutput("x0"),
                   br(),
                   
                   sliderInput('monod.k', 'Monod.K multiplier', value=1, min=1, max = 10, step = 1),
                   
                   br(), 
                   sliderInput('t_end', 'simulation time', value=1000, min=10, max = 10000, step = 10)),
      mainPanel(tabsetPanel(tabPanel("Growth rates",
                                     htmlOutput("betaDist"),
                                     br(),
                                     tags$div(sliderInput("alpha", "alpha", min = 0, max=10, value=0.5, step = 0.1),  style="display:inline-block"),
                                     tags$div(sliderInput("beta", "beta", min = 0, max=10, value=0.5, step = 0.1),  style="display:inline-block"),
                                     textOutput("gr.rates1"),
                                     verbatimTextOutput("gr.rates2"),
                                     textOutput("gr.rates3"),
                                     plotOutput("gr.expect", width = 150, height = 150),
                                     ),
                            tabPanel("Resources",
                                     tags$div(sliderInput("cons.w", "consumption weight", min = 0, max=1, value=1, step = 0.1),  style="display:inline-block"),
                                     tags$div(sliderInput("prod.w", "production weight", min = 0, max=1, value=1, step = 0.1),  style="display:inline-block"),
                                     tags$div(sliderInput("maint.w", "maintenance", min = 0, max=1, value=.2, step = 0.1),  style="display:inline-block"),
                                     br(),
                                     br(),
                                     plotOutput("matrixHmap", width = 600, height = 400),
                                     br(),
                                     br(),
                                     tags$div(sliderInput("r.conc", "resource concentration", min = 0, max=1000, value=100, step = 1),  style="display:inline-block"),
                                     tags$div(sliderInput("r.even", "resource evenness", min = 1, max=1000, value=1, step = 1),  style="display:inline-block"),
                                     plotOutput("resourcePlot", width = 600, height = 400)),
                              
                              tabPanel("Simulations",
                                       br(),
                                       plotOutput("communityPlot", width = 600, height = 400),
                                       br(),
                                       plotOutput("resourceDynPlot", width = 600, height = 400)
    

  
)))))

# Define server logic required to draw a histogram
server <- function(input, output) {
    nSpecies <- reactive(input$n.sp)
    nResources <- reactive(input$n.res)
    alpha <- reactive(input$alpha)
    beta <- reactive(input$beta)
    consumption.w <- reactive(input$cons.w)
    production.w <- reactive(input$prod.w)
    maintenance <- reactive(input$maint.w)
    evenness  <- reactive(input$r.even)
    concentration <- reactive(input$r.conc)
    monod <- reactive(input$monod.k)
    t_end <- reactive(input$t_end)
    mean.prod <- reactive(production.w() * nResources())
    mean.con <- reactive(consumption.w() * nResources())
    
    resource_dist <- reactive(rdirichlet(1, rep(.5, nResources())*evenness()))
    
    resources = reactive(resource_dist()*concentration()*nResources())
    
    output$resourcePlot <- renderPlot(makePiePlot(resources()[,], label = 'concentration', title='resources'))
    
    output$x0 <- renderText(paste('Initial populations are set to constant and equal to ', sp.x0))
    g.rates <- reactive({rbeta(nSpecies(), shape1 = alpha(), shape2 = beta())})
    output$betaDist <- renderText('<b>Distribution of growth rates</b> (
<a href="https://en.wikipedia.org/wiki/Beta_distribution">beta distribution</a>)')
    output$gr.rates1 <- renderText("Growth rates are: \n")
    output$gr.rates2 <- renderPrint(sprintf('%.2f', g.rates())) 
    output$gr.rates3 <- renderText("The expected distribution is:\n")
    output$gr.expect <- renderPlot({ggplot(data.frame(x = c(-4, 4)), aes(x = x)) + stat_function(fun = dbeta, args = list(shape1=alpha(), shape2=beta())) + xlim(0,1)}, res=60)
    matrix.E <- reactive(randomE(n_species = nSpecies(),
                                 n_resources = nResources(),
                                 mean.con = mean.con(),
                                 mean.prod = mean.prod(),
                                 maintenance = maintenance()))
    
    monod.matrix <- reactive(matrix(rgamma(n=nSpecies()*nResources(), shape = monod()*max(resources()),
                                             rate = 1), nrow = nSpecies()))
    
    output$matrixHmap <- renderPlot(makeHeatmap(matrix.E(), 'Consumption/production\nmatrix'))
    
    simul <- reactive(simulateConsumerResource(n_species = nSpecies(),
                                               n_resources = nResources(), 
                                               eff = matrix.E(),
                                               x0 =rep(sp.x0, nSpecies()),
                                               resources = resources(),
                                               growth_rates = g.rates(),
                                               monod.k = monod.matrix(),
                                               norm = FALSE,
                                               t_end = t_end()))
    
    output$communityPlot <- renderPlot(makePlot(simul()$matrix))
    
    output$resourceDynPlot <- renderPlot(makePlotRes(simul()$resources))

    
}

# Run the application 
shinyApp(ui = ui, server = server)
