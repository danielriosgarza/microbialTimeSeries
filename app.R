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
library(reshape2)
library(gtools)

Rfiles = gsub(" ", "", paste("./R/", list.files("./R")))
sapply(Rfiles, source)

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

makeHeatmap <-function(out.matrix, midpoint_color){
    out.matrix = t(out.matrix[,colnames(out.matrix)!="time"])
    df = melt(out.matrix)
    names(df)<- c("x", "y", "abundance")
    df$y <- factor(df$y, levels=rev(unique(sort(df$y))))
    fig <- ggplot(df, aes(x,y,fill=abundance)) + 
        geom_tile(colour="#edfaf9",size=0.005) +
        theme(axis.title = element_blank()) + 
        scale_fill_gradient2('abundance', low = "white", high = "magenta3", midpoint = midpoint_color) +  
        theme_void() +
        scale_y_discrete(expand=c(0,0))
    
    if (ncol(out.matrix)<21){
        fig <- fig + geom_text(aes(label = round(abundance, 1)))
    }
    fig
}

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Hubbell neutral model"),
    br(),
    sidebarLayout(
        sidebarPanel(width = 2,
    sliderInput('n.sp', 'Number of species in metacommunity', value=5, min=1, max = 100, step = 1),
    
    sliderInput("multiplier", "Metacommunity eveness", min = 1, max=1000, value=1.0, step = 1),
    
    sliderInput('n.individuals', 'Number of individuals in local community', value=1000, min=100, max = 10000, step = 10),
    
    br(),
    
    
    
    br(),
    
    
    
    sliderInput("migration.r", "Migration rate", min = 0, max=1, value=.1, step = 0.01),
    
    sliderInput("events", "Simulation events", min = 1, max=100, value=1, step = 1),
    
    sliderInput("steps", "Simulation steps", min = 1000, max=100000, value=1000, step = 10),
    br(),
    br()
        ),
    
    
    
    mainPanel(tabsetPanel(tabPanel("metacommunity",
                     
    
    plotOutput("metacommunityPlot", width = 600, height = 400)),
    
    tabPanel("Initial local community",
    plotOutput("localcommunityPlot", width = 600, height = 400)),
    
    tabPanel("Simulation",
    plotOutput("simulationPlot", width = 600, height = 1000)),
    tabPanel("growth rates",
             radioButtons("gr.bool", "Use growth rates?", c('yes', 'no'), selected ='no' ),
             
             br(),
             htmlOutput("betaDist1"),
             htmlOutput("betaDist"),
             br(),
             tags$div(sliderInput("alpha", "alpha", min = 0, max=10, value=0.5, step = 0.1),  style="display:inline-block"),
             tags$div(sliderInput("beta", "beta", min = 0, max=10, value=0.5, step = 0.1),  style="display:inline-block"),
             textOutput("gr.rates1"),
             verbatimTextOutput("gr.rates2"),
             textOutput("gr.rates3"),
             plotOutput("gr.expect", width = 150, height = 150)))
    
    
    )))

# Define server logic required to draw a histogram
server <- function(input, output) {
    nSpecies <- reactive(input$n.sp)
    nIndividuals <- reactive(input$n.individuals)
    
    eveness <- reactive(input$multiplier)
    
    migration.rate <- reactive(input$migration.r)
    
    events <- reactive(input$events)
    
    steps <- reactive(input$steps)
    
    alpha <- reactive(input$alpha)
    beta <- reactive(input$beta)
    
    
    metacommunity.p <- reactive(rdirichlet(1, rep(.5, nSpecies())*eveness()))
    
    initialComposition <- reactive(rmultinom(1, nIndividuals(), rdirichlet(1, rep(1,nSpecies()))))
    
    
    output$metacommunityPlot <- renderPlot(makePiePlot(metacommunity.p()[,]))
    output$localcommunityPlot <- renderPlot(makePiePlot(initialComposition()[,], "Local\ncommunity", "Initial composition\nlocal community"))
    
    g.rates <- reactive({rbeta(nSpecies(), shape1 = alpha(), shape2 = beta())})
    output$betaDist <- renderText('<b>Distribution of growth rates</b> (
<a href="https://en.wikipedia.org/wiki/Beta_distribution">beta distribution</a>) if used')
    output$gr.rates1 <- renderText("Growth rates are: \n")
    output$gr.rates2 <- renderPrint(sprintf('%.2f', g.rates())) 
    output$gr.rates3 <- renderText("The expected distribution is:\n")
    output$gr.expect <- renderPlot({ggplot(data.frame(x = c(-4, 4)), aes(x = x)) + stat_function(fun = dbeta, args = list(shape1=alpha(), shape2=beta())) + xlim(0,1)}, res=60)
    
    observeEvent(input$gr.bool, {
        if (input$gr.bool=='yes'){
            simul <- reactive(simulateHubbellRates(community.initial = initialComposition()[,],
                                              migration.p = migration.rate(),
                                              metacommunity.p = metacommunity.p(),
                                              k.events = events(),
                                              growth.rates = g.rates(),
                                              norm = FALSE, 
                                              t.end = steps()))
            output$simulationPlot <- renderPlot(makePlot(simul()$matrix))
            
        }else{
            simul <- reactive(simulateHubbell(community.initial = initialComposition()[,],
                                              migration.p = migration.rate(),
                                              metacommunity.p = metacommunity.p(),
                                              k.events = events(),
                                              norm = FALSE, 
                                              t.end = steps(), 
                                              t.store = 100))
            
            output$simulationPlot <- renderPlot(makeHeatmap(simul()$matrix, max(initialComposition())/8))
            
            
        }
    })
    
    
    
    
    
    

    
}

# Run the application 
shinyApp(ui = ui, server = server)
