library(shiny)
library(ggplot2)
library(deSolve)
library(reshape2)
library(gtools)
library(shinyBS)

#### source all files ####
Rfiles = gsub(" ", "", paste("../R/", list.files("../R")))
sapply(Rfiles, source)

#### converting functions ####
text2char <- function(text){
    if(trimws(text) == "") {
        return(NULL)
    } else {
        return(strsplit(x = trimws(text), split = "\\,+\\s+|\\s+\\,+|\\,+|\\;+\\s+|\\s+\\;+|\\;+")[[1]])
    }
}

#### plotting functions ####
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

makePlotRes <- function(out.matrix){
    df <- as.data.frame(out.matrix)
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


#### ui ####
ui <- navbarPage(
    title ="microSimShiny",
    #### tab1 Consumer-Resource Model ####
    tabPanel(
        title = "Consumer-Resource Model", 
        fluidPage(
            fluidRow(
                "Consumer-Resource Model describes the numerical relationships 
                between microorganisms and the resources they consume.
                The idea is originated from ",
                tags$a("Chesson, P., 1990. MacArthur’s consumer-resource model. T
                       heoretical Population Biology 37, 26–38. ", 
                       href = "https://doi.org/10.1016/0040-5809(90)90025-Q"),
                ", and one of implementations (in Python) referred to is ",
                tags$a("Marsland, R., Cui, W., Goldford, J., Mehta, P., 2020. 
                       The Community Simulator: A Python package for microbial 
                       ecology. PLOS ONE 15, e0230430. ", 
                       href = "https://doi.org/10.1371/journal.pone.0230430")
            ),
            br(),
            fluidRow(
                sidebarLayout(
                    sidebarPanel(
                        tags$h3("input controls"),
                        icon("question-circle"),
                        bsTooltip(id = "someInput", title = "This is an input", 
                                  placement = "left", trigger = "hover"),
                        sliderInput("n.species", "number of species", value = 2, min = 2, max = 20),
                        sliderInput("n.resources", "number of resources", value = 4, min = 2, max = 40),
                        checkboxInput("advancedCRM", strong("show advanced parameters"), value = FALSE),
                        conditionalPanel(condition = "input.advancedCRM",
                            textInput("names.species", "names of species"),
                            textInput("names.resources", "names of resources"),
                            textAreaInput("E", "matrix of efficiency"),
                            textInput("x0", "initial abundances of species"),
                            textInput("resources", "initial concentration of resources"),
                            textInput("growth.rates", "maximum growth rates of species"),
                            textAreaInput("monod.constant", "constant of additive monod growth of n.species consuming n.resources"),
                            sliderInput("error.variance", "variance of measurement error", value = 0, min = 0, max = 10, step = 0.1),
                            checkboxInput("norm", strong("returns normalized abundances"), value = FALSE),
                            numericInput("t.end", "final time of the simulation", value = 1000, min = 100, max = 10000),
                        ),
                        actionButton("buttonSimulateCRM", "Run the model", class = "btn btn-primary")
                    ),
                    mainPanel(
                        br(),
                        tags$h3("output"),
                        plotOutput("CRMSpecies"),
                        plotOutput("CRMResources"),

                    )
                )
            )
        )),
    
    #### tab2 Generalized Lotka-Volterra (GLV) Model ####
    tabPanel(
        title = "Generalized Lotka-Volterra (GLV)",
        fluidPage(
            fluidRow(
                "Generalized Lotka-Volterra (GLV) Model is defined by a set of 
                differential equations describing the interspecies interactions."
                ),
            br(),
            fluidRow(
                sidebarLayout(
                    sidebarPanel(
                        tags$h3("input controls"),
                        tags$h4("1. Generate interspecies interactions"),
                        
                        sliderInput("n.speciesGLV", "number of species", value = 2, min = 2, max = 20),
                        checkboxInput("advancedRandomA", strong("show advanced parameters"), value = FALSE),
                        conditionalPanel(condition = "input.advancedRandomA",
                            textInput("names.speciesGLV", "names of species"),
                            sliderInput("diagonalGLV", "diagonal values of matrix A", value = -0.5, min = -2, max = 0, step = 0.1),
                            sliderInput("connectanceGLV", "connectance of matrix A", value = 0.2, min = 0, max = 1),
                            sliderInput("scaleGLV", "scale of matrix A", value = 0.1, min = 0, max = 1),
                            numericInput("mutualismGLV", "relative proportion of mutualism in matrix A", value = 1, min = 0, max = 10,step = 0.05),
                            numericInput("commensalismGLV", "relative proportion of commensalism in matrix A", value = 1, min = 0, max = 10,step = 0.05),
                            numericInput("parasitismGLV", "relative proportion of parasitism in matrix A", value = 1, min = 0, max = 10,step = 0.05),
                            numericInput("amensalismGLV", "relative proportion of amensalism in matrix A", value = 1, min = 0, max = 10,step = 0.05),
                            numericInput("competitionGLV", "relative proportion of competition in matrix A", value = 1, min = 0, max = 10,step = 0.05),
                            textInput("interactionsGLV", "interactions between species"),
                            checkboxInput("symmetricGLV", strong("whether matrix A is symmetric"), value = FALSE),
                            textInput("listAGLV", "a list of previous generated matrix"),
                        ),
                        actionButton("buttonRandomA", "generate random matrix A of interspecies interactions", class = "btn btn-primary"),
                        
                        tags$h4("2. Calculate generalized Lotka-Volterra (GLV) Model"),
                        checkboxInput("advancedGLV", strong("show advanced parameters"), value = FALSE),
                        conditionalPanel(condition = "input.advancedGLV",
                            textInput("x0GLV", "initial abundances of species"),
                            textInput("growth.ratesGLV", "maximum growth rates of species"),
                            checkboxInput("stochasticGLV", "apply stochasticity in the simulation", value = TRUE),
                            conditionalPanel(condition = "input.stochasticGLV",
                                numericInput("sigma.driftGLV", "sigma.driftGLV", value = 0.001, min = 0, max = 1, step = 0.001),
                                numericInput("sigma.epochGLV", "sigma.epochGLV", value = 0.1, min = 0, max = 1, step = 0.001),
                                numericInput("sigma.externalGLV", "sigma.externalGLV", value = 0.3, min = 0, max = 1, step = 0.001),
                                numericInput("epoch.pGLV", "epoch.pGLV", value = 0.001, min = 0, max = 1, step = 0.001),
                                
                            ),
                            numericInput("sigma.migrationGLV", "sigma.migrationGLV", value = 0.01, min = 0, max = 1, step = 0.001),
                            textInput("t.external_eventsGLV", "timepoints of external events"),
                            textInput("t.external_durationsGLV", "time durations of external events"),
                            numericInput("migration.pGLV", "migration.p", value = 0.01, min = 0, max = 1, step = 0.001),
                            textInput("metacommunity.probabilityGLV", "metacommunity.probability"),
                            sliderInput("error.varianceGLV", "variance of measurement error", value = 0, min = 0, max = 10, step = 0.1),
                            checkboxInput("normGLV", strong("returns normalized abundances"), value = FALSE),
                            numericInput("t.endGLV", "final time of the simulation", value = 1000, min = 100, max = 10000),
                        ),
                        actionButton("buttonSimulateGLV", "Run the GLV Model", class = "btn btn-primary"),
                    ),
                    mainPanel(
                        br(),
                        tags$h3("output"),
                        tableOutput("TableA"),
                        plotOutput("GLVSpecies"),
                    )
                )
            )
        )),
    
    #### tab3 Hubbell Model ####
    tabPanel(
        title = "Hubbell Model",
        fluidPage(
            
        )),
    
    #### tab4 Hubbell Model with death rates ####
    tabPanel(
        title = "Hubbell Model with death rates",
        fluidPage(
            
        )),
    
    #### tab5 Logistic Model (with stochasticity) ####
    tabPanel(
        title = "Logistic Model (with stochasticity)",
        fluidPage(
            
        )),
    
)

#### server ####
server <- function(input, output, session) {
    
    #### model1 simulate consumer resource model ####
    addTooltip(session, id = "someInput", title = "This is an input.",
               placement = "left", trigger = "hover")
    
    n.species <- reactive(input$n.species)
    n.resources <- reactive(input$n.resources)
    names.species <- reactive(text2char(input$names.species))
    names.resources <- reactive(text2char(input$names.resources))
    E <- reactive(
        ifelse(input$E == "", 
            return(NULL), 
            return(NULL))) # TODO: convert webpage multiline input to matrix
    x0 <- reactive(text2char(input$x0))
    resources <- reactive(text2char(input$resources))
    growth.rates <- reactive(text2char(input$growth.rates))
    monod.constant <- reactive(text2char(input$monod.constant)) # TODO: convert
    error.variance <- reactive(input$error.variance)
    norm <- reactive(input$norm)
    t.end <- reactive(input$t.end)

    runCRM <- eventReactive(input$buttonSimulateCRM, {
        simulateConsumerResource(
            n.species = n.species(),
            n.resources = n.resources(),
            names.species = names.species(), 
            names.resources = names.resources(), 
            E = E(), 
            x0 = x0(), 
            resources = resources(), 
            growth.rates = growth.rates(), 
            monod.constant = monod.constant(), 
            error.variance = error.variance(), 
            norm = norm(), 
            t.end = t.end()
        )
    })
    
    output$CRMSpecies <- renderPlot(makePlot(runCRM()$matrix))
    output$CRMResources <- renderPlot(makePlot(runCRM()$resources))
    
    #### model2 simulate generalized Lotka-Volterra Model ####
    n.speciesGLV <- reactive(input$n.speciesGLV)
    names.speciesGLV <- reactive(text2char(input$names.speciesGLV))
    diagonalGLV <- reactive(input$diagonalGLV)
    connectanceGLV <- reactive(input$connectanceGLV)
    scaleGLV <- reactive(input$scaleGLV)
    mutualismGLV <- reactive(input$mutualismGLV)
    commensalismGLV <- reactive(input$commensalismGLV)
    parasitismGLV <- reactive(input$parasitismGLV)
    amensalismGLV <- reactive(input$amensalismGLV)
    competitionGLV <- reactive(input$competitionGLV)
    interactionsGLV <- reactive(text2char(input$interactionsGLV))
    symmetricGLV <- reactive(input$symmetricGLV)
    listAGLV <- reactive(text2char(input$listAGLV)) # TODO: convert
    
    generateA <- eventReactive(input$buttonRandomA, {
        randomA(
            n.species = n.speciesGLV(),
            names.species = names.speciesGLV(), 
            diagonal = diagonalGLV(),
            connectance = connectanceGLV(),
            scale = scaleGLV(),
            mutualism = mutualismGLV(),
            commensalism = commensalismGLV(),
            parasitism = parasitismGLV(),
            amensalism = amensalismGLV(),
            competition = competitionGLV(),
            interactions = interactionsGLV(),
            symmetric = symmetricGLV(),
            listA = listAGLV()
        )
    })
    
    output$TableA <- renderTable(generateA())
    
    x0GLV <- reactive(text2char(input$x0GLV))
    growth.ratesGLV <- reactive(text2char(input$growth.ratesGLV))
    stochasticGLV <- reactive(input$stochasticGLV)
    sigma.driftGLV <- reactive(input$sigma.driftGLV)
    sigma.epochGLV <- reactive(input$sigma.epochGLV)
    sigma.externalGLV <- reactive(input$sigma.externalGLV)
    epoch.pGLV <- reactive(input$epoch.pGLV)
    sigma.migrationGLV <- reactive(input$sigma.migrationGLV)
    t.external_eventsGLV <- reactive(as.numeric(text2char(input$t.external_eventsGLV)))
    t.external_durationsGLV <- reactive(as.numeric(text2char(input$t.external_durationsGLV)))
    migration.pGLV <- reactive(input$migration.pGLV)
    metacommunity.probabilityGLV <- reactive(text2char(input$metacommunity.probabilityGLV))
    error.varianceGLV <- reactive(input$error.varianceGLV)
    normGLV <- reactive(input$normGLV)
    t.endGLV <- reactive(input$t.endGLV)
    
    runGLV <- eventReactive(input$buttonSimulateGLV, {
        simulateGLV(
            n.species = n.speciesGLV(),
            names.species = names.speciesGLV(), 
            A = generateA(),
            x0 = x0GLV(), 
            growth.rates = growth.ratesGLV(), 
            sigma.drift = sigma.driftGLV(),
            sigma.epoch = sigma.epochGLV(),
            sigma.external = sigma.externalGLV(),
            sigma.migration = sigma.migrationGLV(),
            epoch.p = epoch.pGLV(),
            t.external_events = t.external_eventsGLV(),
            t.external_durations = t.external_durationsGLV(),
            stochastic = stochasticGLV(),
            migration.p = migration.pGLV(),
            metacommunity.probability = metacommunity.probabilityGLV(),
            error.variance = error.varianceGLV(),
            norm = normGLV(), 
            t.end = t.endGLV()
        )
    })
    
    output$GLVSpecies <- renderPlot(makePlot(runGLV()$matrix))

}

shinyApp(ui, server)