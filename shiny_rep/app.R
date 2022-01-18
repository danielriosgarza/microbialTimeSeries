library(shiny)
library(shinyBS) # for tooltips on mouse floating
library(ggplot2)
library(deSolve)
library(reshape2)
library(gtools)
library(DT) # for formatting data tables in HTML pages
library(formattable)
library(dplyr)

#### source all files ####
Rfiles = gsub(" ", "", paste("../R/", list.files("../R")))
sapply(Rfiles, source)

#### converting functions ####
text2char <- function(text){
    # split words or numbers by separators(',' and ';')
    if(trimws(text) == "") {
        return(NULL)
    } else {
        return(strsplit(x = trimws(text), split = "\\,+\\s+|\\s+\\,+|\\,+|\\;+\\s+|\\s+\\;+|\\;+")[[1]])
    }
}

text2chars <- function(text, len, prefix = NULL, expr = NULL){
    # split words or numbers by separators(',' and ';')
    # if length not enough, generate new values using expr.
    # used for automatic names of species/compounds(resources) with prefix
    # used for random numbers of initial abundances / growth rates with expr
    text <- text2char(text)
    if (length(text) < len){
        if (is.null(prefix) & is.null(expr)){
            stop("'prefix' or 'expr' not provided to function 'text2chars'")
        } else if (is.null(expr)){
            return(c(text, paste0(prefix, seq_len(len))[(length(text)+1):len]))
        } else if (is.null(prefix)){
            return(c(text, eval(parse(text = expr)))[1:len])
        }
    } else if (length(text) > len){
        warning("length of text provided to 'textchars' more than needed")
        return(text[1:len])
    } else {
        return(text)
    }
}

#### plotting functions ####
makePlot <- function(out.matrix, title){
    df <- as.data.frame(out.matrix)
    dft <-  melt(df, id="time")
    names(dft)[2] = "species"
    names(dft)[3] = "x.t"
    lgd = ncol(df)<= 20
    ggplot(dft, aes(time, x.t, col = species)) +
        geom_line(show.legend = lgd, lwd=0.5) +
        ggtitle(title) + 
        theme_linedraw() +
        theme(plot.title = element_text(hjust = 0.5, size = 14))
}

makePlotRes <- function(out.matrix, title){
    df <- as.data.frame(out.matrix)
    dft <-  melt(df, id="time")
    names(dft)[2] = "resources"
    names(dft)[3] = "S.t"
    lgd = ncol(df)<= 20
    ggplot(dft, aes(time, S.t, col = resources)) + 
        geom_line(show.legend = lgd, lwd=0.5) + 
        ggtitle(title) + 
        theme_linedraw() + 
        theme(plot.title = element_text(hjust = 0.5, size = 14))
}

makePiePlot <- function(multinomdist, label = 'Meta\ncommunity', title = "Metacommunity \nspecies abundance\n"){
    df <- data.frame(group = seq(length(multinomdist)), probability = multinomdist)
    fig <- ggplot(df, aes(x=group,y=1,fill=probability, )) + 
        geom_tile(colour="#edfaf9",size=0.005) +
        theme(axis.title = element_blank()) + 
        scale_fill_gradient2(label, low = "white", high = "magenta3", midpoint = max(multinomdist)/8) +  
        theme_void() +
        coord_fixed(ratio = length(multinomdist)/4) +
        ggtitle(title) + theme_linedraw()
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
    
    fig <- fig + labs(x = "species", y = "compounds")+
        theme_linedraw() + 
        theme(
            # axis.title.x = element_text(size = 14, face = "bold.italic"),
            # axis.title.y = element_text(size = 14, face = "bold.italic"),
            plot.title = element_text(hjust = 0.5, size = 14)
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
            titlePanel("Consumer-Resource Model"),
            bsCollapse(id = "contents", 
                       open = "Model", 
                       bsCollapsePanel(title = strong("Model"), value = "Model",
                                       fluidRow(
                                           column(width = 5,
                                                  wellPanel(
                                                      tabsetPanel(
                                                          header = tags$style(HTML("
                                        /* add a border for tabpanes */
                                        .tabbable > .tab-content > .tab-pane {
                                            border-left: 1px solid #ddd;
                                            border-right: 1px solid #ddd;
                                            border-bottom: 1px solid #ddd;
                                            border-radius: 0px 0px 5px 5px;
                                            padding: 10px;
                                        }
                                        .nav-tabs {
                                            margin-bottom: 0;
                                        }
                                    ")),
                                                          tabPanel("Basic",
                                                                   sliderInput("nSpecies", "number of species", value = 2, min = 2, max = 100),
                                                                   bsTooltip("nSpecies", "Number of species in the simulation", "right", options = list(container = "body")),
                                                                   sliderInput("nResources", "number of compounds", value = 4, min = 2, max = 100),
                                                                   bsTooltip("nResources", "Number of compounds in the simulation", "right", options = list(container = "body")),
                                                                   hr(),
                                                                   checkboxInput("changeNamesCRM", strong("custom names of species/compounds"), value = FALSE),
                                                                   conditionalPanel(condition = "input.changeNamesCRM",
                                                                                    helpText("Custom names separate by ',' or ';' (and spaces) will replace default names."),
                                                                                    textInput("namesSpecies", "names of species"),
                                                                                    textInput("namesResources", "names of compounds"),
                                                                   ),
                                                          ),
                                                          
                                                          tabPanel("Compounds",
                                                                   sliderInput("meanConsumption", "consumption weight", value = 0.4, min = 0, max = 1),
                                                                   bsTooltip("meanConsumption", "Mean proportion of compounds consumed by each species.", "right", options = list(container = "body")),
                                                                   sliderInput("meanProduction", "production weight", value = 0.2, min = 0, max = 0.5),
                                                                   bsTooltip("meanProduction", "Mean proportion of compounds produced by each species.", "right", options = list(container = "body")),
                                                                   sliderInput("maintenance", "maintenance weight", value = 0.5, min = 0, max = 1),
                                                                   bsTooltip("maintenance", "How much compounds were used to maintain the microbial community (not involved in further calculation of flux).", "right", options = list(container = "body")),
                                                                   
                                                                   tags$label("Compounds Stochiometry"),
                                                                   dataTableOutput("tableE", width = "100%"),
                                                                   bsTooltip("tableE", "Stochiometric values of consumption and production of compounds by each cell. Positive efficiencies indicate the consumption of resources, whilst negatives indicate that the species would produce the resource.", "right", options = list(container = "body")),
                                                          ),
                                                          
                                                          tabPanel("Growth rates",
                                                                   textInput("x0", "initial abundances of species"),
                                                                   bsTooltip("x0", "If the given initial abundances of species is not enough, random initial abundances will be added.", "right", options = list(container = "body")),
                                                                   verbatimTextOutput("x0Output"),
                                                                   textInput("resources", "initial concentration of compounds"),
                                                                   bsTooltip("resources", "If the given initial concentrations of compounds are not enough, random values will be added.", "right", options = list(container = "body")),
                                                                   verbatimTextOutput("resourcesOutput"),
                                                                   tags$label("Distribution of Growth Rates"),
                                                                   br(),
                                                                   actionButton("buttonBetaEven", "-", class = "btn btn-primary"),
                                                                   bsTooltip("buttonBetaEven", "even distribution"),
                                                                   actionButton("buttonBetaRidge", "^", class = "btn btn-primary"),
                                                                   bsTooltip("buttonBetaRidge", "normal-alike distribution"),
                                                                   actionButton("buttonBetaValley", "˅", class = "btn btn-primary"),
                                                                   bsTooltip("buttonBetaValley", "U shape distribution"),
                                                                   actionButton("buttonBetaLeft", "◟", class = "btn btn-primary"),
                                                                   bsTooltip("buttonBetaLeft", "left skewed distribution"),
                                                                   actionButton("buttonBetaRight", "◞", class = "btn btn-primary"),
                                                                   bsTooltip("buttonBetaRight", "right skewed distribution"),
                                                                   actionButton("buttonBetaLeftTriangle", "\\", class = "btn btn-primary"),
                                                                   bsTooltip("buttonBetaLeftTriangle", "left-triangle distribution"),
                                                                   actionButton("buttonBetaRightTriangle", "/", class = "btn btn-primary"),
                                                                   bsTooltip("buttonBetaRightTriangle", "right-triangle distribution"),
                                                                   
                                                                   sliderInput("alpha", "alpha", value = 1, min = 0, max = 10, step = 0.1),
                                                                   bsTooltip("alpha", "first parameter of beta distribution", "right", options = list(container = "body")),
                                                                   sliderInput("beta", "beta", value = 1, min = 0, max = 10, step = 0.1),
                                                                   bsTooltip("beta", "second parameter of beta distribution", "right", options = list(container = "body")),
                                                                   
                                                                   textInput("growthRates", "maximum growth rates of species"),
                                                                   bsTooltip("growthRates", "If the given growth rates are not enough, random values will be added.", "right", options = list(container = "body")),
                                                                   verbatimTextOutput("growthRatesOutput"),
                                                                   
                                                                   plotOutput("growthRatesDist"),
                                                                   textAreaInput("monodConstant", "constant of additive monod growth of species consuming compounds"),
                                                          ),
                                                          tabPanel("Pertubations",
                                                                   sliderInput("errorVariance", "variance of measurement error", value = 0, min = 0, max = 10, step = 0.1),
                                                                   bsTooltip("errorVariance", "The variance of measurement error. By default it equals to 0, indicating that the result won't contain any measurement error.", "right", options = list(container = "body")),
                                                                   checkboxInput("norm", strong("returns normalized abundances"), value = FALSE),
                                                                   numericInput("tEnd", "final time of the simulation", value = 1000, min = 100, max = 10000, step = 100),
                                                                   bsTooltip("tEnd", "The end time of the simulation.", "right", options = list(container = "body")),
                                                                   # actionButton("buttonSimulateCRM", "Run the model", class = "btn btn-primary")
                                                          ),
                                                      ),
                                                  )
                                           ),
                                           column(width = 7, 
                                                  br(),
                                                  plotOutput("CRMSpecies"),
                                                  plotOutput("CRMResources"),
                                                  plotOutput("CRMPlotE"),
                                                  bsTooltip("CRMPlotE", "Positive values indicate consumption, and negative values indicate production", "left", options = list(container = "body")),
                                                  
                                           ),
                                       )
                       ),
                       bsCollapsePanel(strong("Description"), value = "Description",
                                       fluidRow(
                                           column(width = 12,
                                                  "Consumer-Resource Model describes the numerical relationships between microorganisms and the resources they consume. The idea is originated from ",
                                                  tags$a("Chesson, P., 1990. MacArthur’s consumer-resource model. Theoretical Population Biology 37, 26–38. ", 
                                                         href = "https://doi.org/10.1016/0040-5809(90)90025-Q"),
                                                  ", and one of implementations (in Python) referred to is ",
                                                  tags$a("Marsland, R., Cui, W., Goldford, J., Mehta, P., 2020. The Community Simulator: A Python package for microbial ecology. PLOS ONE 15, e0230430. ", 
                                                         href = "https://doi.org/10.1371/journal.pone.0230430"),
                                           ),
                                       )
                       ),
                       bsCollapsePanel(strong("Inputs"), value = "Inputs",
                                       "This is a panel of input table."),
                       bsCollapsePanel(strong("References"), value = "References",
                                       "Panel of refs.")
            ),
            
        )
    ),
    
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
                        tags$h2("input controls"),
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
                        tags$h2("output"),
                        tableOutput("TableA"),
                        plotOutput("GLVSpecies"),
                    )
                )
            )
        )),
    
    #### tab3 Hubbell Model ####
    tabPanel(
        title = "Hubbell Model (with death rates)",
        fluidPage(
            
        )),
    
    #### tab4 Logistic Model (with stochasticity) ####
    tabPanel(
        title = "Logistic Model (with stochasticity)",
        fluidPage(
            
        )),
)

#### server ####
server <- function(input, output, session) {
    
    #### model1 simulate consumer resource model ####
    ## basic
    n.species <- reactive(input$nSpecies)
    n.resources <- reactive(input$nResources)
    names.species <- reactive(text2chars(input$namesSpecies, len = n.species(), prefix = "sp"))
    names.resources <- reactive(text2chars(input$namesResources, len = n.resources(), prefix = "res"))
    
    ## compounds stochiometry
    mean.consumption <- reactive(input$meanConsumption)
    mean.production <- reactive(input$meanProduction)
    observeEvent(input$meanConsumption | input$meanProduction, {
        updateSliderInput(inputId = "meanProduction", max = 1 - input$meanConsumption)
        updateSliderInput(inputId = "meanConsumption", max = 1 - input$meanProduction)
    })
    maintenance <- reactive(input$maintenance)
    
    # editable table
    E <- reactiveValues(df = NULL)
    observe({
        roundE <- round(randomE(n.species(), 
                                n.resources(), 
                                names.species(), 
                                names.resources(),
                                mean.consumption = as.integer(mean.consumption() * n.resources()),
                                mean.production = as.integer(mean.production() * n.resources()),
                                maintenance = maintenance()
        ), digits = 3)
        E$df <- roundE
    })
    
    # editable table (without color heatmap)
    output$tableE <- renderDataTable(E$df, editable = 'cell', selection = 'none', server = TRUE, options = list(scrollX = TRUE))
    
    output$CRMPlotE <- renderPlot(makeHeatmap(t(E$df), 'Consumption/production matrix'), res = 96)
    
    observeEvent(input$tableE_cell_edit, {
        E$df <<- editData(E$df, input$tableE_cell_edit, 'tableE')
    })
    
    
    ## growth rates
    x0 <- reactive(as.numeric(text2chars(input$x0, len = n.species(), expr = paste0("runif(n = ", n.species() ,", min = 0.1, max = 10)"))))
    output$x0Output <- renderPrint(x0())
    resources <- reactive(as.numeric(text2chars(input$resources, len = n.resources(), expr = paste0("runif(n = ", n.resources() ,", min = 1, max = 100)"))))
    output$resourcesOutput <- renderPrint(resources())
    
    alpha <- reactive(input$alpha)
    beta <- reactive(input$beta)
    # listening to the changes in n.species, alpha, and beta
    observeEvent(input$nSpecies | input$alpha | input$beta, {
        growth.rates <- reactive(as.numeric(text2chars(input$growthRates, len = n.species(), expr = paste0('round(rbeta(',n.species(), ',' ,alpha(), ',' , beta(),'), digits = 3)'))))
    })
    observeEvent(input$nSpecies, {
        x0 <- reactive(as.numeric(text2chars(input$x0, len = n.species(), expr = paste0("runif(n = ", n.species() ,", min = 0.1, max = 10)"))))
    })
    observeEvent(input$nResources, {
        resources <- reactive(as.numeric(text2chars(input$resources, len = n.resources(), expr = paste0("runif(n = ", n.resources() ,", min = 1, max = 100)"))))
    })
    observeEvent(input$buttonBetaEven, {
        updateSliderInput(inputId = "alpha", value = 1)
        updateSliderInput(inputId = "beta", value = 1)
    })
    observeEvent(input$buttonBetaRidge, {
        updateSliderInput(inputId = "alpha", value = 4)
        updateSliderInput(inputId = "beta", value = 4)
    })
    observeEvent(input$buttonBetaValley, {
        updateSliderInput(inputId = "alpha", value = 0.5)
        updateSliderInput(inputId = "beta", value = 0.5)
    })
    observeEvent(input$buttonBetaLeft, {
        updateSliderInput(inputId = "alpha", value = 0.5)
        updateSliderInput(inputId = "beta", value = 1)
    })
    observeEvent(input$buttonBetaRight, {
        updateSliderInput(inputId = "alpha", value = 1)
        updateSliderInput(inputId = "beta", value = 0.5)
    })
    observeEvent(input$buttonBetaLeftTriangle, {
        updateSliderInput(inputId = "alpha", value = 1)
        updateSliderInput(inputId = "beta", value = 2)
    })
    
    observeEvent(input$buttonBetaRightTriangle, {
        updateSliderInput(inputId = "alpha", value = 2)
        updateSliderInput(inputId = "beta", value = 1)
    })
    
    growth.rates <- reactive(as.numeric(text2chars(input$growthRates, len = n.species(), expr = paste0('round(rbeta(',n.species(), ',' ,alpha(), ',' , beta(),'), digits = 3)'))))
    output$growthRatesOutput <- renderPrint(growth.rates())
    
    output$growthRatesDist <- renderPlot({
        ggplot(data.frame(x = c(-4, 4)), aes(x = x)) + 
            stat_function(fun = dbeta, args = list(shape1=alpha(), shape2=beta())) + 
            xlim(0,1) + theme_linedraw()
    }, res = 96) 
    monod.constant <- reactive(text2char(input$monodConstant)) # TODO: edit table
    
    ## miscellaneous
    error.variance <- reactive(input$errorVariance)
    norm <- reactive(input$norm)
    t.end <- reactive(input$tEnd)
    
    runCRM <- reactive(
        simulateConsumerResource(
            n.species = n.species(),
            n.resources = n.resources(),
            names.species = names.species(), 
            names.resources = names.resources(), 
            #E = E(), 
            E = E$df,
            x0 = x0(), 
            resources = resources(), 
            growth.rates = growth.rates(), 
            monod.constant = monod.constant(), 
            error.variance = error.variance(), 
            norm = norm(), 
            t.end = t.end()
        )
    )
    
    output$CRMSpecies <- renderPlot(makePlot(runCRM()$matrix, "abundance of species by time"), res = 96)
    output$CRMResources <- renderPlot(makePlotRes(runCRM()$resources, "quantity of compounds by time"),  res = 96)
    
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