ui <- navbarPage(
    title ="microSimShiny",
    # tab1 Consumer-Resource Model ####
    tabPanel(
        title = "Consumer-Resource Model",
        fluidPage(
            titlePanel("Consumer-Resource Model"),
            bsCollapse(id = "contents", 
                       open = "Model", 
                       ## Model ####
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
                                                          ### Basic ####
                                                          tabPanel("Basic",
                                                                   sliderInput("nSpecies", "number of species", value = 2, min = 1, max = 100),
                                                                   bsTooltip("nSpecies", "Number of species in the simulation", "right", options = list(container = "body")),
                                                                   sliderInput("nResources", "number of compounds", value = 4, min = 1, max = 100),
                                                                   bsTooltip("nResources", "Number of compounds in the simulation", "right", options = list(container = "body")),
                                                                   hr(),
                                                                   checkboxInput("CustomCRM", strong("custom names of species/compounds or simulating time points"), value = FALSE),
                                                                   conditionalPanel(condition = "input.CustomCRM",
                                                                                    helpText("Custom names separate by ',' or ';' (and spaces) will replace default names."),
                                                                                    textInput("namesSpecies", "names of species"),
                                                                                    textInput("namesResources", "names of compounds"),
                                                                   ),
                                                                   conditionalPanel(condition = "input.CustomCRM",
                                                                                    hr(),
                                                                                    numericInput("tStart", "start time of the simulation", value = 0, min = 0, max = 10000, step = 100),
                                                                                    bsTooltip("tStart", "The start time of the simulation.", "right", options = list(container = "body")),
                                                                   ),
                                                                   
                                                                   numericInput("tEnd", "final time of the simulation", value = 1000, min = 100, max = 10000, step = 100),
                                                                   bsTooltip("tEnd", "The end time of the simulation.", "right", options = list(container = "body")),
                                                                   
                                                                   conditionalPanel(condition = "input.CustomCRM",
                                                                                    numericInput("tStep", "time step of the simulation", value = 0.1, min = 0.01, max = 10, step = 0.01),
                                                                                    bsTooltip("tStep", "The time step of the simulation.", "right", options = list(container = "body")),
                                                                                    numericInput("tStore", "stored time points of the simulation", value = 1000, min = 100, max = 10000, step = 100),
                                                                                    bsTooltip("tStore", "The stored time points of the simulation.", "right", options = list(container = "body")),
                                                                   ),
                                                                   
                                                                   
                                                                   
                                                          ),
                                                          ### Compounds ####
                                                          tabPanel("Compounds",
                                                                   sliderInput("resourcesConcentration", "mean initial concentration of compounds", min = 0, max = 1000, value = 100, step = 1),
                                                                   bsTooltip("resourcesConcentration", "mean initial concentration of compounds", "right", options = list(container = "body")),
                                                                   sliderInput("resourcesEvenness", "evenness of compounds", min = 0.1, max = 100, value = 10, step = 0.1),
                                                                   bsTooltip("resourcesEvenness", "higher evenness leads to similar concentrations of initial compounds", "right", options = list(container = "body")),
                                                                   textInput("resourcesCustom", "initial concentration of compounds"),
                                                                   bsTooltip("resourcesCustom", "If the given initial concentrations of compounds are not enough, random values will be added.", "right", options = list(container = "body")),
                                                                   verbatimTextOutput("resourcesOutput"),
                                                                   
                                                                   sliderInput("dilutionRate", "dilution rate", min = 0, max = 1, value = 0, step = 0.001),
                                                                   bsTooltip("dilutionRate", "rate of nutrient exchange as in a Chemostat bioreactor. Flow rate / culture volume", "right", options = list(container = "body")),
                                                                   
                                                                   conditionalPanel(condition = "input.dilutionRate > 0",
                                                                            textInput("resourcesDilution", "resources dilution"),
                                                                            verbatimTextOutput("resourcesDilutionOutput"),
                                                                    ),
                                                                   
                                                                   plotOutput("resourcesPlot"),
                                                                   
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
                                                          ### Growth rates ####
                                                          tabPanel("Growth rates",
                                                                   textInput("x0", "initial abundances of species"),
                                                                   bsTooltip("x0", "If the given initial abundances of species is not enough, random initial abundances will be added.", "right", options = list(container = "body")),
                                                                   verbatimTextOutput("x0Output"),
                                                                   tags$label("Distribution of Growth Rates"),
                                                                   hr(),
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
                                                                   tags$label("Monod Constant"),
                                                                   dataTableOutput("tableMonodConstant", width = "100%"),
                                                          ),
                                                          ### Perturbations ####
                                                          tabPanel("Perturbations",
                                                                   sliderInput("errorVariance", "variance of measurement error", value = 0, min = 0, max = 10, step = 0.1),
                                                                   bsTooltip("errorVariance", "The variance of measurement error. By default it equals to 0, indicating that the result won't contain any measurement error.", "right", options = list(container = "body")),
                                                                   checkboxInput("stochastic", "use stochasitic", value = TRUE),
                                                                   conditionalPanel(condition = "input.stochastic",
                                                                            sliderInput("sigmaDrift", "strength of drift", value = 0, min = 0, max = 1, step = 0.001),
                                                                            sliderInput("sigmaEpoch", "strength of microbial epoch perturbation", value = 0.001, min = 0, max = 1, step = 0.001),
                                                                            sliderInput("sigmaExternal", "strength of external perturbations", value = 0.3, min = 0, max = 1, step = 0.001),
                                                                            sliderInput("sigmaMigration", "intensity of migration", value = 0.01, min = 0, max = 1, step = 0.001),
                                                                            sliderInput("epochP", "epoch.p", value = 0.001, min = 0, max = 1, step = 0.001),
                                                                            textInput("tExternalEvents", "starting time of external events"),
                                                                            verbatimTextOutput("tExternalEventsOutput"),
                                                                            textInput("tExternalDurations", "durations of external events"),
                                                                            verbatimTextOutput("tExternalDurationsOutput"),
                                                                            ),
                                                                   sliderInput("migrationP", "probability/frequency of migration from metacommunity", value = 0.01, min = 0, max = 1),
                                                                   textInput("metacommunityProbability", "metacommunity"),
                                                                   bsTooltip("metacommunityProbability", "Normalized probability distribution of the likelihood that species from the metacommunity can enter the community during the simulation", "right", options = list(container = "body")),
                                                                   verbatimTextOutput("metacommunityProbability"),
                                                                   checkboxInput("norm", strong("returns normalized abundances"), value = FALSE),
                                                                   # actionButton("buttonSimulateCRM", "Run the model", class = "btn btn-primary")
                                                          ),
                                                      ),
                                                  )
                                           ),
                                           ### DisplayPanel(Right) ####
                                           column(width = 7, 
                                                  fluidRow(
                                                      actionButton("CRMEX1", "Ex.1", class = "btn btn-primary", width = "12%"),
                                                      actionButton("CRMEX2", "Ex.2", class = "btn btn-primary", width = "12%"),
                                                      actionButton("CRMEX3", "Ex.3", class = "btn btn-primary", width = "12%"),
                                                      bsButton("CRMEX4pre", "Pre.4", style = "info", width = "7%"),
                                                      bsButton("CRMEX4", "Ex.4", style = "primary", disabled = TRUE, width = "7%"),
                                                      bsTooltip("CRMEX4", "Please run the Pre.4 First", "top", options = list(container = "body")),
                                                      bsButton("CRMEX5pre", "Pre.5", style = "info", width = "7%"),
                                                      bsButton("CRMEX5", "Ex.5", style = "primary", disabled = TRUE, width = "7%"),
                                                      bsTooltip("CRMEX5", "Please run the Pre.5 First", "top", options = list(container = "body")),
                                                      bsButton("CRMEX6pre", "Pre.6", style = "info", width = "7%"),
                                                      bsButton("CRMEX6", "Ex.6", style = "primary", disabled = TRUE, width = "7%"),
                                                      bsTooltip("CRMEX6", "Please run the Pre.6 First", "top", options = list(container = "body")),
                                                      bsButton("CRMEX7pre", "Pre.7", style = "info", width = "7%"),
                                                      bsButton("CRMEX7", "Ex.7", style = "primary", disabled = TRUE, width = "7%"),
                                                      bsTooltip("CRMEX7", "Please run the Pre.7 First", "top", options = list(container = "body")),
                                                  ),
                                                  fluidRow(
                                                      plotOutput("CRMSpecies"),
                                                      plotOutput("CRMResources"),
                                                      plotOutput("CRMPlotE"),
                                                      bsTooltip("CRMPlotE", "Positive values indicate consumption, and negative values indicate production", "left", options = list(container = "body")),
                                                  ),
                                                  
                                                  
                                           ),
                                       )
                       ),
                       ## Description ####
                       bsCollapsePanel(strong("Description"), value = "Description",
                                       fluidRow(
                                           column(width = 12,
                                                  withMathJax(includeMarkdown("crm.md")),
                                           ),
                                       )
                       ),
                       ## Inputs ####
                       bsCollapsePanel(strong("Inputs"), value = "Inputs",
                                       fluidRow(
                                           column(width = 12,
                                                  withMathJax(includeMarkdown("crm_parms.Rmd")),
                                           ),
                                       )
                       ),
                       
                       bsCollapsePanel(strong("References"), value = "References",
                                       "Panel of refs.")
            ),
            
        )
    ),
    
    # tab2 Generalized Lotka-Volterra (GLV) Model ####
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
    
    # tab3 Hubbell Model ####
    tabPanel(
        title = "Hubbell Model (with death rates)",
        fluidPage(
            
        )),
    
    # tab4 Logistic Model (with stochasticity) ####
    tabPanel(
        title = "Logistic Model (with stochasticity)",
        fluidPage(
            
        )),
)
