ui <- navbarPage(
    title ="microSimShiny",
    # tab1 Consumer-Resource Model ####
    tabPanel(
        title = "Consumer-Resource Model",
        titlePanel("Consumer-Resource Model"),
        # enable bsplus tooltips and popovers
        use_bs_tooltip(),
        use_bs_popover(),
        # enable Shinyjs
        useShinyjs(),
        ## CRM Model ####
        bs_accordion(id = "CRMcontents") %>%
            bs_append(
                title = "Model",
                content = 
                    fluidRow(
                        column(
                            width = 5,
                            wellPanel(
                                tabsetPanel(
                                    header = tags$style(
                                        HTML(
                                            "
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
                                            "
                                        )
                                    ),
                                    ### Basic ####
                                    tabPanel(
                                        "Basic",
                                        sliderInput(
                                            "nSpeciesCRM",
                                            "number of species", 
                                            value = 2, 
                                            min = 1, 
                                            max = 100) %>%
                                        shinyInput_label_embed(
                                            shiny_iconlink() %>% 
                                            bs_embed_tooltip(
                                                title =  "Number of species in the simulation", 
                                                placement = "right", 
                                                container = "body"
                                            )
                                        ),
                                        sliderInput(
                                            "nResourcesCRM", 
                                            "number of compounds",
                                            value = 4, 
                                            min = 1, 
                                            max = 100) %>%
                                        shinyInput_label_embed(
                                            shiny_iconlink() %>% 
                                            bs_embed_tooltip(
                                                title =  "Number of compounds in the simulation", 
                                                placement = "right", 
                                                container = "body"
                                            )
                                        ),
                                        tags$hr(),
                                        switchInput(
                                            "CustomCRM",
                                            strong("custom names/times in simulation"),
                                            value = FALSE,
                                            labelWidth = "100%"
                                        ),
                                        conditionalPanel(
                                            condition = "input.CustomCRM",
                                            helpText("Custom names separate by ',' or ';' (and spaces) will replace default names."),
                                            textInput("namesSpeciesCRM", "names of species"),
                                            textInput("namesResourcesCRM", "names of compounds"),
                                        ),
                                        conditionalPanel(
                                            condition = "input.CustomCRM",
                                            tags$hr(),
                                            numericInput("tStartCRM", "start time of the simulation", value = 0, min = 0, max = 10000, step = 100),
                                        ),
                                        numericInput("tEndCRM", "final time of the simulation", value = 1000, min = 100, max = 10000, step = 100),
                                        conditionalPanel(
                                            condition = "input.CustomCRM",
                                            numericInput("tStepCRM", "time step of the simulation", value = 0.1, min = 0.01, max = 10, step = 0.01),
                                            numericInput("tStoreCRM", "stored time points of the simulation", value = 1000, min = 100, max = 10000, step = 100),
                                        ),
                                    ), 
                                    ### Compounds ####
                                    tabPanel(
                                        "Compounds",
                                        sliderInput(
                                            "resourcesConcentrationCRM", 
                                            "mean initial concentration of compounds", 
                                            min = 0, 
                                            max = 1000, 
                                            value = 100, 
                                            step = 1) %>%
                                        shinyInput_label_embed(
                                            shiny_iconlink() %>% 
                                            bs_embed_tooltip(
                                                title =  "initial average concentration of each compound", 
                                                placement = "right", 
                                                container = "body"
                                            )
                                        ),
                                        sliderInput(
                                            "resourcesEvennessCRM",
                                            "evenness of compounds",
                                            min = 0.1,
                                            max = 100,
                                            value = 10,
                                            step = 0.1) %>%
                                        shinyInput_label_embed(
                                            shiny_iconlink() %>% 
                                            bs_embed_tooltip(
                                                title =  "higher evenness leads to similar concentrations of initial compounds", 
                                                placement = "right", 
                                                container = "body"
                                            )
                                        ),
                                        textInput(
                                            "resourcesCustomCRM", 
                                            "initial concentrations of compounds") %>%
                                        shinyInput_label_embed(
                                            shiny_iconlink() %>% 
                                            bs_embed_tooltip(
                                                title =  "If the given initial concentrations of compounds are not enough, random values will be added.", 
                                                placement = "right", 
                                                container = "body"
                                            )
                                        ),
                                        verbatimTextOutput("resourcesOutputCRM"),
                                        sliderInput(
                                            "dilutionRateCRM", 
                                            "dilution rate",
                                            min = 0,
                                            max = 1, 
                                            value = 0,
                                            step = 0.001) %>%
                                        shinyInput_label_embed(
                                            shiny_iconlink() %>% 
                                            bs_embed_tooltip(
                                                title =  "rate of nutrient exchange as in a Chemostat bioreactor. Flow rate / culture volume", 
                                                placement = "right", 
                                                container = "body"
                                            )
                                        ),
                                        
                                        conditionalPanel(
                                            condition = "input.dilutionRateCRM > 0",
                                            textInput("resourcesDilutionCRM", "resources concentration in dilution")  %>%
                                            shinyInput_label_embed(
                                                shiny_iconlink() %>% 
                                                bs_embed_tooltip(
                                                    title =  "concentrations of resources in continuous flow, by default equal to initial concentrations of compounds", 
                                                    placement = "right", 
                                                    container = "body"
                                                )
                                            ),
                                            
                                            verbatimTextOutput("resourcesDilutionCRMOutput"),
                                        ),
                                        
                                        plotOutput("resourcesCRMPlot"),
                                        
                                        sliderInput(
                                            "meanConsumptionCRM", 
                                            "consumption weight", 
                                            value = 0.4,
                                            min = 0, 
                                            max = 1) %>%
                                            shinyInput_label_embed(
                                                shiny_iconlink() %>% 
                                                bs_embed_tooltip(
                                                    title =  "Mean proportion of compounds consumed by each species.", 
                                                    placement = "right", 
                                                    container = "body"
                                                )
                                            ),
                                        sliderInput(
                                            "meanProductionCRM", 
                                            "production weight",
                                            value = 0.2,
                                            min = 0, 
                                            max = 0.5) %>%
                                        shinyInput_label_embed(
                                            shiny_iconlink() %>% 
                                            bs_embed_tooltip(
                                                title =  "Mean proportion of compounds produced by each species.", 
                                                placement = "right", 
                                                container = "body"
                                            )
                                        ),
                                        sliderInput("maintenanceCRM", "maintenance weight", value = 0.5, min = 0, max = 1) %>%
                                        shinyInput_label_embed(
                                            shiny_iconlink() %>% 
                                            bs_embed_tooltip(
                                                title =  "How much compounds were used to maintain the microbial community (not involved in flux of compounds).", 
                                                placement = "right", 
                                                container = "body"
                                            )
                                        ),
                                        tags$div(
                                            tags$label("Compounds Stochiometry"),
                                            tags$div(
                                                class = "pull-right",
                                                shiny_iconlink() %>% 
                                                bs_embed_tooltip(
                                                    title =  "Stochiometric values of consumption and production of compounds by each cell. Positive efficiencies indicate the consumption of resources, whilst negatives indicate that the species would produce the resource.", 
                                                    placement = "right", 
                                                    container = "body"
                                                ),
                                            )
                                        ),
                                        dataTableOutput("tableECRM", width = "100%"),
                                    ),
                                    ### Growth rates ####
                                    tabPanel(
                                        "Growth rates",
                                        textInput("x0CRM", "initial abundances of species")  %>%
                                        shinyInput_label_embed(
                                            shiny_iconlink() %>% 
                                            bs_embed_tooltip(
                                                title =  "If the given initial abundances of species is not enough, random initial abundances will be added.", 
                                                placement = "right", 
                                                container = "body"
                                            )
                                        ),
                                        verbatimTextOutput("x0CRMOutput"),
                                        tags$hr(),
                                        tags$div(
                                            tags$label("Distribution of Growth Rates"),
                                        ),
                                        bs_button(
                                            label = "-", 
                                            button_type = "primary", 
                                            id = "buttonBetaEven", 
                                            style = "width:13%",
                                            class = "btn-default action-button shiny-bound-input"
                                        ) %>% bs_embed_tooltip(title = "even distribution"),
                                        bs_button(
                                            label = "^", 
                                            button_type = "primary", 
                                            id = "buttonBetaRidge", 
                                            style = "width:13%",
                                            class = "btn-default action-button shiny-bound-input"
                                        ) %>% bs_embed_tooltip(title = "normal-alike distribution"),
                                        bs_button(
                                            label = "˅", 
                                            button_type = "primary", 
                                            id = "buttonBetaValley", 
                                            style = "width:13%",
                                            class = "btn-default action-button shiny-bound-input"
                                        ) %>% bs_embed_tooltip(title = "U shape distribution"),
                                        bs_button(
                                            label = "◟",
                                            button_type = "primary", 
                                            id = "buttonBetaLeft", 
                                            style = "width:13%",
                                            class = "btn-default action-button shiny-bound-input"
                                        ) %>% bs_embed_tooltip(title = "left skewed distribution"),
                                        bs_button(
                                            label = "◞", 
                                            button_type = "primary", 
                                            id = "buttonBetaRight", 
                                            style = "width:13%",
                                            class = "btn-default action-button shiny-bound-input"
                                        ) %>% bs_embed_tooltip(title = "right skewed distribution"),
                                        bs_button(
                                            label = "\\", 
                                            button_type = "primary", 
                                            id = "buttonBetaLeftTriangle", 
                                            style = "width:13%",
                                            class = "btn-default action-button shiny-bound-input"
                                        ) %>% bs_embed_tooltip(title = "left-triangle distribution"),
                                        bs_button(
                                            label = "/", 
                                            button_type = "primary", 
                                            id = "buttonBetaRightTriangle", 
                                            style = "width:13%",
                                            class = "btn-default action-button shiny-bound-input"
                                        ) %>% bs_embed_tooltip(title = "right-triangle distribution"),
                                        
                                        sliderInput(
                                            "alphaCRM",
                                            "alpha",
                                            value = 1,
                                            min = 0,
                                            max = 10,
                                            step = 0.1
                                        )  %>%
                                        shinyInput_label_embed(
                                            shiny_iconlink() %>% 
                                            bs_embed_tooltip(
                                                title =  "first parameter of beta distribution", 
                                                placement = "right", 
                                                container = "body"
                                            )
                                        ),
                                        sliderInput(
                                            "betaCRM",
                                            "beta",
                                            value = 1,
                                            min = 0,
                                            max = 10,
                                            step = 0.1
                                        )  %>%
                                        shinyInput_label_embed(
                                            shiny_iconlink() %>% 
                                            bs_embed_tooltip(
                                                title =  "second parameter of beta distribution", 
                                                placement = "right", 
                                                container = "body"
                                            )
                                        ),
                                        textInput(
                                            "growthRatesCRM", 
                                            "maximum growth rates of species"
                                        ) %>%
                                        shinyInput_label_embed(
                                            shiny_iconlink() %>% 
                                            bs_embed_tooltip(
                                                title =  "If the given growth rates are not enough, random values will be added.", 
                                                placement = "right", 
                                                container = "body"
                                            )
                                        ),
                                        verbatimTextOutput("growthRatesCRMOutput"),
                                        
                                        plotOutput("growthRatesCRMDist"),
                                        tags$label("Monod Constant"),
                                        dataTableOutput("tableMonodCRM", width = "100%"),
                                    ),
                                    ### Perturbations ####
                                    tabPanel(
                                        "Perturbations",
                                        sliderInput(
                                            "errorVarianceCRM",
                                            "variance of measurement error",
                                            value = 0,
                                            min = 0,
                                            max = 1,
                                            step = 0.01)  %>%
                                            shinyInput_label_embed(
                                                shiny_iconlink() %>% 
                                                    bs_embed_tooltip(
                                                        title =  "The variance of measurement error. By default it equals to 0, indicating that the result won't contain any measurement error.", 
                                                        placement = "right", 
                                                        container = "body"
                                                    )
                                            ),
                                        switchInput(
                                            "stochasticCRM", 
                                            strong("use stochasitic"), 
                                            value = TRUE,
                                            labelWidth = "100%"
                                        ),
                                        conditionalPanel(
                                            condition = "input.stochasticCRM",
                                            sliderInput("sigmaDriftCRM", "strength of drift", value = 0, min = 0, max = 1, step = 0.001)  %>%
                                            shinyInput_label_embed(
                                                shiny_iconlink() %>% 
                                                bs_embed_tooltip(
                                                    title =  "drift happens on each step of simulation", 
                                                    placement = "right", 
                                                    container = "body"
                                                )
                                            ),
                                            tags$hr(),
                                            sliderInput("epochPCRM", "probability of random periodic (epoch) changes", value = 0.001, min = 0, max = 1, step = 0.001) %>%
                                                shinyInput_label_embed(
                                                    shiny_iconlink() %>% 
                                                        bs_embed_tooltip(
                                                            title =  "microbial epoch perturbations happens by chance", 
                                                            placement = "right", 
                                                            container = "body"
                                                        )
                                                ),
                                            conditionalPanel(
                                                condition = "input.epochPCRM >0",
                                                sliderInput("sigmaEpochCRM", "strength of microbial epoch perturbation", value = 0.001, min = 0, max = 1, step = 0.001), 
                                                
                                            ),
                                            tags$hr(),
                                            sliderInput("sigmaExternalCRM", "strength of external perturbations", value = 0.3, min = 0, max = 1, step = 0.001),
                                            conditionalPanel(
                                                condition = "input.sigmaExternalCRM >0",
                                                textInput("tExternalEventsCRM", "starting time of external events"),
                                                verbatimTextOutput("tExternalEventsCRMOutput"),
                                                textInput("tExternalDurationsCRM", "durations of external events"),
                                                verbatimTextOutput("tExternalDurationsCRMOutput"),
                                            ),
                                        ),
                                        tags$hr(),
                                        sliderInput("migrationPCRM", "probability/frequency of migration from metacommunity", value = 0.01, min = 0, max = 1),
                                        conditionalPanel(
                                            condition = "input.migrationPCRM >0",
                                            sliderInput("sigmaMigrationCRM", "intensity of migration", value = 0.01, min = 0, max = 1, step = 0.001),
                                        ),
                                        textInput(
                                            "metacommunityProbabilityCRM",
                                            "metacommunity") %>%
                                            shinyInput_label_embed(
                                                shiny_iconlink() %>% 
                                                    bs_embed_tooltip(
                                                        title =  "Normalized probability distribution of the likelihood that species from the metacommunity can enter the community during the simulation.", 
                                                        placement = "right", 
                                                        container = "body"
                                                    )
                                            ),
                                        verbatimTextOutput("metacommunityProbabilityCRM"),
                                        tags$hr(),
                                        switchInput(
                                            "normCRM", 
                                            strong("returns normalized abundances"), 
                                            value = FALSE,
                                            labelWidth = "100%"
                                        ),
                                        # actionButton("buttonSimulateCRM", "Run the model", class = "btn btn-primary")
                                    ),
                                ),
                            )
                        ),
                        ### Display Panel ####
                        column(
                            width = 7, 
                            #### example buttons ####
                            fluidRow(
                                style = "padding-left: 15px; padding-right: 15px;",
                                tags$div(
                                    class = "panel panel-default",
                                    tags$div(
                                        class = "panel-heading",
                                        tags$h3(
                                            class = "panel-title",
                                            "Examples",
                                        ),
                                    ),
                                    tags$div(
                                        class = "panel-body",
                                        bs_button(
                                            label = "Ex.1", 
                                            button_type = "primary", 
                                            id = "CRMEX1", 
                                            style = "width:12%",
                                            class = "btn-default action-button shiny-bound-input"
                                        ),
                                        bs_button(
                                            label = "Ex.2", 
                                            button_type = "primary", 
                                            id = "CRMEX2", 
                                            style = "width:12%",
                                            class = "btn-default action-button shiny-bound-input"
                                        ),
                                        bs_button(
                                            label = "Ex.3", 
                                            button_type = "primary", 
                                            id = "CRMEX3", 
                                            style = "width:12%",
                                            class = "btn-default action-button shiny-bound-input"
                                        ),
                                        bs_button(
                                            label = "Pre.4", 
                                            button_type = "info", 
                                            id = "CRMEX4pre", 
                                            style = "width:7%",
                                            class = "btn-default action-button shiny-bound-input"
                                        ),
                                        bs_button(
                                            label = "Ex.4", 
                                            button_type = "primary", 
                                            id = "CRMEX4", 
                                            style = "width:7%",
                                            disabled = "disabled",
                                            class = "btn-default action-button shiny-bound-input"
                                        ) %>% bs_embed_tooltip(title = "Please run the Pre.4 First"),
                                        bs_button(
                                            label = "Pre.5", 
                                            button_type = "info", 
                                            id = "CRMEX5pre", 
                                            style = "width:7%",
                                            class = "btn-default action-button shiny-bound-input"
                                        ),
                                        bs_button(
                                            label = "Ex.5", 
                                            button_type = "primary", 
                                            id = "CRMEX5", 
                                            style = "width:7%",
                                            disabled = "disabled",
                                            class = "btn-default action-button shiny-bound-input"
                                        ) %>% bs_embed_tooltip(title = "Please run the Pre.5 First"),
                                        bs_button(
                                            label = "Pre.6", 
                                            button_type = "info", 
                                            id = "CRMEX6pre", 
                                            style = "width:7%",
                                            class = "btn-default action-button shiny-bound-input"
                                        ),
                                        bs_button(
                                            label = "Ex.6", 
                                            button_type = "primary", 
                                            id = "CRMEX6", 
                                            style = "width:7%",
                                            disabled = "disabled",
                                            class = "btn-default action-button shiny-bound-input"
                                        ) %>% bs_embed_tooltip(title = "Please run the Pre.6 First"),
                                        bs_button(
                                            label = "Pre.7", 
                                            button_type = "info", 
                                            id = "CRMEX7pre", 
                                            style = "width:7%",
                                            class = "btn-default action-button shiny-bound-input"
                                        ),
                                        bs_button(
                                            label = "Ex.7", 
                                            button_type = "primary", 
                                            id = "CRMEX7", 
                                            style = "width:7%",
                                            disabled = "disabled",
                                            class = "btn-default action-button shiny-bound-input"
                                        ) %>% bs_embed_tooltip(title = "Please run the Pre.7 First"),
                                    ),
                                ),
                            ),
                            #### result plots ####
                            tags$br(),
                            fluidRow(
                                style = "padding-left: 15px; padding-right: 15px;",
                                tags$div(
                                    class = "panel panel-default",
                                    tags$div(
                                        class = "panel-heading",
                                        tags$h3(
                                            class = "panel-title",
                                            "Species Change",
                                        ),
                                    ),
                                    tags$div(
                                        class = "panel-body",
                                        plotOutput("CRMSpecies"),
                                    ),
                                ),
                                tags$br(),
                                tags$div(
                                    class = "panel panel-default",
                                    tags$div(
                                        class = "panel-heading",
                                        tags$h3(
                                            class = "panel-title",
                                            "Compounds Change",
                                        ),
                                    ),
                                    tags$div(
                                        class = "panel-body",
                                        plotOutput("CRMResources"),
                                        
                                    ),
                                ),
                                tags$br(),
                                tags$div(
                                    class = "panel panel-default",
                                    tags$div(
                                        class = "panel-heading",
                                        tags$h3(
                                            class = "panel-title",
                                            "Matrix Efficiency",
                                            tags$div(
                                                class = "pull-right",
                                                shiny_iconlink() %>% 
                                                    bs_embed_tooltip(
                                                        title =  "Positive values indicate consumption, and negative values indicate production", 
                                                        placement = "left", 
                                                        container = "body"
                                                    ),
                                            ),
                                        ),
                                    ),
                                    tags$div(
                                        class = "panel-body",
                                        plotOutput("CRMPlotE"),
                                    ),
                                ),
                            ),
                        ),
                    )
            ) %>%
            bs_append(
                ## Description ####
                title = "Description",
                content = 
                    fluidRow(
                        column(
                            width = 12,
                            withMathJax(includeMarkdown("crm.Rmd")),
                        ),
                    )
            ) %>% 
            bs_append(
                ## Inputs ####
                title = "Inputs",
                content = 
                    fluidRow(
                        column(width = 12,
                               withMathJax(includeMarkdown("crm_parms.Rmd")),
                        ),
                    )
            ) %>%
            bs_append(
                ## References ####
                title = "References",
                content = "Panel of refs."
            )
    ),
    
    # tab2 Generalized Lotka-Volterra (GLV) Model ####
    tabPanel(
        title = "Generalized Lotka-Volterra (GLV)",
        titlePanel("Generalized Lotka-Volterra (GLV)"),
        ## GLV Model ####
        bs_accordion(id = "GLVcontents") %>%
            bs_append(
                title = "Model",
                content = 
                    fluidRow(
                        column(
                            width = 5,
                            wellPanel(
                                tabsetPanel(
                                    ### Interspecies interactions ####
                                    tabPanel(
                                        "Interspecies interactions",
                                        sliderInput(
                                            "nSpeciesGLV",
                                            "number of species",
                                            value = 2,
                                            min = 2,
                                            max = 20),
                                        switchInput(
                                            "advancedRandomA",
                                            strong("show advanced parameters"),
                                            labelWidth = "100%"
                                        ),
                                        conditionalPanel(
                                            condition = "input.advancedRandomA",
                                            helpText("Custom names separate by ',' or ';' (and spaces) will replace default names."),
                                            textInput("namesSpeciesGLV", "names of species"),
                                            sliderInput("diagonalGLV", "diagonal values of matrix A", value = -0.5, min = -2, max = 0, step = 0.1),
                                            sliderInput("connectanceGLV", "connectance of matrix A", value = 0.2, min = 0, max = 1),
                                            sliderInput("scaleGLV", "scale of off-diagonal values", value = 0.1, min = 0, max = 1),
                                            tags$hr(),
                                            numericInput("mutualismGLV", "relative proportion of mutualism in matrix A", value = 1, min = 0, max = 10,step = 0.05),
                                            numericInput("commensalismGLV", "relative proportion of commensalism in matrix A", value = 1, min = 0, max = 10,step = 0.05),
                                            numericInput("parasitismGLV", "relative proportion of parasitism in matrix A", value = 1, min = 0, max = 10,step = 0.05),
                                            numericInput("amensalismGLV", "relative proportion of amensalism in matrix A", value = 1, min = 0, max = 10,step = 0.05),
                                            numericInput("competitionGLV", "relative proportion of competition in matrix A", value = 1, min = 0, max = 10,step = 0.05),
                                            
                                            textInput("interactionsCustomGLV", "user-defined interactions between species") %>%
                                            shinyInput_label_embed(
                                                shiny_iconlink() %>% 
                                                bs_embed_tooltip(
                                                    title =  "if the given interactions between species are not enough, random values will be added", 
                                                    placement = "right", 
                                                    container = "body"
                                                )
                                            ),
                                            verbatimTextOutput("interactionsOutputGLV"),
                                            bs_button(
                                                label = "˅", 
                                                button_type = "primary", 
                                                id = "buttonInteractionsU", 
                                                style = "width:30%",
                                                class = "btn-default action-button shiny-bound-input"
                                            ) %>% bs_embed_tooltip(title = "rbeta(n, shape1 = 0.5, shape2 = 0.5)"),
                                            bs_button(
                                                label = "^", 
                                                button_type = "primary", 
                                                id = "buttonInteractionsN", 
                                                style = "width:30%",
                                                class = "btn-default action-button shiny-bound-input"
                                            ) %>% bs_embed_tooltip(title = "rbeta(n, shape1 = 2, shape2 = 2)"),
                                            bs_button(
                                                label = "-",
                                                button_type = "primary",
                                                id = "buttonInteractionsEven",
                                                style = "width:30%",
                                                class = "btn-default action-button shiny-bound-input"
                                            ) %>% bs_embed_tooltip(title = "runif(n, min = 0, max = 1)"),
                                            
                                            tags$hr(),
                                            switchInput(
                                                "symmetricGLV",
                                                strong("matrix A"),
                                                labelWidth = "100%",
                                                onLabel = "symmetric",
                                                offLabel = "non-symmetric"
                                            ),
                                            # textInput("listAGLV", "a list of previous generated matrix"),
                                        ),
                                        # actionButton("buttonRandomA", "generate random matrix A of interspecies interactions", class = "btn btn-primary"),
                                        conditionalPanel(
                                            condition = "input.advancedRandomA",
                                            tags$hr(),
                                            numericInput("tStartGLV", "start time of the simulation", value = 0, min = 0, max = 10000, step = 100),
                                        ),
                                        numericInput("tEndGLV", "final time of the simulation", value = 1000, min = 100, max = 10000, step = 100),
                                        conditionalPanel(
                                            condition = "input.advancedRandomA",
                                            numericInput("tStepGLV", "time step of the simulation", value = 0.1, min = 0.01, max = 10, step = 0.01),
                                            numericInput("tStoreGLV", "stored time points of the simulation", value = 1000, min = 100, max = 10000, step = 100),
                                        ),
                                    ),
                                    ### Growth rates ####
                                    tabPanel(
                                        "Growth rates",
                                        textInput("x0GLV", "initial abundances of species"),
                                        verbatimTextOutput("x0GLVOutput"),
                                        textInput("growthRatesGLV", "maximum growth rates of species"),
                                        verbatimTextOutput("growthRatesGLVOutput"),
                                    ),
                                    ### Perturbations ####
                                    tabPanel(
                                        "Perturbations",
                                        sliderInput(
                                            "errorVarianceGLV", 
                                            "variance of measurement error", 
                                            value = 0, 
                                            min = 0, 
                                            max = 10, 
                                            step = 0.1) %>%
                                        shinyInput_label_embed(
                                            shiny_iconlink() %>% 
                                                bs_embed_tooltip(
                                                    title =  "The variance of measurement error. By default it equals to 0, indicating that the result won't contain any measurement error.", 
                                                    placement = "right", 
                                                    container = "body"
                                                )
                                        ),
                                        
                                        switchInput(
                                            "stochasticGLV", 
                                            strong("use stochasticity"), 
                                            value = TRUE,
                                            labelWidth = "100%"),
                                        
                                        conditionalPanel(
                                            condition = "input.stochasticGLV",
                                            sliderInput(
                                                "sigmaDriftGLV", 
                                                "strength of drift", 
                                                value = 0.001,
                                                min = 0, 
                                                max = 1, 
                                                step = 0.001
                                            ) %>%
                                            shinyInput_label_embed(
                                                shiny_iconlink() %>% 
                                                    bs_embed_tooltip(
                                                        title =  "drift happens on each step of simulation", 
                                                        placement = "right", 
                                                        container = "body"
                                                    )
                                            ),
                                            tags$hr(),
                                            
                                            sliderInput(
                                                "epochPGLV", 
                                                "probability of random periodic (epoch) changes", 
                                                value = 0.001,
                                                min = 0, 
                                                max = 1, 
                                                step = 0.001
                                            ) %>%
                                            shinyInput_label_embed(
                                                shiny_iconlink() %>% 
                                                    bs_embed_tooltip(
                                                        title =  "microbial epoch perturbations happens by chance", 
                                                        placement = "right", 
                                                        container = "body"
                                                    )
                                            ),
                                            conditionalPanel(
                                                condition = "input.epochPGLV > 0",
                                                sliderInput("sigmaEpochGLV", "strength of microbial epoch perturbation", value = 0.001, min = 0, max = 1, step = 0.001),
                                            ),
                                            tags$hr(),
                                            sliderInput("sigmaExternalGLV", "strength of external perturbations", value = 0.3, min = 0, max = 1, step = 0.001),
                                            conditionalPanel(
                                                condition = "input.sigmaExternalGLV > 0",
                                                textInput("tExternalEventsGLV", "timepoints of external events"),
                                                verbatimTextOutput("tExternalEventsGLVOutput"),
                                                textInput("tExternalDurationsGLV", "time durations of external events"),  
                                                verbatimTextOutput("tExternalDurationsGLVOutput"),
                                            ),
                                            tags$hr()
                                        ),
                                        sliderInput("migrationPGLV", "probability/frequency of migration from metacommunity", value = 0.01, min = 0, max = 1, step = 0.01),
                                        conditionalPanel(
                                            condition = "input.migrationPGLV >0",
                                            sliderInput("sigmaMigrationGLV", "intensity of migration", value = 0.01, min = 0, max = 1, step = 0.001),
                                        ),
                                        textInput(
                                            "metacommunityProbabilityGLV", 
                                            "metacommunity"
                                        ) %>%
                                        shinyInput_label_embed(
                                            shiny_iconlink() %>% 
                                                bs_embed_tooltip(
                                                    title =  "Normalized probability distribution of the likelihood that species from the metacommunity can enter the community during the simulation.", 
                                                    placement = "right", 
                                                    container = "body"
                                                )
                                        ),
                                        verbatimTextOutput("metacommunityProbabilityGLV"),
                                        tags$hr(),
                                        switchInput(
                                            "normGLV", 
                                            strong("returns normalized abundances"), 
                                            value = FALSE,
                                            labelWidth = "100%"
                                        ),
                                        
                                    ),
                                ),
                                tags$hr(),
                                tags$h4(
                                    "Press the following button to run the model",
                                    tags$div(
                                        class = "pull-right",
                                        shiny_iconlink() %>% 
                                            bs_embed_tooltip(
                                                title =  "GLV model was not designed responsive to reduce the calculation", 
                                                placement = "left", 
                                                container = "body"
                                            ),
                                    ),
                                ),
                                actionButton("buttonSimulateGLV", "Run the GLV Model", class = "btn btn-primary", width = "100%"),
                            )
                        ),
                        ### Display Panel ####
                        column(
                            width = 7,
                            
                            #### Matrix A ####
                            fluidRow(
                                style = "padding-left: 15px; padding-right: 15px;",
                                tags$div(
                                    class = "panel panel-default",
                                    tags$div(
                                        class = "panel-heading",
                                        tags$h3(
                                            class = "panel-title",
                                            "Matrix of interspecies interactions",
                                        ),
                                    ),
                                    tags$div(
                                        class = "panel-body",
                                        dataTableOutput("TableAGLV", width = "100%"),
                                        plotOutput("GLVPlotA", height = "600px"),
                                    ),
                                ),
                            ),
                            #### result plots ####
                            tags$br(),
                            fluidRow(
                                style = "padding-left: 15px; padding-right: 15px;",
                                tags$div(
                                    class = "panel panel-default",
                                    tags$div(
                                        class = "panel-heading",
                                        tags$h3(
                                            class = "panel-title",
                                            "Species Change",
                                        ),
                                    ),
                                    tags$div(
                                        class = "panel-body",
                                        plotOutput("GLVSpecies"),
                                    ),
                                ),
                            ),
                        ),
                    ) 
                
            ) %>%
            bs_append(
                ## Description ####
                title = "Description",
                content = 
                    fluidRow(
                        column(
                            width = 12,
                            withMathJax(includeMarkdown("glv.Rmd")),
                        ),
                    )
            ) %>% 
            bs_append(
                ## Inputs ####
                title = "Inputs",
                content = 
                    fluidRow(
                        "Inputs of GLV model"
                    )
            ) %>% 
            bs_append(
                ## References ####
                title = "References",
                content = "Panel of refs."
            )
    ),
    # tab3 Hubbell Model ####
    tabPanel(
        title = "Hubbell Model (with death rates)",
        titlePanel("Hubbell Model (with death rates)"),
        ## Hubbell Model ####
        
    ),
    
    # tab4 Logistic Model (with stochasticity) ####
    tabPanel(
        title = "Logistic Model (with stochasticity)",
        titlePanel("Logistic Model (with stochasticity)"),
        ## Logistic Model ####
    )
)
