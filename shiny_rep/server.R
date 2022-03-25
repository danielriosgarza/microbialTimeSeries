# source all files ####
Rfiles = gsub(" ", "", paste("./R/", list.files("./R/")))
sapply(Rfiles, source)

# converting functions ####
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

# server ####
server <- function(input, output, session) {
    # model1 simulate consumer resource model ####
    ## basic ####
    n.species.crm <- reactive(input$nSpeciesCRM)
    n.resources.crm <- reactive(input$nResourcesCRM)
    names.species.crm <- reactive(text2chars(input$namesSpeciesCRM, len = n.species.crm(), prefix = "sp"))
    names.resources.crm <- reactive(text2chars(input$namesResourcesCRM, len = n.resources.crm(), prefix = "res"))
    t.start.crm <- reactive(input$tStartCRM)
    t.end.crm <- reactive(input$tEndCRM)
    observeEvent(input$tStartCRM | input$tEndCRM, {
        updateNumericInput(inputId = "tStartCRM", max = input$tEndCRM)
        updateNumericInput(inputId = "tEndCRM", min = input$tStartCRM)
    })
    
    t.step.crm <- reactive(input$tStepCRM)
    t.store.crm <- reactive(input$tStoreCRM)
    observeEvent(input$tStartCRM | input$tEndCRM | input$tStepCRM | input$tStoreCRM, {
        updateNumericInput(inputId = "tStepCRM", max = (input$tEndCRM-input$tStartCRM)/input$tStoreCRM)
        updateNumericInput(inputId = "tStoreCRM", max = (input$tEndCRM-input$tStartCRM)/input$tStepCRM)
    })
    
    ## compounds ####
    res.conc.crm <- reactive(input$resourcesConcentrationCRM)
    res.even.crm <- reactive(input$resourcesEvennessCRM)
    resources_dist.crm <- reactive(rdirichlet(1, rep(1, n.resources.crm())*res.even.crm())*res.conc.crm()*n.resources.crm())
    res.custom.crm <- reactive(as.numeric(as.vector(text2char(input$resourcesCustomCRM))))
    resources.crm <- reactive({
        if (length(res.custom.crm()) < n.resources.crm()){
            return(c(res.custom.crm(), as.vector(resources_dist.crm())[(length(res.custom.crm())+1):n.resources.crm()]))
        } else {
            return(head(res.custom.crm(), n.resources.crm()))
        }
    })
    output$resourcesOutputCRM <- renderPrint(resources.crm())
    
    ### dilution/influx and outflux ####
    inflow.rate.crm <- reactive(input$inflowRateCRM)
    outflow.rate.crm <- reactive(input$outflowRateCRM)
    volume.crm <- reactive(input$volumeCRM)
    res.dilu.crm <- reactive(as.numeric(as.vector(text2char(input$resourcesDilutionCRM))))
    resources.dilution.crm <- reactive({
        if (length(res.dilu.crm()) < n.resources.crm()){
            return(c(res.dilu.crm(), resources.crm()[(length(res.dilu.crm())+1):n.resources.crm()]))
        } else {
            return(head(res.dilu.crm(), n.resources.crm()))
        }
    })
    output$resourcesDilutionCRMOutput <- renderPrint(resources.dilution.crm())

    output$resourcesCRMPlot <- renderPlot(makePiePlot(resources.crm(), label = 'concentration', title='compounds'))
    mean.consumption.crm <- reactive(input$meanConsumptionCRM)
    mean.production.crm <- reactive(input$meanProductionCRM)
    observeEvent(input$meanConsumptionCRM | input$meanProductionCRM, {
        updateSliderInput(inputId = "meanProductionCRM", max = 1 - input$meanConsumptionCRM)
        updateSliderInput(inputId = "meanConsumptionCRM", max = 1 - input$meanProductionCRM)
    })
    maintenance.crm <- reactive(input$maintenanceCRM)
    
    ### editable matrixECRM and matrixMonodCRM ####
    RV.crm <- reactiveValues(matrixECRM = NULL, matrixMonodCRM = NULL)
    observe({
        roundECRM <- round(
            randomE(
                n.species = n.species.crm(), 
                n.resources = n.resources.crm(), 
                names.species = names.species.crm(), 
                names.resources = names.resources.crm(), 
                mean.consumption = as.integer(mean.consumption.crm() * n.resources.crm()),
                mean.production = as.integer(mean.production.crm() * n.resources.crm()),
                maintenance = maintenance.crm()
            ),
            digits = 3
        )
        RV.crm$matrixECRM <- roundECRM
    })
    output$tableECRM <- renderDataTable(RV.crm$matrixECRM, editable = 'all', selection = 'none', server = TRUE, options = list(scrollX = TRUE))
    output$CRMPlotE <- renderPlot(makeHeatmap(RV.crm$matrixECRM, 'Consumption/production matrix'), res = 96)
    
    observeEvent(input$tableECRM_cell_edit, {
        RV.crm$matrixECRM <<- editData(RV.crm$matrixECRM, input$tableECRM_cell_edit, 'tableECRM')
    })
    
    
    ## growth rates ####
    x0.crm <- reactive(as.numeric(text2chars(input$x0CRM, len = n.species.crm(), expr = paste0("runif(n = ", n.species.crm() ,", min = 0.1, max = 10)"))))
    output$x0CRMOutput <- renderPrint(x0.crm())
    
    
    alpha.crm <- reactive(input$alphaCRM)
    beta.crm <- reactive(input$betaCRM)
    # listening to the changes in nSpeciesCRM, alphaCRM, and betaCRM
    observeEvent(input$nSpeciesCRM | input$alphaCRM | input$betaCRM, {
        growth.rates.crm <- reactive(as.numeric(text2chars(input$growthRatesCRM, len = n.species.crm(), expr = paste0('round(rbeta(',n.species.crm(), ',' ,alpha.crm(), ',' , beta.crm(),'), digits = 3)'))))
    })
    observeEvent(input$nSpeciesCRM, {
        x0.crm <- reactive(as.numeric(text2chars(input$x0CRM, len = n.species.crm(), expr = paste0("runif(n = ", n.species.crm() ,", min = 0.1, max = 10)"))))
    })
    observeEvent(input$buttonBetaEven, {
        updateSliderInput(inputId = "alphaCRM", value = 1)
        updateSliderInput(inputId = "betaCRM", value = 1)
    })
    observeEvent(input$buttonBetaRidge, {
        updateSliderInput(inputId = "alphaCRM", value = 4)
        updateSliderInput(inputId = "betaCRM", value = 4)
    })
    observeEvent(input$buttonBetaValley, {
        updateSliderInput(inputId = "alphaCRM", value = 0.5)
        updateSliderInput(inputId = "betaCRM", value = 0.5)
    })
    observeEvent(input$buttonBetaLeft, {
        updateSliderInput(inputId = "alphaCRM", value = 0.5)
        updateSliderInput(inputId = "betaCRM", value = 1)
    })
    observeEvent(input$buttonBetaRight, {
        updateSliderInput(inputId = "alphaCRM", value = 1)
        updateSliderInput(inputId = "betaCRM", value = 0.5)
    })
    observeEvent(input$buttonBetaLeftTriangle, {
        updateSliderInput(inputId = "alphaCRM", value = 1)
        updateSliderInput(inputId = "betaCRM", value = 2)
    })
    observeEvent(input$buttonBetaRightTriangle, {
        updateSliderInput(inputId = "alphaCRM", value = 2)
        updateSliderInput(inputId = "betaCRM", value = 1)
    })
    
    growth.rates.crm <- reactive(as.numeric(text2chars(input$growthRatesCRM, len = n.species.crm(), expr = paste0('round(rbeta(',n.species.crm(), ',' ,alpha.crm(), ',' , beta.crm(),'), digits = 3)'))))
    output$growthRatesCRMOutput <- renderPrint(growth.rates.crm())
    
    output$growthRatesCRMDist <- renderPlot({
        ggplot(data.frame(x = c(-4, 4)), aes(x = x)) + 
            stat_function(fun = dbeta, args = list(shape1=alpha.crm(), shape2=beta.crm())) + 
            xlim(0,1) + theme_linedraw()
    }, res = 96) 
    observe({
        roundMonodCRM <- round(
            matrix(
                rgamma(
                    n = n.species.crm()*n.resources.crm(), 
                    shape = 50*max(resources.crm()), 
                    rate = 1), 
                nrow = n.species.crm(),
                dimnames = list(names.species.crm(), names.resources.crm())
                ), 
            digits = 3)
        RV.crm$matrixMonodCRM <- roundMonodCRM
    })
    output$tableMonodCRM <- renderDataTable(RV.crm$matrixMonodCRM, editable = 'all', selection = 'none', server = TRUE, options = list(scrollX = TRUE))
    observeEvent(input$tableMonodCRM_cell_edit, {
        RV.crm$matrixMonodCRM <<- editData(RV.crm$matrixMonodCRM, input$tableMonodCRM_cell_edit, 'tableMonodCRM')
    })
    ## perturbation ####
    error.variance.crm <- reactive(input$errorVarianceCRM)
    stochastic.crm <- reactive(input$stochasticCRM)
    sigma.drift.crm <- reactive(input$sigmaDriftCRM)
    sigma.epoch.crm <- reactive(input$sigmaEpochCRM)
    sigma.external.crm <- reactive(input$sigmaExternalCRM)
    sigma.migration.crm <- reactive(input$sigmaMigrationCRM)
    epoch.p.crm <- reactive(input$epochPCRM)
    t.external_events.crm <- reactive(as.numeric(text2char(input$tExternalEventsCRM)))
    output$tExternalEventsCRMOutput <- renderPrint(t.external_events.crm())
    t.external_durations.crm <- reactive(as.numeric(text2char(input$tExternalDurationsCRM)))
    output$tExternalDurationsCRMOutput <- renderPrint(t.external_durations.crm())
    migration.p.crm <- reactive(input$migrationPCRM)
    metacommunity.probability.crm <- reactive(as.numeric(text2chars(input$metacommunityProbabilityCRM, len = n.species.crm(), expr = paste0("rdirichlet(1, alpha = rep(1,", n.species.crm(), "))"))))
    output$metacommunityProbabilityCRM <- renderPrint(metacommunity.probability.crm())
    
    
    norm.crm <- reactive(input$normCRM)
    
    ## examples CRM ####
    observeEvent(input$CRMEX1, {
        updateSliderInput(inputId = "nSpeciesCRM", value = 5)
        updateSliderInput(inputId = "nResourcesCRM", value = 5)
    })
    observeEvent(input$CRMEX2, {
        updateNumericInput(inputId = "tEndCRM", value = 2000)
        updateNumericInput(inputId = "tStoreCRM", value = 500)
        updateSliderInput(inputId = "migrationPCRM", value = 0)
        updateCheckboxInput(inputId = "stochasticCRM", value = FALSE)
        updateSliderInput(inputId = "dilutionRateCRM", value = 0.001)
    })
    observeEvent(input$CRMEX3, {
        updateSliderInput(inputId = "nSpeciesCRM", value = 1)
        updateSliderInput(inputId = "nResourcesCRM", value = 10)
        updateSliderInput(inputId = "maintenanceCRM", value = 0.1)
    })
    observeEvent(input$CRMEX4pre, {
        updateSliderInput(inputId = "nSpeciesCRM", value = 3)
        updateSliderInput(inputId = "nResourcesCRM", value = 4)
        updateTextInput(inputId = "growthRatesCRM", value = "2, 4.5, 2.6")
        updateTextInput(inputId = "x0CRM", value = "1, 2, 1")
        updateTextInput(inputId = "resourcesCustomCRM", value = "10, 0, 0, 0")
        updateTextInput(inputId = "namesSpeciesCRM", value = "homoacetogenic, homofermentative, butyrateProducer")
        updateTextInput(inputId = "namesResourcesCRM", value = "glucose, acetate, lactate, butyrate")
        # updateButton(session, "CRMEX4", disabled = !input$CRMEX4pre)
        shinyjs::enable("CRMEX4")
    })
    observeEvent(input$CRMEX4, {
        # auto update matrixECRM, then take changes after it.
        matrixExample4 <- matrix(c(1, -3, 0, 0, 1, 0, -2, 0, 0, 0, 4, -3), nrow = 3, byrow = TRUE)*c(4.3/4, 2/4, 1/4)
        RV.crm$matrixECRM <- matrixExample4
    })
    observeEvent(input$CRMEX5pre, {
        updateSliderInput(inputId = "nSpeciesCRM", value = 10)
        updateSliderInput(inputId = "nResourcesCRM", value = 10)
        # updateButton(session, "CRMEX5", disabled = !input$CRMEX5pre)
        shinyjs::enable("CRMEX5")
    })
    observeEvent(input$CRMEX5, {
        RV.crm$matrixECRM <- randomE(n.species = 10, n.resources = 10, mean.consumption = 3, mean.production = 1,
                              maintenance = 0.5, trophic.preferences = list(c(5,3,1,1,1,1,1,1,1,1)))
    })
    observeEvent(input$CRMEX6pre, {
        updateSliderInput(inputId = "nSpeciesCRM", value = 20)
        updateSliderInput(inputId = "nResourcesCRM", value = 20)
        # updateButton(session, "CRMEX6", disabled = !input$CRMEX6pre)
        shinyjs::enable("CRMEX6")
    })
    observeEvent(input$CRMEX6, {
        RV.crm$matrixECRM <- randomE(n.species = 20, n.resources = 20, mean.consumption = 3, mean.production = 2, maintenance = 0.0, trophic.levels = c(7, 13))
    })
    observeEvent(input$CRMEX7pre, {
        updateSliderInput(inputId = "nSpeciesCRM", value = 4)
        updateSliderInput(inputId = "nResourcesCRM", value = 11)
        updateTextInput(inputId = "namesSpeciesCRM", value = "A, B, C, D")
        updateTextInput(inputId = "x0CRM", value = "1, 1, 1, 1")
        updateTextInput(inputId = "resourcesCustomCRM", value = "1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5")
        updateCheckboxInput(inputId = "stochasticCRM", value = FALSE)
        updateSliderInput(inputId = "migrationPCRM", value = 0)
        updateSliderInput(inputId = "dilutionRateCRM", value = 0)
        # updateButton(session, "CRMEX7", disabled = !input$CRMEX7pre)
        shinyjs::enable("CRMEX7")
    })
    observeEvent(input$CRMEX7, {
        #secretion of C
        sec.C <- rdirichlet(1, c(1,1,1))*.5
        #The metabolic preferences of A are set to the secretion products of C
        pref.A.D <- list(c(sec.C*1000, rep(1,8)))
        em.A <- randomE(n.species = 1, n.resources = 11, names.species = 'A', trophic.preferences = pref.A.D, mean.production = 3, mean.consumption = 3)
        #secretion of A
        sec.A <- abs(em.A*(em.A<0))
        #The metabolic preferences of D are set to the secretion products of A
        em.D <- randomE(n.species = 1, n.resources = 11, names.species = 'D', trophic.preferences = pref.A.D, mean.production = 3, mean.consumption = 3)
        #secretion of D
        sec.D <- abs(em.D*(em.D<0))
        pref.B <- 1000*((sec.A + sec.D)/(sum(sec.A)+sum(sec.D)))
        pref.B[pref.B==0] <- 1
        pref.B <- list(pref.B[4:11])
        em.B <- randomE(n.species = 1, n.resources = 8, names.species = 'B', trophic.preferences = pref.B, mean.production = 3, mean.consumption = 3)
        #secretion of B
        sec.B <- abs(em.B*(em.B<0))
        #The metabolic preferences of C are set to the secretion products B
        pref.C <- sec.B*1000
        pref.C[pref.C==0] <- 1
        em.B <-t(as.matrix(c(rep(0,3),em.B)))
        row.names(em.B) = 'B'
        em.C <- randomE(n.species = 1, n.resources = 8, names.species = 'C', trophic.preferences = list(pref.C), mean.production = 0, mean.consumption = 3)
        em.C <- cbind(-sec.C, em.C)
        RV.crm$matrixECRM <- rbind(em.A, em.B, em.C, em.D)
    })
    
    ## runCRM ####
    runCRM <- reactive({
    # runCRM <- eventReactive(input$buttonSimulateCRM, {
        simulateConsumerResource(
            n.species = n.species.crm(),
            n.resources = n.resources.crm(),
            names.species = names.species.crm(), 
            names.resources = names.resources.crm(), 
            E = RV.crm$matrixECRM,
            x0 = x0.crm(), 
            resources = resources.crm(), 
            resources.dilution = resources.dilution.crm(),
            inflow.rate = inflow.rate.crm(),
            outflow.rate = outflow.rate.crm(),
            volume = volume.crm(),
            growth.rates = growth.rates.crm(), 
            monod.constant = RV.crm$matrixMonodCRM, 
            sigma.drift = sigma.drift.crm(),
            sigma.epoch = sigma.epoch.crm(),
            sigma.external = sigma.external.crm(),
            sigma.migration = sigma.migration.crm(), 
            epoch.p = epoch.p.crm(), 
            t.external_events = t.external_events.crm(),
            t.external_durations = t.external_durations.crm(),
            stochastic = stochastic.crm(),
            migration.p = migration.p.crm(),
            metacommunity.probability = metacommunity.probability.crm(),
            error.variance = error.variance.crm(), 
            norm = norm.crm(),
            t.start = t.start.crm(),
            t.end = t.end.crm(),
            t.step = t.step.crm(),
            t.store = t.store.crm()
        )}
    )
    
    output$CRMSpecies <- renderPlot(makePlot(runCRM()$matrix, "abundance of species by time"), res = 96)
    output$CRMResources <- renderPlot(makePlotRes(runCRM()$resources, "quantity of compounds by time"),  res = 96)
    output$CRMVolume <- renderPlot(makePlot(runCRM()$volume, "volume changes of the reactor by time"), res = 96)
    
    # model2 simulate generalized Lotka-Volterra Model ####
    ## interspecies interactions ####
    n.species.glv <- reactive(input$nSpeciesGLV)
    names.species.glv <- reactive(text2chars(input$namesSpeciesGLV, len = n.species.glv(), prefix = "sp"))
    diagonal.glv <- reactive(input$diagonalGLV)
    connectance.glv <- reactive(input$connectanceGLV)
    scale.glv <- reactive(input$scaleGLV)
    
    mutualism.glv <- reactive(input$mutualismGLV)
    commensalism.glv <- reactive(input$commensalismGLV)
    parasitism.glv <- reactive(input$parasitismGLV)
    amensalism.glv <- reactive(input$amensalismGLV)
    competition.glv <- reactive(input$competitionGLV)
    
    observeEvent(input$buttonInteractionsU, {
        updateTextInput(inputId = "interactionsCustomGLV", value = round(rbeta(n.species.glv()^2, 0.5, 0.5), digits = 3))
    })
    observeEvent(input$buttonInteractionsN, {
        updateTextInput(inputId = "interactionsCustomGLV", value = round(rbeta(n.species.glv()^2, 2, 2), digits = 3))
    })
    observeEvent(input$buttonInteractionsEven, {
        updateTextInput(inputId = "interactionsCustomGLV", value = round(runif(n.species.glv()^2, 0, 1), digits = 3))
    })
    interactions_dist.glv <- reactive(round(runif(n.species.glv()^2, 0, 1), digits = 3))
    inter.custom.glv <- reactive(as.numeric(as.vector(text2char(input$interactionsCustomGLV))))
    interactions.glv <- reactive({
        if (length(inter.custom.glv()) < n.species.glv()^2){
            return(
                c(
                    inter.custom.glv(), 
                    as.vector(
                        interactions_dist.glv()[(length(inter.custom.glv())+1):(n.species.glv()^2)]
                    )
                )
            )
        } else {
            return(head(inter.custom.glv(), n.species.glv()^2))
        }
    })
    output$interactionsOutputGLV <- renderPrint(interactions.glv())
    symmetric.glv <- reactive(input$symmetricGLV)
    listA.glv <- reactive(text2char(input$listAGLV)) # TODO: convert listA
    
    RV.glv <- reactiveValues(matrixAGLV = NULL)
    # replace generateA() by RV.glv$matrixAGLV
    
    #### editable matrixA ####
    observe({
        roundACRM <- round(
            randomA(
                n.species = n.species.glv(),
                names.species = names.species.glv(), 
                diagonal = diagonal.glv(),
                connectance = connectance.glv(),
                scale = scale.glv(),
                mutualism = mutualism.glv(),
                commensalism = commensalism.glv(),
                parasitism = parasitism.glv(),
                amensalism = amensalism.glv(),
                competition = competition.glv(),
                interactions = interactions.glv(),
                symmetric = symmetric.glv()
            ),
            digits = 3
        )
        RV.glv$matrixAGLV <- roundACRM
    })
    output$TableAGLV <- renderDataTable(RV.glv$matrixAGLV, editable = 'all', selection = 'none', server = TRUE, options = list(scrollX = TRUE))
    output$GLVPlotA <- renderPlot(makeHeatmap(RV.glv$matrixAGLV, "interspecies interaction matrix"), res = 96)
    
    observeEvent(input$tableAGLV_cell_edit, {
        RV.glv$matrixAGLV <<- editData(RV.glv$matrixAGLV, input$tableAGLV_cell_edit, 'tableAGLV')
    })
    
    ## growth rates ####
    x0.glv <- reactive(as.numeric(text2chars(input$x0GLV, len = n.species.glv(), expr = paste0("runif(n =", n.species.glv(), ", min = 0, max = 1)"))))
    output$x0GLVOutput <- renderPrint(x0.glv())
    growth.rates.glv <- reactive(as.numeric(text2chars(input$growthRatesGLV, len = n.species.glv(), expr = paste0("runif(n =", n.species.glv(), ", min = 0, max = 1)"))))
    output$growthRatesGLVOutput <- renderPrint(growth.rates.glv())
    ## perturbations ####
    stochastic.glv <- reactive(input$stochasticGLV)
    sigma.drift.glv <- reactive(input$sigmaDriftGLV)
    sigma.epoch.glv <- reactive(input$sigmaEpochGLV)
    sigma.external.glv <- reactive(input$sigmaExternalGLV)
    epoch.p.glv <- reactive(input$epochPGLV)
    sigma.migration.glv <- reactive(input$sigmaMigrationGLV)
    t.external_events.glv <- reactive(as.numeric(text2char(input$tExternalEventsGLV)))
    output$tExternalEventsGLVOutput <- renderPrint(t.external_events.glv())
    t.external_durations.glv <- reactive(as.numeric(text2char(input$tExternalDurationsGLV)))
    output$tExternalDurationsGLVOutput <- renderPrint(t.external_durations.glv())
    migration.p.glv <- reactive(input$migrationPGLV)
    metacommunity.probability.glv <- reactive(as.numeric(text2chars(input$metacommunityProbabilityGLV, len = n.species.crm(), expr = paste0("rdirichlet(1, alpha = rep(1,", n.species.crm(), "))"))))
    output$metacommunityProbabilityGLV <- renderPrint(metacommunity.probability.glv())
    
    error.variance.glv <- reactive(input$errorVarianceGLV)
    t.start.glv <- reactive(input$tStartGLV)
    t.end.glv <- reactive(input$tEndGLV)
    t.step.glv <- reactive(input$tStepGLV)
    t.store.glv <- reactive(input$tStoreGLV)
    
    norm.glv <- reactive(input$normGLV)
    
    ## examples GLV ####
    observeEvent(input$GLVEX1, {
        updateSliderInput(inputId = "nSpeciesGLV", value = 5)
    })
    observeEvent(input$GLVEX2, {
        updateSliderInput(inputId = "nSpeciesGLV", value = 4)
        updateSliderInput(inputId = "diagonalGLV", value = -1)
        updateSliderInput(inputId = "connectanceGLV", value = 0.5)
        updateSliderInput(inputId = "scaleGLV", value = 0.5)
        updateSwitchInput(inputId = "symmetricGLV", value = TRUE)
        updateSwitchInput(inputId = "stochasticGLV", value = FALSE)
    })
    observeEvent(input$GLVEX3, {
        updateSliderInput(inputId = "nSpeciesGLV", value = 4)
        updateSliderInput(inputId = "diagonalGLV", value = -1)
        updateSliderInput(inputId = "connectanceGLV", value = 0.5)
        updateSliderInput(inputId = "scaleGLV", value = 0.5)
        updateSwitchInput(inputId = "symmetricGLV", value = TRUE)
        updateSwitchInput(inputId = "stochasticGLV", value = FALSE)
        updateSliderInput(inputId = "migrationPCRM", value = 0)
    })
    observeEvent(input$GLVEX4, {
        updateSliderInput(inputId = "nSpeciesGLV", value = 4)
        updateSliderInput(inputId = "diagonalGLV", value = -1)
        updateSliderInput(inputId = "connectanceGLV", value = 0.5)
        updateSliderInput(inputId = "scaleGLV", value = 0.5)
        updateSwitchInput(inputId = "symmetricGLV", value = TRUE)
        updateSwitchInput(inputId = "stochasticGLV", value = FALSE)
        updateSliderInput(inputId = "migrationPCRM", value = 0)
        updateSliderInput(inputId = "errorVarianceGLV", value = 0.001)
    })
    
    ## runGLV ####
    
    # use button or not?
    runGLV <- reactive({
    # runGLV <- eventReactive(input$buttonSimulateGLV, {
        simulateGLV(
            n.species = n.species.glv(),
            names.species = names.species.glv(), 
            A = RV.glv$matrixAGLV,
            x0 = x0.glv(), 
            growth.rates = growth.rates.glv(), 
            sigma.drift = sigma.drift.glv(),
            sigma.epoch = sigma.epoch.glv(),
            sigma.external = sigma.external.glv(),
            sigma.migration = sigma.migration.glv(),
            epoch.p = epoch.p.glv(),
            t.external_events = t.external_events.glv(),
            t.external_durations = t.external_durations.glv(),
            stochastic = stochastic.glv(),
            migration.p = migration.p.glv(),
            metacommunity.probability = metacommunity.probability.glv(),
            error.variance = error.variance.glv(),
            norm = norm.glv(), 
            t.start = t.start.glv(),
            t.end = t.end.glv(),
            t.step = t.step.glv(),
            t.store = t.store.glv()
        )
    })
    output$GLVSpecies <- renderPlot(makePlot(runGLV()$matrix))
    
    # model3 simulate Hubbell neutral model with growth rates ####
    ## basic ####
    n.species.hub <- reactive(input$nSpeciesHUB)
    x0.hub <- reactive(as.numeric(text2chars(input$x0HUB, len = n.species.hub(), expr = paste0("rep(100, ", n.species.hub() ,")"))))
    output$x0HUBOutput <- renderPrint(x0.hub())
    names.species.hub <- reactive(text2chars(input$namesSpeciesHUB, len = n.species.hub(), prefix = "sp"))
    growth.rates.hub <- reactive(as.numeric(text2chars(input$growthRatesHUB, len = n.species.hub(), expr = paste0('rep(1, ',n.species.hub(), ')'))))
    output$growthRatesHUBOutput <- renderPrint(growth.rates.hub())
    
    t.start.hub <- reactive(input$tStartHUB)
    t.end.hub <- reactive(input$tEndHUB)
    observeEvent(input$tStartHUB | input$tEndHUB, {
        updateNumericInput(inputId = "tStartHUB", max = input$tEndHUB)
        updateNumericInput(inputId = "tEndHUB", min = input$tStartHUB)
    })
    
    t.step.hub <- reactive(input$tStepHUB)
    t.store.hub <- reactive(input$tStoreHUB)
    observeEvent(input$tStartHUB | input$tEndHUB | input$tStepHUB | input$tStoreHUB, {
        updateNumericInput(inputId = "tStepHUB", max = (input$tEndHUB-input$tStartHUB)/input$tStoreHUB)
        updateNumericInput(inputId = "tStoreHUB", max = (input$tEndHUB-input$tStartHUB)/input$tStepHUB)
    })
    
    ## perturbations ####
    error.variance.hub <- reactive(input$errorVarianceHUB)
    k.events.hub <- reactive(input$kEventsHUB)
    migration.p.hub <- reactive(input$migrationPHUB)
    metacommunity.probability.hub <- reactive(as.numeric(text2chars(input$metacommunityProbabilityHUB, len = n.species.hub(), expr = paste0("rdirichlet(1, alpha = rep(1,", n.species.hub(), "))"))))
    output$metacommunityProbabilityHUB <- renderPrint(metacommunity.probability.hub())
    
    norm.hub <- reactive(input$normHUB)
    
    ## examples HUB ####
    observeEvent(input$HUBEX1, {
        updateSliderInput(inputId = "nSpeciesHUB", value = 5)
    })
    observeEvent(input$HUBEX2, {
        updateSliderInput(inputId = "nSpeciesHUB", value = 5)
        updateSliderInput(inputId = "migrationPHUB", value = 0)
    })
    observeEvent(input$HUBEX3, {
        updateSliderInput(inputId = "nSpeciesHUB", value = 5)
        updateSliderInput(inputId = "migrationPHUB", value = 1)
        updateTextInput(inputId = "metacommunityProbabilityHUB", value = "0.1, 0.15, 0.2, 0.25, 0.3")
        updateTextInput(inputId = "tEndHUB", value = 20)
        updateTextInput(inputId = "tStoreHUB", value = 200)
    })
    observeEvent(input$HUBEX4, {
        updateSliderInput(inputId = "nSpeciesHUB", value = 5)
        updateSliderInput(inputId = "migrationPHUB", value = 1)
        updateTextInput(inputId = "metacommunityProbabilityHUB", value = "0.1, 0.15, 0.2, 0.25, 0.3")
        updateTextInput(inputId = "tEndHUB", value = 20)
        updateTextInput(inputId = "tStoreHUB", value = 200)
        updateSliderInput(inputId = "errorVarianceHUB", value = 100)
    })
    observeEvent(input$HUBEX5, {
        updateSliderInput(inputId = "nSpeciesHUB", value = 5)
        updateSliderInput(inputId = "migrationPHUB", value = 0.1)
        updateTextInput(inputId = "metacommunityProbabilityHUB", value = "0.1, 0.15, 0.2, 0.25, 0.3")
        updateTextInput(inputId = "tEndHUB", value = 20)
        updateTextInput(inputId = "tStoreHUB", value = 1000)
        updateSliderInput(inputId = "kEventsHUB", value = 5)
        updateTextInput(inputId = "growthRatesHUB", value = "1.1, 1.05, 1, 0.95, 0.9")
    })
    
    ## runHUB ####
    runHUB <- reactive({
    # runHUB <- eventReactive(input$buttonSimulateHUB, {
        simulateHubbellRates(
            n.species = n.species.hub(),
            x0 = x0.hub(),
            names.species = names.species.hub(),
            migration.p = migration.p.hub(),
            metacommunity.probability = metacommunity.probability.hub(),
            k.events = k.events.hub(),
            growth.rates = growth.rates.hub(),
            error.variance = error.variance.hub(),
            norm = norm.hub(),
            t.start = t.start.hub(),
            t.end = t.end.hub(),
            t.step = t.step.hub(),
            t.store = t.store.hub()
        )
    })
    output$HUBSpecies <- renderPlot(makePlot(runHUB()$matrix, "abundance of species by time"), res = 96)
    
    # model4 simulate stochastic logistic model ####
    ## basic ####
    n.species.log <- reactive(input$nSpeciesLOG)
    names.species.log <- reactive(text2chars(input$namesSpeciesLOG, len = n.species.log(), prefix = "sp"))
    t.start.log <- reactive(input$tStartLOG)
    t.end.log <- reactive(input$tEndLOG)
    observeEvent(input$tStartLOG | input$tEndLOG, {
        updateNumericInput(inputId = "tStartLOG", max = input$tEndLOG)
        updateNumericInput(inputId = "tEndLOG", min = input$tStartLOG)
    })
    
    t.step.log <- reactive(input$tStepLOG)
    t.store.log <- reactive(input$tStoreLOG)
    observeEvent(input$tStartLOG | input$tEndLOG | input$tStepLOG | input$tStoreLOG, {
        updateNumericInput(inputId = "tStepLOG", max = (input$tEndLOG-input$tStartLOG)/input$tStoreLOG)
        updateNumericInput(inputId = "tStoreLOG", max = (input$tEndLOG-input$tStartLOG)/input$tStepLOG)
    })
    
    ## growth and death rates ####
    x0.log <- reactive(as.numeric(text2chars(input$x0LOG, len = n.species.log(), expr = paste0("rep(100, ", n.species.log() ,")"))))
    output$x0LOGOutput <- renderPrint(x0.log())
    growth.rates.log <- reactive(as.numeric(text2chars(input$growthRatesLOG, len = n.species.log(), expr = paste0("runif(n = ", n.species.log() ,", min = 0.1, max = 0.2)"))))
    output$growthRatesLOGOutput <- renderPrint(growth.rates.log())
    death.rates.log <- reactive(as.numeric(text2chars(input$deathRatesLOG, len = n.species.log(), expr = paste0("runif(n = ", n.species.log() ,", min = 0.0005, max = 0.0025)"))))
    output$deathRatesLOGOutput <- renderPrint(death.rates.log())
    carrying.capacities.log <- reactive(as.numeric(text2chars(input$carryingCapacitiesLOG, len = n.species.log(), expr = paste0("runif(n = ", n.species.log() ,", min = 1000, max = 2000)"))))
    output$carryingCapacitiesLOGOutput <- renderPrint(carrying.capacities.log())
    
    ## perturbations ####
    error.variance.log <- reactive(input$errorVarianceLOG)
    stochastic.log <- reactive(input$stochasticLOG)
    sigma.drift.log <- reactive(input$sigmaDriftLOG)
    epoch.p.log <- reactive(input$epochPLOG)
    sigma.epoch.log <- reactive(input$sigmaEpochLOG)
    sigma.external.log <- reactive(input$sigmaExternalLOG)
    t.external_events.log <- reactive(as.numeric(text2char(input$tExternalEventsLOG)))
    output$tExternalEventsLOGOutput <- renderPrint(t.external_events.log())
    t.external_durations.log <- reactive(as.numeric(text2char(input$tExternalDurationsLOG)))
    output$tExternalDurationsLOGOutput <- renderPrint(t.external_durations.log())
    migration.p.log <- reactive(input$migrationPLOG)
    sigma.migration.log <- reactive(input$sigmaMigrationLOG)
    metacommunity.probability.log <- reactive(as.numeric(text2chars(input$metacommunityProbabilityLOG, len = n.species.log(), expr = paste0("rdirichlet(1, alpha = rep(1,", n.species.log(), "))"))))
    output$metacommunityProbabilityLOG <- renderPrint(metacommunity.probability.log())
    
    norm.log <- reactive(input$normLOG)
    
    
    ## examples LOG ####
    observeEvent(input$LOGEX1, {
        updateSliderInput(inputId = "nSpeciesLOG", value = 5)
        updateSwitchInput(inputId = "stochasticLOG", value = FALSE)
    })
    
    observeEvent(input$LOGEX2, {
        updateSliderInput(inputId = "nSpeciesLOG", value = 5)
        updateSwitchInput(inputId = "stochasticLOG", value = FALSE)
        updateTextInput(inputId = "deathRatesLOG", value = "0,0,0,0,0")
    })
    
    observeEvent(input$LOGEX3, {
        updateSliderInput(inputId = "nSpeciesLOG", value = 5)
        updateSwitchInput(inputId = "stochasticLOG", value = FALSE)
        updateTextInput(inputId = "deathRatesLOG", value = "0,0,0,0,0")
        updateTextInput(inputId = "growthRatesLOG", value = "0.1, 0.2, 0.3, 0.4, 0.5")
        updateTextInput(inputId = "deathRatesLOG", value = "0.001, 0.0008, 0.0006, 0.0004, 0.0002")
        updateTextInput(inputId = "carryingCapacitiesLOG", value = "1000, 1200, 1400, 1600, 1800")
        
    })
    
    observeEvent(input$LOGEX4, {
        updateSliderInput(inputId = "nSpeciesLOG", value = 5)
        updateTextInput(inputId = "deathRatesLOG", value = "0,0,0,0,0")
        updateTextInput(inputId = "growthRatesLOG", value = "0.1, 0.2, 0.3, 0.4, 0.5")
        updateTextInput(inputId = "deathRatesLOG", value = "0.001, 0.0008, 0.0006, 0.0004, 0.0002")
        updateTextInput(inputId = "carryingCapacitiesLOG", value = "1000, 1200, 1400, 1600, 1800")
        updateSliderInput(inputId = "errorVarianceLOG", value = 500)
        updateSwitchInput(inputId = "normLOG", value = TRUE)
    })
    
    
    
    ## runLOG ####
    runLOG <- reactive({
    # runLOG <- eventReactive(input$buttonSimulateLOG, {
        simulateStochasticLogistic(
            n.species = n.species.log(),
            names.species = names.species.log(),
            growth.rates = growth.rates.log(),
            carrying.capacities = carrying.capacities.log(),
            death.rates = death.rates.log(),
            x0 = x0.log(),
            sigma.drift = sigma.drift.log(),
            sigma.epoch = sigma.epoch.log(),
            sigma.external = sigma.external.log(),
            sigma.migration = sigma.migration.log(),
            epoch.p = epoch.p.log(),
            t.external_events = t.external_events.log(),
            t.external_durations = t.external_durations.log(),
            migration.p = migration.p.log(),
            metacommunity.probability = metacommunity.probability.log(),
            stochastic = stochastic.log(),
            error.variance = error.variance.log(),
            norm = norm.log(),
            t.start = t.start.log(),
            t.end = t.end.log(),
            t.step = t.step.log(),
            t.store = t.store.log()
        )
    })
    
    output$LOGSpecies <- renderPlot(makePlot(runLOG()$matrix, "abundance of species by time"), res = 96)
}
