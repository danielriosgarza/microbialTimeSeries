# source all files ####
Rfiles = gsub(" ", "", paste("../R/", list.files("../R")))
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

# plotting functions ####
makePlot <- function(out.matrix, title = "abundance of species by time"){
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

makePlotRes <- function(out.matrix, title = "quantity of compounds by time"){
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
        ggtitle(title) # + theme_linedraw()
    fig
}

makeHeatmap <-function(matrix.A, title = "Consumption/production matrix"){
    df = melt(t(matrix.A))
    names(df)<- c("x", "y", "strength")
    df$y <- factor(df$y, levels=rev(unique(sort(df$y))))
    fig <- ggplot(df, aes(x,y,fill=strength)) + geom_tile() + coord_equal() +
        theme(axis.title = element_blank()) + 
        scale_fill_gradient2('strength', low = "red", mid = "white", high = "blue", midpoint = 0)+
        theme_void() + ggtitle(title)
    
    if (ncol(matrix.A)<=10 & nrow(matrix.A)<=10){
        fig <- fig + geom_text(aes(label = round(strength, 2)))
    } else if (ncol(matrix.A)<=15 & nrow(matrix.A)<=15){
        fig <- fig + geom_text(aes(label = round(strength, 1)))
    } else {
        fig <- fig
    }
    
    fig <- fig + labs(x = "compounds", y = "species")+
        theme_linedraw() + 
        theme(plot.title = element_text(hjust = 0.5, size = 14))
    
    if (nrow(matrix.A) >= 20){
        # too many species 
        fig <- fig + theme(
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
        )
    }
    if (ncol(matrix.A) >= 20){
        # too many resources
        fig <- fig + theme(
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()
        )
    }
    fig
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
    
    ## compounds stochiometry ####
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
    
    ## dilution ####
    res.dilu.crm <- reactive(as.numeric(as.vector(text2char(input$resourcesDilutionCRM))))
    resources.dilution.crm <- reactive({
        if (length(res.dilu.crm()) < n.resources.crm()){
            return(c(res.dilu.crm(), resources.crm()[(length(res.dilu.crm())+1):n.resources.crm()]))
        } else {
            return(head(res.dilu.crm(), n.resources.crm()))
        }
    })
    output$resourcesDilutionCRMOutput <- renderPrint(resources.dilution.crm())
    dilution.rate.crm <- reactive(input$dilutionRateCRM)
    
    output$resourcesCRMPlot <- renderPlot(makePiePlot(resources.crm(), label = 'concentration', title='compounds'))
    mean.consumption.crm <- reactive(input$meanConsumptionCRM)
    mean.production.crm <- reactive(input$meanProductionCRM)
    observeEvent(input$meanConsumptionCRM | input$meanProductionCRM, {
        updateSliderInput(inputId = "meanProductionCRM", max = 1 - input$meanConsumptionCRM)
        updateSliderInput(inputId = "meanConsumptionCRM", max = 1 - input$meanProductionCRM)
    })
    maintenance.crm <- reactive(input$maintenanceCRM)
    
    ## editable matrixECRM and matrixMonodCRM ####
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
                maintenance = maintenance.crm()),
            digits = 3)
        RV.crm$matrixECRM <- roundECRM
    })
    output$tableECRM <- renderDataTable(RV.crm$matrixECRM, editable = 'cell', selection = 'none', server = TRUE, options = list(scrollX = TRUE))
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
                nrow = n.species.crm()), 
            digits = 3)
        RV.crm$matrixMonodCRM <- roundMonodCRM
    })
    output$tableMonodCRM <- renderDataTable(RV.crm$matrixMonodCRM, editable = 'cell', selection = 'none', server = TRUE, options = list(scrollX = TRUE))
    observeEvent(input$tableMonodCRM_cell_edit, {
        RV.crm$matrixMonodCRM <<- editData(RV.crm$matrixMonodCRM, input$tableMonodCRM_cell_edit, 'tableMonodCRM')
    })
    ## pertubation ####
    error.variance.crm <- reactive(input$errorVarianceCRM)
    norm.crm <- reactive(input$normCRM)
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
    
    ## examples ####
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
        updateButton(session, "CRMEX4", disabled = !input$CRMEX4pre)
    })
    observeEvent(input$CRMEX4, {
        # auto update matrixECRM, then take changes after it.
        matrixExample4 <- matrix(c(1, -3, 0, 0, 1, 0, -2, 0, 0, 0, 4, -3), nrow = 3, byrow = TRUE)*c(4.3/4, 2/4, 1/4)
        RV.crm$matrixECRM <- matrixExample4
    })
    observeEvent(input$CRMEX5pre, {
        updateSliderInput(inputId = "nSpeciesCRM", value = 10)
        updateSliderInput(inputId = "nResourcesCRM", value = 10)
        updateButton(session, "CRMEX5", disabled = !input$CRMEX5pre)
        
    })
    observeEvent(input$CRMEX5, {
        RV.crm$matrixECRM <- randomE(n.species = 10, n.resources = 10, mean.consumption = 3, mean.production = 1,
                              maintenance = 0.5, trophic.preferences = list(c(5,3,1,1,1,1,1,1,1,1)))
    })
    observeEvent(input$CRMEX6pre, {
        updateSliderInput(inputId = "nSpeciesCRM", value = 20)
        updateSliderInput(inputId = "nResourcesCRM", value = 20)
        updateButton(session, "CRMEX6", disabled = !input$CRMEX6pre)
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
        updateButton(session, "CRMEX7", disabled = !input$CRMEX7pre)
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
    runCRM <- reactive(
        simulateConsumerResource(
            n.species = n.species.crm(),
            n.resources = n.resources.crm(),
            names.species = names.species.crm(), 
            names.resources = names.resources.crm(), 
            E = RV.crm$matrixECRM,
            x0 = x0.crm(), 
            resources = resources.crm(), 
            resources.dilution = resources.dilution.crm(),
            growth.rates = growth.rates.crm(), 
            monod.constant = RV.crm$matrixMonodCRM, 
            dilution.rate = dilution.rate.crm(),
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
            t.end = t.end.crm(),
            t.start = t.start.crm(),
            t.step = t.step.crm(),
            t.store = t.store.crm()
        )
    )
    
    output$CRMSpecies <- renderPlot(makePlot(runCRM()$matrix, "abundance of species by time"), res = 96)
    output$CRMResources <- renderPlot(makePlotRes(runCRM()$resources, "quantity of compounds by time"),  res = 96)
    
    # model2 simulate generalized Lotka-Volterra Model ####
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
    
    observeEvent(input$buttonInteractionsU, {
        updateTextInput(inputId = "interactionsGLV", value = rbeta(n.speciesGLV()^2, 0.5, 0.5))
    })
    observeEvent(input$buttonInteractionsN, {
        updateTextInput(inputId = "interactionsGLV", value = rbeta(n.speciesGLV()^2, 2, 2))
    })
    observeEvent(input$buttonInteractionsEven, {
        updateTextInput(inputId = "interactionsGLV", value = runif(n.speciesGLV()^2, 0, 1))
    })
    
    interactionsGLV <- reactive(as.numeric(text2char(input$interactionsGLV)))
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
            symmetric = symmetricGLV()
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
