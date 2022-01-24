library(shiny)
library(shinyBS)
library(ggplot2)
library(deSolve)
library(reshape2)
library(gtools)
library(DT)
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
    
    if (ncol(matrix.A)<=12 & nrow(matrix.A)<=12){
        fig <- fig + geom_text(aes(label = round(strength, 2)))
    } else if (ncol(matrix.A)<=20 & nrow(matrix.A)<=20){
        fig <- fig + geom_text(aes(label = round(strength, 1)))
    }
    
    fig <- fig + labs(x = "compounds", y = "species")+
        theme_linedraw() + 
        theme(plot.title = element_text(hjust = 0.5, size = 14))
    
    if (nrow(matrix.A) > 20){
        # too many species 
        fig <- fig + theme(
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
        )
    }
    if (ncol(matrix.A) > 20){
        # too many resources
        fig <- fig + theme(
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()
        )
    }
    fig
}

#### ui ####
source("ui.R")

#### server ####
server <- function(input, output, session) {
    
    #### model1 simulate consumer resource model ####
    ## basic
    n.species <- reactive(input$nSpecies)
    n.resources <- reactive(input$nResources)
    names.species <- reactive(text2chars(input$namesSpecies, len = n.species(), prefix = "sp"))
    names.resources <- reactive(text2chars(input$namesResources, len = n.resources(), prefix = "res"))
    t.end <- reactive(input$tEnd)
    ## compounds stochiometry
    # observeEvent(input$nResources, {
    #     resources <- reactive(as.numeric(text2chars(input$resources, len = n.resources(), expr = paste0("runif(n = ", n.resources() ,", min = 1, max = 100)"))))
    # })
    # resources_custom <- reactive(as.numeric(text2chars(input$resources_custom, len = n.resources(), expr = paste0("runif(n = ", n.resources() ,", min = 1, max = 100)"))))
    res.conc <- reactive(input$resourcesConcentration)
    res.even <- reactive(input$resourcesEvenness)
    resources_dist <- reactive(rdirichlet(1, rep(1, n.resources())*res.even())*res.conc()*n.resources())
    res.custom <- reactive(as.numeric(as.vector(text2char(input$resourcesCustom))))
    resources <- reactive({
        if (length(res.custom()) < n.resources()){
            return(c(res.custom(), as.vector(resources_dist())[(length(res.custom())+1):n.resources()]))
        } else {
            return(head(res.custom(), n.resources()))
        }
    })
    output$resourcesOutput <- renderPrint(resources())
    
    # dilution
    res.dilu <- reactive(as.numeric(as.vector(text2char(input$resourcesDilution))))
    resources.dilution <- reactive({
        if (length(res.dilu()) < n.resources()){
            return(c(res.dilu(), resources()[(length(res.dilu())+1):n.resources()]))
        } else {
            return(head(res.dilu(), n.resources()))
        }
    })
    output$resourcesDilutionOutput <- renderPrint(resources.dilution())
    dilution.rate <- reactive(input$dilutionRate)
    
    output$resourcesPlot <- renderPlot(makePiePlot(resources(), label = 'concentration', title='compounds'))
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
    
    output$CRMPlotE <- renderPlot(makeHeatmap(E$df, 'Consumption/production matrix'), res = 96)
    
    observeEvent(input$tableE_cell_edit, {
        E$df <<- editData(E$df, input$tableE_cell_edit, 'tableE')
    })
    
    
    ## growth rates
    x0 <- reactive(as.numeric(text2chars(input$x0, len = n.species(), expr = paste0("runif(n = ", n.species() ,", min = 0.1, max = 10)"))))
    output$x0Output <- renderPrint(x0())
    
    
    alpha <- reactive(input$alpha)
    beta <- reactive(input$beta)
    # listening to the changes in n.species, alpha, and beta
    observeEvent(input$nSpecies | input$alpha | input$beta, {
        growth.rates <- reactive(as.numeric(text2chars(input$growthRates, len = n.species(), expr = paste0('round(rbeta(',n.species(), ',' ,alpha(), ',' , beta(),'), digits = 3)'))))
    })
    observeEvent(input$nSpecies, {
        x0 <- reactive(as.numeric(text2chars(input$x0, len = n.species(), expr = paste0("runif(n = ", n.species() ,", min = 0.1, max = 10)"))))
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
    
    ## pertubation
    error.variance <- reactive(input$errorVariance)
    norm <- reactive(input$norm)
    stochastic <- reactive(input$stochastic)
    sigma.drift <- reactive(input$sigmaDrift)
    sigma.epoch <- reactive(input$sigmaEpoch)
    sigma.external <- reactive(input$sigmaExternal)
    sigma.migration <- reactive(input$sigmaMigration)
    epoch.p <- reactive(input$epochP)
    t.external_events <- reactive(as.numeric(text2char(input$tExternalEvents)))
    output$tExternalEventsOutput <- renderPrint(t.external_events())
    t.external_durations <- reactive(as.numeric(text2char(input$tExternalDurations)))
    output$tExternalDurationsOutput <- renderPrint(t.external_durations())
    migration.p <- reactive(input$migrationP)
    metacommunity.probability <- reactive(as.numeric(text2chars(input$metacommunityProbability, len = n.species(), expr = paste0("rdirichlet(1, alpha = rep(1,", n.species(), "))"))))
    output$metacommunityProbability <- renderPrint(metacommunity.probability())
    runCRM <- reactive(
        simulateConsumerResource(
            n.species = n.species(),
            n.resources = n.resources(),
            names.species = names.species(), 
            names.resources = names.resources(), 
            E = E$df,
            x0 = x0(), 
            resources = resources(), 
            resources.dilution = resources.dilution(),
            growth.rates = growth.rates(), 
            monod.constant = monod.constant(), 
            dilution.rate = dilution.rate(),
            sigma.drift = sigma.drift(),
            sigma.epoch = sigma.epoch(),
            sigma.external = sigma.external(),
            sigma.migration = sigma.migration(), 
            epoch.p = epoch.p(), 
            t.external_events = t.external_events(),
            t.external_durations = t.external_durations(),
            stochastic = stochastic(),
            migration.p = migration.p(),
            metacommunity.probability = metacommunity.probability(),
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