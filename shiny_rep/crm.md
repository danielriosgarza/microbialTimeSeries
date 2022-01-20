---
title: "Consumer Resource Model"
author: "Daniel Garza"
date: "1/19/2022"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## simulateConsumerResource

Function to simulate McArthur's consumer-resource model and extensions for microbiomes. The model describes the explicit interaction between microbes in a community via the consumption and production of metabolites.

#### **Model Description**

Each microbial species from a pool of $n$ different species has preferences for specific subset of the $k$ different substrates that are available to the community. Their growth rates depend on the availability of these resources, through the following differential equation:

$\frac{dX_i}{dt} = \mu_{i}X_i (\sum_{j=1}^{k} e_{i,j} f_{i,j} S_j - \delta_i)$

where, where $\mu_{i}$ is the maximum growth rate of species $i$, while $X_i$ is its abundance at time $t$.

$e_{i,j}$ is the yield of substrate $j$ for species $i$, $f_{i,j}$ is the feeding form, which constraints the quantity of the substrate that can be consumed at each time step, while $S_j$ is the concentration of the substrate $j$ at time $t$, and $\delta_i$ is a dilution term.

The form of the feeding term is the Monod equation $f_{i,j}=\frac{S_j}{k_{i,j}+S_j}$, where $k_{i,j}$ is the Monod constant.

Resources change according to the following differential equation:

$\frac{dS_j}{dt}=\phi_j - \sum_{i=1}^{n} e_{i,j}f_{i,j}X_i + \sum_{i=1}^{n}p_{i,j} X_{i} (\sum_{j=1}^{k}e_{i,j}f_{i,j}-\delta_i)$

where $\phi_j$ is a dilution term and $p_{i,j}$ is the yield of production of substrate $j$ that results from the growth species $i$.

In practice, a species either consumes or produces a metabolite (or is indifferent to its presence). Allowing us to summarize $e_{i,j}$ and $p_{i,j}$ in a single matrix ($E$, containing $n$ rows and $k$ columns). To distinguish them, the entries of $E$ are, respectively, positive and negative for the consumed and produced substrates.

There are many possibilities for structuring $E$ according to specific assumptions. We will later summarize the options that are built into miaSim. We will first jump to some quick examples of how to simulate a consumer resource model with miaSim.

### **Examples**

To illustrate the basic parameters, we begin with a model of a single species that grows in an environment with five independent substrates. Some of the examples below are compatible with the Shiny app, click on the button to explore them with the app.

If not provided by the user, all parameters have defaults except for the the number of species and the number of resources. Check the list of parameters below for a complete description of the parameters and their defaults

```{n.species <- 1}
n.resources <- 5


#simulate the model
CRMsimul <- simulateConsumerResource(n.species = n.species, n.resources = n.resources)

#visualize the result
makePlot(CRMsimul$matrix) #species plot
makePlot(CRMsimul$resources) #resources plot
```

MiaSim provides a helper function to generate matrix the matrix $E$. Below is a simulation of the same model where the user has more control over the parameters.

```{#generate the matrix E}
E = randomE(n.species = n.species, n.resources = n.resources)
print(E)

#positive entries are consumed, negative produced, zero has no influence on the species

#define some simulation parameters
t.end = 2000 #when to stop
t.store = 500 #how many samples of the simulation to store (evenly spaced)
migration.p = 0 #whether to allow migration from a metacommunity
stochastic = 0 # whether to use noise in the simulation
dilution.rate = 0.01 #adding a dilution rate

#simulate the model
CRMsimul <- simulateConsumerResource(n.species = n.species, n.resources = n.resources, stochastic = stochastic, migration.p = migration.p, E=E, t.end = t.end, t.store = t.store, dilution.rate=dilution.rate)

#visualize the result
makePlot(CRMsimul$matrix) #species plot
makePlot(CRMsimul$resources) #resources plot



```

Click on the example button to

## .
